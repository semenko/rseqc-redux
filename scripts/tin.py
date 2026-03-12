#!/usr/bin/env python
"""Calculate transcript integrity number (TIN) for each transcript in a BED file."""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections.abc import Generator

import pysam
from bx.intervals import Intersecter, Interval
from numpy import mean, median, std

from rseqc import getBamFiles
from rseqc.cli_common import printlog
from rseqc.SAM import _pysam_iter


def uniqify(seq: list) -> list:
    """
    duplicated members only keep one copy. [1,2,2,3,3,4] => [1,2,3,4].
    """
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]


def shannon_entropy(arg: list[float]) -> float:
    """
    calculate shannon's H = -sum(P*log(P)). arg is a list of float numbers. Note we used
    natural log here.
    """
    lst_sum = sum(arg)
    entropy = 0.0
    for i in arg:
        entropy += (i / lst_sum) * math.log(i / lst_sum)
    if entropy == 0:
        return 0
    else:
        return -entropy


def build_bitsets(list: list) -> dict:
    """
    build intevalTree from list
    """
    ranges = {}
    for entry in list:
        chrom = entry[0]
        st = entry[1]
        end = entry[2]
        if chrom not in ranges:
            ranges[chrom] = Intersecter()
        ranges[chrom].add_interval(Interval(st, end))
    return ranges


def union_exons(refbed: str) -> dict:
    """
    take the union of all exons defined in refbed file and build bitset
    """
    from rseqc import BED

    tmp = BED.ParseBED(refbed)
    all_exons = tmp.getExon()
    unioned_exons = BED.unionBed3(all_exons)
    exon_ranges = build_bitsets(unioned_exons)
    return exon_ranges


def estimate_bg_noise(chrom: str, tx_st: int, tx_end: int, samfile: pysam.AlignmentFile, e_ranges: dict) -> float:
    """
    estimate background noise level for a particular transcript
    """
    intron_sig = 0.0  # reads_num * reads_len
    alignedReads = samfile.fetch(chrom, tx_st, tx_end)
    for aligned_read in _pysam_iter(alignedReads):
        if aligned_read.is_qcfail:
            continue
        if aligned_read.is_unmapped:
            continue
        if aligned_read.is_secondary:
            continue
        read_start = aligned_read.pos
        if read_start < tx_st:
            continue
        if read_start >= tx_end:
            continue
        read_len = aligned_read.qlen
        if len(e_ranges[chrom].find(read_start, read_start + read_len)) > 0:
            continue
        intron_sig += read_len
    return intron_sig


def genomic_positions(refbed: str, sample_size: int) -> Generator:
    """
    return genomic positions of each nucleotide in mRNA. sample_size: Number of nucleotide
    positions sampled from mRNA.
    """
    if refbed is None:
        print("You must specify a bed file representing gene model\n", file=sys.stderr)
        sys.exit(1)

    with open(refbed, "r") as _fh:
        for line in _fh:
            try:
                if line.startswith(("#", "track", "browser")):
                    continue
                # Parse fields from gene tabls
                fields = line.split()
                chrom = fields[0]
                tx_start = int(fields[1])
                tx_end = int(fields[2])
                geneName = fields[3]
                mRNA_size = sum(int(i) for i in fields[10].strip(",").split(","))

                exon_starts = [int(x) for x in fields[11].rstrip(",\n").split(",")]
                exon_starts = [x + tx_start for x in exon_starts]
                exon_ends = [int(x) for x in fields[10].rstrip(",\n").split(",")]
                exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
                intron_size = tx_end - tx_start - mRNA_size
                if intron_size < 0:
                    intron_size = 0
            except (IndexError, ValueError):
                print("[NOTE:input bed must be 12-column] skipped this line: " + line, end=" ", file=sys.stderr)
                continue

            chose_bases = [tx_start + 1, tx_end]
            exon_bounds = []
            gene_all_base = []
            if mRNA_size <= sample_size:  # return all bases of mRNA
                for st, end in zip(exon_starts, exon_ends):
                    chose_bases.extend(
                        list(range(st + 1, end + 1))
                    )  # 1-based coordinates on genome, include exon boundaries
                yield (geneName, chrom, tx_start, tx_end, intron_size, chose_bases)
            elif mRNA_size > sample_size:
                step_size = int(mRNA_size / sample_size)
                for st, end in zip(exon_starts, exon_ends):
                    gene_all_base.extend(list(range(st + 1, end + 1)))
                    exon_bounds.append(st + 1)
                    exon_bounds.append(end)
                indx = list(range(0, len(gene_all_base), step_size))
                chose_bases = [gene_all_base[i] for i in indx]
                yield (geneName, chrom, tx_start, tx_end, intron_size, uniqify(exon_bounds + chose_bases))


def check_min_reads(samfile: pysam.AlignmentFile, chrom: str, tx_st: int, tx_end: int, cutoff: int) -> bool:
    """
    make sure the gene has minimum reads coverage. if cutoff = 10, each gene must have
    10 *different* reads.
    """
    tmp = False
    read_count = set()
    try:
        alignedReads = samfile.fetch(chrom, tx_st, tx_end)
        for aligned_read in _pysam_iter(alignedReads):
            if aligned_read.is_qcfail:
                continue
            if aligned_read.is_unmapped:
                continue
            if aligned_read.is_secondary:
                continue
            read_start = aligned_read.pos
            if read_start < tx_st:
                continue
            if read_start >= tx_end:
                continue
            read_count.add(read_start)
            if len(read_count) > cutoff:  # no need to loop anymore
                tmp = True
                break
        return tmp
    except (KeyError, ValueError):
        return False


def genebody_coverage(
    samfile: pysam.AlignmentFile, chrom: str, positions: list[int], bg_level: float = 0
) -> list[float]:
    """
    calculate coverage for each nucleotide in *positions*. some times len(cvg) < len(positions)
    because positions where there is no mapped reads were ignored.
    """
    cvg = []
    start = positions[0] - 1
    end = positions[-1]

    try:
        for pileupcolumn in _pysam_iter(samfile.pileup(chrom, start, end, truncate=True)):
            ref_pos = pileupcolumn.pos + 1
            if ref_pos not in positions:
                continue
            if pileupcolumn.n == 0:
                cvg.append(0.0)
                continue
            cover_read = 0.0
            for pileupread in pileupcolumn.pileups:
                if pileupread.is_del:
                    continue
                if pileupread.alignment.is_qcfail:
                    continue
                if pileupread.alignment.is_secondary:
                    continue
                if pileupread.alignment.is_unmapped:
                    continue
                cover_read += 1.0
            cvg.append(cover_read)
    except (KeyError, ValueError):
        cvg = []

    if bg_level <= 0:
        return cvg
    else:
        tmp = []
        for i in cvg:
            subtracted_sig = int(i - bg_level)
            if subtracted_sig > 0:
                tmp.append(subtracted_sig)
            else:
                tmp.append(0)
        return tmp


def tin_score(cvg: list[float], length: int) -> float:
    """calcualte TIN score"""

    if len(cvg) == 0:
        tin = 0
        return tin

    cvg_eff = [float(i) for i in cvg if float(i) > 0]  # remove positions with 0 read coverage
    entropy = shannon_entropy(cvg_eff)

    tin = 100 * (math.exp(entropy)) / length
    return tin


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input",
        dest="input_files",
        help=(
            'Input BAM file(s). "-i" takes these input:'
            " 1) a single BAM file."
            ' 2) "," separated BAM files (no spaces allowed).'
            " 3) directory containing one or more bam files."
            " 4) plain text file containing the path of one or more"
            " bam files (Each row is a BAM file path). All BAM files"
            " should be sorted and indexed using samtools. [required]"
        ),
    )
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]",
    )
    parser.add_argument(
        "-c",
        "--minCov",
        type=int,
        dest="minimum_coverage",
        default=10,
        help="Minimum number of read mapped to a transcript. default=%(default)s",
    )
    parser.add_argument(
        "-n",
        "--sample-size",
        type=int,
        dest="sample_size",
        default=100,
        help=(
            "Number of equal-spaced nucleotide positions picked from"
            " mRNA. Note: if this number is larger than the length of"
            " mRNA (L), it will be halved until it's smaller than L."
            " default=%(default)s"
        ),
    )
    parser.add_argument(
        "-s",
        "--subtract-background",
        action="store_true",
        dest="subtract_bg",
        help=(
            "Subtract background noise (estimated from intronic"
            " reads). Only use this option if there are substantial"
            " intronic reads."
        ),
    )
    args = parser.parse_args()

    # if '-s' was set
    if args.subtract_bg:
        exon_ranges = union_exons(args.ref_gene_model)

    if args.sample_size < 0:
        print("Number of nucleotide can't be negative", file=sys.stderr)
        sys.exit(1)
    elif args.sample_size > 1000:
        print("Warning: '-n' is too large! Please try smaller '-n' valeu if program is running slow.", file=sys.stderr)

    if not (args.input_files and args.ref_gene_model):
        parser.print_help()
        sys.exit(1)

    if not os.path.exists(args.ref_gene_model):
        print("\n\n" + args.ref_gene_model + " does NOT exists" + "\n", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    printlog("Get BAM file(s) ...")
    bamfiles = sorted(getBamFiles.get_bam_files(args.input_files))

    if len(bamfiles) <= 0:
        print("No BAM file found, exit.", file=sys.stderr)
        sys.exit(1)
    else:
        print("Total %d BAM file(s):" % len(bamfiles), file=sys.stderr)
        for f in bamfiles:
            print("\t" + f, file=sys.stderr)

    for f in bamfiles:
        printlog("Processing " + f)

        with (
            open(os.path.basename(f).replace("bam", "") + "summary.txt", "w") as SUM,
            open(os.path.basename(f).replace("bam", "") + "tin.xls", "w") as OUT,
        ):
            print("\t".join(["Bam_file", "TIN(mean)", "TIN(median)", "TIN(stdev)"]), file=SUM)
            print("\t".join(["geneID", "chrom", "tx_start", "tx_end", "TIN"]), file=OUT)

            samfile = pysam.AlignmentFile(f, "rb")
            sample_TINs = []  # sample level TIN, values are from different genes
            finish = 0
            noise_level = 0.0
            for gname, i_chr, i_tx_start, i_tx_end, intron_size, pick_positions in genomic_positions(
                refbed=args.ref_gene_model, sample_size=args.sample_size
            ):
                finish += 1

                # check minimum reads coverage
                if check_min_reads(samfile, i_chr, i_tx_start, i_tx_end, args.minimum_coverage) is not True:
                    print("\t".join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, 0.0)]), file=OUT)
                    continue

                # estimate background noise if '-s' was specified
                if args.subtract_bg:
                    intron_signals = estimate_bg_noise(i_chr, i_tx_start, i_tx_end, samfile, exon_ranges)
                    if intron_size > 0:
                        noise_level = intron_signals / intron_size

                coverage = genebody_coverage(samfile, i_chr, sorted(pick_positions), noise_level)

                tin1 = tin_score(cvg=coverage, length=len(pick_positions))
                sample_TINs.append(tin1)
                print("\t".join([str(i) for i in (gname, i_chr, i_tx_start, i_tx_end, tin1)]), file=OUT)
                print(" %d transcripts finished\r" % (finish), end=" ", file=sys.stderr)

            print(
                "\t".join(
                    [str(i) for i in (os.path.basename(f), mean(sample_TINs), median(sample_TINs), std(sample_TINs))]
                ),
                file=SUM,
            )
            samfile.close()


if __name__ == "__main__":
    main()
