#!/usr/bin/env python
"""
calculate fragment size for each gene/transcript. For each transcript/gene, it Will report:
1) # of fragment that was used.
2) mean of fragment size
3) median of fragment size
4) stdev of fragment size
"""

import argparse
import os
import sys

import pysam
from numpy import mean, median, std


def overlap_length2(lst1, lst2):
    overlap_len = 0
    for x in lst1:
        for y in lst2:
            overlap_len += len(list(range(max(x[0], y[0]), min(x[-1], y[-1]) + 1)))
    return overlap_len


def fragment_size(bedfile, samfile, qcut=30, ncut=5):
    """calculate the fragment size for each gene"""
    with open(bedfile, "r") as _fh:
        for line in _fh:
            exon_range = []
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            chrom = fields[0]
            tx_start = int(fields[1])
            tx_end = int(fields[2])
            geneName = fields[3]
            exon_starts = list(map(int, fields[11].rstrip(",\n").split(",")))
            exon_starts = [x + tx_start for x in exon_starts]
            exon_ends = list(map(int, fields[10].rstrip(",\n").split(",")))
            exon_ends = [x + y for x, y in zip(exon_starts, exon_ends)]
            geneID = "\t".join([str(i) for i in (chrom, tx_start, tx_end, geneName)])

            for st, end in zip(exon_starts, exon_ends):
                exon_range.append([st + 1, end + 1])
                # exon_range.append([chrom, st,end])

            try:
                alignedReads = samfile.fetch(chrom, tx_start, tx_end)
            except (KeyError, ValueError):
                yield "\t".join([str(i) for i in (geneID, 0, 0, 0)])
                continue

            frag_sizes = []
            for aligned_read in alignedReads:
                if not aligned_read.is_paired:  # skip single sequencing
                    continue
                if aligned_read.is_read2:
                    continue
                if aligned_read.mate_is_unmapped:
                    continue
                if aligned_read.is_qcfail:
                    continue  # skip low quanlity
                if aligned_read.is_duplicate:
                    continue  # skip duplicate read
                if aligned_read.is_secondary:
                    continue  # skip non primary hit
                if aligned_read.mapq < qcut:
                    continue

                read_st = aligned_read.pos
                mate_st = aligned_read.pnext
                if read_st > mate_st:
                    (read_st, mate_st) = (mate_st, read_st)
                if read_st < tx_start or mate_st > tx_end:
                    continue
                read_len = aligned_read.qlen
                map_range = [[read_st + 1, mate_st]]
                # map_range = [[chrom, read_st, mate_st]]
                frag_len = overlap_length2(exon_range, map_range) + read_len
                frag_sizes.append(frag_len)
            if len(frag_sizes) < ncut:
                yield "\t".join([str(i) for i in (geneID, len(frag_sizes), 0, 0, 0)])
            else:
                yield "\t".join(
                    [str(i) for i in (geneID, len(frag_sizes), mean(frag_sizes), median(frag_sizes), std(frag_sizes))]
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument("-i", "--input", dest="input_file", help="Input BAM file")
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        help="Reference gene model in BED format. Must be strandard 12-column BED file. [required]",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help=(
            "Minimum mapping quality (phred scaled) for an alignment"
            ' to be called "uniquely mapped". default=%(default)s'
        ),
    )
    parser.add_argument(
        "-n",
        "--frag-num",
        type=int,
        dest="fragment_num",
        default=3,
        help="Minimum number of fragment. default=%(default)s",
    )

    args = parser.parse_args()

    if not (args.input_file and args.refgene_bed):
        parser.print_help()
        sys.exit(1)
    if not os.path.exists(args.input_file + ".bai"):
        print("cannot find index file of input BAM file", file=sys.stderr)
        print(args.input_file + ".bai" + " does not exists", file=sys.stderr)
        sys.exit(1)

    for file in (args.input_file, args.refgene_bed):
        if not os.path.exists(file):
            print(file + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(1)

    print(
        "\t".join(
            [
                str(i)
                for i in ("chrom", "tx_start", "tx_end", "symbol", "frag_count", "frag_mean", "frag_median", "frag_std")
            ]
        )
    )
    for tmp in fragment_size(args.refgene_bed, pysam.Samfile(args.input_file), args.map_qual, args.fragment_num):
        print(tmp)


if __name__ == "__main__":
    main()
