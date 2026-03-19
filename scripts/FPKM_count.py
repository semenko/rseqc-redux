#!/usr/bin/env python

"""
calculate raw read count, FPM (fragment per million) and FPKM (fragment per million mapped
reads per kilobase exon) for each gene in BED file.
"""

import sys

from bx.intervals import Intersecter, Interval

from rseqc import SAM
from rseqc.cli_common import (
    add_input_bam_arg,
    add_mapq_arg,
    add_output_prefix_arg,
    add_refgene_arg,
    create_parser,
    iter_bed12,
    validate_bam_index,
    validate_files_exist,
)
from rseqc.SAM import _parse_strand_rule, _pysam_iter


def build_range(refgene: str) -> dict:
    """build ranges for exonic region"""
    ranges = {}
    for record in iter_bed12(refgene):
        chrom = record.chrom.upper()
        for st, end in zip(record.exon_starts, record.exon_ends):
            if chrom not in ranges:
                ranges[chrom] = Intersecter()
            ranges[chrom].add_interval(Interval(st, end))
    return ranges


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser, help="Alignment file in BAM format (SAM is not supported). [required]")
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    add_refgene_arg(parser, dest="refgene_bed")
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "How read(s) were stranded during sequencing. For example:"
            " --strand='1++,1--,2+-,2-+' means that this is a pair-end,"
            " strand-specific RNA-seq, and the strand rule is: read1"
            " mapped to '+' => parental gene on '+'; read1 mapped to"
            " '-' => parental gene on '-'; read2 mapped to '+' =>"
            " parental gene on '-'; read2 mapped to '-' => parental"
            " gene on '+'. If you are not sure about the strand rule,"
            " run 'infer_experiment.py' default=%(default)s"
            " (Not a strand specific RNA-seq data)"
        ),
    )
    parser.add_argument(
        "-u",
        "--skip-multi-hits",
        action="store_true",
        dest="skip_multi",
        help=("How to deal with multiple hit reads. Presence this option renders program to skip multiple hits reads."),
    )
    parser.add_argument(
        "-e",
        "--only-exonic",
        action="store_true",
        dest="only_exon",
        help=(
            "How to count total reads. Presence of this option renders"
            " program only used exonic (UTR exons and CDS exons) reads,"
            " otherwise use all reads."
        ),
    )
    add_mapq_arg(parser)
    parser.add_argument(
        "-s",
        "--single-read",
        type=float,
        dest="single_read",
        default=1,
        help=(
            "How to count read-pairs that only have one end mapped."
            " 0: ignore it. 0.5: treat it as half fragment."
            " 1: treat it as whole fragment. default=%(default)s"
        ),
    )

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file and args.refgene_bed):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file, args.refgene_bed)
    validate_bam_index(args.input_file)

    obj = SAM.ParseBAM(args.input_file)
    with open(args.output_prefix + ".FPKM.xls", "w") as OUT:
        strandRule = _parse_strand_rule(args.strand_rule)

        # ++++++++++++++++++++++++++++++++++++counting fragments
        print("Extract exon regions from  " + args.refgene_bed + "...", file=sys.stderr)
        gene_ranges = build_range(args.refgene_bed)
        print("Counting total fragment ... ", end=" ", file=sys.stderr)

        total_frags = 0.0
        exonic_frags = 0.0

        for aligned_read in _pysam_iter(obj.samfile):
            if aligned_read.is_qcfail:
                continue  # skip low quanlity
            if aligned_read.is_duplicate:
                continue  # skip duplicate read
            if aligned_read.is_secondary:
                continue  # skip non primary hit
            if args.skip_multi:
                if aligned_read.mapq < args.map_qual:
                    continue
            try:
                chrom = obj.samfile.getrname(aligned_read.tid).upper()
            except (KeyError, ValueError):
                continue
            read_st = aligned_read.pos
            read_end = read_st + aligned_read.rlen  # not exactly the end position in case of splicing, insertion,etc

            if not aligned_read.is_paired:  # if read is NOT paired in sequencing (single-end sequencing)
                total_frags += 1
                if (chrom in gene_ranges) and gene_ranges[chrom].find(read_st, read_end):
                    exonic_frags += 1
            else:  # for pair-end sequencing
                if aligned_read.is_read2:
                    continue  # only count read1
                mate_st = aligned_read.pnext
                mate_end = mate_st + aligned_read.rlen

                if aligned_read.is_unmapped:  # read1 unmapped
                    if aligned_read.mate_is_unmapped:
                        continue  # both unmap
                    else:  # read2 is mapped
                        total_frags += args.single_read
                        if (chrom in gene_ranges) and gene_ranges[chrom].find(mate_st, mate_end):
                            exonic_frags += args.single_read
                else:
                    if aligned_read.mate_is_unmapped:
                        total_frags += args.single_read
                        if (chrom in gene_ranges) and gene_ranges[chrom].find(read_st, read_end):
                            exonic_frags += args.single_read
                    else:
                        total_frags += 1
                        if (
                            (chrom in gene_ranges)
                            and gene_ranges[chrom].find(read_st, read_end)
                            and gene_ranges[chrom].find(mate_st, mate_end)
                        ):
                            exonic_frags += 1

        print("Done", file=sys.stderr)
        print(f"Total fragment = {str(total_frags):<20}", file=sys.stderr)
        print(f"Total exonic fragment = {str(exonic_frags):<20}", file=sys.stderr)

        if total_frags > 0 and exonic_frags > 0:
            if args.only_exon:
                denominator = exonic_frags
            else:
                denominator = total_frags
        else:
            print("Total tags cannot be 0 or negative number", file=sys.stderr)
            sys.exit(1)

        # ++++++++++++++++++++++++++++++++++++++++++++++++
        obj = SAM.ParseBAM(args.input_file)
        print(
            "\t".join(("#chrom", "st", "end", "accession", "mRNA_size", "gene_strand", "Frag_count", "FPM", "FPKM")),
            file=OUT,
        )

        gene_finished = 0

        # calculate raw count, FPM, FPKM for each gene
        for record in iter_bed12(args.refgene_bed):
            frag_count_f = 0.0
            frag_count_r = 0.0
            frag_count_fr = 0.0
            mRNA_size = 0.0
            exon_ranges = Intersecter()
            chrom = record.chrom
            tx_start = record.tx_start
            tx_end = record.tx_end
            geneName = record.gene_name
            gstrand = record.strand.replace(" ", "_")

            for st, end in zip(record.exon_starts, record.exon_ends):
                mRNA_size += end - st
                exon_ranges.add_interval(Interval(st, end))

            # extract reads mapped gene region
            try:
                alignedReads = obj.samfile.fetch(chrom, tx_start, tx_end)
            except (KeyError, ValueError):
                continue
            for aligned_read in _pysam_iter(alignedReads):
                if aligned_read.is_qcfail:
                    continue  # skip low quanlity
                if aligned_read.is_duplicate:
                    continue  # skip duplicate read
                if aligned_read.is_secondary:
                    continue  # skip non primary hit
                if args.skip_multi:
                    if aligned_read.mapq < args.map_qual:
                        continue

                # single end sequencing
                if not aligned_read.is_paired:
                    frag_st = aligned_read.pos
                    frag_end = (
                        frag_st + aligned_read.rlen
                    )  # not exactly the end position in case of splicing, insertion,etc
                    if aligned_read.is_reverse:
                        strand_key = "-"
                    else:
                        strand_key = "+"

                    if exon_ranges.find(frag_st, frag_end):
                        if args.strand_rule is None:
                            frag_count_fr += 1
                        elif strand_key in strandRule and strandRule[strand_key] == "+":
                            frag_count_f += 1
                        elif strand_key in strandRule and strandRule[strand_key] == "-":
                            frag_count_r += 1

                # pair-end sequencing
                else:
                    frag_st = aligned_read.pos
                    frag_end = aligned_read.pnext
                    if not exon_ranges.find(frag_st, frag_st + 1) and not exon_ranges.find(frag_end, frag_end + 1):
                        continue
                    if aligned_read.is_read2:
                        continue
                    if aligned_read.is_reverse:
                        strand_key = "1-"
                    else:
                        strand_key = "1+"

                    if args.strand_rule is None:
                        if aligned_read.is_unmapped:
                            if aligned_read.mate_is_unmapped:  # both unmapped
                                continue
                            else:  # only read2 mapped
                                frag_count_fr += args.single_read
                        else:
                            if aligned_read.mate_is_unmapped:  # only read1 mapped
                                frag_count_fr += args.single_read
                            else:  # both mapped
                                frag_count_fr += 1
                    else:
                        if strand_key in strandRule and strandRule[strand_key] == "+":
                            if aligned_read.is_unmapped:
                                if aligned_read.mate_is_unmapped:  # both unmapped
                                    continue
                                else:  # only read2 mapped
                                    frag_count_f += args.single_read
                            else:
                                if aligned_read.mate_is_unmapped:  # only read1 mapped
                                    frag_count_f += args.single_read
                                else:  # both mapped
                                    frag_count_f += 1
                        if strand_key in strandRule and strandRule[strand_key] == "-":
                            if aligned_read.is_unmapped:
                                if aligned_read.mate_is_unmapped:  # both unmapped
                                    continue
                                else:  # only read2 mapped
                                    frag_count_r += args.single_read
                            else:
                                if aligned_read.mate_is_unmapped:  # only read1 mapped
                                    frag_count_r += args.single_read
                                else:  # both mapped
                                    frag_count_r += 1

            FPM_fr = frag_count_fr * 1000000 / denominator
            FPM_f = frag_count_f * 1000000 / denominator
            FPM_r = frag_count_r * 1000000 / denominator
            FPKM_fr = frag_count_fr * 1000000000 / (denominator * mRNA_size)
            FPKM_f = frag_count_f * 1000000000 / (denominator * mRNA_size)
            FPKM_r = frag_count_r * 1000000000 / (denominator * mRNA_size)

            if args.strand_rule is None:
                frag_count, fpm, fpkm = frag_count_fr, FPM_fr, FPKM_fr
            elif gstrand == "+":
                frag_count, fpm, fpkm = frag_count_f, FPM_f, FPKM_f
            elif gstrand == "-":
                frag_count, fpm, fpkm = frag_count_r, FPM_r, FPKM_r
            else:
                frag_count = None

            if frag_count is not None:
                print(
                    "\t".join(
                        str(i) for i in (chrom, tx_start, tx_end, geneName, mRNA_size, gstrand, frag_count, fpm, fpkm)
                    ),
                    file=OUT,
                )

            gene_finished += 1
            print(f" {gene_finished} transcripts finished\r", end=" ", file=sys.stderr)


if __name__ == "__main__":
    main()
