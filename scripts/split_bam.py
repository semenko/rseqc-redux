#!/usr/bin/env python
"""Split BAM file according to an input gene list (BED format)."""

import sys

import pysam

from rseqc import BED
from rseqc.cli_common import (
    add_input_bam_arg,
    add_output_prefix_arg,
    build_bitsets,
    create_parser,
    validate_files_exist,
)
from rseqc.SAM import _pysam_iter


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser, help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed.")
    parser.add_argument(
        "-r",
        "--genelist",
        dest="gene_list",
        help=(
            "Gene list in bed format. All reads hits to exon regions"
            " (defined by this gene list) will be saved into one BAM"
            " file, the remaining reads will saved into another BAM"
            " file."
        ),
    )
    add_output_prefix_arg(
        parser,
        help=(
            'Prefix of output BAM files. "prefix.in.bam" file'
            " contains reads mapped to the gene list specified by"
            ' "-r", "prefix.ex.bam" contains reads that cannot'
            ' mapped to gene list. "prefix.junk.bam" contains'
            " qcfailed or unmapped reads."
        ),
    )
    args = parser.parse_args()

    if not (args.input_file and args.gene_list):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.gene_list, args.input_file)

    # build bitset for gene list
    print("reading " + args.gene_list + " ... ", end=" ", file=sys.stderr)
    obj = BED.ParseBED(args.gene_list)
    exons = obj.getExon()
    exon_ranges = build_bitsets(exons)
    print("Done", file=sys.stderr)

    samfile = pysam.AlignmentFile(args.input_file, "rb")
    out1 = pysam.AlignmentFile(
        args.output_prefix + ".in.bam", "wb", template=samfile
    )  # bam file containing reads hit to exon region
    out2 = pysam.AlignmentFile(
        args.output_prefix + ".ex.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region
    out3 = pysam.AlignmentFile(
        args.output_prefix + ".junk.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region

    total_alignment = 0
    in_alignment = 0
    ex_alignment = 0
    bad_alignment = 0
    print("spliting " + args.input_file + " ...", end=" ", file=sys.stderr)
    for aligned_read in _pysam_iter(samfile):
        total_alignment += 1

        if aligned_read.is_qcfail:
            bad_alignment += 1
            out3.write(aligned_read)
            continue
        if aligned_read.is_unmapped:
            bad_alignment += 1
            out3.write(aligned_read)
            continue

        chrom = samfile.getrname(aligned_read.tid)
        chrom = chrom.upper()
        read_start = aligned_read.pos
        mate_start = aligned_read.mpos

        if aligned_read.mate_is_unmapped:  # only one end mapped
            if chrom not in exon_ranges:
                out2.write(aligned_read)
                ex_alignment += 1
                continue
            else:
                if len(exon_ranges[chrom].find(read_start, read_start + 1)) >= 1:
                    out1.write(aligned_read)
                    in_alignment += 1
                    continue
                elif len(exon_ranges[chrom].find(read_start, read_start + 1)) == 0:
                    out2.write(aligned_read)
                    ex_alignment += 1
                    continue
        else:  # both end mapped
            if chrom not in exon_ranges:
                out2.write(aligned_read)
                ex_alignment += 1
                continue
            else:
                if (len(exon_ranges[chrom].find(read_start, read_start + 1)) >= 1) or (
                    len(exon_ranges[chrom].find(mate_start, mate_start + 1)) >= 1
                ):
                    out1.write(aligned_read)
                    in_alignment += 1
                else:
                    out2.write(aligned_read)
                    ex_alignment += 1

    print("Done", file=sys.stderr)

    print("%-55s%d" % ("Total records:", total_alignment))
    print("%-55s%d" % (args.output_prefix + ".in.bam (Alignments consumed by input gene list):", in_alignment))
    print("%-55s%d" % (args.output_prefix + ".ex.bam (Alignments not consumed by input gene list):", ex_alignment))
    print("%-55s%d" % (args.output_prefix + ".junk.bam (qcfailed, unmapped reads):", bad_alignment))


if __name__ == "__main__":
    main()
