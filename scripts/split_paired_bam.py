#!/usr/bin/env python
"""Split a paired-end BAM file into two single-end BAM files."""

import sys

import pysam

from rseqc.cli_common import add_input_bam_arg, create_parser, validate_files_exist
from rseqc.SAM import _pysam_iter


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(
        parser,
        help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help=(
            'Prefix of output BAM files. "prefix.R1.bam" file'
            ' contains the 1st read, "prefix.R2.bam" file'
            " contains the 2nd read"
        ),
    )
    args = parser.parse_args()

    if not args.input_file:
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)

    samfile = pysam.AlignmentFile(args.input_file, "rb")
    OUT1 = pysam.AlignmentFile(
        args.output_prefix + ".R1.bam", "wb", template=samfile
    )  # bam file containing reads hit to exon region
    OUT2 = pysam.AlignmentFile(
        args.output_prefix + ".R2.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region
    OUT3 = pysam.AlignmentFile(
        args.output_prefix + ".unmap.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region

    total_alignment = 0
    r1_alignment = 0
    r2_alignment = 0
    unmapped = 0

    print("spliting " + args.input_file + " ...", end=" ", file=sys.stderr)
    for old_alignment in _pysam_iter(samfile):
        new_alignment = pysam.AlignedRead()  # create AlignedRead object
        total_alignment += 1

        new_alignment.qname = old_alignment.qname
        new_alignment.tid = old_alignment.tid
        new_alignment.pos = old_alignment.pos
        new_alignment.mapq = old_alignment.mapq
        new_alignment.cigar = old_alignment.cigar
        new_alignment.seq = old_alignment.seq
        new_alignment.qual = old_alignment.qual
        new_alignment.tags = old_alignment.tags
        new_alignment.flag = 0x0000
        if old_alignment.is_unmapped:
            OUT3.write(old_alignment)
            unmapped += 1
            continue
        if old_alignment.is_reverse:
            new_alignment.flag = new_alignment.flag | 0x0010

        if old_alignment.is_secondary:
            new_alignment.flag = new_alignment.flag | 0x0100
        if old_alignment.is_qcfail:
            new_alignment.flag = new_alignment.flag | 0x0200
        if old_alignment.is_duplicate:
            new_alignment.flag = new_alignment.flag | 0x0400
        if old_alignment.is_read1:
            OUT1.write(new_alignment)
            r1_alignment += 1
        else:
            OUT2.write(new_alignment)
            r2_alignment += 1

    print("Done", file=sys.stderr)

    print("%-55s%d" % ("Total records:", total_alignment))
    print("%-55s%d" % (args.output_prefix + "Read 1:", r1_alignment))
    print("%-55s%d" % (args.output_prefix + "Read 2:", r2_alignment))
    print("%-55s%d" % (args.output_prefix + "Unmapped:", unmapped))


if __name__ == "__main__":
    main()
