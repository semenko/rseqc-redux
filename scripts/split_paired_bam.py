#!/usr/bin/env python
"""-------------------------------------------------------------------------------------------------
Split bam file (pair-end) into 2 single-end bam file
-------------------------------------------------------------------------------------------------"""

import argparse
import os
import sys

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format. BAM file should be sorted and indexed",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help='Prefix of output BAM files. "prefix.R1.bam" file contains the 1st read, "prefix.R2.bam" file contains the 2nd read',
    )
    args = parser.parse_args()

    if not (args.input_file):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(args.input_file):
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(0)

    samfile = pysam.Samfile(args.input_file, "rb")
    OUT1 = pysam.Samfile(
        args.output_prefix + ".R1.bam", "wb", template=samfile
    )  # bam file containing reads hit to exon region
    OUT2 = pysam.Samfile(
        args.output_prefix + ".R2.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region
    OUT3 = pysam.Samfile(
        args.output_prefix + ".unmap.bam", "wb", template=samfile
    )  # bam file containing reads not hit to exon region

    total_alignment = 0
    r1_alignment = 0
    r2_alignment = 0
    unmapped = 0

    print("spliting " + args.input_file + " ...", end=" ", file=sys.stderr)
    try:
        while 1:
            new_alignment = pysam.AlignedRead()  # create AlignedRead object
            old_alignment = next(samfile)
            total_alignment += 1

            new_alignment.qname = old_alignment.qname  # 1st column. read name.
            # new_alignment.flag = old_alignment.flag        # 2nd column. subject to change. flag value
            new_alignment.tid = old_alignment.tid  # 3rd column. samfile.getrname(tid) == chrom name
            new_alignment.pos = (
                old_alignment.pos
            )  # 4th column. reference Start position of the aligned part (of read) [0-based]
            new_alignment.mapq = old_alignment.mapq  # 5th column. mapping quality
            new_alignment.cigar = old_alignment.cigar  # 6th column. subject to change.
            # new_alignment.rnext = old_alignment.rnext      # 7th column. tid of the reference (mate read mapped to)
            # new_alignment.pnext = old_alignment.pnext      # 8th column. position of the reference (0 based, mate read mapped to)
            # new_alignment.tlen = old_alignment.tlen        # 9th column. insert size
            new_alignment.seq = old_alignment.seq  # 10th column. read sequence. all bases.
            new_alignment.qual = old_alignment.qual  # 11th column. read sequence quality. all bases.
            new_alignment.tags = old_alignment.tags  # 12 - columns
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

    except StopIteration:
        print("Done", file=sys.stderr)

    print("%-55s%d" % ("Total records:", total_alignment))
    print("%-55s%d" % (args.output_prefix + "Read 1:", r1_alignment))
    print("%-55s%d" % (args.output_prefix + "Read 2:", r2_alignment))
    print("%-55s%d" % (args.output_prefix + "Unmapped:", unmapped))


if __name__ == "__main__":
    main()
