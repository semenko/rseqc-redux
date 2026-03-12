#!/usr/bin/env python
"""
Equally divide BAM file (m alignments) into n parts. Each part contains roughly m/n alignments
that are randomly sampled from total alignments.
"""

import argparse
import os
import sys
from random import randrange

import pysam


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM format. BAM file should be sorted and indexed.",
    )
    parser.add_argument("-n", "--subset-num", type=int, dest="subset_num", help="Number of small BAM files")
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help='Prefix of output BAM files. Output "Prefix_num.bam".',
    )
    parser.add_argument("-s", "--skip-unmap", action="store_true", dest="skip_unmap", help="Skip unmapped reads.")
    args = parser.parse_args()

    if not (args.input_file and args.subset_num and args.output_prefix):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(args.input_file):
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(0)

    samfile = pysam.Samfile(args.input_file, "rb")

    sub_bam = {}
    count_bam = {}
    for i in range(0, args.subset_num):
        sub_bam[i] = pysam.Samfile(args.output_prefix + "_" + str(i) + ".bam", "wb", template=samfile)
        count_bam[i] = 0

    total_alignment = 0
    print("Dividing " + args.input_file + " ...", end=" ", file=sys.stderr)
    try:
        while 1:
            aligned_read = next(samfile)
            if aligned_read.is_unmapped and args.skip_unmap is True:
                continue
            total_alignment += 1
            tmp = randrange(0, args.subset_num)
            sub_bam[tmp].write(aligned_read)
            count_bam[tmp] += 1

    except (StopIteration, ValueError):
        print("Done", file=sys.stderr)

    for i in range(0, args.subset_num):
        print("%-55s%d" % (args.output_prefix + "_" + str(i) + ".bam", count_bam[i]))


if __name__ == "__main__":
    main()
