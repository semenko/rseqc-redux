#!/usr/bin/env python
"""
Equally divide BAM file (m alignments) into n parts. Each part contains roughly m/n alignments
that are randomly sampled from total alignments.
"""

import contextlib
import sys
from random import randrange

import pysam

from rseqc.cli_common import create_parser, validate_files_exist
from rseqc.SAM import _pysam_iter


def main() -> None:
    parser = create_parser(__doc__)
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
        sys.exit(1)
    validate_files_exist(args.input_file)

    with contextlib.ExitStack() as stack:
        samfile = stack.enter_context(pysam.AlignmentFile(args.input_file, "rb"))

        sub_bam = {}
        count_bam = {}
        for i in range(0, args.subset_num):
            sub_bam[i] = stack.enter_context(
                pysam.AlignmentFile(args.output_prefix + "_" + str(i) + ".bam", "wb", template=samfile)
            )
            count_bam[i] = 0

        total_alignment = 0
        print("Dividing " + args.input_file + " ...", end=" ", file=sys.stderr)
        for aligned_read in _pysam_iter(samfile):
            if aligned_read.is_unmapped and args.skip_unmap is True:
                continue
            total_alignment += 1
            tmp = randrange(0, args.subset_num)
            sub_bam[tmp].write(aligned_read)
            count_bam[tmp] += 1

        print("Done", file=sys.stderr)

    for i in range(0, args.subset_num):
        print(f"{args.output_prefix + '_' + str(i) + '.bam':<55}{count_bam[i]}")


if __name__ == "__main__":
    main()
