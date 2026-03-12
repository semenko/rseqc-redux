#!/usr/bin/env python
"""
Calculate the distribution of mismatches across reads. Note that the "MD" tag must exist in BAM file.
"""

import sys

from rseqc import SAM
from rseqc.cli_common import add_mapq_arg, add_output_prefix_arg, create_parser, run_rscript, validate_files_exist


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument("-i", "--input", dest="input_bam", help="Input BAM file. [required]")
    parser.add_argument(
        "-l",
        "--read-align-length",
        type=int,
        dest="read_alignment_length",
        help=(
            "Alignment length of read. It is usually set to the"
            " original read length. For example, all these cigar"
            ' strings ("101M", "68M140N33M", "53M1D48M") suggest'
            " the read alignment length is 101. [required]"
        ),
    )
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-n",
        "--read-num",
        type=int,
        default=1000000,
        dest="read_number",
        help="Number of aligned reads with mismatches used to calculate the mismatch profile. default=%(default)s",
    )
    add_mapq_arg(parser, help="Minimum mapping quality. default=%(default)s")
    args = parser.parse_args()

    if not args.input_bam:
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_bam)

    if not args.output_prefix:
        print("\n\n You must specify the output prefix", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    if not args.read_alignment_length:
        print("\n\n You must specify read alignment length. It is usually the read length.", file=sys.stderr)
        parser.print_help()
        sys.exit(1)

    obj = SAM.ParseBAM(args.input_bam)
    obj.mismatchProfile(
        read_length=args.read_alignment_length,
        read_num=args.read_number,
        q_cut=args.map_qual,
        outfile=args.output_prefix,
    )
    run_rscript(args.output_prefix + ".mismatch_profile.r")


if __name__ == "__main__":
    main()
