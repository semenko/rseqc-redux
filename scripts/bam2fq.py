#!/usr/bin/env python
"""
Description: Convert alignments in BAM or SAM format into fastq format.
"""

import subprocess
import sys

from rseqc import SAM
from rseqc.cli_common import add_input_bam_arg, add_output_prefix_arg, create_parser, validate_files_exist


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser)
    add_output_prefix_arg(parser, help="Prefix of output fastq files(s).")
    parser.add_argument(
        "-s",
        "--single-end",
        action="store_true",
        dest="single",
        help="Specificy '-s' or '--single-end' for single-end sequencing.",
    )
    parser.add_argument(
        "-c",
        "--compress",
        action="store_true",
        dest="gzip",
        help="Specificy '-c' or '--compress' to compress output fastq file(s) using 'gzip' command.",
    )
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    if args.single is True:
        obj.bam2fq(prefix=args.output_prefix, paired=False)
        if args.gzip is True:
            try:
                print("run gzip ... ", end=" ", file=sys.stderr)
                subprocess.run(["gzip", args.output_prefix + ".fastq"], check=False)
                print("Done.", file=sys.stderr)
            except OSError:
                pass
    else:
        obj.bam2fq(prefix=args.output_prefix, paired=True)
        if args.gzip is True:
            try:
                print("run gzip ...", file=sys.stderr)
                subprocess.run(["gzip", args.output_prefix + ".R1.fastq"], check=False)
                subprocess.run(["gzip", args.output_prefix + ".R2.fastq"], check=False)
                print("Done.", file=sys.stderr)
            except OSError:
                pass


if __name__ == "__main__":
    main()
