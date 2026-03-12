#!/usr/bin/env python
"""
Summarizing mapping statistics of a BAM or SAM file.
"""

import argparse
import os
import sys

from rseqc import SAM


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help='Minimum mapping quality (phred scaled) to determine "uniquely mapped" reads. default=%(default)s',
    )
    args = parser.parse_args()

    if not (args.input_file):
        parser.print_help()
        sys.exit(1)
    if not os.path.exists(args.input_file):
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(1)

    obj = SAM.ParseBAM(args.input_file)
    obj.stat(q_cut=args.map_qual)


if __name__ == "__main__":
    main()
