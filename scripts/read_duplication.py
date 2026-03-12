#!/usr/bin/env python
"""
Calculte reads' duplication rate.
Sequence-based: Reads with identical sequence are considered as "duplicate reads".
Mapping-based: Reads mapped to the exact same location are considered as "duplicate reads".
"""

import argparse
import os
import sys

from rseqc import SAM
from rseqc.cli_common import run_rscript


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_argument("-o", "--out-prefix", dest="output_prefix", help="Prefix of output files(s).")
    parser.add_argument(
        "-u",
        "--up-limit",
        type=int,
        dest="upper_limit",
        default=500,
        help="Upper limit of reads' occurrence. Only used for plotting, default=%(default)s (times)",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help=(
            "Minimum mapping quality (phred scaled) for an alignment"
            ' to be considered as "uniquely mapped".'
            " default=%(default)s"
        ),
    )
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    if os.path.exists(args.input_file):
        obj = SAM.ParseBAM(args.input_file)
        obj.readDupRate(outfile=args.output_prefix, up_bound=args.upper_limit, q_cut=args.map_qual)
        run_rscript(args.output_prefix + ".DupRate_plot.r")
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
