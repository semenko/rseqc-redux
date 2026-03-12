#!/usr/bin/env python
"""
Calculate the distributions of inserted nucleotides across reads
Note CIGAR strings within SAM/BAM file should have 'I' operation
"""

import argparse
import os
import subprocess
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
    parser.add_argument("-o", "--out-prefix", dest="output_prefix", help="Prefix of output files(s).")
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
    parser.add_argument(
        "-s",
        "--sequencing",
        dest="layout",
        help='Sequencing layout. "SE"(single-end) or "PE"(pair-end). ',
    )
    args = parser.parse_args()

    if not (args.input_file and args.output_prefix and args.layout):
        parser.print_help()
        sys.exit(1)
    for input_file in [args.input_file]:
        if not os.path.exists(input_file):
            print("\n\n" + input_file + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(1)

    obj = SAM.ParseBAM(args.input_file)
    if args.layout == "SE":
        obj.insertion_profile(outfile=args.output_prefix, q_cut=args.map_qual, PE=False)
    elif args.layout == "PE":
        obj.insertion_profile(outfile=args.output_prefix, q_cut=args.map_qual, PE=True)
    else:
        print('unknow sequencing layout. Must be "SE" or "PE"', file=sys.stderr)
    try:
        subprocess.run(["Rscript", args.output_prefix + ".insertion_profile.r"], check=False)
    except OSError:
        print("Cannot generate pdf file from " + args.output_prefix + ".insertion_profile.r", file=sys.stderr)
        pass


if __name__ == "__main__":
    main()
