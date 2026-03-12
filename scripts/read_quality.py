#!/usr/bin/env python
"""Calculate Phred quality score distribution for each position on the read."""

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
        help="Alignment file in BAM or SAM format. [required]",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output files(s). [required]",
    )
    parser.add_argument(
        "-r",
        "--reduce",
        type=int,
        dest="reduce_fold",
        default=1,
        help=(
            "To avoid making huge vector in R, nucleotide with"
            " particular phred score less frequent than this number"
            " will be ignored. Increase this number save more memory"
            " while reduce precision. Set to 1 achieves maximum"
            " precision (i.e. every nucleotide will be considered)."
            " This option only applies to the 'boxplot'."
            " default=%(default)s"
        ),
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help=(
            "Minimum mapping quality (phred scaled) for an alignment"
            ' to be called "uniquely mapped". default=%(default)s'
        ),
    )
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    if os.path.exists(args.input_file):
        obj = SAM.ParseBAM(args.input_file)
        obj.readsQual_boxplot(outfile=args.output_prefix, q_cut=args.map_qual, shrink=args.reduce_fold)
        run_rscript(args.output_prefix + ".qual.r")
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
