#!/usr/bin/env python
"""Check nucleotide frequency at each read position (5'->3') and generate an NVC plot."""

import argparse
import os
import subprocess
import sys
from time import strftime

from rseqc import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)
    with open("class.log", "a") as LOG:
        print(mesg, file=LOG)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Input file in BAM or SAM format.[required]",
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output files(s). [required]",
    )
    parser.add_argument(
        "-x",
        "--nx",
        action="store_true",
        dest="unknown_nucleotide",
        help="Flag option. Presense of this flag tells program to include N,X in output NVC plot [required]",
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
        obj.readsNVC(outfile=args.output_prefix, nx=args.unknown_nucleotide, q_cut=args.map_qual)
        try:
            subprocess.run(["Rscript", args.output_prefix + ".NVC_plot.r"], check=False)
        except OSError:
            pass
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
