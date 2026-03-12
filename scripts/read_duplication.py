#!/usr/bin/env python
"""
Calculte reads' duplication rate.
Sequence-based: Reads with identical sequence are considered as "duplicate reads".
Mapping-based: Reads mapped to the exact same location are considered as "duplicate reads".
"""

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
        try:
            subprocess.call("Rscript " + args.output_prefix + ".DupRate_plot.r", shell=True)
        except Exception:
            pass
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
