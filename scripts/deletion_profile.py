#!/usr/bin/env python
"""
Calculate the distributions of deleted nucleotides across reads.
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument("-i", "--input", dest="input_bam", help="Input BAM file. [required]")
    parser.add_argument(
        "-l",
        "--read-align-length",
        type=int,
        dest="read_alignment_length",
        help='Alignment length of read. It is usually set to the orignial read length. For example, all these cigar strings ("101M", "68M140N33M", "53M1D48M") suggest the read alignment length is 101. [required]',
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output files(s). [required]",
    )
    parser.add_argument(
        "-n",
        "--read-num",
        type=int,
        default=1000000,
        dest="read_number",
        help="Number of aligned reads with deletions used to calculate the deletion profile. default=%(default)s",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help="Minimum mapping quality. default=%(default)s",
    )
    args = parser.parse_args()

    if not (args.input_bam):
        parser.print_help()
        sys.exit(0)
    for f in [args.input_bam]:
        if not os.path.exists(f):
            print("\n\n" + f + " does NOT exists" + "\n", file=sys.stderr)
            parser.print_help()
            sys.exit(0)

    if not (args.output_prefix):
        print("\n\n You must specify the output prefix", file=sys.stderr)
        parser.print_help()
        sys.exit(0)

    if not (args.read_alignment_length):
        print("\n\n You must specify read alignment length. It is usually the read length.", file=sys.stderr)
        parser.print_help()
        sys.exit(0)

    obj = SAM.ParseBAM(args.input_bam)
    obj.deletionProfile(
        read_length=args.read_alignment_length,
        read_num=args.read_number,
        q_cut=args.map_qual,
        outfile=args.output_prefix,
    )

    try:
        subprocess.call("Rscript " + args.output_prefix + ".deletion_profile.r", shell=True)
    except Exception:
        print("Cannot generate pdf file from " + args.output_prefix + ".deletion_profile.r", file=sys.stderr)
        pass


if __name__ == "__main__":
    main()
