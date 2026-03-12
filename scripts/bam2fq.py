#!/usr/bin/env python
"""
Description: Convert alignments in BAM or SAM format into fastq format.
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
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output fastq files(s).",
    )
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

    # print args.single
    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    if os.path.exists(args.input_file):
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
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
