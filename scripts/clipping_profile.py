#!/usr/bin/env python
"""
This program is to estimate clipping profile of RNA-seq reads from BAM or SAM file
Note that to use this funciton, CIGAR strings within SAM/BAM file should have 'S' operation
(This means your reads mapper should support clipped mapping).
"""

import os
import sys

if sys.version_info[0] != 3:
    print(
        "\nYou are using python"
        + str(sys.version_info[0])
        + "."
        + str(sys.version_info[1])
        + " This verion of RSeQC needs python3!\n",
        file=sys.stderr,
    )
    sys.exit()

import argparse
import subprocess

from qcmodule import SAM


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
        help='Minimum mapping quality (phred scaled) for an alignment to be considered as "uniquely mapped". default=%(default)s',
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
        sys.exit(0)
    for input_file in [args.input_file]:
        if not os.path.exists(input_file):
            print("\n\n" + input_file + " does NOT exists" + "\n", file=sys.stderr)
            sys.exit(0)

    obj = SAM.ParseBAM(args.input_file)
    if args.layout == "SE":
        obj.clipping_profile(outfile=args.output_prefix, q_cut=args.map_qual, type="S", PE=False)
    elif args.layout == "PE":
        obj.clipping_profile(outfile=args.output_prefix, q_cut=args.map_qual, type="S", PE=True)
    else:
        print('unknow sequencing layout. Must be "SE" or "PE"', file=sys.stderr)
    try:
        subprocess.call("Rscript " + args.output_prefix + ".clipping_profile.r", shell=True)
    except Exception:
        print("Cannot generate pdf file from " + args.output_prefix + ".clipping_profile.r", file=sys.stderr)
        pass


if __name__ == "__main__":
    main()
