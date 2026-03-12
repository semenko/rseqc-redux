#!/usr/bin/env python
"""-------------------------------------------------------------------------------------------------
Calculate distribution of reads' GC content
-------------------------------------------------------------------------------------------------"""

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
from time import strftime

from qcmodule import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG = open("class.log", "a")
    print(mesg, file=sys.stderr)
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
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
        help='Minimum mapping quality (phred scaled) for an alignment to be called "uniquely mapped". default=%(default)s',
    )
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(0)
    if os.path.exists(args.input_file):
        obj = SAM.ParseBAM(args.input_file)
        obj.readGC(outfile=args.output_prefix, q_cut=args.map_qual)
        try:
            subprocess.call("Rscript " + args.output_prefix + ".GC_plot.r", shell=True)
        except Exception:
            pass
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    main()
