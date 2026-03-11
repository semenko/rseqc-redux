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

import subprocess
from optparse import OptionParser
from time import strftime

from qcmodule import SAM


def printlog(mesg):
    """print progress into stderr and log file"""
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    LOG = open("class.log", "a")
    print(mesg, file=sys.stderr)
    print(mesg, file=LOG)


def main():
    usage = "%prog [options]" + "\n" + __doc__ + "\n"
    parser = OptionParser(usage, version="%prog 5.0.2")
    parser.add_option(
        "-i",
        "--input-file",
        action="store",
        type="string",
        dest="input_file",
        help="Alignment file in BAM or SAM format.",
    )
    parser.add_option(
        "-o", "--out-prefix", action="store", type="string", dest="output_prefix", help="Prefix of output files(s)."
    )
    parser.add_option(
        "-q",
        "--mapq",
        action="store",
        type="int",
        dest="map_qual",
        default=30,
        help='Minimum mapping quality (phred scaled) for an alignment to be called "uniquely mapped". default=%default',
    )
    (options, args) = parser.parse_args()

    if not (options.output_prefix and options.input_file):
        parser.print_help()
        sys.exit(0)
    if os.path.exists(options.input_file):
        obj = SAM.ParseBAM(options.input_file)
        obj.readGC(outfile=options.output_prefix, q_cut=options.map_qual)
        try:
            subprocess.call("Rscript " + options.output_prefix + ".GC_plot.r", shell=True)
        except Exception:
            pass
    else:
        print("\n\n" + options.input_file + " does NOT exists" + "\n", file=sys.stderr)
        # parser.print_help()
        sys.exit(0)


if __name__ == "__main__":
    main()
