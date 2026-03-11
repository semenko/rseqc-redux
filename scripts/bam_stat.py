#!/usr/bin/env python
"""
Summarizing mapping statistics of a BAM or SAM file.
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

from optparse import OptionParser

from qcmodule import SAM


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
        "-q",
        "--mapq",
        action="store",
        type="int",
        dest="map_qual",
        default=30,
        help='Minimum mapping quality (phred scaled) to determine "uniquely mapped" reads. default=%default',
    )
    (options, args) = parser.parse_args()

    if not (options.input_file):
        parser.print_help()
        sys.exit(0)
    if not os.path.exists(options.input_file):
        print("\n\n" + input_file + " does NOT exists" + "\n", file=sys.stderr)  # noqa: F821 — known bug: should be options.input_file
        sys.exit(0)

    obj = SAM.ParseBAM(options.input_file)
    obj.stat(q_cut=options.map_qual)


if __name__ == "__main__":
    main()
