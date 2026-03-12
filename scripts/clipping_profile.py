#!/usr/bin/env python
"""
This program is to estimate clipping profile of RNA-seq reads from BAM or SAM file
Note that to use this funciton, CIGAR strings within SAM/BAM file should have 'S' operation
(This means your reads mapper should support clipped mapping).
"""

import sys

from rseqc import SAM
from rseqc.cli_common import (
    add_input_bam_arg,
    add_mapq_arg,
    add_output_prefix_arg,
    create_parser,
    run_rscript,
    validate_files_exist,
)


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser)
    add_output_prefix_arg(parser)
    add_mapq_arg(
        parser,
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
    validate_files_exist(args.input_file)

    obj = SAM.ParseBAM(args.input_file)
    if args.layout == "SE":
        obj.clipping_profile(outfile=args.output_prefix, q_cut=args.map_qual, type="S", PE=False)
    elif args.layout == "PE":
        obj.clipping_profile(outfile=args.output_prefix, q_cut=args.map_qual, type="S", PE=True)
    else:
        print('unknow sequencing layout. Must be "SE" or "PE"', file=sys.stderr)
    run_rscript(args.output_prefix + ".clipping_profile.r")


if __name__ == "__main__":
    main()
