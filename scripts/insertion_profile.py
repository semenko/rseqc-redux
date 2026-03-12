#!/usr/bin/env python
"""
Calculate the distributions of inserted nucleotides across reads
Note CIGAR strings within SAM/BAM file should have 'I' operation
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
        obj.insertion_profile(outfile=args.output_prefix, q_cut=args.map_qual, PE=False)
    elif args.layout == "PE":
        obj.insertion_profile(outfile=args.output_prefix, q_cut=args.map_qual, PE=True)
    else:
        print('unknow sequencing layout. Must be "SE" or "PE"', file=sys.stderr)
    run_rscript(args.output_prefix + ".insertion_profile.r")


if __name__ == "__main__":
    main()
