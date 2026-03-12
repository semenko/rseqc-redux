#!/usr/bin/env python
"""
Calculte reads' duplication rate.
Sequence-based: Reads with identical sequence are considered as "duplicate reads".
Mapping-based: Reads mapped to the exact same location are considered as "duplicate reads".
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
    parser.add_argument(
        "-u",
        "--up-limit",
        type=int,
        dest="upper_limit",
        default=500,
        help="Upper limit of reads' occurrence. Only used for plotting, default=%(default)s (times)",
    )
    add_mapq_arg(
        parser,
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
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.readDupRate(outfile=args.output_prefix, up_bound=args.upper_limit, q_cut=args.map_qual)
    run_rscript(args.output_prefix + ".DupRate_plot.r")


if __name__ == "__main__":
    main()
