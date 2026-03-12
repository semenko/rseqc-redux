#!/usr/bin/env python
"""Check nucleotide frequency at each read position (5'->3') and generate an NVC plot."""

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
    add_input_bam_arg(parser, help="Input file in BAM or SAM format.[required]")
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-x",
        "--nx",
        action="store_true",
        dest="unknown_nucleotide",
        help="Flag option. Presense of this flag tells program to include N,X in output NVC plot [required]",
    )
    add_mapq_arg(parser)
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.readsNVC(outfile=args.output_prefix, nx=args.unknown_nucleotide, q_cut=args.map_qual)
    run_rscript(args.output_prefix + ".NVC_plot.r")


if __name__ == "__main__":
    main()
