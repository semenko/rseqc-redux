#!/usr/bin/env python
"""Calculate the GC content distribution of mapped reads."""

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
    add_mapq_arg(parser)
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.readGC(outfile=args.output_prefix, q_cut=args.map_qual)
    run_rscript(args.output_prefix + ".GC_plot.r")


if __name__ == "__main__":
    main()
