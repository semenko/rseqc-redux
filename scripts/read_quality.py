#!/usr/bin/env python
"""Calculate Phred quality score distribution for each position on the read."""

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
    add_input_bam_arg(parser, help="Alignment file in BAM or SAM format. [required]")
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-r",
        "--reduce",
        type=int,
        dest="reduce_fold",
        default=1,
        help=(
            "To avoid making huge vector in R, nucleotide with"
            " particular phred score less frequent than this number"
            " will be ignored. Increase this number save more memory"
            " while reduce precision. Set to 1 achieves maximum"
            " precision (i.e. every nucleotide will be considered)."
            " This option only applies to the 'boxplot'."
            " default=%(default)s"
        ),
    )
    add_mapq_arg(parser)
    args = parser.parse_args()

    if not (args.output_prefix and args.input_file):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.readsQual_boxplot(outfile=args.output_prefix, q_cut=args.map_qual, shrink=args.reduce_fold)
    run_rscript(args.output_prefix + ".qual.r")


if __name__ == "__main__":
    main()
