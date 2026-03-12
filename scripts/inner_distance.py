#!/usr/bin/env python
"""
Calculate the inner distance (insert size)  of RNA-seq fragments.

               RNA fragment
 _________________||_________________
|                                    |
|                                    |
||||||||||------------------||||||||||
  read_1      insert_size     read_2

fragment size = read_1 + insert_size + read_2
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
    add_output_prefix_arg(parser, help="Prefix of output files(s)")
    parser.add_argument("-r", "--refgene", dest="ref_gene", help="Reference gene model in BED format.")
    parser.add_argument(
        "-k",
        "--sample-size",
        type=int,
        dest="sampleSize",
        default=1000000,
        help="Number of read-pairs used to estimate inner distance. default=%(default)s",
    )
    parser.add_argument(
        "-l",
        "--lower-bound",
        type=int,
        dest="lower_bound_size",
        default=-250,
        help="Lower bound of inner distance (bp). This option is used for ploting histograme. default=%(default)s",
    )
    parser.add_argument(
        "-u",
        "--upper-bound",
        type=int,
        dest="upper_bound_size",
        default=250,
        help="Upper bound of inner distance (bp). This option is used for plotting histogram. default=%(default)s",
    )
    parser.add_argument(
        "-s",
        "--step",
        type=int,
        dest="step_size",
        default=5,
        help="Step size (bp) of histograme. This option is used for plotting histogram. default=%(default)s",
    )
    add_mapq_arg(parser)

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file and args.ref_gene):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file, args.ref_gene)
    if args.step_size <= 0:
        print("step size is a positive interger", file=sys.stderr)
        sys.exit(1)
    obj = SAM.ParseBAM(args.input_file)
    obj.mRNA_inner_distance(
        outfile=args.output_prefix,
        low_bound=args.lower_bound_size,
        up_bound=args.upper_bound_size,
        step=args.step_size,
        refbed=args.ref_gene,
        sample_size=args.sampleSize,
        q_cut=args.map_qual,
    )
    run_rscript(args.output_prefix + ".inner_distance_plot.r")


if __name__ == "__main__":
    main()
