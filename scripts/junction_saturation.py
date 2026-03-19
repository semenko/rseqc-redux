#!/usr/bin/env python
"""
Check if sequencing depth is saturated or not, based on the concept that when sequencing depth is
approaching saturation, less NEW junctions will be discovered.
See http://rseqc.sourceforge.net/ for details.
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
    add_input_bam_arg(parser, help="Alignment file in BAM or SAM format.[required]")
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-r",
        "--refgene",
        dest="refgene_bed",
        help=(
            "Reference gene model in bed format. This gene model is"
            " used to determine known splicing junctions. [required]"
        ),
    )
    parser.add_argument(
        "-l",
        "--percentile-floor",
        type=int,
        dest="percentile_low_bound",
        default=5,
        help="Sampling starts from this percentile. A integer between 0 and 100. default=%(default)s",
    )
    parser.add_argument(
        "-u",
        "--percentile-ceiling",
        type=int,
        dest="percentile_up_bound",
        default=100,
        help="Sampling ends at this percentile. A integer between 0 and 100. default=%(default)s",
    )
    parser.add_argument(
        "-s",
        "--percentile-step",
        type=int,
        dest="percentile_step",
        default=5,
        help=(
            "Sampling frequency. Smaller value means more sampling"
            " times. A integer between 0 and 100."
            " default=%(default)s"
        ),
    )
    parser.add_argument(
        "-m",
        "--min-intron",
        type=int,
        dest="minimum_intron_size",
        default=50,
        help="Minimum intron size (bp). default=%(default)s",
    )
    parser.add_argument(
        "-v",
        "--min-coverage",
        type=int,
        dest="minimum_splice_read",
        default=1,
        help="Minimum number of supportting reads to call a junction. default=%(default)s",
    )
    add_mapq_arg(parser)

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file and args.refgene_bed):
        parser.print_help()
        sys.exit(1)
    if args.percentile_low_bound < 0 or args.percentile_low_bound > 100:
        print("percentile_low_bound must be larger than 0 and smaller than 100", file=sys.stderr)
        sys.exit(1)
    if args.percentile_up_bound < 0 or args.percentile_up_bound > 100:
        print("percentile_up_bound must be larger than 0 and smaller than 100", file=sys.stderr)
        sys.exit(1)
    if args.percentile_up_bound < args.percentile_low_bound:
        print("percentile_up_bound must be larger than percentile_low_bound", file=sys.stderr)
        sys.exit(1)
    if args.percentile_step < 0 or args.percentile_step > args.percentile_up_bound:
        print("percentile_step must be larger than 0 and smaller than percentile_up_bound", file=sys.stderr)
        sys.exit(1)
    validate_files_exist(args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.saturation_junction(
        outfile=args.output_prefix,
        refgene=args.refgene_bed,
        sample_start=args.percentile_low_bound,
        sample_end=args.percentile_up_bound,
        sample_step=args.percentile_step,
        min_intron=args.minimum_intron_size,
        recur=args.minimum_splice_read,
        q_cut=args.map_qual,
    )
    run_rscript(args.output_prefix + ".junctionSaturation_plot.r")


if __name__ == "__main__":
    main()
