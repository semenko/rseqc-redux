#!/usr/bin/env python
"""Infer RNA-seq experiment design (strandedness and layout) from a SAM/BAM file."""

import sys

from rseqc import SAM
from rseqc.cli_common import add_input_bam_arg, add_mapq_arg, add_refgene_arg, create_parser, validate_files_exist


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser, help="Input alignment file in SAM or BAM format")
    add_refgene_arg(parser, dest="refgene_bed")
    parser.add_argument(
        "-s",
        "--sample-size",
        type=int,
        dest="sample_size",
        default=200000,
        help="Number of reads sampled from SAM/BAM file. default=%(default)s",
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

    if not (args.input_file and args.refgene_bed):
        parser.print_help()
        print("\n\n" + __doc__, file=sys.stderr)
        sys.exit(1)
    validate_files_exist(args.input_file, args.refgene_bed)
    if args.sample_size < 1000:
        print("Warn: Sample Size too small to give a accurate estimation", file=sys.stderr)
    obj = SAM.ParseBAM(args.input_file)
    (protocol, sp1, sp2, other) = obj.configure_experiment(
        refbed=args.refgene_bed, sample_size=args.sample_size, q_cut=args.map_qual
    )
    if other < 0:
        other = 0.0
    if protocol == "PairEnd":
        print("\n\nThis is PairEnd Data")
        print("Fraction of reads failed to determine: %.4f" % other)
        print('Fraction of reads explained by "1++,1--,2+-,2-+": %.4f' % sp1)
        print('Fraction of reads explained by "1+-,1-+,2++,2--": %.4f' % sp2)

    elif protocol == "SingleEnd":
        print("\n\nThis is SingleEnd Data")
        print("Fraction of reads failed to determine: %.4f" % other)
        print('Fraction of reads explained by "++,--": %.4f' % sp1)
        print('Fraction of reads explained by "+-,-+": %.4f' % sp2)

    else:
        print("Unknown Data type")


if __name__ == "__main__":
    main()
