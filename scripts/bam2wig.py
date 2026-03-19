#!/usr/bin/env python
"""
Convert BAM file into wig file. BAM file must be sorted and indexed using SAMtools.
Note: SAM format file is not supported.
"""

import sys

from rseqc import SAM
from rseqc.cli_common import add_mapq_arg, create_parser, load_chromsize, validate_bam_index, validate_files_exist


def main() -> None:
    parser = create_parser(__doc__)
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help=(
            "Alignment file in BAM format. BAM file must be sorted and"
            " indexed using samTools. .bam and .bai files should be"
            " placed in the same directory."
        ),
    )
    parser.add_argument(
        "-s",
        "--chromSize",
        dest="chromSize",
        help=(
            "Chromosome size file. Tab or space separated text file"
            " with 2 columns: first column is chromosome name/ID,"
            " second column is chromosome size. Chromosome name"
            ' (such as "chr1") should be consistent between this'
            " file and the BAM file."
        ),
    )
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help=(
            "Prefix of output wiggle files(s). One wiggle file will"
            " be generated for non strand-specific data, two wiggle"
            ' files ("Prefix_Forward.wig" and "Prefix_Reverse.wig")'
            " will be generated for strand-specific RNA-seq data."
        ),
    )
    parser.add_argument(
        "-t",
        "--wigsum",
        type=int,
        dest="total_wigsum",
        help=(
            "Specified wigsum. Eg: 1,000,000,000 equals to coverage"
            " of 10 million 100nt reads. Ignore this option to"
            " disable normalization"
        ),
    )
    parser.add_argument(
        "-u", "--skip-multi-hits", action="store_true", dest="skip_multi", help="Skip non-unique hit reads."
    )
    parser.add_argument(
        "-d",
        "--strand",
        dest="strand_rule",
        default=None,
        help=(
            "How read(s) were stranded during sequencing. For example:"
            " --strand='1++,1--,2+-,2-+' means that this is a pair-end,"
            " strand-specific RNA-seq data, and the strand rule is:"
            " read1 mapped to '+' => parental gene on '+';"
            " read1 mapped to '-' => parental gene on '-';"
            " read2 mapped to '+' => parental gene on '-';"
            " read2 mapped to '-' => parental gene on '+'."
            " If you are not sure about the strand rule, run"
            " 'infer_experiment.py' default=%(default)s"
            " (Not a strand specific RNA-seq data)."
        ),
    )
    add_mapq_arg(parser, help='Minimum mapping quality to determine "uniquely mapped". default=%(default)s')

    args = parser.parse_args()

    if not (args.output_prefix and args.input_file and args.chromSize):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.input_file, args.chromSize)
    validate_bam_index(args.input_file)

    if args.skip_multi:
        print("Skip multi-hits:True")
    else:
        print("Skip multi-hits:False")

    chromSizes = load_chromsize(args.chromSize)

    norm_factor = None
    if args.total_wigsum:
        obj = SAM.ParseBAM(args.input_file)
        wig_sum = obj.calWigSum(chrom_sizes=chromSizes, skip_multi=args.skip_multi, q_cut=args.map_qual)
        print("\n\ntotal wigsum is:" + str(wig_sum) + "\n", file=sys.stderr)
        try:
            norm_factor = args.total_wigsum / wig_sum
        except ZeroDivisionError:
            norm_factor = None

    obj = SAM.ParseBAM(args.input_file)
    obj.bamTowig(
        outfile=args.output_prefix,
        chrom_sizes=chromSizes,
        chrom_file=args.chromSize,
        q_cut=args.map_qual,
        skip_multi=args.skip_multi,
        strand_rule=args.strand_rule,
        WigSumFactor=norm_factor,
    )


if __name__ == "__main__":
    main()
