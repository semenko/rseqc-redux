#!/usr/bin/env python
"""
Annotate splicing reads against gene model in two levels: reads level and  juncion level.
Note:
1) A read, especially long read, can be spliced 2 or more times
2) Multiple splicing reads spanning the same intron can be consolidated into one splicing junction.
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


def generate_bed12(infile: str, size: int = 1) -> None:
    """
    infile: input file. eg: chrX    66766604        66788677        348     partial_novel
    size: the block size representing exons
    """

    outfile = infile.replace(".xls", ".bed")
    with open(outfile, "w") as OUT:
        with open(infile, "r") as _fh:
            for line in _fh:
                if line.startswith("chrom"):
                    continue
                line = line.strip()
                f = line.split()
                if len(f) != 5:
                    continue
                chrom = f[0]
                start = int(f[1]) - size
                end = int(f[2]) + size
                score = int(f[3])
                strand = "."
                name = f[4]
                thick_st = start
                thick_end = end
                if name == "annotated":
                    color = "205,0,0"
                elif name == "partial_novel":
                    color = "0,205,0"
                elif name == "complete_novel":
                    color = "0,0,205"
                else:
                    color = "0,0,0"
                blockCount = 2
                blockSizes = ",".join((str(size), str(size)))
                blockStarts = "0," + str(end - size - start)
                print(
                    "\t".join(
                        [
                            str(i)
                            for i in [
                                chrom,
                                start,
                                end,
                                name,
                                score,
                                strand,
                                thick_st,
                                thick_end,
                                color,
                                blockCount,
                                blockSizes,
                                blockStarts,
                            ]
                        ]
                    ),
                    file=OUT,
                )


def generate_interact(infile: str, bam_file: str, size: int = 1) -> None:
    """
    infile: input file. eg: chrX    66766604        66788677        348     partial_novel
    size: the block size representing exons
    """

    outfile = infile.replace(".xls", ".Interact.bed")
    with open(outfile, "w") as OUT:
        print(
            'track type=interact name="Splice junctions"'
            f' description="Splice junctions detected from {bam_file}"'
            " maxHeightPixels=200:200:50 visibility=full",
            file=OUT,
        )
        with open(infile, "r") as _fh:
            for line in _fh:
                if line.startswith("chrom"):
                    continue
                line = line.strip()
                f = line.split()
                if len(f) != 5:
                    continue
                chrom = f[0]
                chromStart = int(f[1]) - size
                chromEnd = int(f[2]) + size
                name1 = f[4]
                if name1 == "annotated":
                    color = "205,0,0"
                elif name1 == "partial_novel":
                    color = "0,205,0"
                elif name1 == "complete_novel":
                    color = "0,0,205"
                else:
                    color = "0,0,0"

                name = chrom + ":" + str(chromStart) + "-" + str(chromEnd) + "_" + name1
                score = int(f[3])
                value = float(score)
                exp = "RNAseq_junction"

                sourceChrom = chrom
                sourceStart = chromStart
                sourceEnd = chromStart + size
                sourceName = sourceChrom + ":" + str(sourceStart) + "-" + str(sourceEnd)
                sourceStrand = "."

                targetChrom = chrom
                targetStart = chromEnd - size
                targetEnd = chromEnd
                targetName = targetChrom + ":" + str(targetStart) + "-" + str(targetEnd)
                targetStrand = "."
                print(
                    "\t".join(
                        [
                            str(i)
                            for i in [
                                chrom,
                                chromStart,
                                chromEnd,
                                name,
                                score,
                                value,
                                exp,
                                color,
                                sourceChrom,
                                sourceStart,
                                sourceEnd,
                                sourceName,
                                sourceStrand,
                                targetChrom,
                                targetStart,
                                targetEnd,
                                targetName,
                                targetStrand,
                            ]
                        ]
                    ),
                    file=OUT,
                )


def main() -> None:
    parser = create_parser(__doc__)
    add_input_bam_arg(parser)
    parser.add_argument(
        "-r",
        "--refgene",
        dest="ref_gene_model",
        help=(
            "Reference gene model in bed format. This file is better"
            " to be a pooled gene model as it will be used to"
            " annotate splicing junctions [required]"
        ),
    )
    add_output_prefix_arg(parser, help="Prefix of output files(s). [required]")
    parser.add_argument(
        "-m",
        "--min-intron",
        type=int,
        dest="min_intron",
        default=50,
        help="Minimum intron length (bp). default=%(default)s [optional]",
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

    if not (args.output_prefix and args.input_file and args.ref_gene_model):
        parser.print_help()
        sys.exit(1)
    validate_files_exist(args.ref_gene_model, args.input_file)
    obj = SAM.ParseBAM(args.input_file)
    obj.annotate_junction(
        outfile=args.output_prefix,
        refgene=args.ref_gene_model,
        min_intron=args.min_intron,
        q_cut=args.map_qual,
    )
    run_rscript(args.output_prefix + ".junction_plot.r")
    try:
        print("Create BED file ...", file=sys.stderr)
        generate_bed12(args.output_prefix + ".junction.xls")
    except (OSError, IndexError, ValueError):
        pass
    try:
        print("Create Interact file ...", file=sys.stderr)
        generate_interact(args.output_prefix + ".junction.xls", args.input_file)
    except (OSError, IndexError, ValueError):
        pass


if __name__ == "__main__":
    main()
