#!/usr/bin/env python
"""
Annotate splicing reads against gene model in two levels: reads level and  juncion level.
Note:
1) A read, especially long read, can be spliced 2 or more times
2) Multiple splicing reads spanning the same intron can be consolidated into one splicing junction.
"""

import argparse
import os
import sys

from rseqc import SAM
from rseqc.cli_common import run_rscript


def generate_bed12(infile, size=1):
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


def generate_interact(infile, bam_file, size=1):
    """
    infile: input file. eg: chrX    66766604        66788677        348     partial_novel
    size: the block size representing exons
    """

    outfile = infile.replace(".xls", ".Interact.bed")
    with open(outfile, "w") as OUT:
        print(
            'track type=interact name="Splice junctions"'
            ' description="Splice junctions detected from %s"'
            " maxHeightPixels=200:200:50 visibility=full" % bam_file,
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--input-file",
        dest="input_file",
        help="Alignment file in BAM or SAM format.",
    )
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
    parser.add_argument(
        "-o",
        "--out-prefix",
        dest="output_prefix",
        help="Prefix of output files(s). [required]",
    )
    parser.add_argument(
        "-m",
        "--min-intron",
        type=int,
        dest="min_intron",
        default=50,
        help="Minimum intron length (bp). default=%(default)s [optional]",
    )
    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        dest="map_qual",
        default=30,
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
    if not os.path.exists(args.ref_gene_model):
        print("\n\n" + args.ref_gene_model + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(1)
    if os.path.exists(args.input_file):
        obj = SAM.ParseBAM(args.input_file)
        obj.annotate_junction(
            outfile=args.output_prefix,
            refgene=args.ref_gene_model,
            min_intron=args.min_intron,
            q_cut=args.map_qual,
        )
        run_rscript(args.output_prefix + ".junction_plot.r")
    else:
        print("\n\n" + args.input_file + " does NOT exists" + "\n", file=sys.stderr)
        sys.exit(1)
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
