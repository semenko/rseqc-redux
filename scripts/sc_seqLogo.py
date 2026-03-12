#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program generates a DNA sequence logo from fasta or fastq format file.
It is useful to visualize the nucleotide compositions of sample barcodes,
cell barcodes and molecular barcodes.
"""

import argparse
import logging
import sys

from rseqc import fastq

__contributor__ = "Liguo Wang"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--infile",
        dest="in_file",
        help='Input DNA sequence file in FASTQ (https://en.wikipedia.org/wiki/FASTQ_format#),\
FASTA (https://en.wikipedia.org/wiki/FASTA_format) or pure sequence format. All sequences must be\
 the same length. This file can be plain text or compressed format (".gz", ".Z",\
".z",".bz", ".bz2", ".bzip2").',
    )
    parser.add_argument("-o", "--outfile", dest="out_file", help="The prefix of output files.")
    parser.add_argument(
        "--iformat",
        dest="in_format",
        default="fq",
        help="The format of input file. Must be 'fq' or 'fa'. default='%(default)s'",
    )
    parser.add_argument(
        "--oformat",
        dest="out_format",
        default="pdf",
        help="The format of output logo file. Must be 'pdf', 'png' or 'svg'. default='%(default)s'",
    )
    parser.add_argument(
        "-n",
        "--nseq-limit",
        type=int,
        dest="max_seq",
        default=None,
        help="Only process this many sequences and stop. default=%(default)s (generate logo from\
ALL sequences).",
    )
    parser.add_argument(
        "--font-name",
        dest="font_name",
        default="sans",
        help="The font of logo characters. For a list of valid font names, run logomaker.list_font_names().\
default='%(default)s'",
    )
    parser.add_argument(
        "--stack-order",
        dest="stack_order",
        default="big_on_top",
        help="Must be 'big_on_top', 'small_on_top', or 'fixed'. 'big_on_top' : nucleotide with the highest\
frequency will be on the top; 'small_on_top' : nucleotide with the lowest frequency will be on the\
top; 'fixed' : nucleotides from top to bottom are in the same order as characters appear in the\
data frame. default='%(default)s'",
    )
    parser.add_argument(
        "--flip-below",
        action="store_true",
        dest="flip_below",
        default=False,
        help="If set, characters below the X-axis (which correspond to negative values in the matrix)\
will be flipped upside down. default=%(default)s",
    )
    parser.add_argument(
        "--shade-below",
        type=float,
        dest="shade_below",
        default=0.0,
        help="The amount of shading to use for characters drawn below the X-axis. 0 <= shade_below <= 1.\
Larger values correspond to more shading. default=%(default)s",
    )
    parser.add_argument(
        "--fade-below",
        type=float,
        dest="fade_below",
        default=0.0,
        help="The amount of fading to use for characters drawn below the X-axis. 0 <= shade_below <= 1.\
 Larger values correspond to more fading. default=%(default)s",
    )
    parser.add_argument(
        "--excludeN",
        action="store_true",
        dest="exclude_N",
        default=False,
        help='If set, exclude all DNA sequences containing "N".',
    )
    parser.add_argument(
        "--highlight-start",
        type=int,
        dest="hi_start",
        default=None,
        help="Highlight logo from this position. Must be within [0, sequence_length-1].\
default=%(default)s (no highlight)",
    )
    parser.add_argument(
        "--highlight-end",
        type=int,
        dest="hi_end",
        default=None,
        help="Highlight logo to this position. Must be within [0, len(logo)-1].\
default=%(default)s (no highlight)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="debug",
        default=False,
        help="If set, print detailed information for debugging.",
    )

    args = parser.parse_args()

    # DEGUB->INFO->WARNING->ERROR->CRITICAL
    if args.debug:
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s]  %(message)s", datefmt="%Y-%m-%d %I:%M:%S", level=logging.DEBUG
        )
    else:
        logging.basicConfig(
            format="%(asctime)s [%(levelname)s]  %(message)s", datefmt="%Y-%m-%d %I:%M:%S", level=logging.INFO
        )

    for file in [args.in_file, args.out_file]:
        if not (file):
            parser.print_help()
            sys.exit(1)
    if args.in_format not in ["fa", "fq"]:
        logging.error("--format takes 'fq' or 'fa' as argument.")
        parser.print_help()
        sys.exit(1)
    if args.shade_below < 0 or args.shade_below > 1:
        logging.error("The shade_below value must be between 0 and 1")
        sys.exit(1)
    if args.fade_below < 0 or args.fade_below > 1:
        logging.error("The fade_below value must be between 0 and 1")
        sys.exit(1)

    if args.in_format.lower() == "fq":
        file_iter = fastq.fastq_iter(args.in_file, mode="seq")
    else:
        file_iter = fastq.fasta_iter(args.in_file)
    mat = fastq.seq2countMat(file_iter, step_size=10000, exclude_N=args.exclude_N, limit=args.max_seq)
    fastq.write_matrix_csv(mat, args.out_file + ".count_matrix.csv", index_label="Index")
    fastq.make_logo(
        mat,
        outfile=args.out_file,
        exclude_N=args.exclude_N,
        font_name=args.font_name,
        stack_order=args.stack_order,
        flip_below=args.flip_below,
        shade_below=args.shade_below,
        fade_below=args.fade_below,
        highlight_start=args.hi_start,
        highlight_end=args.hi_end,
        oformat=args.out_format,
    )


if __name__ == "__main__":
    main()
