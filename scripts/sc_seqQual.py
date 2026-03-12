#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This program generates heatmap from a FASTQ file to visualize the sequencing quality.
"""

import argparse
import logging
import sys

from rseqc import fastq, heatmap

__contributor__ = "Liguo Wang"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument(
        "-i",
        "--infile",
        dest="in_file",
        help="Input file in FASTQ (https://en.wikipedia.org/wiki/FASTQ_format#) format.",
    )
    parser.add_argument("-o", "--outfile", dest="out_file", help="The prefix of output files.")
    parser.add_argument(
        "-n",
        "--nseq-limit",
        type=int,
        dest="max_seq",
        default=None,
        help="Only process this many sequences and stop. default=%(default)s (generate logo from ALL sequences).",
    )
    parser.add_argument(
        "--cell-width",
        type=int,
        dest="cell_width",
        default=12,
        help="Cell width (in points) of the heatmap. default=%(default)s",
    )
    parser.add_argument(
        "--cell-height",
        type=int,
        dest="cell_height",
        default=10,
        help="Cell height (in points) of the heatmap. default=%(default)s",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        dest="font_size",
        default=6,
        help="Font size in points. If --display-num was set, fontsize_number = 0.8 * font_size. default=%(default)s",
    )
    parser.add_argument(
        "--angle",
        type=int,
        dest="col_angle",
        default=45,
        help="The angle (must be 0, 45, 90, 270, 315) of column text lables under the heatmap. default=%(default)s",
    )
    parser.add_argument(
        "--text-color",
        dest="text_color",
        default="black",
        help="The color of numbers in each cell. default=%(default)s",
    )
    parser.add_argument(
        "--file-type",
        dest="file_type",
        default="pdf",
        help="The file type of heatmap. Choose one of 'pdf', 'png', 'tiff', 'bmp', 'jpeg'. default=%(default)s",
    )
    parser.add_argument(
        "--no-num",
        action="store_true",
        dest="no_num",
        default=False,
        help="if set, will not print numerical values to cells. default=%(default)s",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        dest="debug",
        default=False,
        help="If set, will produce detailed information for debugging.",
    )

    args = parser.parse_args()
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

    file_iter = fastq.fastq_iter(args.in_file, mode="qual")
    qual_mat = fastq.qual2countMat(file_iter, limit=args.max_seq)

    # Write count CSV (rows=quality_scores descending, cols=positions)
    fastq.write_matrix_csv(
        qual_mat,
        args.out_file + ".qual_count.csv",
        index_label="Index",
        transpose=True,
        sort_index_descending=True,
    )
    # Write percent CSV (same layout, normalized per column)
    fastq.write_matrix_csv(
        qual_mat,
        args.out_file + ".qual_percent.csv",
        index_label="Index",
        transpose=True,
        sort_index_descending=True,
        normalize=True,
    )

    heatmap.make_heatmap(
        infile=args.out_file + ".qual_percent.csv",
        outfile=args.out_file + ".qual_heatmap",
        filetype=args.file_type,
        cell_width=args.cell_width,
        cell_height=args.cell_height,
        col_angle=args.col_angle,
        font_size=args.font_size,
        text_color=args.text_color,
        no_numbers=args.no_num,
    )


if __name__ == "__main__":
    main()
