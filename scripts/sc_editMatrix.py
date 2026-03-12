#!/usr/bin/env python3
"""
This program generates heatmaps to visualize the **positions** (X-axis),
type of edits (Y-axis, such as "C" to "T") and  **frequencies** (color)
of error-corrected nucleotides in cell barcodes and UMIs.

"""

import argparse
import logging
import sys

from rseqc import heatmap, scbam


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version="5.0.2")
    parser.add_argument("-i", "--infile", dest="in_file", help="Input file in BAM foramt.")
    parser.add_argument("-o", "--outfile", dest="out_file", help="The prefix of output files.")
    parser.add_argument(
        "--limit",
        type=int,
        dest="reads_num",
        default=None,
        help="Number of alignments to process. default=%(default)s",
    )
    parser.add_argument(
        "--cr-tag",
        dest="CR_tag",
        default="CR",
        help="Tag of cellular barcode reported by the sequencer in BAM file. default='%(default)s'",
    )
    parser.add_argument(
        "--cb-tag",
        dest="CB_tag",
        default="CB",
        help="Tag of error-corrected cellular barcode in BAM file. default='%(default)s'",
    )
    parser.add_argument(
        "--ur-tag",
        dest="UR_tag",
        default="UR",
        help="Tag of UMI reported by the sequencer in BAM file. default='%(default)s'",
    )
    parser.add_argument(
        "--ub-tag",
        dest="UB_tag",
        default="UB",
        help="Tag of error-corrected UMI in BAM file. default='%(default)s'",
    )
    parser.add_argument(
        "--cell-width",
        type=int,
        dest="cell_width",
        default=15,
        help="Points of cell width in the heatmap. default=%(default)s",
    )
    parser.add_argument(
        "--cell-height",
        type=int,
        dest="cell_height",
        default=10,
        help="Points of cell height in the heatmap. default=%(default)s",
    )
    parser.add_argument(
        "--font-size",
        type=int,
        dest="font_size",
        default=8,
        help="Font size. If --display-num was set, fontsize_number = 0.8 * font_size. default=%(default)s",
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
        "--verbose",
        action="store_true",
        dest="debug",
        default=False,
        help="If set, detailed running information is printed to screen.",
    )
    parser.add_argument(
        "--no-num",
        action="store_true",
        dest="no_num",
        default=False,
        help="If set, will not print numerical values to cells. default=%(default)s",
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

    scbam.barcode_edits(
        infile=args.in_file,
        outfile=args.out_file,
        limit=args.reads_num,
        CR_tag=args.CR_tag,
        CB_tag=args.CB_tag,
        UR_tag=args.UR_tag,
        UB_tag=args.UB_tag,
    )

    CB_mat_file = args.out_file + ".CB_edits_count.csv"
    CB_heatmap_file = args.out_file + ".CB_edits_heatmap"
    heatmap.make_heatmap(
        infile=CB_mat_file,
        outfile=CB_heatmap_file,
        filetype=args.file_type,
        cell_width=args.cell_width,
        cell_height=args.cell_height,
        col_angle=args.col_angle,
        font_size=args.font_size,
        text_color=args.text_color,
        no_numbers=args.no_num,
        log2_scale=True,
    )

    UMI_mat_file = args.out_file + ".UMI_edits_count.csv"
    UMI_heatmap_file = args.out_file + ".UMI_edits_heatmap"
    heatmap.make_heatmap(
        infile=UMI_mat_file,
        outfile=UMI_heatmap_file,
        filetype=args.file_type,
        cell_width=args.cell_width,
        cell_height=args.cell_height,
        col_angle=args.col_angle,
        font_size=args.font_size,
        text_color=args.text_color,
        no_numbers=args.no_num,
        log2_scale=True,
    )


if __name__ == "__main__":
    main()
