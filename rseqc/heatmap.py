"""
Created on Wed Aug  5 15:24:22 2020

@author: m102324
"""

import logging
import subprocess


def make_heatmap(
    infile: str,
    outfile: str,
    filetype: str,
    cell_width: int,
    cell_height: int,
    col_angle: int,
    font_size: int,
    text_color: str,
    no_numbers: bool = False,
    log2_scale: bool = False,
) -> None:
    """
    Use pheatmap to generate heatmap.

    Parameters
    ----------
    infile : str
            Input file
    outfile : str
            Output file prefix.
    filetype : str
            Filetype is decided by the extension. Currently following formats are supported: png, pdf, tiff, bmp, jpeg.
    cell_width : int
            Points of cell width in the heatmap.
    cell_height : int
            Points of cell height in the heatmap.
    cell_angle : int
            The angle (must be 0, 45, 90, 270, 315) of column text lables under the heatmap.
    font_size : int
            Font size. If --display-num was set, the font size of number = 0.8 * font_size.
    text_color : str
            The color of numbers in each cell.
    dispaly_num : bool
            If set, numbers will be displayed on each cell
    """
    logging.info(f'Writing R code to "{outfile}.r"')
    with open(outfile + ".r", "w") as ROUT:
        print('if(!require(pheatmap)){install.packages("pheatmap")}', file=ROUT)
        print("library(pheatmap)", file=ROUT)
        print(f"dat = read.table(file = '{infile}',sep=',',header=T,row.names=1,check.names=FALSE)", file=ROUT)

        color_ramp = "color = colorRampPalette(c('#f1eef6','#d7b5d8','#df65b0', '#ce1256'))(50)"
        if no_numbers:
            logging.info("Does not displayed numerical values on heatmap")
            if log2_scale:
                print(
                    f"pheatmap(log2(as.matrix(dat)+1), filename='{outfile}.{filetype}', "
                    f"cellwidth = {cell_width}, cellheight = {cell_height}, "
                    f"display_numbers = FALSE, angle_col={col_angle}, fontsize={font_size},"
                    f"cluster_rows=F, cluster_cols=F, scale='none',{color_ramp})",
                    file=ROUT,
                )
            else:
                print(
                    f"pheatmap(as.matrix(dat), filename='{outfile}.{filetype}', "
                    f"cellwidth = {cell_width}, cellheight = {cell_height}, "
                    f"display_numbers = FALSE, angle_col={col_angle}, fontsize={font_size},"
                    f"cluster_rows=F, cluster_cols=F, scale='none',{color_ramp})",
                    file=ROUT,
                )
        else:
            logging.info("Displayed numerical values on heatmap")
            if log2_scale:
                logging.info("Numbers will be displayed on log2 scale")
                print(
                    f"pheatmap(log2(as.matrix(dat)+1), filename='{outfile}.{filetype}', "
                    f"cellwidth = {cell_width}, cellheight = {cell_height}, "
                    f"display_numbers = TRUE, angle_col={col_angle}, fontsize={font_size}, "
                    f"number_color='{text_color}',cluster_rows=F, cluster_cols=F, "
                    f"scale='none',{color_ramp})",
                    file=ROUT,
                )
            else:
                print(
                    f"pheatmap(as.matrix(dat), filename='{outfile}.{filetype}', "
                    f"cellwidth = {cell_width}, cellheight = {cell_height}, "
                    f"display_numbers = TRUE, angle_col={col_angle}, fontsize={font_size}, "
                    f"number_color='{text_color}',cluster_rows=F, cluster_cols=F, "
                    f"scale='none',{color_ramp})",
                    file=ROUT,
                )

    logging.info(f'Running R script file "{outfile}.r"')
    try:
        subprocess.run(["Rscript", outfile + ".r"], check=False)
    except OSError:
        logging.error('Failed to run Rscript file "%s"', outfile + ".r")
