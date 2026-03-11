"""Tests for rseqc.heatmap."""

import importlib

from rseqc import heatmap


def test_import():
    """Verify that rseqc.heatmap can be imported."""
    mod = importlib.import_module("rseqc.heatmap")
    assert mod is not None


def _read_r_script(outfile):
    """Helper to read generated R script content."""
    with open(outfile + ".r") as f:
        return f.read()


def test_make_heatmap_basic(tmp_path):
    """make_heatmap generates an R script with pheatmap call."""
    outfile = str(tmp_path / "test_heatmap")
    heatmap.make_heatmap(
        infile="input.csv",
        outfile=outfile,
        filetype="png",
        cell_width=20,
        cell_height=20,
        col_angle=45,
        font_size=10,
        text_color="black",
    )
    content = _read_r_script(outfile)
    assert "library(pheatmap)" in content
    assert "input.csv" in content
    assert "display_numbers = TRUE" in content
    assert "number_color='black'" in content


def test_make_heatmap_no_numbers(tmp_path):
    """no_numbers=True produces display_numbers = FALSE."""
    outfile = str(tmp_path / "test_heatmap")
    heatmap.make_heatmap(
        infile="input.csv",
        outfile=outfile,
        filetype="pdf",
        cell_width=15,
        cell_height=15,
        col_angle=90,
        font_size=12,
        text_color="red",
        no_numbers=True,
    )
    content = _read_r_script(outfile)
    assert "display_numbers = FALSE" in content
    assert "number_color" not in content


def test_make_heatmap_log2_scale(tmp_path):
    """log2_scale=True produces log2 transformation in R code."""
    outfile = str(tmp_path / "test_heatmap")
    heatmap.make_heatmap(
        infile="input.csv",
        outfile=outfile,
        filetype="png",
        cell_width=20,
        cell_height=20,
        col_angle=45,
        font_size=10,
        text_color="black",
        log2_scale=True,
    )
    content = _read_r_script(outfile)
    assert "log2(as.matrix(dat)+1)" in content


def test_make_heatmap_filetype_in_output(tmp_path):
    """Output filename includes the correct filetype."""
    outfile = str(tmp_path / "test_heatmap")
    heatmap.make_heatmap(
        infile="input.csv",
        outfile=outfile,
        filetype="tiff",
        cell_width=20,
        cell_height=20,
        col_angle=45,
        font_size=10,
        text_color="black",
    )
    content = _read_r_script(outfile)
    assert "test_heatmap.tiff" in content


def test_make_heatmap_cell_dimensions(tmp_path):
    """Cell dimensions appear in R script."""
    outfile = str(tmp_path / "test_heatmap")
    heatmap.make_heatmap(
        infile="input.csv",
        outfile=outfile,
        filetype="png",
        cell_width=30,
        cell_height=25,
        col_angle=270,
        font_size=14,
        text_color="blue",
    )
    content = _read_r_script(outfile)
    assert "cellwidth = 30" in content
    assert "cellheight = 25" in content
    assert "angle_col=270" in content
    assert "fontsize=14" in content
