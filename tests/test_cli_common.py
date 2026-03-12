"""Tests for rseqc.cli_common."""

import importlib
import os
import tempfile

from rseqc.cli_common import build_bitsets, load_chromsize, printlog


def test_import():
    """Verify that rseqc.cli_common can be imported."""
    mod = importlib.import_module("rseqc.cli_common")
    assert mod is not None


def test_printlog(capsys):
    """printlog writes timestamped message to stderr."""
    printlog("test message")
    captured = capsys.readouterr()
    assert "test message" in captured.err
    assert "@" in captured.err


def test_build_bitsets_basic():
    """build_bitsets creates interval trees from entries."""
    entries = [["chr1", 100, 200], ["chr1", 300, 400], ["chr2", 500, 600]]
    result = build_bitsets(entries)
    assert "CHR1" in result
    assert "CHR2" in result
    assert len(result["CHR1"].find(150, 160)) == 1
    assert len(result["CHR1"].find(250, 260)) == 0
    assert len(result["CHR2"].find(550, 560)) == 1


def test_build_bitsets_uppercases_chrom():
    """build_bitsets uppercases chromosome names."""
    entries = [["chrX", 0, 100]]
    result = build_bitsets(entries)
    assert "CHRX" in result
    assert len(result["CHRX"].find(50, 60)) == 1


def test_build_bitsets_converts_to_int():
    """build_bitsets handles string coordinates."""
    entries = [["chr1", "100", "200"]]
    result = build_bitsets(entries)
    assert len(result["CHR1"].find(150, 160)) == 1


def test_load_chromsize():
    """load_chromsize reads tab-separated chrom sizes."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
        f.write("chr1\t248956422\n")
        f.write("chr2\t242193529\n")
        f.write("# comment line\n")
        f.write("\n")  # blank line
        f.name
    try:
        result = load_chromsize(f.name)
        assert result == {"chr1": 248956422, "chr2": 242193529}
    finally:
        os.unlink(f.name)
