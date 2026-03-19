"""Tests for rseqc.cli_common."""

import importlib
import os
import tempfile

from rseqc.cli_common import build_bitsets, iter_bed12, load_chromsize, printlog


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


# --- iter_bed12 ---


def test_iter_bed12_basic(tmp_path):
    """iter_bed12 parses a standard BED12 line."""
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n")
    records = list(iter_bed12(str(bed)))
    assert len(records) == 1
    r = records[0]
    assert r.chrom == "chr1"
    assert r.tx_start == 1000
    assert r.tx_end == 5000
    assert r.gene_name == "gene1"
    assert r.strand == "+"
    assert r.exon_starts == [1000, 2500, 4500]
    assert r.exon_ends == [1500, 3100, 5000]


def test_iter_bed12_skips_comments(tmp_path):
    """iter_bed12 skips #, track, and browser lines."""
    bed = tmp_path / "test.bed"
    bed.write_text(
        "# comment\n"
        "track name=test\n"
        "browser position chr1:1-100\n"
        "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n"
    )
    records = list(iter_bed12(str(bed)))
    assert len(records) == 1
    assert records[0].gene_name == "gene1"


def test_iter_bed12_skips_malformed(tmp_path, capsys):
    """iter_bed12 skips malformed lines with a warning."""
    bed = tmp_path / "test.bed"
    bed.write_text("bad line\nchr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t1\t4000,\t0,\n")
    records = list(iter_bed12(str(bed)))
    assert len(records) == 1
    captured = capsys.readouterr()
    assert "skipped this line" in captured.err


def test_iter_bed12_empty_file(tmp_path):
    """iter_bed12 yields nothing for an empty file."""
    bed = tmp_path / "test.bed"
    bed.write_text("")
    records = list(iter_bed12(str(bed)))
    assert records == []


def test_iter_bed12_fields_accessible(tmp_path):
    """iter_bed12 exposes raw fields for callers needing extra columns."""
    bed = tmp_path / "test.bed"
    bed.write_text("chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500\n")
    records = list(iter_bed12(str(bed)))
    r = records[0]
    assert r.fields[0] == "chr1"
    assert r.fields[3] == "gene1"
    # fields[10] contains block sizes
    assert r.fields[10] == "500,600,500"
