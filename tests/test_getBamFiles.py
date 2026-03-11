"""Tests for rseqc.getBamFiles."""

import importlib

from rseqc import getBamFiles


def test_import():
    """Verify that rseqc.getBamFiles can be imported."""
    mod = importlib.import_module("rseqc.getBamFiles")
    assert mod is not None


# --- isbamfile ---


def test_isbamfile_nonexistent():
    assert getBamFiles.isbamfile("/nonexistent/file.bam") is False


def test_isbamfile_not_bam(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("hello")
    assert getBamFiles.isbamfile(str(f)) is False


def test_isbamfile_empty_bam(tmp_path):
    f = tmp_path / "test.bam"
    f.write_text("")
    assert getBamFiles.isbamfile(str(f)) is False


def test_isbamfile_no_index(tmp_path):
    f = tmp_path / "test.bam"
    f.write_text("data")
    assert getBamFiles.isbamfile(str(f)) is False


def test_isbamfile_with_index(tmp_path):
    f = tmp_path / "test.bam"
    f.write_text("data")
    idx = tmp_path / "test.bam.bai"
    idx.write_text("index")
    assert getBamFiles.isbamfile(str(f)) is True


# --- get_bam_files ---


def test_get_bam_files_empty_dir(tmp_path):
    result = getBamFiles.get_bam_files(str(tmp_path))
    assert result == []


def test_get_bam_files_with_bams(tmp_path):
    for name in ["a.bam", "b.bam"]:
        (tmp_path / name).write_text("data")
        (tmp_path / (name + ".bai")).write_text("index")
    result = getBamFiles.get_bam_files(str(tmp_path))
    assert len(result) == 2


def test_get_bam_files_from_text_file(tmp_path):
    bam = tmp_path / "test.bam"
    bam.write_text("data")
    (tmp_path / "test.bam.bai").write_text("index")
    listing = tmp_path / "bams.txt"
    listing.write_text(str(bam) + "\n")
    result = getBamFiles.get_bam_files(str(listing))
    assert len(result) == 1
