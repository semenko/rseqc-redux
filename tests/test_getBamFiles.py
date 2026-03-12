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


def test_get_bam_files_single_bam(tmp_path):
    """Single BAM file passed directly."""
    bam = tmp_path / "test.bam"
    bam.write_text("data")
    (tmp_path / "test.bam.bai").write_text("index")
    result = getBamFiles.get_bam_files(str(bam))
    assert len(result) == 1
    assert result[0] == str(bam)


def test_get_bam_files_text_file_with_comments(tmp_path):
    """Text file with comment lines should skip them."""
    bam = tmp_path / "test.bam"
    bam.write_text("data")
    (tmp_path / "test.bam.bai").write_text("index")
    listing = tmp_path / "bams.txt"
    listing.write_text(f"# this is a comment\n{bam}\n# another comment\n")
    result = getBamFiles.get_bam_files(str(listing))
    assert len(result) == 1


def test_get_bam_files_comma_separated(tmp_path):
    """Comma-separated BAM files."""
    bam1 = tmp_path / "a.bam"
    bam1.write_text("data")
    (tmp_path / "a.bam.bai").write_text("index")
    bam2 = tmp_path / "b.bam"
    bam2.write_text("data")
    (tmp_path / "b.bam.bai").write_text("index")
    result = getBamFiles.get_bam_files(f"{bam1},{bam2}")
    assert len(result) == 2


def test_get_bam_files_nonexistent():
    """Non-existent path that isn't a dir, file, or bam."""
    result = getBamFiles.get_bam_files("/nonexistent/path/to/nothing")
    assert result == []


def test_get_bam_files_printit(tmp_path, capsys):
    """printit=True should print file paths."""
    bam = tmp_path / "test.bam"
    bam.write_text("data")
    (tmp_path / "test.bam.bai").write_text("index")
    getBamFiles.get_bam_files(str(bam), printit=True)
    captured = capsys.readouterr()
    assert str(bam) in captured.out
