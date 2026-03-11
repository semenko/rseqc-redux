"""Tests for rseqc.ireader."""

import gzip
import importlib

from rseqc import ireader


def test_import():
    """Verify that rseqc.ireader can be imported."""
    mod = importlib.import_module("rseqc.ireader")
    assert mod is not None


# --- reader with plain text ---


def test_reader_plain_text(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("line1\nline2\nline3\n")
    lines = list(ireader.reader(str(f)))
    assert lines == ["line1", "line2", "line3"]


def test_reader_plain_text_empty(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("")
    lines = list(ireader.reader(str(f)))
    assert lines == []


def test_reader_strips_carriage_return(tmp_path):
    f = tmp_path / "crlf.txt"
    f.write_bytes(b"line1\r\nline2\r\n")
    lines = list(ireader.reader(str(f)))
    assert lines == ["line1", "line2"]


# --- reader with gzip ---


def test_reader_gzip(tmp_path):
    f = tmp_path / "test.gz"
    with gzip.open(str(f), "wb") as gz:
        gz.write(b"gzline1\ngzline2\n")
    lines = list(ireader.reader(str(f)))
    assert lines == ["gzline1", "gzline2"]


# --- nopen ---


def test_nopen_plain(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("hello\n")
    fh = ireader.nopen(str(f), "rb")
    content = fh.read()
    fh.close()
    assert b"hello" in content


def test_nopen_passthrough():
    # passing a non-string should return it as-is
    import io

    buf = io.BytesIO(b"test")
    assert ireader.nopen(buf) is buf
