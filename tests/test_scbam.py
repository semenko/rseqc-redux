"""Tests for rseqc.scbam — pure-logic functions only (no BAM files needed)."""

import importlib

from rseqc import scbam


def test_import():
    """Verify that rseqc.scbam can be imported."""
    mod = importlib.import_module("rseqc.scbam")
    assert mod is not None


# --- diff_str ---


def test_diff_str_identical():
    assert scbam.diff_str("ACGT", "ACGT") == []


def test_diff_str_single_diff():
    result = scbam.diff_str("ACGT", "ACTT")
    assert len(result) == 1
    assert result[0] == [2, "G", "T"]


def test_diff_str_multiple_diffs():
    result = scbam.diff_str("AAAA", "TTTT")
    assert len(result) == 4


def test_diff_str_different_lengths():
    assert scbam.diff_str("AC", "ACG") == []


def test_diff_str_empty():
    assert scbam.diff_str("", "") == []


# --- read_match_type ---


def test_read_match_type_consecutive():
    assert scbam.read_match_type("50M") == "Map_consecutively"


def test_read_match_type_splicing():
    assert scbam.read_match_type("25M1000N25M") == "Map_with_splicing"


def test_read_match_type_clip_left():
    assert scbam.read_match_type("5S45M") == "Map_with_clipping"


def test_read_match_type_clip_right():
    assert scbam.read_match_type("45M5S") == "Map_with_clipping"


def test_read_match_type_splice_clip_right():
    assert scbam.read_match_type("20M500N25M5S") == "Map_with_splicing_and_clipping"


def test_read_match_type_splice_clip_left():
    assert scbam.read_match_type("5S20M500N25M") == "Map_with_splicing_and_clipping"


def test_read_match_type_others():
    assert scbam.read_match_type("10M5I35M") == "Others"


def test_read_match_type_empty():
    assert scbam.read_match_type("") == "Others"


# --- list2str ---


def test_list2str_simple():
    assert scbam.list2str([(0, 50)]) == "50M"


def test_list2str_splice():
    assert scbam.list2str([(0, 25), (3, 1000), (0, 25)]) == "25M1000N25M"


def test_list2str_clip():
    assert scbam.list2str([(4, 5), (0, 45)]) == "5S45M"


def test_list2str_insertion():
    assert scbam.list2str([(0, 30), (1, 5), (0, 20)]) == "30M5I20M"


def test_list2str_deletion():
    assert scbam.list2str([(0, 30), (2, 5), (0, 20)]) == "30M5D20M"


def test_list2str_all_ops():
    # Test all CIGAR operation codes
    ops = [(0, 10), (1, 2), (2, 3), (3, 100), (4, 5), (5, 1), (6, 1), (7, 10), (8, 2)]
    result = scbam.list2str(ops)
    assert result == "10M2I3D100N5S1H1P10=2X"


def test_list2str_empty():
    assert scbam.list2str([]) == ""
