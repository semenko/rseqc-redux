"""Tests for rseqc.bam_cigar — pure functions operating on CIGAR tuples."""

import importlib

from rseqc import bam_cigar


def test_import():
    """Verify that rseqc.bam_cigar can be imported."""
    mod = importlib.import_module("rseqc.bam_cigar")
    assert mod is not None


# --- map_bounds ---


def test_map_bounds_simple_match():
    assert bam_cigar.map_bounds(100, [(0, 50)]) == (100, 150)


def test_map_bounds_with_intron():
    assert bam_cigar.map_bounds(100, [(0, 25), (3, 1000), (0, 25)]) == (100, 1150)


def test_map_bounds_insertion_ignored():
    assert bam_cigar.map_bounds(100, [(0, 10), (1, 5), (0, 10)]) == (100, 120)


def test_map_bounds_deletion():
    assert bam_cigar.map_bounds(100, [(0, 10), (2, 3), (0, 10)]) == (100, 123)


def test_map_bounds_soft_clip_ignored():
    assert bam_cigar.map_bounds(100, [(4, 5), (0, 45)]) == (100, 145)


# --- fetch_exon ---


def test_fetch_exon_simple():
    result = bam_cigar.fetch_exon("chr1", 100, [(0, 50)])
    assert result == [("chr1", 100, 150)]


def test_fetch_exon_spliced():
    result = bam_cigar.fetch_exon("chr1", 100, [(0, 25), (3, 1000), (0, 25)])
    assert result == [("chr1", 100, 125), ("chr1", 1125, 1150)]


def test_fetch_exon_with_insertion():
    result = bam_cigar.fetch_exon("chr1", 100, [(0, 10), (1, 5), (0, 10)])
    assert result == [("chr1", 100, 110), ("chr1", 110, 120)]


def test_fetch_exon_with_deletion():
    result = bam_cigar.fetch_exon("chr1", 100, [(0, 10), (2, 3), (0, 10)])
    assert result == [("chr1", 100, 110), ("chr1", 113, 123)]


def test_fetch_exon_with_soft_clip():
    result = bam_cigar.fetch_exon("chr1", 100, [(4, 5), (0, 45)])
    assert result == [("chr1", 105, 150)]


# --- fetch_intron ---


def test_fetch_intron_none():
    assert bam_cigar.fetch_intron("chr1", 100, [(0, 50)]) == []


def test_fetch_intron_single():
    result = bam_cigar.fetch_intron("chr1", 100, [(0, 25), (3, 1000), (0, 25)])
    assert result == [("chr1", 125, 1125)]


def test_fetch_intron_multiple():
    result = bam_cigar.fetch_intron("chr1", 100, [(0, 10), (3, 500), (0, 10), (3, 300), (0, 10)])
    assert result == [("chr1", 110, 610), ("chr1", 620, 920)]


# --- fetch_clip ---


def test_fetch_clip_none():
    assert bam_cigar.fetch_clip("chr1", 100, [(0, 50)]) == []


def test_fetch_clip_head():
    result = bam_cigar.fetch_clip("chr1", 100, [(4, 5), (0, 45)])
    assert result == [("chr1", 100, 105)]


def test_fetch_clip_tail():
    result = bam_cigar.fetch_clip("chr1", 100, [(0, 45), (4, 5)])
    assert result == [("chr1", 145, 150)]


def test_fetch_clip_both():
    result = bam_cigar.fetch_clip("chr1", 100, [(4, 3), (0, 44), (4, 3)])
    assert len(result) == 2


# --- fetch_deletion ---


def test_fetch_deletion_none():
    assert bam_cigar.fetch_deletion("chr1", 100, [(0, 50)]) == []


def test_fetch_deletion_single():
    result = bam_cigar.fetch_deletion("chr1", 100, [(0, 25), (2, 3), (0, 25)])
    assert result == [("chr1", 125, 128)]


# --- fetch_insertion ---


def test_fetch_insertion_none():
    assert bam_cigar.fetch_insertion("chr1", 100, [(0, 50)]) == []


def test_fetch_insertion_single():
    result = bam_cigar.fetch_insertion("chr1", 100, [(0, 25), (1, 5), (0, 25)])
    assert result == [("chr1", 125, 5)]


# --- fetch_deletion_range ---


def test_fetch_deletion_range_none():
    assert bam_cigar.fetch_deletion_range([(0, 50)]) == []


def test_fetch_deletion_range_single():
    result = bam_cigar.fetch_deletion_range([(0, 25), (2, 3), (0, 25)])
    assert result == [(25, 3)]


# --- fetch_insertion_range ---


def test_fetch_insertion_range_none():
    assert bam_cigar.fetch_insertion_range([(0, 50)]) == []


def test_fetch_insertion_range_single():
    result = bam_cigar.fetch_insertion_range([(0, 25), (1, 5), (0, 25)])
    assert result == [(25, 5)]


# --- list2str ---


def test_list2str():
    assert bam_cigar.list2str([(4, 1), (0, 9)]) == "1S9M"


def test_list2str_complex():
    assert bam_cigar.list2str([(0, 25), (3, 1000), (0, 25)]) == "25M1000N25M"


def test_list2str_empty():
    assert bam_cigar.list2str([]) == ""


# --- list2longstr ---


def test_list2longstr_simple():
    assert bam_cigar.list2longstr([(0, 3)]) == "MMM"


def test_list2longstr_with_clip():
    assert bam_cigar.list2longstr([(4, 2), (0, 3)]) == "SSMMM"


def test_list2longstr_skips_gaps():
    assert bam_cigar.list2longstr([(0, 2), (3, 100), (0, 2)]) == "MMMM"


def test_list2longstr_insertion():
    assert bam_cigar.list2longstr([(0, 2), (1, 3), (0, 2)]) == "MMIIIMM"
