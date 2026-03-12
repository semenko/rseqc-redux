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


# --- Edge cases: CIGAR with only insertions ---


def test_fetch_exon_insertion_only():
    """Insertion-only CIGAR produces no exon bounds."""
    result = bam_cigar.fetch_exon("chr1", 100, [(1, 10)])
    assert result == []


def test_fetch_exon_hard_clip():
    """Hard clipping (op 5) falls through to continue."""
    result = bam_cigar.fetch_exon("chr1", 100, [(5, 3), (0, 10), (5, 3)])
    assert result == [("chr1", 100, 110)]


def test_fetch_intron_with_insertion():
    """Insertion in intron-containing CIGAR is skipped."""
    result = bam_cigar.fetch_intron("chr1", 100, [(0, 10), (1, 5), (3, 500), (0, 10)])
    assert result == [("chr1", 110, 610)]


def test_fetch_intron_with_deletion():
    """Deletion advances position but produces no intron."""
    result = bam_cigar.fetch_intron("chr1", 100, [(0, 10), (2, 5), (0, 10)])
    assert result == []


def test_fetch_intron_with_soft_clip():
    """Soft clip does not advance position in intron context."""
    result = bam_cigar.fetch_intron("chr1", 100, [(4, 5), (0, 10), (3, 500), (0, 10)])
    assert result == [("chr1", 110, 610)]


def test_fetch_intron_hard_clip():
    """Hard clip (op 5) falls through to continue."""
    result = bam_cigar.fetch_intron("chr1", 100, [(5, 3), (0, 10), (3, 500), (0, 10)])
    assert result == [("chr1", 110, 610)]


def test_fetch_clip_with_insertion():
    """Insertion doesn't produce clip bounds."""
    result = bam_cigar.fetch_clip("chr1", 100, [(1, 5), (0, 10)])
    assert result == []


def test_fetch_clip_with_deletion():
    """Deletion advances position but no clip."""
    result = bam_cigar.fetch_clip("chr1", 100, [(0, 10), (2, 3), (4, 5)])
    assert result == [("chr1", 113, 118)]


def test_fetch_clip_with_gap():
    """Gap (intron) advances position."""
    result = bam_cigar.fetch_clip("chr1", 100, [(0, 10), (3, 500), (0, 10), (4, 5)])
    assert result == [("chr1", 620, 625)]


def test_fetch_clip_hard_clip():
    """Hard clip (op 5) falls through."""
    result = bam_cigar.fetch_clip("chr1", 100, [(5, 3), (4, 5), (0, 10)])
    assert result == [("chr1", 100, 105)]


def test_fetch_deletion_with_insertion():
    """Insertion doesn't affect deletion detection."""
    result = bam_cigar.fetch_deletion("chr1", 100, [(0, 10), (1, 5), (2, 3), (0, 10)])
    assert result == [("chr1", 110, 113)]


def test_fetch_deletion_with_gap():
    """Gap advances position."""
    result = bam_cigar.fetch_deletion("chr1", 100, [(0, 10), (3, 500), (2, 3), (0, 10)])
    assert result == [("chr1", 610, 613)]


def test_fetch_deletion_with_soft_clip():
    """Soft clip advances position."""
    result = bam_cigar.fetch_deletion("chr1", 100, [(4, 5), (0, 10), (2, 3), (0, 10)])
    assert result == [("chr1", 115, 118)]


def test_fetch_deletion_hard_clip():
    """Hard clip falls through."""
    result = bam_cigar.fetch_deletion("chr1", 100, [(5, 3), (0, 10), (2, 3)])
    assert result == [("chr1", 110, 113)]


# --- fetch_deletion_range edge cases ---


def test_fetch_deletion_range_with_soft_clip():
    result = bam_cigar.fetch_deletion_range([(4, 5), (0, 10), (2, 3), (0, 10)])
    assert result == [(15, 3)]


def test_fetch_deletion_range_with_insertion():
    result = bam_cigar.fetch_deletion_range([(0, 10), (1, 5), (2, 3), (0, 10)])
    assert result == [(15, 3)]


def test_fetch_deletion_range_with_gap():
    """Gap (op 3) is skipped in range calculation."""
    result = bam_cigar.fetch_deletion_range([(0, 10), (3, 500), (2, 3)])
    assert result == [(10, 3)]


def test_fetch_deletion_range_hard_clip():
    result = bam_cigar.fetch_deletion_range([(5, 3), (0, 10), (2, 3)])
    assert result == [(10, 3)]


# --- fetch_insertion_range edge cases ---


def test_fetch_insertion_range_with_soft_clip():
    result = bam_cigar.fetch_insertion_range([(4, 5), (0, 10), (1, 3)])
    assert result == [(15, 3)]


def test_fetch_insertion_range_with_deletion():
    """Deletion (op 2) is skipped."""
    result = bam_cigar.fetch_insertion_range([(0, 10), (2, 5), (1, 3)])
    assert result == [(10, 3)]


def test_fetch_insertion_range_with_gap():
    """Gap (op 3) is skipped."""
    result = bam_cigar.fetch_insertion_range([(0, 10), (3, 500), (1, 3)])
    assert result == [(10, 3)]


def test_fetch_insertion_range_hard_clip():
    result = bam_cigar.fetch_insertion_range([(5, 3), (0, 10), (1, 3)])
    assert result == [(10, 3)]


# --- fetch_insertion edge cases ---


def test_fetch_insertion_with_deletion():
    result = bam_cigar.fetch_insertion("chr1", 100, [(0, 10), (2, 5), (1, 3), (0, 10)])
    assert result == [("chr1", 115, 3)]


def test_fetch_insertion_with_gap():
    result = bam_cigar.fetch_insertion("chr1", 100, [(0, 10), (3, 500), (1, 3), (0, 10)])
    assert result == [("chr1", 610, 3)]


def test_fetch_insertion_with_soft_clip():
    result = bam_cigar.fetch_insertion("chr1", 100, [(4, 5), (0, 10), (1, 3)])
    assert result == [("chr1", 115, 3)]


def test_fetch_insertion_hard_clip():
    result = bam_cigar.fetch_insertion("chr1", 100, [(5, 3), (0, 10), (1, 3)])
    assert result == [("chr1", 110, 3)]


# --- map_bounds edge cases ---


def test_map_bounds_empty_cigar():
    assert bam_cigar.map_bounds(100, []) == (100, 100)


def test_map_bounds_hard_clip_only():
    """Hard clip (op 5) is ignored."""
    assert bam_cigar.map_bounds(100, [(5, 10)]) == (100, 100)


def test_map_bounds_complex():
    """S+M+I+D+N+M+S — only M/D/N contribute to span."""
    cigar = [(4, 3), (0, 10), (1, 5), (2, 3), (3, 1000), (0, 10), (4, 3)]
    assert bam_cigar.map_bounds(100, cigar) == (100, 1123)


# --- list2longstr edge cases ---


def test_list2longstr_with_equals_and_mismatch():
    """Op 7 (=) and 8 (X) should be included."""
    assert bam_cigar.list2longstr([(7, 3), (8, 2)]) == "===XX"


def test_list2longstr_deletion_skipped():
    """Op 2 (D) should be excluded."""
    assert bam_cigar.list2longstr([(0, 2), (2, 5), (0, 2)]) == "MMMM"
