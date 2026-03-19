"""Tests for rseqc.bam_cigar — pure functions operating on CIGAR tuples."""

import importlib

from rseqc import bam_cigar


def test_import():
    """Verify that rseqc.bam_cigar can be imported."""
    mod = importlib.import_module("rseqc.bam_cigar")
    assert mod is not None


# --- fetch_intron ---


def test_fetch_intron_none():
    assert bam_cigar.fetch_intron(100, [(0, 50)]) == []


def test_fetch_intron_single():
    result = bam_cigar.fetch_intron(100, [(0, 25), (3, 1000), (0, 25)])
    assert result == [(125, 1125)]


def test_fetch_intron_multiple():
    result = bam_cigar.fetch_intron(100, [(0, 10), (3, 500), (0, 10), (3, 300), (0, 10)])
    assert result == [(110, 610), (620, 920)]


def test_fetch_intron_with_insertion():
    """Insertion in intron-containing CIGAR is skipped."""
    result = bam_cigar.fetch_intron(100, [(0, 10), (1, 5), (3, 500), (0, 10)])
    assert result == [(110, 610)]


def test_fetch_intron_with_deletion():
    """Deletion advances position but produces no intron."""
    result = bam_cigar.fetch_intron(100, [(0, 10), (2, 5), (0, 10)])
    assert result == []


def test_fetch_intron_with_soft_clip():
    """Soft clip does not advance position in intron context."""
    result = bam_cigar.fetch_intron(100, [(4, 5), (0, 10), (3, 500), (0, 10)])
    assert result == [(110, 610)]


def test_fetch_intron_hard_clip():
    """Hard clip (op 5) falls through to continue."""
    result = bam_cigar.fetch_intron(100, [(5, 3), (0, 10), (3, 500), (0, 10)])
    assert result == [(110, 610)]


# --- fetch_deletion_range ---


def test_fetch_deletion_range_none():
    assert bam_cigar.fetch_deletion_range([(0, 50)]) == []


def test_fetch_deletion_range_single():
    result = bam_cigar.fetch_deletion_range([(0, 25), (2, 3), (0, 25)])
    assert result == [(25, 3)]


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
