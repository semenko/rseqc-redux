"""Tests for rseqc.cigar."""

import importlib

from rseqc import cigar


def test_import():
    """Verify that rseqc.cigar can be imported."""
    mod = importlib.import_module("rseqc.cigar")
    assert mod is not None


# --- list2str ---


def test_list2str_simple_match():
    assert cigar.list2str([(0, 50)]) == "50M"


def test_list2str_complex():
    # pysam returns tuples of (op_code, length)
    assert cigar.list2str([(4, 5), (0, 45)]) == "5S45M"


def test_list2str_all_ops():
    ops = [(0, 10), (1, 2), (2, 3), (3, 1000), (4, 5), (5, 6), (6, 7), (7, 8), (8, 9)]
    result = cigar.list2str(ops)
    assert result == "10M2I3D1000N5S6H7P8=9X"


def test_list2str_empty():
    assert cigar.list2str([]) == ""


# --- fetch_head_clip ---


def test_fetch_head_clip_with_clip():
    # 5S45M: head clip of 5 bases before position 100
    result = cigar.fetch_head_clip("chr1", 100, "5S45M")
    assert result == [["chr1", 95, 100]]


def test_fetch_head_clip_no_clip():
    result = cigar.fetch_head_clip("chr1", 100, "50M")
    assert result == []


# --- fetch_tail_clip ---


def test_fetch_tail_clip_with_clip():
    # 45M5S: tail clip of 5 bases after the match
    result = cigar.fetch_tail_clip("chr1", 100, "45M5S")
    # ref_part matches M, I, S, N, D, =, X
    # ref_part for "45M5S" => [45, 5], sum=50
    # h_len=0, chrom_end = 100 + (50 - 0) = 150, chrom_st = 150 - 5 = 145
    assert result == [["chr1", 145, 150]]


def test_fetch_tail_clip_no_clip():
    result = cigar.fetch_tail_clip("chr1", 100, "50M")
    assert result == []


def test_fetch_tail_clip_both_clips():
    # 5S40M5S
    result = cigar.fetch_tail_clip("chr1", 100, "5S40M5S")
    # ref_part for "5S40M5S" => [5, 40, 5], sum=50
    # h_len=5, chrom_end = 100 + (50 - 5) = 145, chrom_st = 145 - 5 = 140
    assert result == [["chr1", 140, 145]]


# --- fetch_exon ---


def test_fetch_exon_simple():
    # 50M at position 100
    result = cigar.fetch_exon("chr1", 100, "50M")
    assert result == [["chr1", 100, 150]]


def test_fetch_exon_spliced():
    # 25M1000N25M: two exons separated by an intron
    result = cigar.fetch_exon("chr1", 100, "25M1000N25M")
    assert len(result) == 2
    assert result[0] == ["chr1", 100, 125]
    assert result[1] == ["chr1", 1125, 1150]


def test_fetch_exon_with_head_clip():
    # 5S45M at position 100
    result = cigar.fetch_exon("chr1", 100, "5S45M")
    assert result == [["chr1", 100, 145]]


# --- fetch_intron ---


def test_fetch_intron_no_intron():
    result = cigar.fetch_intron("chr1", 100, "50M")
    assert result == []


def test_fetch_intron_single():
    # 25M1000N25M
    result = cigar.fetch_intron("chr1", 100, "25M1000N25M")
    assert result == [["chr1", 125, 1125]]


def test_fetch_intron_multiple():
    # 10M500N10M300N10M
    result = cigar.fetch_intron("chr1", 100, "10M500N10M300N10M")
    assert len(result) == 2
    assert result[0] == ["chr1", 110, 610]
    assert result[1] == ["chr1", 620, 920]


# --- fetch_insertion ---


def test_fetch_insertion_none():
    result = cigar.fetch_insertion("chr1", 100, "50M")
    assert result == []


def test_fetch_insertion_single():
    # 25M5I20M — insertion after 25 matched bases
    result = cigar.fetch_insertion("chr1", 100, "25M5I20M")
    assert len(result) == 1
    assert result[0][0] == "chr1"
    assert result[0][1] == 125  # position after 25M
    assert result[0][2] == 130  # insertion length 5


# --- fetch_deletion ---


def test_fetch_deletion_none():
    result = cigar.fetch_deletion("chr1", 100, "50M")
    assert result == []


def test_fetch_deletion_single():
    # 25M3D25M — deletion of 3 ref bases after 25 matched
    result = cigar.fetch_deletion("chr1", 100, "25M3D25M")
    assert len(result) == 1
    assert result[0][0] == "chr1"
    assert result[0][1] == 125
    assert result[0][2] == 128
