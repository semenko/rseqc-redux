"""Tests for rseqc.FrameKmer."""

import importlib

from rseqc import FrameKmer


def test_import():
    """Verify that rseqc.FrameKmer can be imported."""
    mod = importlib.import_module("rseqc.FrameKmer")
    assert mod is not None


# --- word_generator ---


def test_word_generator_basic():
    words = list(FrameKmer.word_generator("ACGTACGT", word_size=3, step_size=1))
    assert words[0] == "ACG"
    assert words[1] == "CGT"
    assert len(words) == 6  # 8-3+1


def test_word_generator_step3():
    words = list(FrameKmer.word_generator("ACGTACGTAC", word_size=3, step_size=3))
    assert words == ["ACG", "TAC", "GTA"]


def test_word_generator_frame1():
    words = list(FrameKmer.word_generator("ACGTACGT", word_size=3, step_size=3, frame=1))
    assert words[0] == "CGT"


def test_word_generator_short_seq():
    words = list(FrameKmer.word_generator("AC", word_size=3, step_size=1))
    assert words == []


def test_word_generator_exact():
    words = list(FrameKmer.word_generator("ACG", word_size=3, step_size=1))
    assert words == ["ACG"]


# --- all_possible_kmer ---


def test_all_possible_kmer_1():
    kmers = list(FrameKmer.all_possible_kmer(1))
    assert set(kmers) == {"A", "C", "G", "T", "N"}


def test_all_possible_kmer_2():
    kmers = list(FrameKmer.all_possible_kmer(2))
    assert len(kmers) == 25  # 5^2


def test_all_possible_kmer_3():
    kmers = list(FrameKmer.all_possible_kmer(3))
    assert len(kmers) == 125  # 5^3


# --- kmer_ratio ---


def test_kmer_ratio_basic():
    seq = "ATGATGATGATGATG"
    # build coding/noncoding frequency tables
    coding = {"ATG": 10, "TGA": 5, "GAT": 8}
    noncoding = {"ATG": 5, "TGA": 10, "GAT": 8}
    result = FrameKmer.kmer_ratio(seq, word_size=3, step_size=3, coding=coding, noncoding=noncoding)
    assert isinstance(result, float)


def test_kmer_ratio_short_seq():
    result = FrameKmer.kmer_ratio("AT", word_size=3, step_size=3, coding={}, noncoding={})
    assert result == 0
