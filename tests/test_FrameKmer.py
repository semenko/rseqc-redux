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


# --- seq_generator ---


def test_seq_generator_fasta(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 2
    # seq_generator uppercases the entire line including headers
    assert results[0] == ["SEQ1", "ACGTACGT"]
    assert results[1] == ["SEQ2", "TTTTAAAA"]


def test_seq_generator_fastq_header(tmp_path):
    """seq_generator also supports @ headers (FASTQ-style)."""
    fa = tmp_path / "test.fa"
    fa.write_text("@read1\nACGT\n@read2\nTGCA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 2
    assert results[0][0] == "READ1"  # uppercased


def test_seq_generator_skips_comments(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text("# comment\n>seq1\nACGT\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 1
    assert results[0][0] == "SEQ1"  # uppercased


def test_seq_generator_filters_non_dna(tmp_path):
    """Lines with non-ACGTN characters are not appended to sequence."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n+\nIIII\n>seq2\nTTTT\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    # '+' line starts with special char but not '>' or '@' — not a header
    # 'IIII' doesn't match DNA_pat, so not appended
    assert any(r[0] == "SEQ1" for r in results)


def test_seq_generator_multiline(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\nTGCA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert results[0][1] == "ACGTTGCA"


# --- kmer_freq_file ---


def test_kmer_freq_file_basic(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n")
    result = FrameKmer.kmer_freq_file(str(fa), word_size=2, step_size=1)
    assert isinstance(result, dict)
    assert result["AC"] == 2
    assert result["CG"] == 2
    assert result["GT"] == 2
    assert result["TA"] == 1


def test_kmer_freq_file_multiple_seqs(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n>seq2\nACGT\n")
    result = FrameKmer.kmer_freq_file(str(fa), word_size=2, step_size=1)
    assert result["AC"] == 2  # 1 from each seq


def test_kmer_freq_file_min_count(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n")
    result = FrameKmer.kmer_freq_file(str(fa), word_size=2, step_size=1, min_count=2)
    # Only kmers with count >= 2 should be included
    assert "TA" not in result  # TA only appears once
    assert "AC" in result  # AC appears twice


def test_kmer_freq_file_with_frame(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n")
    result = FrameKmer.kmer_freq_file(str(fa), word_size=3, step_size=3, frame=0)
    assert isinstance(result, dict)
    assert "ACG" in result


# --- kmer_ratio edge cases ---


def test_kmer_ratio_coding_only():
    """coding > 0, noncoding == 0 → adds 1."""
    coding = {"ATG": 10}
    noncoding = {"ATG": 0}
    result = FrameKmer.kmer_ratio("ATGATG", word_size=3, step_size=3, coding=coding, noncoding=noncoding)
    assert result == 1.0


def test_kmer_ratio_noncoding_only():
    """coding == 0, noncoding > 0 → subtracts 1."""
    coding = {"ATG": 0}
    noncoding = {"ATG": 10}
    result = FrameKmer.kmer_ratio("ATGATG", word_size=3, step_size=3, coding=coding, noncoding=noncoding)
    assert result == -1.0


def test_kmer_ratio_both_zero():
    """coding == 0, noncoding == 0 → skipped."""
    coding = {"ATG": 0, "TGA": 5}
    noncoding = {"ATG": 0, "TGA": 5}
    # Only TGA contributes: log(5/5) = 0
    result = FrameKmer.kmer_ratio("ATGTGA", word_size=3, step_size=3, coding=coding, noncoding=noncoding)
    assert result == 0.0


def test_kmer_ratio_missing_kmer():
    """kmer not in coding or noncoding → skipped, causing ZeroDivisionError."""
    import pytest

    coding = {"AAA": 5}
    noncoding = {"AAA": 5}
    # ATG not in either dict → all words skipped → frame0_count stays 0 → ZeroDivisionError
    with pytest.raises(ZeroDivisionError):
        FrameKmer.kmer_ratio("ATGATG", word_size=3, step_size=3, coding=coding, noncoding=noncoding)


def test_kmer_ratio_equal_frequencies():
    """Equal coding and noncoding → log(1) = 0."""
    coding = {"ATG": 5, "TGA": 5}
    noncoding = {"ATG": 5, "TGA": 5}
    result = FrameKmer.kmer_ratio("ATGTGA", word_size=3, step_size=3, coding=coding, noncoding=noncoding)
    assert abs(result) < 1e-10
