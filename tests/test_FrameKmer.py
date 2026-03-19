"""Tests for rseqc.FrameKmer."""

import gzip
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


# --- seq_generator ---


def test_seq_generator_fasta(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 2
    assert results[0] == ["SEQ1", "ACGTACGT"]
    assert results[1] == ["SEQ2", "TTTTAAAA"]


def test_seq_generator_fastq_header(tmp_path):
    """seq_generator also supports @ headers (FASTQ-style)."""
    fa = tmp_path / "test.fa"
    fa.write_text("@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nIIII\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 2
    assert results[0][0] == "READ1"


def test_seq_generator_skips_comments(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text("# comment\n>seq1\nACGT\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 1
    assert results[0][0] == "SEQ1"


def test_seq_generator_filters_non_dna(tmp_path):
    """Records with non-ACGTN characters in the sequence are skipped."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTXYZ\n>seq2\nTTTT\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    # seq1 has non-DNA characters, should be filtered
    assert len(results) == 1
    assert results[0] == ["SEQ2", "TTTT"]


def test_seq_generator_multiline(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\nTGCA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert results[0][1] == "ACGTTGCA"


def test_seq_generator_lowercase(tmp_path):
    """Lowercase sequences are uppercased."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nacgtacgt\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert results[0][1] == "ACGTACGT"


def test_seq_generator_gzip(tmp_path):
    """seq_generator handles gzip-compressed FASTA."""
    fa = tmp_path / "test.fa.gz"
    with gzip.open(str(fa), "wt") as gz:
        gz.write(">seq1\nACGT\n>seq2\nTGCA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 2
    assert results[0][1] == "ACGT"
    assert results[1][1] == "TGCA"


def test_seq_generator_empty_file(tmp_path):
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert results == []


def test_seq_generator_n_bases(tmp_path):
    """Sequences with N bases are valid DNA and should be kept."""
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGNTGCA\n")
    results = list(FrameKmer.seq_generator(str(fa)))
    assert len(results) == 1
    assert results[0][1] == "ACGNTGCA"


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
