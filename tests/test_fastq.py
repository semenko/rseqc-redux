"""Tests for rseqc.fastq."""

import importlib

from rseqc import fastq


def test_import():
    """Verify that rseqc.fastq can be imported."""
    mod = importlib.import_module("rseqc.fastq")
    assert mod is not None


# --- fasta_iter ---


def test_fasta_iter(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGTACGT\n>seq2\nTTTTAAAA\n")
    seqs = list(fastq.fasta_iter(str(fa)))
    assert seqs == ["ACGTACGT", "TTTTAAAA"]


def test_fasta_iter_empty_lines(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\n\nACGT\n\n")
    seqs = list(fastq.fasta_iter(str(fa)))
    assert seqs == ["ACGT"]


# --- fastq_iter ---


def test_fastq_iter_seq(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nIIIIIIII\n")
    seqs = list(fastq.fastq_iter(str(fq), mode="seq"))
    assert seqs == ["ACGTACGT", "TTTTAAAA"]


def test_fastq_iter_qual(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTTTTAAAA\n+\nJJJJJJJJ\n")
    quals = list(fastq.fastq_iter(str(fq), mode="qual"))
    assert quals == ["IIIIIIII", "JJJJJJJJ"]


# --- seq2countMat ---


def test_seq2countmat_basic(tmp_path):
    fq = tmp_path / "test.fq"
    fq.write_text("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\nIIII\n")
    s_obj = fastq.fastq_iter(str(fq), mode="seq")
    mat = fastq.seq2countMat(s_obj, limit=None)
    assert mat.shape[0] == 4  # 4 positions
    assert mat.loc[0, "A"] == 2
    assert mat.loc[1, "C"] == 2
