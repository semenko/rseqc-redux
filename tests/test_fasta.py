"""Tests for rseqc.fasta."""

import importlib

from rseqc.fasta import Fasta


def test_import():
    """Verify that rseqc.fasta can be imported."""
    mod = importlib.import_module("rseqc.fasta")
    assert mod is not None


def test_fasta_from_file(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGT\n>chr2\nTTTTAAAA\n")
    obj = Fasta(str(fa))
    assert "chr1" in obj.seqs
    assert "chr2" in obj.seqs
    assert obj.seqs["chr1"] == "ACGTACGT"
    assert obj.seqs["chr2"] == "TTTTAAAA"


def test_fasta_get_names(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n>seq2\nTGCA\n")
    obj = Fasta(str(fa))
    assert obj.getNames() == ["seq1", "seq2"]


def test_fasta_get_seq(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n")
    obj = Fasta(str(fa))
    assert obj.getSeq("seq1") == "ACGT"


def test_fasta_get_all_seqs(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n>seq2\nTGCA\n")
    obj = Fasta(str(fa))
    seqs = obj.getSeq()
    assert len(seqs) == 2


def test_fasta_add_seq():
    obj = Fasta()
    obj.addSeq("test", "acgt")
    assert obj.seqs["test"] == "ACGT"


def test_fasta_rev_comp(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nACGT\n")
    obj = Fasta(str(fa))
    result = obj.revComp("seq1")
    assert result == "ACGT"  # ACGT is its own reverse complement


def test_fasta_rev_comp_asymmetric(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">seq1\nAAAA\n")
    obj = Fasta(str(fa))
    result = obj.revComp("seq1")
    assert result == "TTTT"


def test_fasta_fetch_seq(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGTACGT\n")
    obj = Fasta(str(fa))
    result = obj.fetchSeq(chr="chr1", st=2, end=6)
    assert result == "GTAC"


def test_fasta_multiline(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGT\nTGCA\n")
    obj = Fasta(str(fa))
    assert obj.seqs["chr1"] == "ACGTTGCA"


def test_fasta_get_seq_len(tmp_path):
    fa = tmp_path / "test.fa"
    fa.write_text(">chr1\nACGTACGT\n")
    obj = Fasta(str(fa))
    lengths = obj.getSeqLen()
    assert lengths["chr1"] == 8


def test_fasta_empty_file(tmp_path):
    """Bug #2 regression: empty FASTA should not raise UnboundLocalError."""
    fa = tmp_path / "empty.fa"
    fa.write_text("")
    obj = Fasta(str(fa))
    assert obj.seqs == {}
    assert obj.IDs == []


def test_fasta_no_header(tmp_path):
    """Bug #2 regression: FASTA with no header lines should handle gracefully."""
    fa = tmp_path / "noheader.fa"
    fa.write_text("ACGTACGT\n")
    obj = Fasta(str(fa))
    assert obj.seqs == {}
    assert obj.IDs == []
