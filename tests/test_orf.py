"""Tests for rseqc.orf."""

import importlib

from rseqc import orf


def test_import():
    """Verify that rseqc.orf can be imported."""
    mod = importlib.import_module("rseqc.orf")
    assert mod is not None


# --- _reverse_comp ---


def test_reverse_comp_basic():
    assert orf._reverse_comp("ACGT") == "ACGT"  # palindrome


def test_reverse_comp_simple():
    assert orf._reverse_comp("AAAA") == "TTTT"


def test_reverse_comp_with_n():
    assert orf._reverse_comp("ACNGT") == "ACNGT"


def test_reverse_comp_single():
    assert orf._reverse_comp("A") == "T"


# --- longest_orf ---
# NOTE: When sc=None, Python >=3.12 raises UnboundLocalError because the
# parameter `sc` shadows the module-level `start_coden` via `for sc in start_coden`.
# We must pass sc/tc explicitly to work around this bug.


def test_longest_orf_simple():
    seq = "AAAATGAAATAACCC"  # ATG at 3, TAA at 9, length=6, in-frame
    result = orf.longest_orf(seq, "+", sc="ATG", tc="TAG,TAA,TGA")
    assert "ATG" in result
    assert len(result) % 3 == 0


def test_longest_orf_no_start():
    seq = "CCCCCCTAACCC"
    result = orf.longest_orf(seq, "+", sc="ATG", tc="TAG,TAA,TGA")
    assert result == ""


def test_longest_orf_no_stop():
    seq = "ATGAAACCC"
    result = orf.longest_orf(seq, "+", sc="ATG", tc="TAG,TAA,TGA")
    assert result == ""


def test_longest_orf_minus_strand():
    # reverse complement of ATGAAATAA is TTATTTCAT
    seq = "TTATTTCAT"
    result = orf.longest_orf(seq, "-", sc="ATG", tc="TAG,TAA,TGA")
    assert "ATG" in result


def test_longest_orf_multiple():
    # two ORFs, should return the longest
    seq = "ATGAAATAACCCATGAAAAAATAACCC"
    result = orf.longest_orf(seq, "+", sc="ATG", tc="TAG,TAA,TGA")
    assert len(result) >= 6


def test_longest_orf_default_codons():
    """Calling without sc/tc should use module-level start/stop codons."""
    result = orf.longest_orf("ATGAAATAACCC", "+")
    assert "ATG" in result
    assert len(result) % 3 == 0


def test_longest_orf_custom_start_codon():
    """Custom start codon CTG."""
    seq = "CTGAAATAACCC"
    result = orf.longest_orf(seq, "+", sc="CTG", tc="TAA")
    assert result.startswith("CTG")


def test_longest_orf_custom_stop_codon():
    """Custom stop codon."""
    seq = "ATGAAAGGGCCC"
    result = orf.longest_orf(seq, "+", sc="ATG", tc="GGG")
    assert "ATG" in result


def test_longest_orf_empty_seq():
    result = orf.longest_orf("", "+", sc="ATG", tc="TAA")
    assert result == ""


def test_longest_orf_all_stop_codons():
    """Test all three stop codons work."""
    for stop in ["TAA", "TAG", "TGA"]:
        seq = f"ATG{'A' * 6}{stop}CCC"
        result = orf.longest_orf(seq, "+", sc="ATG", tc=stop)
        assert len(result) > 0


# --- longest_orf_bed ---


def test_longest_orf_bed_plus_strand():
    """Test longest_orf_bed with a simple single-exon gene on + strand."""
    seq = "ATGAAATAACCC"
    # BED12 line: chrom, start, end, name, score, strand, thickStart, thickEnd, rgb, blockCount, blockSizes, blockStarts
    bedline = "chr1\t100\t112\tgene1\t0\t+\t100\t112\t0\t1\t12,\t0,"
    result = orf.longest_orf_bed(seq, bedline, sc="ATG", tc="TAA")
    assert result is not None
    fields = result.split("\t")
    assert fields[0] == "chr1"
    # CDS start and end should be within gene bounds
    assert int(fields[6]) >= 100
    assert int(fields[7]) <= 112


def test_longest_orf_bed_minus_strand():
    """Test longest_orf_bed on minus strand."""
    # reverse complement of ATGAAATAA is TTATTTCAT
    seq = "TTATTTCAT"
    bedline = "chr1\t200\t209\tgene2\t0\t-\t200\t209\t0\t1\t9,\t0,"
    result = orf.longest_orf_bed(seq, bedline, sc="ATG", tc="TAA")
    assert result is not None
    fields = result.split("\t")
    assert fields[5] == "-"


def test_longest_orf_bed_no_orf():
    """No ORF found returns None."""
    seq = "CCCCCCCCC"
    bedline = "chr1\t100\t109\tgene1\t0\t+\t100\t109\t0\t1\t9,\t0,"
    result = orf.longest_orf_bed(seq, bedline, sc="ATG", tc="TAA")
    assert result is None


# --- _reverse_comp edge cases ---


def test_reverse_comp_empty():
    assert orf._reverse_comp("") == ""


def test_reverse_comp_all_bases():
    assert orf._reverse_comp("ACGTX") == "XACGT"
