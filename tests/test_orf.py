"""Tests for rseqc.orf."""

import importlib

import pytest

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


@pytest.mark.xfail(
    reason="Known bug: sc parameter shadows module-level start_coden, UnboundLocalError when sc=None",
    raises=UnboundLocalError,
    strict=True,
)
def test_longest_orf_default_codons_bug():
    """Document the bug: calling without sc/tc raises UnboundLocalError."""
    orf.longest_orf("ATGAAATAACCC", "+")
