"""Tests for rseqc.annoGene."""

import importlib
from pathlib import Path

from rseqc import annoGene


def test_import():
    """Verify that rseqc.annoGene can be imported."""
    mod = importlib.import_module("rseqc.annoGene")
    assert mod is not None


FIXTURES_DIR = Path(__file__).parent / "fixtures"


# --- getCDSExonFromFile ---


def test_get_cds_exon_from_file():
    result = annoGene.getCDSExonFromFile(str(FIXTURES_DIR / "mini.bed"))
    assert len(result) >= 5
    assert any(":+" in r[0] for r in result)


# --- getUTRExonFromFile ---


def test_get_utr_exon_from_file():
    result = annoGene.getUTRExonFromFile(str(FIXTURES_DIR / "mini.bed"), utr=35)
    assert len(result) >= 2


def test_get_utr_exon_5prime_only():
    result = annoGene.getUTRExonFromFile(str(FIXTURES_DIR / "mini.bed"), utr=5)
    assert len(result) >= 1


# --- getExonFromFile ---


def test_get_exon_from_file():
    result = annoGene.getExonFromFile(str(FIXTURES_DIR / "mini.bed"))
    assert len(result) == 6
    assert ":" in result[0][0]


# --- getUTRExonFromLine (Bug #5 regression) ---


def test_get_utr_exon_from_line():
    """Bug #5 regression: chromm typo was NameError."""
    bedline = "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500"
    result = annoGene.getUTRExonFromLine(bedline, utr=35)
    assert result is not None
    assert len(result) >= 1


# --- getCDSExonFromLine (Bug #5 regression) ---


def test_get_cds_exon_from_line():
    """Bug #5 regression: chromm typo was NameError."""
    bedline = "chr1\t1000\t5000\tgene1\t0\t+\t1200\t4800\t0,0,0\t3\t500,600,500\t0,1500,3500"
    result = annoGene.getCDSExonFromLine(bedline)
    assert result is not None
    assert len(result) >= 1


# --- getExonFromFile2 (Bug #6 regression) ---


def test_get_exon_from_file2():
    """Bug #6 regression: txstart/txEnd undefined, append() wrong args."""
    result = annoGene.getExonFromFile2(str(FIXTURES_DIR / "mini.bed"))
    assert len(result) >= 1
