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
