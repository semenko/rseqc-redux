"""Tests for rseqc.bam_cigar."""

import importlib


def test_import():
    """Verify that rseqc.bam_cigar can be imported."""
    mod = importlib.import_module("rseqc.bam_cigar")
    assert mod is not None
