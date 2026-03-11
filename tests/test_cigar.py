"""Tests for rseqc.cigar."""

import importlib


def test_import():
    """Verify that rseqc.cigar can be imported."""
    mod = importlib.import_module("rseqc.cigar")
    assert mod is not None
