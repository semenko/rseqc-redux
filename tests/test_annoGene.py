"""Tests for rseqc.annoGene."""

import importlib


def test_import():
    """Verify that rseqc.annoGene can be imported."""
    mod = importlib.import_module("rseqc.annoGene")
    assert mod is not None
