"""Tests for rseqc.fasta."""

import importlib


def test_import():
    """Verify that rseqc.fasta can be imported."""
    mod = importlib.import_module("rseqc.fasta")
    assert mod is not None
