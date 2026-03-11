"""Tests for rseqc.fastq."""

import importlib


def test_import():
    """Verify that rseqc.fastq can be imported."""
    mod = importlib.import_module("rseqc.fastq")
    assert mod is not None
