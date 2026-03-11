"""Tests for rseqc.BED."""

import importlib


def test_import():
    """Verify that rseqc.BED can be imported."""
    mod = importlib.import_module("rseqc.BED")
    assert mod is not None
