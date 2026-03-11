"""Tests for rseqc.twoList."""

import importlib


def test_import():
    """Verify that rseqc.twoList can be imported."""
    mod = importlib.import_module("rseqc.twoList")
    assert mod is not None
