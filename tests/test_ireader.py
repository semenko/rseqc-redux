"""Tests for rseqc.ireader."""

import importlib


def test_import():
    """Verify that rseqc.ireader can be imported."""
    mod = importlib.import_module("rseqc.ireader")
    assert mod is not None
