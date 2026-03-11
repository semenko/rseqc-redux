"""Tests for rseqc.dotProduct."""

import importlib


def test_import():
    """Verify that rseqc.dotProduct can be imported."""
    mod = importlib.import_module("rseqc.dotProduct")
    assert mod is not None
