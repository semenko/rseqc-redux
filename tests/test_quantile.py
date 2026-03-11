"""Tests for rseqc.quantile."""

import importlib


def test_import():
    """Verify that rseqc.quantile can be imported."""
    mod = importlib.import_module("rseqc.quantile")
    assert mod is not None
