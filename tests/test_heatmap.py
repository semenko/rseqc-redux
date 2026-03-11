"""Tests for rseqc.heatmap."""

import importlib


def test_import():
    """Verify that rseqc.heatmap can be imported."""
    mod = importlib.import_module("rseqc.heatmap")
    assert mod is not None
