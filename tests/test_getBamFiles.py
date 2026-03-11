"""Tests for rseqc.getBamFiles."""

import importlib


def test_import():
    """Verify that rseqc.getBamFiles can be imported."""
    mod = importlib.import_module("rseqc.getBamFiles")
    assert mod is not None
