"""Tests for rseqc.changePoint."""

import importlib


def test_import():
    """Verify that rseqc.changePoint can be imported."""
    mod = importlib.import_module("rseqc.changePoint")
    assert mod is not None
