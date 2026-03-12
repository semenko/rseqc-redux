"""Tests for rseqc.wiggle."""

import importlib


def test_import():
    """Verify that rseqc.wiggle can be imported."""
    mod = importlib.import_module("rseqc.wiggle")
    assert mod is not None
