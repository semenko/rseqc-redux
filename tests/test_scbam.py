"""Tests for rseqc.scbam."""

import importlib


def test_import():
    """Verify that rseqc.scbam can be imported."""
    mod = importlib.import_module("rseqc.scbam")
    assert mod is not None
