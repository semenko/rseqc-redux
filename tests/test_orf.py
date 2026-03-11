"""Tests for rseqc.orf."""

import importlib


def test_import():
    """Verify that rseqc.orf can be imported."""
    mod = importlib.import_module("rseqc.orf")
    assert mod is not None
