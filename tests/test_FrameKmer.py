"""Tests for rseqc.FrameKmer."""

import importlib


def test_import():
    """Verify that rseqc.FrameKmer can be imported."""
    mod = importlib.import_module("rseqc.FrameKmer")
    assert mod is not None
