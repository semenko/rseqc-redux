"""Tests for rseqc.dotProduct."""

import importlib

import pytest

from rseqc import dotProduct


def test_import():
    """Verify that rseqc.dotProduct can be imported."""
    mod = importlib.import_module("rseqc.dotProduct")
    assert mod is not None


# --- d0 (loop-based dot product) ---


def test_d0_basic():
    assert dotProduct.d0([1, 2, 3], [4, 5, 6]) == 32


def test_d0_zeros():
    assert dotProduct.d0([0, 0, 0], [1, 2, 3]) == 0


def test_d0_single():
    assert dotProduct.d0([5], [3]) == 15


def test_d0_negative():
    assert dotProduct.d0([-1, 2], [3, -4]) == -11


def test_d0_length_mismatch():
    with pytest.raises(ValueError):
        dotProduct.d0([1, 2], [1, 2, 3])


# --- d1 (imap-based) ---


def test_d1_basic():
    assert dotProduct.d1([1, 2, 3], [4, 5, 6]) == 32


def test_d1_zeros():
    assert dotProduct.d1([0, 0, 0], [1, 2, 3]) == 0


# --- d2 (map-based) ---


def test_d2_basic():
    assert dotProduct.d2([1, 2, 3], [4, 5, 6]) == 32


# --- d3 (starmap-based) — note: returns iterator, not sum ---


def test_d3_returns_iterable():
    result = dotProduct.d3([1, 2, 3], [4, 5, 6])
    # d3 returns a starmap iterator, not a sum
    assert sum(result) == 32


# --- check ---


def test_check_same_length():
    # should not raise
    dotProduct.check([1, 2, 3], [4, 5, 6])


def test_check_different_length():
    with pytest.raises(ValueError):
        dotProduct.check([1, 2], [1, 2, 3])
