"""Tests for rseqc.changePoint."""

import importlib

from rseqc import changePoint


def test_import():
    """Verify that rseqc.changePoint can be imported."""
    mod = importlib.import_module("rseqc.changePoint")
    assert mod is not None


# --- S_diff ---


def test_s_diff_constant():
    # constant list: S_diff should be 0
    result = changePoint.S_diff([5, 5, 5, 5])
    index_list, diff = result
    assert diff == 0.0


def test_s_diff_step_change():
    # clear step change: [1,1,1,10,10,10]
    result = changePoint.S_diff([1, 1, 1, 10, 10, 10])
    index_list, diff = result
    assert diff == 13.5
    # S_cum = [0, -4.5, -9.0, -13.5, -9.0, -4.5, 0.0]
    # max is at index 0, which is what nlargest returns
    assert index_list == [0]


def test_s_diff_increasing():
    result = changePoint.S_diff([1, 2, 3, 4, 5])
    index_list, diff = result
    assert diff > 0


def test_s_diff_two_elements():
    result = changePoint.S_diff([0, 10])
    index_list, diff = result
    assert diff > 0
