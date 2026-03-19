"""Tests for rseqc.mystat."""

import importlib

from rseqc import mystat


def test_import():
    """Verify that rseqc.mystat can be imported."""
    mod = importlib.import_module("rseqc.mystat")
    assert mod is not None


# --- percentile_list ---


def test_percentile_list_short():
    # lists shorter than 100 are returned as-is
    data = list(range(50))
    assert mystat.percentile_list(data) == data


def test_percentile_list_empty():
    assert mystat.percentile_list([]) is None


def test_percentile_list_100():
    data = list(range(200))
    result = mystat.percentile_list(data)
    assert len(result) == 100
    # last element should be close to max
    assert result[-1] == 199 or result[-1] == 198


def test_percentile_list_exactly_100_elements():
    data = list(range(100))
    result = mystat.percentile_list(data)
    assert len(result) == 100
    assert result[0] == 1  # 1st percentile
    assert result[-1] == 99  # 100th percentile


def test_percentile_list_101_elements():
    data = list(range(101))
    result = mystat.percentile_list(data)
    assert len(result) == 100
    assert result[-1] == 100


def test_percentile_list_floats():
    data = sorted([0.5 * i for i in range(200)])
    result = mystat.percentile_list(data)
    assert len(result) == 100
    assert all(isinstance(x, int) for x in result)


def test_percentile_list_large():
    data = list(range(10000))
    result = mystat.percentile_list(data)
    assert len(result) == 100
    assert result[0] == 100  # 1st percentile of 0..9999
    assert result[-1] == 9999
