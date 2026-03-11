"""Tests for rseqc.quantile."""

import importlib

from rseqc.quantile import quantile


def test_import():
    """Verify that rseqc.quantile can be imported."""
    mod = importlib.import_module("rseqc.quantile")
    assert mod is not None


def test_quantile_median_odd():
    x = [1, 2, 3, 4, 5]
    result = quantile(x, 0.5)
    assert result == 3


def test_quantile_median_even():
    x = [1, 2, 3, 4]
    result = quantile(x, 0.5)
    assert result == 2.5


def test_quantile_min():
    x = [10, 20, 30]
    result = quantile(x, 0.0)
    assert result == 10


def test_quantile_max():
    x = [10, 20, 30]
    result = quantile(x, 1.0)
    assert result == 30


def test_quantile_quarter():
    x = [11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75]
    result = quantile(x, 0.25)
    assert isinstance(result, (int, float))
    assert 17.0 <= result <= 22.0


def test_quantile_invalid_type():
    x = [1, 2, 3]
    assert quantile(x, 0.5, qtype=0) is None
    assert quantile(x, 0.5, qtype=10) is None


def test_quantile_unsorted_input():
    x = [5, 1, 3, 2, 4]
    result = quantile(x, 0.5)
    assert result == 3


def test_quantile_all_types():
    x = [11.4, 17.3, 21.3, 25.9, 40.1, 50.5, 60.0, 70.0, 75]
    for qtype in range(1, 10):
        result = quantile(x, 0.35, qtype)
        assert isinstance(result, (int, float))
        assert 11.4 <= result <= 75
