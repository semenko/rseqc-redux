"""Tests for rseqc.twoList."""

import importlib

import numpy as np
import pytest

from rseqc import twoList


def test_import():
    """Verify that rseqc.twoList can be imported."""
    mod = importlib.import_module("rseqc.twoList")
    assert mod is not None


# --- check_list ---


def test_check_list_same_size():
    v1 = np.array([1, 2, 3])
    v2 = np.array([4, 5, 6])
    twoList.check_list(v1, v2)  # should not raise


def test_check_list_different_size():
    v1 = np.array([1, 2])
    v2 = np.array([1, 2, 3])
    with pytest.raises(ValueError):
        twoList.check_list(v1, v2)


# --- Add ---


def test_add():
    v1 = np.array([1, 2, 3])
    v2 = np.array([4, 5, 6])
    result = twoList.Add(v1, v2)
    np.testing.assert_array_equal(result, [5, 7, 9])


# --- Subtract ---


def test_subtract():
    v1 = np.array([10, 20, 30])
    v2 = np.array([1, 2, 3])
    result = twoList.Subtract(v1, v2)
    np.testing.assert_array_equal(result, [9, 18, 27])


# --- Product ---


def test_product():
    v1 = np.array([2, 3, 4])
    v2 = np.array([5, 6, 7])
    result = twoList.Product(v1, v2)
    np.testing.assert_array_equal(result, [10, 18, 28])


# --- Average ---


def test_average():
    v1 = np.array([10.0, 20.0])
    v2 = np.array([20.0, 40.0])
    result = twoList.Average(v1, v2)
    np.testing.assert_array_equal(result, [15.0, 30.0])


# --- geometricMean ---


def test_geometric_mean():
    v1 = np.array([4.0, 9.0])
    v2 = np.array([16.0, 1.0])
    result = twoList.geometricMean(v1, v2)
    np.testing.assert_array_almost_equal(result, [8.0, 3.0])


# --- euclidean_distance ---


def test_euclidean_distance_same():
    v1 = np.array([1.0, 2.0, 3.0])
    v2 = np.array([1.0, 2.0, 3.0])
    assert twoList.euclidean_distance(v1, v2) == 0.0


def test_euclidean_distance_basic():
    v1 = np.array([0.0, 0.0])
    v2 = np.array([3.0, 4.0])
    # sqrt((9+16)/2) = sqrt(12.5)
    expected = (25 / 2) ** 0.5
    assert abs(twoList.euclidean_distance(v1, v2) - expected) < 1e-10


# --- Division ---


def test_division():
    v1 = np.array([9.0, 19.0, 29.0])
    v2 = np.array([4.0, 9.0, 14.0])
    result = twoList.Division(v1, v2)
    # (v1+1)/(v2+1) = [10/5, 20/10, 30/15] = [2.0, 2.0, 2.0]
    np.testing.assert_array_almost_equal(result, [2.0, 2.0, 2.0])


def test_division_zeros():
    v1 = np.array([0.0, 0.0])
    v2 = np.array([0.0, 0.0])
    result = twoList.Division(v1, v2)
    # (0+1)/(0+1) = 1.0
    np.testing.assert_array_almost_equal(result, [1.0, 1.0])


def test_division_length_mismatch():
    v1 = np.array([1.0, 2.0])
    v2 = np.array([1.0])
    with pytest.raises(ValueError):
        twoList.Division(v1, v2)


# --- Max ---


def test_max():
    v1 = np.array([1, 5, 3])
    v2 = np.array([4, 2, 6])
    result = list(twoList.Max(v1, v2))
    assert result == [4, 5, 6]


def test_max_equal():
    v1 = np.array([3, 3])
    v2 = np.array([3, 3])
    result = list(twoList.Max(v1, v2))
    assert result == [3, 3]


# --- Min ---


def test_min():
    v1 = np.array([1, 5, 3])
    v2 = np.array([4, 2, 6])
    result = list(twoList.Min(v1, v2))
    assert result == [1, 2, 3]


def test_min_equal():
    v1 = np.array([3, 3])
    v2 = np.array([3, 3])
    result = list(twoList.Min(v1, v2))
    assert result == [3, 3]
