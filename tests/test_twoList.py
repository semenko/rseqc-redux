"""Tests for array manipulation functions in scripts.overlay_bigwig."""

import numpy as np
import pytest

from scripts.overlay_bigwig import _ACTIONS, _apply

# --- _apply (size validation) ---


def test_apply_same_size():
    v1 = np.array([1, 2, 3])
    v2 = np.array([4, 5, 6])
    _apply(lambda a, b: a + b, v1, v2)  # should not raise


def test_apply_different_size():
    v1 = np.array([1, 2])
    v2 = np.array([1, 2, 3])
    with pytest.raises(ValueError):
        _apply(lambda a, b: a + b, v1, v2)


# --- Add ---


def test_add():
    v1 = np.array([1, 2, 3])
    v2 = np.array([4, 5, 6])
    result = _apply(_ACTIONS["Add"], v1, v2)
    np.testing.assert_array_equal(result, [5, 7, 9])


# --- Subtract ---


def test_subtract():
    v1 = np.array([10, 20, 30])
    v2 = np.array([1, 2, 3])
    result = _apply(_ACTIONS["Subtract"], v1, v2)
    np.testing.assert_array_equal(result, [9, 18, 27])


# --- Product ---


def test_product():
    v1 = np.array([2, 3, 4])
    v2 = np.array([5, 6, 7])
    result = _apply(_ACTIONS["Product"], v1, v2)
    np.testing.assert_array_equal(result, [10, 18, 28])


# --- Average ---


def test_average():
    v1 = np.array([10.0, 20.0])
    v2 = np.array([20.0, 40.0])
    result = _apply(_ACTIONS["Average"], v1, v2)
    np.testing.assert_array_equal(result, [15.0, 30.0])


# --- geometricMean ---


def test_geometric_mean():
    v1 = np.array([4.0, 9.0])
    v2 = np.array([16.0, 1.0])
    result = _apply(_ACTIONS["geometricMean"], v1, v2)
    np.testing.assert_array_almost_equal(result, [8.0, 3.0])


# --- Division ---


def test_division():
    v1 = np.array([9.0, 19.0, 29.0])
    v2 = np.array([4.0, 9.0, 14.0])
    result = _apply(_ACTIONS["Division"], v1, v2)
    # (v1+1)/(v2+1) = [10/5, 20/10, 30/15] = [2.0, 2.0, 2.0]
    np.testing.assert_array_almost_equal(result, [2.0, 2.0, 2.0])


def test_division_zeros():
    v1 = np.array([0.0, 0.0])
    v2 = np.array([0.0, 0.0])
    result = _apply(_ACTIONS["Division"], v1, v2)
    # (0+1)/(0+1) = 1.0
    np.testing.assert_array_almost_equal(result, [1.0, 1.0])


def test_division_length_mismatch():
    v1 = np.array([1.0, 2.0])
    v2 = np.array([1.0])
    with pytest.raises(ValueError):
        _apply(_ACTIONS["Division"], v1, v2)


# --- Max ---


def test_max():
    v1 = np.array([1, 5, 3])
    v2 = np.array([4, 2, 6])
    result = _apply(_ACTIONS["Max"], v1, v2)
    np.testing.assert_array_equal(result, [4, 5, 6])


def test_max_equal():
    v1 = np.array([3, 3])
    v2 = np.array([3, 3])
    result = _apply(_ACTIONS["Max"], v1, v2)
    np.testing.assert_array_equal(result, [3, 3])


# --- Min ---


def test_min():
    v1 = np.array([1, 5, 3])
    v2 = np.array([4, 2, 6])
    result = _apply(_ACTIONS["Min"], v1, v2)
    np.testing.assert_array_equal(result, [1, 2, 3])


def test_min_equal():
    v1 = np.array([3, 3])
    v2 = np.array([3, 3])
    result = _apply(_ACTIONS["Min"], v1, v2)
    np.testing.assert_array_equal(result, [3, 3])
