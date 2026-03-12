"""Tests for rseqc.mystat."""

import importlib
import math

from rseqc import mystat


def test_import():
    """Verify that rseqc.mystat can be imported."""
    mod = importlib.import_module("rseqc.mystat")
    assert mod is not None


# --- RSS ---


def test_rss_single():
    assert mystat.RSS("3") == 3.0


def test_rss_multiple():
    # sqrt(3^2 + 4^2) = sqrt(25) = 5
    assert mystat.RSS("3,4") == 5.0


def test_rss_pythagorean():
    # sqrt(5^2 + 12^2) = sqrt(169) = 13
    assert mystat.RSS("5,12") == 13.0


# --- H_mean ---


def test_hmean_basic():
    # harmonic mean of 1,2,4 = 3 / (1/1 + 1/2 + 1/4) = 3/1.75
    result = mystat.H_mean("1,2,4")
    assert abs(result - 3 / 1.75) < 1e-10


def test_hmean_with_zero():
    # zeros are skipped; harmonic mean of just [2, 4] = 2/0.75
    result = mystat.H_mean("0,2,4")
    assert abs(result - 2 / 0.75) < 1e-10


def test_hmean_all_zeros():
    result = mystat.H_mean("0,0,0")
    assert result == "NA"


# --- shannon_entropy ---


def test_shannon_entropy_uniform():
    # uniform distribution [1,1,1,1] => entropy = log(4)
    result = mystat.shannon_entropy([1, 1, 1, 1])
    assert abs(result - math.log(4)) < 1e-10


def test_shannon_entropy_single():
    # single element => entropy = 0
    result = mystat.shannon_entropy([5])
    assert result == 0


def test_shannon_entropy_zeros_filtered():
    # zeros are filtered out
    result = mystat.shannon_entropy([0, 0, 5])
    assert result == 0


# --- simpson_index ---


def test_simpson_index_uniform():
    # uniform [1,1,1,1]: simpson = 1 - 4*(1/4)^2 = 1 - 0.25 = 0.75
    result = mystat.simpson_index("1,1,1,1")
    assert abs(result - 0.75) < 1e-10


def test_simpson_index_single():
    # single element: simpson = 1 - 1 = 0
    result = mystat.simpson_index("5")
    assert abs(result - 0.0) < 1e-10


# --- percentile ---


def test_percentile_median():
    data = [1, 2, 3, 4, 5]
    result = mystat.percentile(data, 50)
    assert result == 3


def test_percentile_min():
    data = [1, 2, 3, 4, 5]
    result = mystat.percentile(data, 0)
    assert result == 1


def test_percentile_max():
    data = [1, 2, 3, 4, 5]
    result = mystat.percentile(data, 100)
    assert result == 5


def test_percentile_interpolation():
    data = [10, 20, 30, 40]
    result = mystat.percentile(data, 25)
    # k = 3 * 0.25 = 0.75, f=0, c=1
    # d0 = 10 * (1 - 0.75) = 2.5, d1 = 20 * 0.75 = 15
    assert abs(result - 17.5) < 1e-10


def test_percentile_empty():
    assert mystat.percentile([], 50) is None


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


# --- Hill_number ---


def test_hill_number_uniform():
    # Hill number with q=2 for uniform should be N (number of species)
    result = mystat.Hill_number("1,1,1,1", qvalue=2)
    assert abs(result - 4.0) < 1e-10


# --- shannon_entropy_es ---


def test_shannon_entropy_es_single():
    result = mystat.shannon_entropy_es([5])
    assert result == 0


def test_shannon_entropy_es_empty():
    result = mystat.shannon_entropy_es([])
    assert result == "NA"


def test_shannon_entropy_es_multi():
    """Chao-Shen estimator with multi-element input should be positive."""
    result = mystat.shannon_entropy_es([1, 2, 3, 4])
    assert isinstance(result, float)
    assert result > 0


# --- shannon_entropy_ht ---


def test_shannon_entropy_ht_uniform():
    # HT estimator for uniform distribution, should be positive
    result = mystat.shannon_entropy_ht("1,1,1,1")
    assert isinstance(result, float)
    assert result > 0


def test_shannon_entropy_ht_single():
    result = mystat.shannon_entropy_ht("5")
    assert result == 0


def test_shannon_entropy_ht_na():
    # all zeros filtered out => empty => NA
    result = mystat.shannon_entropy_ht("0,0,0")
    assert result == "NA"


def test_shannon_entropy_ht_two_values():
    result = mystat.shannon_entropy_ht("1,3")
    assert isinstance(result, float)
    assert result > 0


# --- simpson_index_es ---


def test_simpson_index_es_uniform():
    result = mystat.simpson_index_es("1,1,1,1")
    # For n_i = 1 each, sum(n_i * (n_i - 1)) = 0, so index = 1 - 0 = 1.0
    # But sum = 4, so 1 - 0 / (4*3) = 1.0
    assert abs(result - 1.0) < 1e-10


def test_simpson_index_es_single():
    result = mystat.simpson_index_es("5")
    # 1 - 5*4/(5*4) = 0
    assert abs(result - 0.0) < 1e-10


def test_simpson_index_es_two_values():
    result = mystat.simpson_index_es("3,7")
    assert isinstance(result, float)
    assert 0 <= result <= 1


# --- Hill_number extended ---


def test_hill_number_q_half():
    # q=0.5, should be >= q=2 result (favors rare species)
    result_half = mystat.Hill_number("1,1,1,1", qvalue=0.5)
    result_two = mystat.Hill_number("1,1,1,1", qvalue=2)
    # For uniform distribution, all q values give same result = 4
    assert abs(result_half - 4.0) < 1e-10
    assert abs(result_two - 4.0) < 1e-10


def test_hill_number_q_half_nonuniform():
    result_half = mystat.Hill_number("1,2,3,4", qvalue=0.5)
    result_two = mystat.Hill_number("1,2,3,4", qvalue=2)
    # q=0.5 favors rare species, so result should be >= q=2
    assert result_half >= result_two


def test_hill_number_q1():
    """q=1 should return exp(H) where H is Shannon entropy."""
    # For uniform distribution [1,1,1,1], exp(log(4)) = 4.0
    result = mystat.Hill_number("1,1,1,1", qvalue=1)
    assert abs(result - 4.0) < 1e-10


def test_hill_number_q1_nonuniform():
    """q=1 for non-uniform distribution."""
    result = mystat.Hill_number("1,2,3,4", qvalue=1)
    assert isinstance(result, float)
    assert result > 0
