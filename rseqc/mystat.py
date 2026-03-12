import math
from collections.abc import Callable
from typing import Any


def RSS(arg: str) -> Any:
    """calculate Square root of sum of square. Input is ',' separated numbers"""
    lst = arg.split(",")
    return sum(int(i) ** 2 for i in lst) ** 0.5


def H_mean(arg: str) -> float | str:
    """calculate harmornic mean. Input is ',' separated numbers"""
    lst = [1 / float(i) for i in arg.split(",") if float(i) != 0]
    if len(lst) == 0:
        return "NA"
    else:
        return len(lst) / (sum(lst))


def shannon_entropy(arg: list[int | float]) -> float | int:
    """calculate shannon's entropy (or Shannon-Wiener index)."""
    lst = arg
    lst = [float(i) for i in lst if float(i) > 0]
    total = sum(lst)
    entropy = 0.0
    for i in lst:
        freq = i / total
        entropy += freq * math.log(freq)
    if entropy == 0:
        return 0
    else:
        return -entropy


def shannon_entropy_es(arg: list[int | float]) -> float | str | int:
    """calculate estimator of shannon's entropy (Chao & Shen, 2003)"""
    lst = arg
    lst = [float(i) for i in lst if float(i) > 0]
    if sum(lst) <= 0 or min(lst) < 0:
        return "NA"  # if there is no fragmental splicing
    if len(lst) == 1:
        return 0  # if there is only 1 fragmental splicing
    lst.append(2)

    # estimate C_bar
    singleton = 0
    entropy = 0.0
    for i in lst:
        if i == 1:
            singleton += 1

    total = sum(lst)
    C_bar = 1 - (singleton / total)
    for i in lst:
        freq = C_bar * i / total
        entropy += (freq * math.log(freq)) / (1 - (1 - freq) ** total)
    if entropy == 0:
        return 0
    else:
        return -entropy


def shannon_entropy_ht(arg: str) -> float | str | int:
    """calculate estimator of shannon's entropy based on Horzitz-Thompson"""
    lst: list[float] = [float(i) for i in arg.split(",") if float(i) > 0]
    if sum(lst) <= 0 or min(lst) < 0:
        return "NA"  # if there is no fragmental splicing
    if len(lst) == 1:
        return 0  # if there is only 1 fragmental splicing

    # estimate C_bar
    total = sum(lst)
    entropy = 0.0
    for i in lst:
        freq = i / total
        entropy += (freq * math.log(freq)) / (1 - (1 - freq) ** total)
    return -entropy


def simpson_index(arg: str) -> float | int:
    """calculate Gini-Simpson's index. Input is ',' separated numbers"""
    lst: list[float] = [float(i) for i in arg.split(",") if float(i) > 0]
    simpson = 0.0

    try:
        total = sum(lst)
        for i in lst:
            simpson = simpson + (i / total) ** 2
        return 1 - simpson
    except ZeroDivisionError:
        return 0


def simpson_index_es(arg: str) -> float | int:
    """calculate estimator Gini-Simpson's index. Input is ',' separated numbers"""
    lst: list[float] = [float(i) for i in arg.split(",") if float(i) > 0]
    simpson = 0.0

    try:
        total = sum(lst)
        for i in lst:
            simpson = simpson + i * (i - 1)
        return 1 - (simpson / (total * (total - 1)))
    except ZeroDivisionError:
        return 0


def Hill_number(arg: str, qvalue: int | float = 1) -> float:
    """Calculate real diversity (Hill's number). Input is ',' separated numbers. qvalue is the only
    parameter for Hill's function. When q=1, it return exp(H) which is the effective number of junctions
    calculated by Shannon's entropy. When q<1, Hill's function was favors low frequency junctions.
    When q>1, Hill's function was favors high frequency junctions (common junctions). Simpon's Index
    is particular case of Hill's function as q=2"""

    lst: list[float] = [float(i) for i in arg.split(",") if float(i) > 0]
    total = sum(lst)
    freq = [(i / total) ** qvalue for i in lst]
    try:
        return (sum(freq)) ** (1 / (1 - qvalue))  # type: ignore[no-any-return]
    except ZeroDivisionError:
        return math.exp(shannon_entropy([float(i) for i in arg.split(",") if float(i) > 0]))


def percentile(N: list[Any], percent: float, key: Callable[[Any], float] = lambda x: x) -> float | None:
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0 to 100.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N) - 1) * percent / 100.0
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c - k)
    d1 = key(N[int(c)]) * (k - f)
    return d0 + d1


def percentile_list(N: list[int | float]) -> list[int | float] | None:
    """
    Find the percentile of a list of values.
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @return - the list of percentile of the values
    """
    if not N:
        return None
    if len(N) < 100:
        return N
    per_list: list[int | float] = []
    for i in range(1, 101):
        k = (len(N) - 1) * i / 100.0
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            per_list.append(int(N[int(k)]))
        else:
            d0 = N[int(f)] * (c - k)
            d1 = N[int(c)] * (k - f)
            per_list.append(int(round(d0 + d1)))
    return per_list
