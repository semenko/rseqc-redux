"""manipulate ndarray list"""

from typing import Any

from numpy.typing import NDArray


def check_list(v1: NDArray[Any], v2: NDArray[Any]) -> None:
    """check if the length of two list is same"""
    if v1.size != v2.size:
        raise ValueError("the lenght of both arrays must be the same")
    pass


def Add(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """add two list"""
    check_list(v1, v2)
    return v1.__add__(v2)


def Subtract(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """subtract v2 from v1"""
    check_list(v1, v2)
    return v1.__sub__(v2)


def Product(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """return product of two list"""
    check_list(v1, v2)
    return v1.__mul__(v2)


def Division(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """return divide v1 by v2. add 1 to both v1 and v2"""
    check_list(v1, v2)
    return (v1 + 1).__div__(v2 + 1)  # type: ignore[attr-defined]


def Average(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """return arithmetic mean of two list"""
    check_list(v1, v2)
    return v1.__add__(v2) / 2


def geometricMean(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """return geometric mean of two list"""
    check_list(v1, v2)
    return (v1.__mul__(v2)) ** 0.5


def Max(v1: NDArray[Any], v2: NDArray[Any]) -> map:
    """pairwise comparison two list. return  the max one between two paried number"""
    check_list(v1, v2)
    return map(max, zip(v1, v2))


def Min(v1: NDArray[Any], v2: NDArray[Any]) -> map:
    """pairwise comparison two list. return  the max one between two paried number"""
    check_list(v1, v2)
    return map(min, zip(v1, v2))


def euclidean_distance(v1: NDArray[Any], v2: NDArray[Any]) -> Any:
    """return euclidean distance"""
    check_list(v1, v2)
    return (sum((v1.__sub__(v2)) ** 2) / v1.size) ** 0.5
