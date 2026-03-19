import numpy as np


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
    result = np.percentile(N, range(1, 101), method="linear")
    return [int(round(x)) for x in result]
