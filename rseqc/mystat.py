import math


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
