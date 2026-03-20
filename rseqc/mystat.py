import math


def percentile_list(N: list[int | float]) -> list[int | float] | None:
    """
    Return 100 percentile values (1st–100th) from a sorted list using linear interpolation.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @return - 100 interpolated values, or N itself if len < 100.
    """
    if not N:
        return None
    if len(N) < 100:
        return N
    result: list[int | float] = []
    for i in range(1, 101):
        k = (len(N) - 1) * i / 100.0
        f = math.floor(k)
        c = math.ceil(k)
        if f == c:
            result.append(int(N[int(k)]))
        else:
            d0 = N[int(f)] * (c - k)
            d1 = N[int(c)] * (k - f)
            result.append(int(round(d0 + d1)))
    return result
