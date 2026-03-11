from heapq import nlargest
from random import shuffle


def S_diff(lst: list[int | float]) -> list:  # type: ignore[type-arg]
    """Given a list of int or float, calculate S_diff and S_point"""

    S_avg = sum(lst) / len(lst)
    S_dist = [i - S_avg for i in lst]  # distance to average
    S_cum: list[float] = []  # list of cumulative sum
    S_cum.append(0.0)
    for i in range(0, len(S_dist)):
        S_cum.append(S_cum[i] + S_dist[i])
    return [nlargest(1, list(range(0, len(S_cum))), key=lambda i: S_cum[i]), (max(S_cum) - min(S_cum))]
    # return the index of maximum_diff index, and maximum_diff


def bootstrap(lst: list[int | float], obs: float, rep: int = 1000) -> float:
    """Given a list of int or float (lst) and an observation value(obs). calcualte the chance (pvalue)
    of getting this observation through bootstrapping."""

    shuffled_diff: list[list] = []  # type: ignore[type-arg]
    count = 0
    tmp = lst
    for i in range(0, rep):
        shuffle(tmp)
        shuffled_diff.append(S_diff(tmp))

    for i in sorted(shuffled_diff):  # type: ignore[type-var,assignment]
        if i >= obs:  # type: ignore[operator]
            count += 1
    if count / rep < 0.5:
        return count / rep
    else:
        return 1 - count / rep
