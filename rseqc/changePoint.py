from heapq import nlargest


def S_diff(lst: list[int | float]) -> list[list[int] | float]:
    """Given a list of int or float, calculate S_diff and S_point"""

    S_avg = sum(lst) / len(lst)
    S_dist = [i - S_avg for i in lst]  # distance to average
    S_cum: list[float] = []  # list of cumulative sum
    S_cum.append(0.0)
    for i in range(0, len(S_dist)):
        S_cum.append(S_cum[i] + S_dist[i])
    return [nlargest(1, list(range(0, len(S_cum))), key=lambda i: S_cum[i]), (max(S_cum) - min(S_cum))]
    # return the index of maximum_diff index, and maximum_diff
