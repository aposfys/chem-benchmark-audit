"""Metrics, bootstrap intervals, and the cross-split comparison table.

The deliverable of this module is one table: model x split regime, with intervals. A
difference between two models is only reported when their intervals do not overlap --
the whole argument of the repo is that small differences on a leaky split are noise.
"""

from __future__ import annotations

import math
import random
from collections.abc import Sequence


def rmse(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """Root mean squared error."""
    _check_pair(y_true, y_pred)
    total = sum((a - b) ** 2 for a, b in zip(y_true, y_pred, strict=True))
    return math.sqrt(total / len(y_true))


def spearman(y_true: Sequence[float], y_pred: Sequence[float]) -> float:
    """Spearman rank correlation, with average ranks for ties."""
    _check_pair(y_true, y_pred)
    rank_true = _ranks(y_true)
    rank_pred = _ranks(y_pred)
    n = len(y_true)
    mean = (n - 1) / 2.0
    cov = sum((a - mean) * (b - mean) for a, b in zip(rank_true, rank_pred, strict=True))
    var_true = sum((a - mean) ** 2 for a in rank_true)
    var_pred = sum((b - mean) ** 2 for b in rank_pred)
    if var_true == 0.0 or var_pred == 0.0:
        return 0.0
    return cov / math.sqrt(var_true * var_pred)


def bootstrap_ci(
    y_true: Sequence[float],
    y_pred: Sequence[float],
    *,
    statistic: str = "rmse",
    n_resamples: int = 1000,
    alpha: float = 0.05,
    seed: int = 0,
) -> tuple[float, float]:
    """Percentile bootstrap interval for ``rmse`` or ``spearman``.

    Resampling is over test compounds, which is the source of variation that matters: two
    models whose intervals overlap have not been shown to differ. Every comparison in this
    repo is read through these intervals rather than through the point estimates.

    Degenerate resamples -- ones where every drawn label is identical, which makes Spearman
    undefined -- are discarded rather than scored as zero correlation, which would drag the
    lower bound down for a reason that has nothing to do with the model.
    """
    _check_pair(y_true, y_pred)
    func = {"rmse": rmse, "spearman": spearman}[statistic]
    rng = random.Random(seed)
    n = len(y_true)
    values = []
    for _ in range(n_resamples):
        idx = [rng.randrange(n) for _ in range(n)]
        resampled_true = [y_true[i] for i in idx]
        if len(set(resampled_true)) < 2:
            continue
        values.append(func(resampled_true, [y_pred[i] for i in idx]))
    if not values:
        point = func(y_true, y_pred)
        return point, point
    values.sort()
    lo = values[int((alpha / 2) * len(values))]
    hi = values[min(int((1 - alpha / 2) * len(values)), len(values) - 1)]
    return lo, hi


def intervals_overlap(a: tuple[float, float], b: tuple[float, float]) -> bool:
    """Whether two confidence intervals overlap.

    The repo's central claim is that small differences on a leaky split are noise, so a
    difference is only ever reported as a difference when this returns ``False``.
    """
    return not (a[1] < b[0] or b[1] < a[0])


def _ranks(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        average = (i + j) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    return ranks


def _check_pair(y_true: Sequence[float], y_pred: Sequence[float]) -> None:
    if len(y_true) != len(y_pred):
        raise ValueError(f"length mismatch: {len(y_true)} true vs {len(y_pred)} predicted")
    if not y_true:
        raise ValueError("cannot score an empty set")
