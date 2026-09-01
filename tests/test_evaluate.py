"""Metric invariants. A metric that is wrong in a flattering direction is worse than none."""

from __future__ import annotations

import pytest

from chembench import evaluate


def test_rmse_is_zero_for_perfect_prediction() -> None:
    assert evaluate.rmse([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]) == 0.0


def test_rmse_matches_hand_computation() -> None:
    """Errors of 1 and 2 give sqrt((1 + 4) / 2)."""
    assert evaluate.rmse([1.0, 2.0], [2.0, 4.0]) == pytest.approx(2.5**0.5)


def test_spearman_is_one_for_monotone_prediction() -> None:
    """Spearman ranks, so a monotone but badly scaled prediction still scores 1."""
    assert evaluate.spearman([1.0, 2.0, 3.0], [10.0, 200.0, 3000.0]) == pytest.approx(1.0)


def test_spearman_is_minus_one_when_reversed() -> None:
    assert evaluate.spearman([1.0, 2.0, 3.0], [3.0, 2.0, 1.0]) == pytest.approx(-1.0)


def test_spearman_handles_ties_without_dividing_by_zero() -> None:
    assert evaluate.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]) == 0.0


def test_bootstrap_interval_brackets_the_point_estimate() -> None:
    y_true = [float(i) for i in range(40)]
    y_pred = [float(i) + (0.5 if i % 2 else -0.5) for i in range(40)]
    point = evaluate.rmse(y_true, y_pred)
    lo, hi = evaluate.bootstrap_ci(y_true, y_pred, n_resamples=200, seed=0)
    assert lo <= point <= hi


def test_length_mismatch_is_an_error_not_a_truncation() -> None:
    with pytest.raises(ValueError):
        evaluate.rmse([1.0, 2.0], [1.0])


def test_bootstrap_survives_a_constant_truth_vector() -> None:
    """Every resample is degenerate here; Spearman is undefined and must not be faked."""
    low, high = evaluate.bootstrap_ci(
        [1.0] * 8,
        [0.3, 0.9, 0.2, 0.5, 0.7, 0.1, 0.4, 0.6],
        statistic="spearman",
        n_resamples=50,
    )
    assert low == high  # falls back to the point estimate rather than inventing spread


def test_intervals_overlap_detects_both_orders() -> None:
    assert evaluate.intervals_overlap((0.1, 0.5), (0.4, 0.9))
    assert evaluate.intervals_overlap((0.4, 0.9), (0.1, 0.5))
    assert not evaluate.intervals_overlap((0.1, 0.3), (0.4, 0.9))
    assert not evaluate.intervals_overlap((0.4, 0.9), (0.1, 0.3))


def test_intervals_that_touch_count_as_overlapping() -> None:
    assert evaluate.intervals_overlap((0.1, 0.5), (0.5, 0.9))
