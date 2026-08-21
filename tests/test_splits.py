"""The splitting invariants the whole repo depends on.

These are not parsing tests. Each one asserts a property that, if it broke silently,
would make every downstream number wrong in the direction that flatters the model.
"""

from __future__ import annotations

import pytest

from chembench import splits

# Twelve compounds over four scaffolds, deliberately unbalanced so that a
# largest-group-first assignment has something to decide.
SCAFFOLDS = {
    "c01": "c1ccccc1",
    "c02": "c1ccccc1",
    "c03": "c1ccccc1",
    "c04": "c1ccccc1",
    "c05": "c1ccncc1",
    "c06": "c1ccncc1",
    "c07": "c1ccncc1",
    "c08": "C1CCCCC1",
    "c09": "C1CCCCC1",
    "c10": "C1CCNCC1",
    "c11": "C1CCNCC1",
    "c12": "C1CCNCC1",
}


def test_scaffold_split_leaks_nothing() -> None:
    """A scaffold split that leaks at all is not a scaffold split."""
    train, test = splits.scaffold_split(SCAFFOLDS, test_frac=0.25)
    assert splits.scaffold_leakage(train, test, SCAFFOLDS) == 0.0


def test_random_split_leaks_and_that_is_the_point() -> None:
    """The baseline regime must demonstrably leak, or the repo has no subject."""
    train, test = splits.random_split(list(SCAFFOLDS), test_frac=0.25, seed=0)
    assert splits.scaffold_leakage(train, test, SCAFFOLDS) > 0.0


def test_splits_partition_the_input() -> None:
    """No compound may be dropped or duplicated by either regime."""
    for train, test in (
        splits.scaffold_split(SCAFFOLDS, test_frac=0.25),
        splits.random_split(list(SCAFFOLDS), test_frac=0.25, seed=0),
    ):
        assert sorted(train + test) == sorted(SCAFFOLDS)
        assert not set(train) & set(test)


def test_scaffold_split_keeps_groups_whole() -> None:
    """A scaffold group split across the boundary is leakage by another name."""
    train, test = splits.scaffold_split(SCAFFOLDS, test_frac=0.25)
    train_scaffolds = {SCAFFOLDS[key] for key in train}
    test_scaffolds = {SCAFFOLDS[key] for key in test}
    assert not train_scaffolds & test_scaffolds


def test_scaffold_split_is_deterministic() -> None:
    assert splits.scaffold_split(SCAFFOLDS, test_frac=0.25) == splits.scaffold_split(
        SCAFFOLDS, test_frac=0.25
    )


def test_duplicate_groups_finds_repeated_structures() -> None:
    """Duplicates are the quiet form of leakage; they must never be silently tolerated."""
    identifiers = {"a": "KEY-1", "b": "KEY-1", "c": "KEY-2"}
    assert splits.duplicate_groups(identifiers) == {"KEY-1": ["a", "b"]}


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.1, 1.5])
def test_rejects_impossible_fractions(bad: float) -> None:
    with pytest.raises(ValueError):
        splits.random_split(list(SCAFFOLDS), test_frac=bad)
