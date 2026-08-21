"""Splitting regimes, and the leakage report that justifies using them.

Three regimes are supported, in increasing order of honesty:

``random``
    What most published numbers use. Near-identical molecules land on both sides.

``scaffold``
    Bemis-Murcko scaffold groups are kept whole and assigned entirely to one side.
    Note that "scaffold split" is not one algorithm -- generic and non-generic atom
    typing give materially different difficulty -- so the variant used is recorded
    alongside every result.

``activity_cliff``
    Test set enriched for pairs that are structurally similar but differ sharply in
    potency, in the spirit of MoleculeACE. This is where descriptor baselines have
    repeatedly matched or beaten deep models.

Nothing here imports RDKit. Scaffolds arrive as a precomputed mapping so that the
splitting logic stays pure, fast and testable without a chemistry toolkit installed.
"""

from __future__ import annotations

import random
from collections import defaultdict
from collections.abc import Mapping, Sequence

Split = tuple[list[str], list[str]]


def random_split(keys: Sequence[str], test_frac: float = 0.2, seed: int = 0) -> Split:
    """Shuffle and cut. Provided as the baseline to be argued against, not with."""
    _check_frac(test_frac)
    shuffled = list(keys)
    random.Random(seed).shuffle(shuffled)
    n_test = round(len(shuffled) * test_frac)
    return shuffled[n_test:], shuffled[:n_test]


def scaffold_split(
    scaffolds: Mapping[str, str], test_frac: float = 0.2, seed: int = 0
) -> Split:
    """Assign whole scaffold groups to one side, largest group first.

    Largest-first is deterministic and is what most published scaffold splits do, but it
    biases the test set towards singleton scaffolds. ``seed`` only breaks ties between
    equally sized groups, so the split is reproducible without being an accident of
    dictionary ordering.
    """
    _check_frac(test_frac)
    groups: dict[str, list[str]] = defaultdict(list)
    for key, scaffold in scaffolds.items():
        groups[scaffold].append(key)

    rng = random.Random(seed)
    ordered = sorted(groups.values(), key=lambda members: (-len(members), rng.random()))

    n_test_target = round(len(scaffolds) * test_frac)
    train: list[str] = []
    test: list[str] = []
    for members in ordered:
        if len(test) < n_test_target:
            test.extend(members)
        else:
            train.extend(members)
    return train, test


def activity_cliff_split(
    scaffolds: Mapping[str, str],
    activities: Mapping[str, float],
    test_frac: float = 0.2,
    seed: int = 0,
) -> Split:
    """Enrich the test set for activity cliffs.

    Requires a similarity matrix, so it belongs downstream of fingerprinting rather than
    here. Implemented in the modelling milestone; see README milestone 3.
    """
    raise NotImplementedError("activity-cliff splitting lands with the fingerprint backend")


def scaffold_leakage(
    train: Sequence[str], test: Sequence[str], scaffolds: Mapping[str, str]
) -> float:
    """Fraction of test compounds whose scaffold is also present in train.

    This is the number the whole repo turns on. A random split typically returns
    something far from zero; a correct scaffold split must return exactly zero, and the
    test suite asserts that rather than trusting it.
    """
    if not test:
        return 0.0
    train_scaffolds = {scaffolds[key] for key in train}
    leaked = sum(1 for key in test if scaffolds[key] in train_scaffolds)
    return leaked / len(test)


def duplicate_groups(identifiers: Mapping[str, str]) -> dict[str, list[str]]:
    """Group record ids by structural identifier (InChIKey, or its skeleton block).

    Duplicates are the quiet form of leakage: the same parent compound recorded twice
    under different salts or assay ids can be split across train and test, which reads
    as generalisation and is memorisation.
    """
    groups: dict[str, list[str]] = defaultdict(list)
    for record_id, identifier in identifiers.items():
        groups[identifier].append(record_id)
    return {key: members for key, members in groups.items() if len(members) > 1}


def _check_frac(test_frac: float) -> None:
    if not 0.0 < test_frac < 1.0:
        raise ValueError(f"test_frac must lie in (0, 1), got {test_frac}")
