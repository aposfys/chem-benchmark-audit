"""The experiment: every model, every split regime, one curation.

The design constraint that makes the result mean anything is that curation happens once,
per target, before any splitting -- so the only thing that varies between cells of the
table is the split and the model. If a number moves, there is exactly one reason it could
have.
"""

from __future__ import annotations

import json
import time
from collections.abc import Sequence
from dataclasses import asdict, dataclass
from pathlib import Path

from chembench.curate import (
    Record,
    find_activity_cliffs,
    murcko_scaffold,
    read_curated,
)
from chembench.evaluate import bootstrap_ci, rmse, spearman
from chembench.models import available, build
from chembench.splits import (
    activity_cliff_split,
    cliff_enrichment,
    random_split,
    scaffold_leakage,
    scaffold_split,
)

SPLIT_NAMES = ("random", "scaffold", "activity_cliff")


@dataclass
class Cell:
    """One (target, split, model) result."""

    target_id: str
    split: str
    model: str
    n_train: int
    n_test: int
    rmse: float
    rmse_low: float
    rmse_high: float
    spearman: float
    spearman_low: float
    spearman_high: float
    seconds: float


@dataclass
class SplitDescription:
    """What a split regime actually produced, so the label can be checked."""

    name: str
    n_train: int
    n_test: int
    scaffold_leakage: float
    cliff_enrichment: float


def prepare_target(
    path: Path, *, generic_scaffolds: bool = False, seed: int = 0
) -> tuple[list[Record], dict[str, str], set[str], int]:
    """Load a curated target and compute everything the splits need.

    Scaffolds and cliffs are computed once per target and shared across every split and
    model, so no cell can differ because it saw a different definition.
    """
    records, _ = read_curated(path)
    scaffolds = {
        record.record_id: murcko_scaffold(record.smiles, generic=generic_scaffolds)
        for record in records
    }
    cliff_members, cliff_pairs = find_activity_cliffs(records)
    return records, scaffolds, cliff_members, cliff_pairs


def make_splits(
    records: Sequence[Record],
    scaffolds: dict[str, str],
    cliff_members: set[str],
    *,
    test_frac: float = 0.2,
    seed: int = 0,
) -> dict[str, tuple[list[str], list[str], SplitDescription]]:
    """Build all three splits and describe each one."""
    keys = [record.record_id for record in records]
    baseline = len(cliff_members) / len(keys) if keys else 0.0

    built = {
        "random": random_split(keys, test_frac=test_frac, seed=seed),
        "scaffold": scaffold_split(scaffolds, test_frac=test_frac, seed=seed),
        "activity_cliff": activity_cliff_split(
            keys, cliff_members, test_frac=test_frac, seed=seed
        ),
    }

    described = {}
    for name, (train, test) in built.items():
        described[name] = (
            train,
            test,
            SplitDescription(
                name=name,
                n_train=len(train),
                n_test=len(test),
                scaffold_leakage=scaffold_leakage(train, test, scaffolds),
                cliff_enrichment=cliff_enrichment(test, cliff_members, baseline),
            ),
        )
    return described


def run_cell(
    records: Sequence[Record],
    train_keys: Sequence[str],
    test_keys: Sequence[str],
    model_name: str,
    *,
    target_id: str,
    split_name: str,
    seed: int = 0,
) -> Cell:
    """Fit one model on one split and score it."""
    by_id = {record.record_id: record for record in records}
    train_smiles = [by_id[key].smiles for key in train_keys]
    train_y = [by_id[key].activity for key in train_keys]
    test_smiles = [by_id[key].smiles for key in test_keys]
    test_y = [by_id[key].activity for key in test_keys]

    started = time.time()
    model = build(model_name, seed=seed)
    model.fit(train_smiles, train_y)
    predictions = model.predict(test_smiles)
    elapsed = time.time() - started

    rmse_low, rmse_high = bootstrap_ci(test_y, predictions, statistic="rmse", seed=seed)
    rho_low, rho_high = bootstrap_ci(test_y, predictions, statistic="spearman", seed=seed)

    return Cell(
        target_id=target_id,
        split=split_name,
        model=model_name,
        n_train=len(train_keys),
        n_test=len(test_keys),
        rmse=rmse(test_y, predictions),
        rmse_low=rmse_low,
        rmse_high=rmse_high,
        spearman=spearman(test_y, predictions),
        spearman_low=rho_low,
        spearman_high=rho_high,
        seconds=elapsed,
    )


def run(
    curated_dir: Path,
    results_dir: Path,
    *,
    models: Sequence[str],
    splits: Sequence[str],
    generic_scaffolds: bool = False,
    test_frac: float = 0.2,
    seed: int = 0,
    max_targets: int = 0,
) -> dict:
    """Run the whole grid and write ``findings.json``."""
    usable = []
    unavailable = {}
    for name in models:
        ok, why = available(name)
        (usable.append(name) if ok else unavailable.update({name: why}))
    if not usable:
        raise RuntimeError(f"no model backends available: {unavailable}")

    paths = sorted(curated_dir.glob("CHEMBL*.json"))
    if max_targets:
        paths = paths[:max_targets]
    if not paths:
        raise RuntimeError(f"no curated targets in {curated_dir}; run `make data` first")

    cells: list[Cell] = []
    descriptions: list[dict] = []
    target_summaries: list[dict] = []

    for path in paths:
        records, scaffolds, cliff_members, cliff_pairs = prepare_target(
            path, generic_scaffolds=generic_scaffolds, seed=seed
        )
        target_id = records[0].target_id
        built = make_splits(records, scaffolds, cliff_members, test_frac=test_frac, seed=seed)
        target_summaries.append(
            {
                "target_id": target_id,
                "n_compounds": len(records),
                "n_scaffolds": len(set(scaffolds.values())),
                "cliff_compounds": len(cliff_members),
                "cliff_pairs": cliff_pairs,
                "cliff_fraction": round(len(cliff_members) / len(records), 4),
            }
        )
        print(
            f"{target_id}: {len(records)} compounds, {len(set(scaffolds.values()))} scaffolds, "
            f"{len(cliff_members)} cliff compounds ({cliff_pairs} pairs)",
            flush=True,
        )

        for split_name in splits:
            train, test, description = built[split_name]
            descriptions.append({"target_id": target_id, **asdict(description)})
            for model_name in usable:
                cell = run_cell(
                    records,
                    train,
                    test,
                    model_name,
                    target_id=target_id,
                    split_name=split_name,
                    seed=seed,
                )
                cells.append(cell)
                print(
                    f"  {split_name:<15} {model_name:<11} "
                    f"RMSE {cell.rmse:.3f} [{cell.rmse_low:.3f}, {cell.rmse_high:.3f}]  "
                    f"rho {cell.spearman:.3f}  ({cell.seconds:.0f}s)",
                    flush=True,
                )

    findings = {
        "configuration": {
            "models": usable,
            "models_unavailable": unavailable,
            "splits": list(splits),
            "test_frac": test_frac,
            "seed": seed,
            # "scaffold split" alone is not a reproducible statement.
            "scaffold_variant": "generic" if generic_scaffolds else "bemis-murcko",
            "cliff_definition": {
                "similarity": "ECFP4 Tanimoto >= 0.9",
                "activity": "|delta pChEMBL| >= 1.0",
            },
        },
        "targets": target_summaries,
        "splits": descriptions,
        "cells": [asdict(cell) for cell in cells],
    }
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
    return findings
