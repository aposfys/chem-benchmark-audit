"""The rendered report, and the per-cell counts it prints."""

from __future__ import annotations

from chembench import report


def _cell(target: str, split: str, model: str, rmse: float, half: float) -> dict:
    return {
        "target_id": target,
        "split": split,
        "model": model,
        "rmse": rmse,
        "rmse_low": rmse - half,
        "rmse_high": rmse + half,
    }


CELLS = [
    _cell("T1", "random", "ecfp_svm", 0.70, 0.05),
    _cell("T1", "random", "chemprop", 0.72, 0.05),
    _cell("T1", "scaffold", "ecfp_svm", 0.75, 0.03),
    _cell("T1", "scaffold", "chemprop", 0.90, 0.03),
    _cell("T2", "random", "ecfp_svm", 0.60, 0.05),
    _cell("T2", "random", "chemprop", 0.61, 0.05),
    _cell("T2", "scaffold", "ecfp_svm", 0.80, 0.05),
    _cell("T2", "scaffold", "chemprop", 0.79, 0.05),
]


def test_cell_overlaps_counts_and_names_the_separated_cell() -> None:
    overlapping, total, separated = report.cell_overlaps(CELLS, "ecfp_svm", "chemprop")
    assert (overlapping, total) == (3, 4)
    assert [(c["target_id"], c["split"]) for c in separated] == [("T1", "scaffold")]


def test_render_prints_the_overlap_count_and_the_enrichment_ceiling() -> None:
    findings = {
        "configuration": {
            "models": ["ecfp_svm", "chemprop"],
            "splits": ["random", "scaffold", "activity_cliff"],
            "scaffold_variant": "bemis-murcko",
            "cliff_definition": {"similarity": "s", "activity": "a"},
            "test_frac": 0.2,
            "seed": 0,
        },
        "targets": [
            {
                "target_id": t,
                "n_compounds": 100,
                "n_scaffolds": 40,
                "cliff_compounds": 5,
                "cliff_fraction": 0.05,
                "cliff_pairs": 3,
            }
            for t in ("T1", "T2")
        ],
        "splits": [
            {"target_id": t, "name": s, "scaffold_leakage": 0.5, "cliff_enrichment": 1.0}
            for t in ("T1", "T2")
            for s in ("random", "scaffold", "activity_cliff")
        ],
        "cells": CELLS
        + [
            _cell(t, "activity_cliff", m, 0.8, 0.05)
            for t in ("T1", "T2")
            for m in ("ecfp_svm", "chemprop")
        ],
    }
    text = report.render(findings)
    assert "intervals overlap in 5 of 6 cells" in text
    assert "T1 Scaffold (0.750 against 0.900)" in text
    assert "1/test_frac (5.00x here)" in text
