"""Turn ``findings.json`` into the table the repository exists to produce.

Every comparison is read through the bootstrap intervals. Where two intervals overlap, the
difference is reported as not established. That is the repo's own argument applied to its
own results, and it applies just as much when the answer is inconvenient.
"""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from statistics import mean

from chembench.evaluate import intervals_overlap

SPLIT_LABELS = {
    "random": "Random",
    "scaffold": "Scaffold",
    "activity_cliff": "Activity cliff",
}
MODEL_LABELS = {
    "ecfp_svm": "ECFP4 + SVM",
    "chemprop": "chemprop (D-MPNN)",
    "foundation": "ChemBERTa + ridge",
}


def _by(cells: list[dict], *keys: str) -> dict[tuple, list[dict]]:
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for cell in cells:
        grouped[tuple(cell[key] for key in keys)].append(cell)
    return grouped


def cell_overlaps(
    cells: list[dict], model_a: str, model_b: str
) -> tuple[int, int, list[dict]]:
    """Count target x split cells where two models' RMSE intervals overlap.

    Returns ``(overlapping, total, separated)``, where ``separated`` lists the cells whose
    intervals do not overlap. This is an unpaired comparison of two percentile intervals,
    which is conservative for two models scored on the same test compounds.
    """
    by_cell = _by(cells, "target_id", "split")
    overlapping = 0
    total = 0
    separated: list[dict] = []
    for (target_id, split), rows in by_cell.items():
        found = {row["model"]: row for row in rows}
        if model_a not in found or model_b not in found:
            continue
        a, b = found[model_a], found[model_b]
        total += 1
        if intervals_overlap((a["rmse_low"], a["rmse_high"]), (b["rmse_low"], b["rmse_high"])):
            overlapping += 1
        else:
            separated.append(
                {
                    "target_id": target_id,
                    "split": split,
                    model_a: a["rmse"],
                    model_b: b["rmse"],
                }
            )
    return overlapping, total, separated


def render(findings: dict) -> str:
    """Render the findings as Markdown."""
    config = findings["configuration"]
    cells = findings["cells"]
    models = config["models"]
    splits = config["splits"]
    lines: list[str] = []

    lines.append("# Results\n")
    lines.append(
        f"{len(findings['targets'])} ChEMBL targets, "
        f"{sum(t['n_compounds'] for t in findings['targets']):,} curated compounds, "
        f"{len(models)} model families, {len(splits)} split regimes.\n"
    )

    # ---- what the splits actually did -------------------------------------------------
    lines.append("## The splits are what they claim to be\n")
    lines.append(
        "Measured per target and averaged. Scaffold leakage is the fraction of test "
        "compounds whose Murcko scaffold also appears in training; cliff enrichment is how "
        "much denser the test set is in activity-cliff compounds than the dataset overall.\n"
    )
    lines.append("| Split | Scaffold leakage | Cliff enrichment |")
    lines.append("| --- | ---: | ---: |")
    by_split = _by(findings["splits"], "name")
    for split in splits:
        rows = by_split[(split,)]
        lines.append(
            f"| {SPLIT_LABELS.get(split, split)} "
            f"| {mean(r['scaffold_leakage'] for r in rows):.1%} "
            f"| {mean(r['cliff_enrichment'] for r in rows):.2f}x |"
        )
    lines.append("")
    targets = findings["targets"]
    if "activity_cliff" in splits and targets:
        ceiling = 1 / config["test_frac"]
        low = min(t["cliff_fraction"] for t in targets)
        high = max(t["cliff_fraction"] for t in targets)
        text = (
            f"Cliff enrichment cannot exceed 1/test_frac ({ceiling:.2f}x here). It reaches "
            "that ceiling whenever the cliff compounds fit inside the test set"
        )
        if high <= config["test_frac"]:
            text += (
                f", which they do for every target ({low:.1%} to {high:.1%} of compounds "
                f"are cliff members against a test fraction of {config['test_frac']:.0%}). "
                "Every cliff compound is then in test and no cliff pair is split across "
                "train and test."
            )
        else:
            text += "."
        lines.append(text + "\n")
    random_leak = mean(r["scaffold_leakage"] for r in by_split[("random",)])
    lines.append(
        f"**{random_leak:.0%} of a random split's test compounds share a scaffold with "
        "something the model trained on.** That is the leak the rest of this table prices.\n"
    )

    # ---- the headline table -----------------------------------------------------------
    lines.append("## RMSE by model and split\n")
    lines.append(
        "pChEMBL units, averaged over targets, with the mean of the per-target 95% "
        "bootstrap intervals.\n"
    )
    header = (
        "| Model | "
        + " | ".join(SPLIT_LABELS.get(s, s) for s in splits)
        + " | Random to scaffold |"
    )
    lines.append(header)
    lines.append("| --- | " + " | ".join("---:" for _ in splits) + " | ---: |")

    grouped = _by(cells, "model", "split")
    for model in models:
        row = [f"| {MODEL_LABELS.get(model, model)} "]
        values = {}
        for split in splits:
            rows = grouped[(model, split)]
            point = mean(r["rmse"] for r in rows)
            low = mean(r["rmse_low"] for r in rows)
            high = mean(r["rmse_high"] for r in rows)
            values[split] = point
            row.append(f"| {point:.3f} [{low:.3f}, {high:.3f}] ")
        if "random" in values and "scaffold" in values:
            delta = values["scaffold"] - values["random"]
            row.append(f"| +{delta:.3f} |" if delta >= 0 else f"| {delta:.3f} |")
        else:
            row.append("| — |")
        lines.append("".join(row))
    lines.append("")

    # ---- model comparison, read through the intervals ---------------------------------
    lines.append("## Which model wins, and whether that is established\n")
    for split in splits:
        lines.append(f"**{SPLIT_LABELS.get(split, split)}.**")
        per_model = {}
        for model in models:
            rows = grouped[(model, split)]
            per_model[model] = (
                mean(r["rmse"] for r in rows),
                mean(r["rmse_low"] for r in rows),
                mean(r["rmse_high"] for r in rows),
            )
        ranked = sorted(per_model.items(), key=lambda kv: kv[1][0])
        best, best_stats = ranked[0]
        text = f" Lowest RMSE: {MODEL_LABELS.get(best, best)} at {best_stats[0]:.3f}."
        if len(ranked) > 1:
            runner, runner_stats = ranked[1]
            if intervals_overlap(
                (best_stats[1], best_stats[2]), (runner_stats[1], runner_stats[2])
            ):
                text += (
                    f" Its interval overlaps {MODEL_LABELS.get(runner, runner)} "
                    f"({runner_stats[0]:.3f}), so the difference is **not established**."
                )
            else:
                text += (
                    f" The interval is separated from {MODEL_LABELS.get(runner, runner)} "
                    f"({runner_stats[0]:.3f}), so the difference **is** established."
                )
        lines.append(text + "\n")

    lines.append(
        "Per cell, read through unpaired overlap of the two 95% bootstrap intervals "
        "(resampled over test compounds, one split seed):\n"
    )
    for i, model_a in enumerate(models):
        for model_b in models[i + 1 :]:
            overlapping, total, separated = cell_overlaps(cells, model_a, model_b)
            if not total:
                continue
            text = (
                f"- {MODEL_LABELS.get(model_a, model_a)} against "
                f"{MODEL_LABELS.get(model_b, model_b)}: intervals overlap in "
                f"{overlapping} of {total} cells."
            )
            if separated and len(separated) <= 3:
                text += " Separated: " + "; ".join(
                    f"{cell['target_id']} {SPLIT_LABELS.get(cell['split'], cell['split'])} "
                    f"({cell[model_a]:.3f} against {cell[model_b]:.3f})"
                    for cell in separated
                )
                text += "."
            lines.append(text)
    lines.append("")

    # ---- per target -------------------------------------------------------------------
    lines.append("## Per target\n")
    lines.append("| Target | Compounds | Scaffolds | Cliff compounds | Cliff pairs |")
    lines.append("| --- | ---: | ---: | ---: | ---: |")
    for target in findings["targets"]:
        lines.append(
            f"| {target['target_id']} | {target['n_compounds']:,} | {target['n_scaffolds']:,} "
            f"| {target['cliff_compounds']:,} ({target['cliff_fraction']:.1%}) "
            f"| {target['cliff_pairs']:,} |"
        )
    lines.append("")

    lines.append(
        "| Target | Split | " + " | ".join(MODEL_LABELS.get(m, m) for m in models) + " |"
    )
    lines.append("| --- | --- | " + " | ".join("---:" for _ in models) + " |")
    per_target = _by(cells, "target_id", "split", "model")
    for target in findings["targets"]:
        for split in splits:
            row = [f"| {target['target_id']} | {SPLIT_LABELS.get(split, split)} "]
            for model in models:
                found = per_target.get((target["target_id"], split, model))
                row.append(f"| {found[0]['rmse']:.3f} " if found else "| — ")
            lines.append("".join(row) + "|")
    lines.append("")

    # ---- configuration ----------------------------------------------------------------
    lines.append("## Configuration\n")
    lines.append(f"- Scaffold variant: `{config['scaffold_variant']}`")
    lines.append(
        f"- Activity cliff: {config['cliff_definition']['similarity']} and "
        f"{config['cliff_definition']['activity']}"
    )
    lines.append(f"- Test fraction: {config['test_frac']}, seed {config['seed']}")
    if config.get("models_unavailable"):
        lines.append(f"- **Not run** (backend missing): {config['models_unavailable']}")
    lines.append("")
    return "\n".join(lines)


def write(findings_path: Path, out_path: Path) -> Path:
    findings = json.loads(findings_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render(findings))
    return out_path
