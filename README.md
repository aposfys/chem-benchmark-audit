# chem-benchmark-audit — how much of molecular ML survives an honest split?

[![CI](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Ruff](https://img.shields.io/badge/lint-ruff-261230)](https://docs.astral.sh/ruff/)

> **Status: skeleton.** Structure, interfaces and tests are in place; the analysis is not
> run yet. The headline below is the claim this repo exists to earn or refute — it is not
> a result. No number here is real until `make analysis` produces it.

Molecular property prediction is reported almost everywhere on random splits. Random splits
put near-identical molecules in both train and test, so the model is rewarded for
memorising a scaffold rather than generalising from it. This repo measures the size of that
effect on curated ChEMBL targets, under three splitting regimes, across three model families
of increasing complexity.

**The question:** when the split stops leaking, how much of the deep model's advantage over
a 2010-era baseline is left?

| | |
| --- | --- |
| **Source** | ChEMBL, curated in-repo (Checker → Standardizer → GetParent) |
| **Targets** | 3–5, selected for assay consistency and activity-cliff density |
| **Splits** | random · Bemis–Murcko scaffold · activity-cliff (MoleculeACE-style) |
| **Models** | ECFP4 + SVM · chemprop (D-MPNN) · a molecular foundation model |
| **Metric** | RMSE / Spearman, reported per split regime, with confidence intervals |

## Why this is worth doing

The field is currently saying this out loud, in print:

- **LIT-PCBA** was constructed specifically as the unbiased answer to DUD-E, and a 2025
  audit found it carries data leakage, duplication and pervasive analog redundancy.
- **DUD-E and MUV** carry analog bias and decoy artifacts — models separate actives from
  decoys on crude physicochemical properties rather than on binding.
- **MoleculeACE** (30 ChEMBL targets, 48.7K molecules) repeatedly shows descriptor + SVM
  matching or beating deep models once activity cliffs are in the test set.
- **Polaris** exists as the community's response to exactly this problem.

What is missing is a single, reproducible, end-to-end repo where the *same* curation feeds
the *same* models under all three splits, so the size of the effect is directly readable.

## Traps this pipeline is built to avoid

Filled in as they are hit. Two are already known and are handled in `curate.py`:

- **ChEMBL does not canonicalise tautomers; PubChem does.** Merging bioactivity data from
  both sources without a common tautomer convention lets the model learn which *database* a
  compound came from. Standardisation happens before any merge, and the test suite asserts
  a single tautomer convention across the merged set.
- **Salts and solvates inflate apparent diversity.** Two records of the same parent compound
  as different salts look like two molecules to a fingerprint and to a scaffold splitter, so
  they can land on opposite sides of a split. `GetParent` runs before splitting, not after.
- **Scaffold splitting is not one algorithm.** Bemis–Murcko with and without generic atom
  types produce materially different difficulty. The variant is recorded in `findings.json`,
  because "scaffold split" alone is not a reproducible statement.

## Layout

```
src/chembench/
  curate.py     ChEMBL fetch, validity checking, standardisation, parent extraction
  splits.py     random / scaffold / activity-cliff splitting + a leakage report
  models.py     ECFP+SVM baseline, chemprop, foundation model — one interface
  evaluate.py   metrics, bootstrap CIs, and the cross-split comparison table
  cli.py        `python -m chembench.cli`
```

## Running it

```bash
make install     # pip install -e ".[dev]"
make data        # fetch and curate ChEMBL targets (cached)
make analysis    # all splits x all models -> results/
make test
```

RDKit and the model backends are optional extras, imported lazily, so the base install and
CI stay light.

## References

- LIT-PCBA (*J Chem Inf Model* 2020) and its 2025 leakage audit
- MoleculeACE — activity cliff benchmark, 30 ChEMBL targets
- An open source chemical structure curation pipeline using RDKit (*J Cheminform* 2020)
- Polaris — curated benchmarks for drug discovery ML
