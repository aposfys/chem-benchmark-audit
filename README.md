# chem-benchmark-audit
How much of molecular ML survives an honest split?

[![CI](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** Structure, interfaces and tests are in place; the analysis has not been run. Nothing in this repo is a result yet.

Molecular property prediction is reported almost everywhere on random splits, which put near-identical molecules in both train and test. This repo feeds the same curated ChEMBL targets to the same models under three splitting regimes — random, Bemis–Murcko scaffold, and activity-cliff — so the size of that effect is directly readable.

### Running it
```
make install     # pip install -e ".[dev]"
make data        # fetch and curate ChEMBL targets (cached)
make analysis    # all splits x all models -> results/
make test
```
RDKit and the model backends are optional extras, imported lazily, so the base install and CI stay light.

### Layout
```
src/chembench/
  curate.py     ChEMBL fetch, validity checking, standardisation, parent extraction
  splits.py     random / scaffold / activity-cliff splitting + a leakage report
  models.py     ECFP+SVM baseline, chemprop, foundation model — one interface
  evaluate.py   metrics, bootstrap CIs, and the cross-split comparison table
  cli.py        `python -m chembench.cli`
```

### Design notes
[Why this framing, the traps the pipeline is built to avoid, and the references](docs/DESIGN.md)
