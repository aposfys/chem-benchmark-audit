# chem-benchmark-audit
How much of molecular ML survives an honest split?

[![CI](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

5 ChEMBL targets, **20,255 curated compounds**, 3 model families, 3 split regimes — 45 cells, one curation, so the only things that vary between cells are the split and the model.

### The splits are what they claim to be

| Split | Scaffold leakage | Cliff enrichment |
| --- | ---: | ---: |
| Random | **67.1%** | 1.06× |
| Scaffold | **0.0%** | 0.99× |
| Activity cliff | 64.0% | **5.00×** |

**Two-thirds of a random split's test compounds share a Murcko scaffold with something the model trained on.** The scaffold split leaks exactly nothing, which a test asserts rather than assumes.

### The split costs about 0.13 RMSE, and it costs every model the same

| Model | Random | Scaffold | Activity cliff | Random → scaffold |
| --- | ---: | ---: | ---: | ---: |
| ECFP4 + SVM | 0.683 [0.632, 0.737] | 0.813 [0.767, 0.864] | 0.762 [0.706, 0.818] | **+0.130** |
| chemprop (D-MPNN) | 0.700 [0.653, 0.747] | 0.832 [0.786, 0.882] | 0.805 [0.748, 0.863] | **+0.132** |
| ChemBERTa + ridge | 0.941 [0.888, 0.995] | 1.058 [1.003, 1.119] | 1.013 [0.954, 1.073] | **+0.117** |

That is the robust finding: closing the leak costs roughly 0.13 pChEMBL units, and it costs the descriptor baseline, the graph network and the transformer almost identically. **The leak is not a deep-learning problem. It is an evaluation problem, and it inflates everyone equally.**

### The deep model's advantage is not there to survive

The question this repo was built to ask was how much of the deep model's advantage over a 2010-era baseline is left once the split stops leaking. The answer is that **there was no established advantage to begin with**.

Read through the bootstrap intervals, ECFP4 + SVM and chemprop are **indistinguishable in 14 of 15 target × split cells**. The single cell where the intervals separate is thrombin under a scaffold split — and it favours the SVM (0.750 [0.701, 0.801] against 0.872 [0.826, 0.918]).

**This comparison is biased toward the SVM and says so.** The baseline gets a small grid search over `C` and `gamma`, cross-validated inside the training fold; chemprop and the ChemBERTa head run at defaults. That is the cheapest advantage to give a baseline, and it is given deliberately — but it means the honest reading is "no difference was established here", not "the deep model cannot win". A tuned D-MPNN might.

ChemBERTa is the one clear loser, at 0.941 against 0.683 on random splits. It is used the way the "just use a foundation model" claim usually means it — frozen, mean-pooled, with a ridge head — and it was not fine-tuned. That is a statement about frozen embeddings, not about the checkpoint.

### Running it

```
make install
make data        # fetch and curate 5 ChEMBL targets (~9 min, cached)
make analysis    # 45 cells: 3 models x 3 splits x 5 targets (~3 h on CPU)
make test
```

`python -m chembench.cli evaluate --report-only` re-renders the tables from an existing `findings.json` without re-running anything.

Curation keeps 20,255 of 24,565 fetched measurements: records without a pChEMBL value, censored relations and ChEMBL's own validity-flagged rows are refused, and 4,266 duplicate measurements of the same parent structure are collapsed to their median. Every refusal is counted by reason and written out beside the data.

### Layout

```
src/chembench/
  curate.py     ChEMBL fetch, validity checking, standardisation, parent extraction, cliffs
  splits.py     random / scaffold / activity-cliff splitting + the leakage report
  models.py     ECFP+SVM, chemprop, ChemBERTa — one interface, lazy backends
  evaluate.py   RMSE, Spearman, bootstrap intervals, overlap test
  experiment.py the 45-cell grid
  report.py     the tables
  cli.py        `python -m chembench.cli`
```

24 tests, none needing a chemistry toolkit.

### More

- [Full results, per target and per cell](results/RESULTS.md)
- [Why this framing, and the traps the pipeline is built to avoid](docs/DESIGN.md)
