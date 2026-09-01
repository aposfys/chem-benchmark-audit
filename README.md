# chem-benchmark-audit
How much of molecular ML survives an honest split?

[![CI](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

5 ChEMBL targets, 20,255 curated compounds, 3 model families, 3 split regimes —
45 cells, one curation, so the only things varying between cells are the split
and the model.

```
make install
make data        # fetch and curate 5 ChEMBL targets (~9 min, cached)
make analysis    # 45 cells: 3 models x 3 splits x 5 targets (~3 h on CPU)
make test        # 27 tests, none needing a chemistry toolkit
```

`python -m chembench.cli evaluate --report-only` re-renders the tables from an
existing `findings.json` without re-running anything.

### The splits are what they claim to be

| Split | Scaffold leakage | Cliff enrichment |
| --- | ---: | ---: |
| Random | **67.1%** | 1.06× |
| Scaffold | **0.0%** | 0.99× |
| Activity cliff | 64.0% | **5.00×** |

Two-thirds of a random split's test compounds share a Murcko scaffold with
something the model trained on. The scaffold split leaks exactly nothing, which a
test asserts rather than assumes.

### The split costs ~0.13 RMSE, and it costs every model the same

| Model | Random | Scaffold | Random → scaffold |
| --- | ---: | ---: | ---: |
| ECFP4 + SVM | 0.683 [0.632, 0.737] | 0.813 [0.767, 0.864] | **+0.130** |
| chemprop (D-MPNN) | 0.700 [0.653, 0.747] | 0.832 [0.786, 0.882] | **+0.132** |
| ChemBERTa + ridge | 0.941 [0.888, 0.995] | 1.058 [1.003, 1.119] | **+0.117** |

Closing the leak costs the descriptor baseline, the graph network and the
transformer almost identically. **The leak is not a deep-learning problem. It is
an evaluation problem, and it inflates everyone equally.**

### The deep model's advantage is not there to survive

ECFP4 + SVM and chemprop are **indistinguishable in 14 of 15 target × split
cells**. The one cell where the intervals separate favours the SVM.

**This comparison is biased toward the SVM and says so.** The baseline gets a
small grid search cross-validated inside the training fold; chemprop and the
ChemBERTa head run at defaults. That is the cheapest advantage to give a
baseline, and it is given deliberately — but the honest reading is "no difference
was established here", not "the deep model cannot win". ChemBERTa is used frozen
and mean-pooled with a ridge head, so its 0.941 is a statement about frozen
embeddings, not about the checkpoint.

### More

- [Analysis](ANALYSIS.md) — what was done and why, including curation refusals
- [Results](results/RESULTS.md) — full results, per target and per cell
- [Design](docs/DESIGN.md) — why this framing, the layout, and the traps it avoids
