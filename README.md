# chem-benchmark-audit
How much of molecular ML survives an honest split?

[![CI](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/chem-benchmark-audit/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

5 ChEMBL Ki targets, 20,255 curated compounds, 3 models and 3 splits, 45 cells from one
curation, so the only things that vary between cells are the split and the model.

```
make install     # package, RDKit, chemprop, transformers and dev tools
make data        # fetch and curate 5 ChEMBL targets (~9 min, cached)
make analysis    # 45 cells: 3 models x 3 splits x 5 targets (~3 h on CPU)
make test        # unit tests, none needing a chemistry toolkit
chembench evaluate --report-only   # re-render RESULTS.md from findings.json
```

### The splits do what they say

| Split | Scaffold leakage | Cliff enrichment |
| --- | ---: | ---: |
| Random | **67.1%** | 1.06× |
| Scaffold | **0.0%** | 0.99× |
| Activity cliff | 64.0% | **5.00×** |

Two thirds of a random split's test compounds share a Murcko scaffold with a training
compound. The scaffold split leaks nothing, and a test asserts it. The 5.00× is the ceiling
of 1/test_frac, reached because every cliff compound fits in the test set.

### The scaffold split costs every model about 0.13 RMSE

| Model | Random | Scaffold | Random → scaffold |
| --- | ---: | ---: | ---: |
| ECFP4 + SVM | 0.683 [0.632, 0.737] | 0.813 [0.767, 0.864] | **+0.130** |
| chemprop (D-MPNN) | 0.700 [0.653, 0.747] | 0.832 [0.786, 0.882] | **+0.132** |
| ChemBERTa + ridge | 0.941 [0.888, 0.995] | 1.058 [1.003, 1.119] | **+0.117** |

pChEMBL units, averaged over targets. Brackets are the mean of the per-target 95% bootstrap
intervals, not an interval on the mean. The descriptor baseline, the graph network and the
frozen transformer lose almost the same amount, so the leak inflates all three alike.

No deep-model advantage was established. ECFP4 + SVM and chemprop are not separated by
unpaired interval overlap in 14 of 15 target × split cells, and the one separated cell
favours the SVM. The comparison leans towards the SVM, which gets a small grid search inside
the training fold, while chemprop runs 40 epochs untuned and ChemBERTa is frozen with a
RidgeCV head.

### More

- [Analysis](ANALYSIS.md) covers what was done and why, what the numbers do not support,
  and how this sits against published work
- [Results](results/RESULTS.md) has every target and cell, with the overlap counts
- [Curation report](results/curation_report.json) has per-target counts and the ChEMBL release
- [Design](docs/DESIGN.md) has the framing, the layout and the traps it avoids
