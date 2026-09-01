# chem-benchmark-audit — design notes

Written before the analysis was run, so the evaluation cannot be adjusted to fit the
results.

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

## References

- LIT-PCBA (*J Chem Inf Model* 2020) and its 2025 leakage audit
- MoleculeACE — activity cliff benchmark, 30 ChEMBL targets
- An open source chemical structure curation pipeline using RDKit (*J Cheminform* 2020)
- Polaris — curated benchmarks for drug discovery ML

## Layout

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

27 tests, none needing a chemistry toolkit.

Curation keeps 20,255 of 24,565 fetched measurements: records without a pChEMBL
value, censored relations and ChEMBL's own validity-flagged rows are refused, and
4,266 duplicate measurements of the same parent structure are collapsed to their
median. Every refusal is counted by reason and written out beside the data.
