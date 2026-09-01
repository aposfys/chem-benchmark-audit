# Analysis

What was built, why it was built that way, and what the numbers do and do not support.

## The question, and how it changed

The repository was set up to ask: *when the split stops leaking, how much of the deep
model's advantage over a 2010-era baseline is left?*

That question presupposes an advantage. The run does not find one. Read through bootstrap
intervals, ECFP4 + SVM and chemprop are indistinguishable in **14 of 15** target × split
cells, and the single separated cell favours the SVM. So the honest restatement is: *there
was no established advantage to erode, and the thing that actually moves the numbers is the
split.*

Leaving the original framing in place and reporting a null against it would have been the
easy write-up. The framing is restated instead, because a question whose premise failed is
a result, not a disappointment.

## Design decisions, and the reasoning

**One curation, computed once per target, before any splitting.** Scaffolds and activity
cliffs are computed in `prepare_target` and shared across every split and model. If a number
moves between two cells, there is exactly one thing that could have caused it. Recomputing
per cell would have made a curation difference indistinguishable from a model difference.

**`GetParent` runs before splitting, never after.** Two records of the same parent compound
as different salts look like two molecules to a fingerprint and to a scaffold splitter, so
they can land on opposite sides of a split and be scored as generalisation. Duplicate
measurements are collapsed to their **median** rather than dropped or kept: keeping them is
leakage, dropping all but one is arbitrary, and the median uses the replicates without
letting a heavily-measured compound appear on both sides.

**The baseline is tuned and the deep models are not.** This is deliberate and it runs
against the repository's own thesis. A small grid over `C` and `gamma`, cross-validated
*inside the training fold*, is what an SVM needs to be competitive. chemprop and the
ChemBERTa head run at defaults.

The direction of that bias matters for how the result reads. It favours the SVM — and the
SVM still does not clearly win. That makes "no difference was established" robust to the
bias, and it is why the conclusion is phrased that way rather than as "the baseline wins".
It does **not** rule out that a tuned D-MPNN would win; that experiment was not run.

**Every comparison is read through intervals, including the inconvenient ones.**
`report.py` calls `intervals_overlap` and prints "not established" when they overlap. The
repository's whole argument is that small differences on a leaky split are noise, and that
argument applies to its own results.

**The scaffold variant is recorded.** Bemis–Murcko with and without generic atom typing
produce materially different difficulty. "Scaffold split" alone is not a reproducible
statement, so the variant goes into `findings.json`.

## What the numbers mean

Random splits leak **67.1%** of test scaffolds into training. Closing that leak costs
**+0.130** RMSE for the SVM, **+0.132** for chemprop and **+0.117** for ChemBERTa. The near
identity of those three is the finding: the leak is not a deep-learning artefact, it is an
evaluation artefact, and it inflates every model family by about the same amount.

ChemBERTa is the one clear loser (0.941 against 0.683 on random splits). It is used frozen,
mean-pooled, with a ridge head — the way "just use a foundation model" usually means. That
is a statement about frozen embeddings, not about the checkpoint, and fine-tuning was not
attempted.

## What is not established

- That deep models cannot beat the baseline here. Only that they did not, untuned.
- Anything about targets outside this panel of five.
- Anything about classification; every metric here is regression on pChEMBL.

## What would change the conclusion

A tuned chemprop. If a hyperparameter search inside the training fold moved chemprop's
intervals clear of the SVM's, the original framing would be back and the split effect would
still stand. That is the obvious next run and it is not cheap: chemprop is roughly four
minutes per fit on CPU, so the 45-cell grid takes about three hours before any tuning.
