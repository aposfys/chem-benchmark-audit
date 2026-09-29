# Analysis

What was built, why it was built that way, and what the numbers do and do not support.

## The question, and how it changed

The repository was set up to ask: *when the split stops leaking, how much of the deep
model's advantage over a 2010-era baseline is left?*

That question presupposes an advantage. The run does not find one. Read through bootstrap
intervals, ECFP4 + SVM and chemprop are **not separated in 14 of 15** target × split cells,
and the single separated cell (CHEMBL204, scaffold) favours the SVM. So the honest restatement is: *there
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

**The baseline is tuned and the deep models are barely tuned.** This is deliberate and it
runs against the repository's own thesis. A small grid over `C` and `gamma`, cross-validated
*inside the training fold*, is what an SVM needs to be competitive. chemprop runs a fixed 40
epochs with no architecture or learning-rate search, and the ChemBERTa head is a RidgeCV over
four alphas on frozen embeddings.

The direction of that bias matters for how the result reads. It favours the SVM, and the
SVM still does not clearly win. That makes "no difference was established" robust to the
bias, and it is why the conclusion is phrased that way rather than as "the baseline wins".
It does **not** rule out that a tuned D-MPNN would win. That experiment was not run.

**Every comparison is read through intervals, including the inconvenient ones.**
`report.py` calls `intervals_overlap` and prints "not established" when they overlap, and it
prints the per-cell overlap count so the "14 of 15" is rendered rather than asserted. The
repository's whole argument is that small differences on a leaky split are noise, and that
argument applies to its own results.

**What the two non-random splits actually do.** The scaffold split fills the test set with
the largest scaffold groups, so the test set is a handful of large congeneric series and
every singleton scaffold is in training. That is the reverse of DeepChem's convention and it
makes the split hard, but it also means a bootstrap over test compounds resamples inside a
dozen or so series and understates the uncertainty. The cliff split fills the test set with
cliff compounds first, and because cliff members are 3.8% to 7.8% of each target against a
20% test fraction, all of them land in test. Both members of every cliff pair are therefore
in test, no pair straddles the split, and the 5.00x enrichment is the 1/test_frac ceiling
rather than a measured degree of enrichment. Both behaviours are stated in the docstrings in
`splits.py` and pinned by tests.

**The scaffold variant is recorded.** Bemis–Murcko with and without generic atom typing
produce materially different difficulty. "Scaffold split" alone is not a reproducible
statement, so the variant goes into `findings.json`.

## What the numbers mean

Random splits leak **67.1%** of test scaffolds into training. Closing that leak costs
**+0.130** RMSE for the SVM, **+0.132** for chemprop and **+0.117** for ChemBERTa. The near
identity of those three is the finding: the leak is not a deep-learning artefact, it is an
evaluation artefact, and it inflates every model family by about the same amount.

ChemBERTa is the one clear loser (0.941 against 0.683 on random splits). It is used frozen,
mean-pooled, with a ridge head, the way "just use a foundation model" usually means. That
is a statement about frozen embeddings, not about the checkpoint, and fine-tuning was not
attempted.

## Where this sits against published work

Both findings above are consistent with published work rather than novel, and the
scaffold-split question is largely settled in the literature.

- Fooladi et al., *JCIM* 2025, [10.1021/acs.jcim.5c00475](https://doi.org/10.1021/acs.jcim.5c00475),
  14 models over 8 datasets and 10 splitting strategies. They report that classical ML and
  graph neural networks perform not substantially different from random splitting under
  Bemis-Murcko scaffold splits, and that chemical-similarity clustering (UMAP over ECFP4) is
  the hard split for both families.
- Guo et al. 2024, *Scaffold Splits Overestimate Virtual Screening Performance*. Molecules
  with different scaffolds are often similar, so a scaffold split still leaves
  unrealistically high train-test similarity.

The classical-versus-GNN result here agrees with Fooladi et al. The scaffold penalty does
not. They find Bemis-Murcko splits about as easy as random, while here the scaffold split
costs about 0.13 RMSE for every model and the averaged intervals do not overlap. The likely
reason is the splitter, not the target panel. This scaffold split sends the largest scaffold
series to test and keeps every singleton in training, which is harder than the usual
fill-training-first convention, and Ki regression on five deeply measured targets is a
narrower task than the eight mixed bioactivity and ADMET datasets they use.

Neither paper's harder split is implemented here. Similarity clustering over ECFP4 is the
obvious next measurement, because 0.0% scaffold leakage is a statement about Murcko
scaffolds and not about chemical similarity. Read this repository as a controlled and
reproducible demonstration of a known effect on five ChEMBL targets, with the leakage
assertion under test rather than assumed.

## What is not established

- That deep models cannot beat the baseline here. Only that they did not, barely tuned.
- Anything a paired test would settle. The comparison is an unpaired overlap of two
  percentile intervals on the same test compounds, at one split seed, and per-compound
  predictions are not committed, so it cannot be re-tested without rerunning the grid.
- Anything about targets outside this panel of five.
- Anything about classification. Every metric here is regression on pChEMBL.

## What would change the conclusion

A tuned chemprop. If a hyperparameter search inside the training fold moved chemprop's
intervals clear of the SVM's, the original framing would be back and the split effect would
still stand. A paired, scaffold-clustered bootstrap would also be a sharper test than the
unpaired overlap used here. That is the obvious next run and it is not cheap: chemprop is roughly four
minutes per fit on CPU, so the 45-cell grid takes about three hours before any tuning.
