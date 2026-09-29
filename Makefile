.PHONY: install data analysis quick test clean clean-data all

PYTHON ?= python3

all: analysis

## Install the package with RDKit, the three model backends and dev tooling.
## CI installs only ".[dev]", which is enough for the tests.
install:
	$(PYTHON) -m pip install -e ".[dev,chem,models]"

## Fetch the selected ChEMBL targets and run the curation pipeline
## (Checker -> Standardizer -> GetParent). Cached in data/, so later runs skip it.
data:
	$(PYTHON) -m chembench.cli curate

## Every split regime against every available model, plus the comparison table
analysis: data
	$(PYTHON) -m chembench.cli evaluate

## One target, ECFP+SVM only, enough to see the split effect
quick: data
	$(PYTHON) -m chembench.cli evaluate --models ecfp_svm --targets 1

test:
	$(PYTHON) -m pytest -q

## Remove build and run debris. The committed results/*.json and RESULTS.md stay.
clean:
	rm -rf results/models results/_smoke
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete the cached ChEMBL downloads and the curated sets
clean-data: clean
	rm -f data/*.csv data/*.parquet
