.PHONY: install data analysis quick test clean clean-data all

PYTHON ?= python3

all: analysis

## Install the package plus dev tooling. RDKit and the model backends are
## optional extras:  pip install -e ".[chem,models]"
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Fetch the selected ChEMBL targets and run the curation pipeline
## (Checker -> Standardizer -> GetParent). Cached in data/, so later runs skip it.
data:
	$(PYTHON) -m chembench.cli curate

## Every split regime against every available model, plus the comparison table
analysis: data
	$(PYTHON) -m chembench.cli evaluate

## One target, ECFP+SVM only — enough to see the split effect
quick: data
	$(PYTHON) -m chembench.cli evaluate --models ecfp_svm --targets 1

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete the cached ChEMBL downloads and the curated sets
clean-data: clean
	rm -f data/*.csv data/*.parquet
