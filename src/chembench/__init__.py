"""Leakage-aware evaluation of molecular property prediction on curated ChEMBL targets.

The package is deliberately split so that curation, splitting and modelling can be
audited independently: a result is only meaningful if the curation that produced it and
the split that scored it are both stated.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
