"""Three model families behind one interface, so the split is the only thing that varies.

The comparison is only meaningful if every model sees byte-identical train and test sets,
so models never touch curation or splitting -- they receive keys and return predictions.

``ecfp_svm``
    ECFP4 + support vector regression. The 2010-era baseline that keeps refusing to lose.
``chemprop``
    Directed message-passing neural network over the molecular graph.
``foundation``
    A pretrained molecular representation with a light head on top.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Protocol

MODEL_NAMES = ("ecfp_svm", "chemprop", "foundation")


class Model(Protocol):
    """Minimal interface every model family implements."""

    name: str

    def fit(self, smiles: Sequence[str], y: Sequence[float]) -> None: ...

    def predict(self, smiles: Sequence[str]) -> list[float]: ...


def build(name: str, **params: object) -> Model:
    """Construct a model by name. Backends are imported lazily inside each branch."""
    if name not in MODEL_NAMES:
        raise ValueError(f"unknown model {name!r}; expected one of {MODEL_NAMES}")
    raise NotImplementedError("milestone 2: model backends")
