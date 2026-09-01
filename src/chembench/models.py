"""Three model families behind one interface, so the split is the only thing that varies.

The comparison is only meaningful if every model sees byte-identical train and test sets,
so models never touch curation or splitting -- they receive SMILES and return predictions.

``ecfp_svm``
    ECFP4 + support vector regression. The 2010-era baseline that keeps refusing to lose.
``chemprop``
    Directed message-passing neural network over the molecular graph.
``foundation``
    ChemBERTa-77M, a transformer pretrained on 77M PubChem molecules, frozen, with a
    ridge head on its mean-pooled embeddings.

**The baseline is tuned and the deep models are not.** That asymmetry is deliberate and it
runs *against* the repo's own thesis: a small grid search over ``C`` and ``gamma`` is what
an SVM needs to be competitive, while chemprop and the foundation head run at defaults. If
the classical baseline still wins under an honest split, it wins having been given the
advantage that is cheapest to give. If it loses, the loss is real. Either way the direction
of the bias is stated rather than left for a reader to guess.

Backends are imported lazily inside each class, so the splitting and metric code -- and CI
-- run without RDKit, torch or chemprop installed.
"""

from __future__ import annotations

import os
import warnings
from collections.abc import Sequence
from typing import Any, Protocol

MODEL_NAMES = ("ecfp_svm", "chemprop", "foundation")

#: Pretrained checkpoint for the ``foundation`` family. Pinned, because "a foundation
#: model" is not a reproducible statement.
CHEMBERTA_CHECKPOINT = "DeepChem/ChemBERTa-77M-MLM"


class Model(Protocol):
    """Minimal interface every model family implements."""

    name: str

    def fit(self, smiles: Sequence[str], y: Sequence[float]) -> None: ...

    def predict(self, smiles: Sequence[str]) -> list[float]: ...


# ---------------------------------------------------------------------------------------


class EcfpSvm:
    """ECFP4 fingerprints into an RBF support vector regressor.

    The grid is small on purpose: it is what a practitioner would actually try, not an
    exhaustive search. Tuning happens by cross-validation *within the training fold*, so no
    test compound influences hyperparameter choice -- which is the usual quiet way a
    baseline gets an unfair advantage.
    """

    name = "ecfp_svm"

    def __init__(self, n_bits: int = 2048, radius: int = 2, seed: int = 0, tune: bool = True):
        self.n_bits = n_bits
        self.radius = radius
        self.seed = seed
        self.tune = tune
        self._model: Any = None

    def _featurise(self, smiles: Sequence[str]):
        import numpy as np
        from rdkit import Chem, DataStructs, RDLogger
        from rdkit.Chem import rdFingerprintGenerator

        RDLogger.DisableLog("rdApp.*")
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=self.radius, fpSize=self.n_bits
        )
        matrix = np.zeros((len(smiles), self.n_bits), dtype=np.uint8)
        for row, text in enumerate(smiles):
            mol = Chem.MolFromSmiles(text)
            if mol is None:
                raise ValueError(f"cannot parse {text!r}")
            arr = np.zeros((self.n_bits,), dtype=np.uint8)
            DataStructs.ConvertToNumpyArray(generator.GetFingerprint(mol), arr)
            matrix[row] = arr
        return matrix

    def fit(self, smiles: Sequence[str], y: Sequence[float]) -> None:
        from sklearn.model_selection import GridSearchCV
        from sklearn.svm import SVR

        features = self._featurise(smiles)
        if self.tune and len(smiles) >= 50:
            search = GridSearchCV(
                SVR(kernel="rbf"),
                {"C": [1.0, 10.0, 100.0], "gamma": ["scale", 0.001, 0.01]},
                cv=3,
                scoring="neg_root_mean_squared_error",
                n_jobs=-1,
            )
            search.fit(features, list(y))
            self._model = search.best_estimator_
        else:
            self._model = SVR(kernel="rbf", C=10.0).fit(features, list(y))

    def predict(self, smiles: Sequence[str]) -> list[float]:
        if self._model is None:
            raise RuntimeError("fit() before predict()")
        return [float(value) for value in self._model.predict(self._featurise(smiles))]


# ---------------------------------------------------------------------------------------


class Chemprop:
    """Directed message-passing neural network (chemprop v2), on CPU."""

    name = "chemprop"

    def __init__(self, epochs: int = 40, batch_size: int = 64, seed: int = 0):
        self.epochs = epochs
        self.batch_size = batch_size
        self.seed = seed
        self._model: Any = None
        self._trainer: Any = None
        self._featurizer: Any = None

    def _datapoints(self, smiles: Sequence[str], y: Sequence[float] | None):
        import numpy as np
        from chemprop import data

        if y is None:
            return [data.MoleculeDatapoint.from_smi(text, None) for text in smiles]
        return [
            data.MoleculeDatapoint.from_smi(text, np.array([value]))
            for text, value in zip(smiles, y, strict=True)
        ]

    def fit(self, smiles: Sequence[str], y: Sequence[float]) -> None:
        import lightning.pytorch as pl
        import torch
        from chemprop import data, featurizers, models, nn

        warnings.filterwarnings("ignore")
        torch.manual_seed(self.seed)

        self._featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        train = data.MoleculeDataset(self._datapoints(smiles, y), self._featurizer)
        # Targets are scaled on the training fold only, and the inverse transform is
        # attached to the model, so predictions come back in pChEMBL units without the
        # test fold ever contributing to the scaler.
        scaler = train.normalize_targets()
        loader = data.build_dataloader(
            train, batch_size=self.batch_size, num_workers=0, seed=self.seed
        )

        ffn = nn.RegressionFFN(
            output_transform=nn.UnscaleTransform.from_standard_scaler(scaler)
        )
        self._model = models.MPNN(
            nn.BondMessagePassing(),
            nn.MeanAggregation(),
            ffn,
            batch_norm=True,
            metrics=[nn.metrics.RMSE()],
        )
        self._trainer = pl.Trainer(
            max_epochs=self.epochs,
            accelerator="cpu",
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=False,
            enable_model_summary=False,
            deterministic=True,
        )
        self._trainer.fit(self._model, loader)

    def predict(self, smiles: Sequence[str]) -> list[float]:
        import numpy as np
        from chemprop import data

        if self._model is None or self._trainer is None:
            raise RuntimeError("fit() before predict()")
        dataset = data.MoleculeDataset(self._datapoints(smiles, None), self._featurizer)
        loader = data.build_dataloader(
            dataset, batch_size=self.batch_size, num_workers=0, shuffle=False
        )
        batches = self._trainer.predict(self._model, loader)
        return [float(v) for v in np.concatenate([b.numpy() for b in batches]).ravel()]


# ---------------------------------------------------------------------------------------


class Foundation:
    """Frozen ChemBERTa embeddings with a ridge head.

    The transformer is not fine-tuned. That is the cheap, standard way a pretrained
    molecular representation actually gets used, and it is the configuration the "just use
    a foundation model" claim usually refers to. Fine-tuning would be a different and more
    expensive experiment; that it was not done is recorded rather than glossed.
    """

    name = "foundation"

    def __init__(
        self, checkpoint: str = CHEMBERTA_CHECKPOINT, seed: int = 0, batch_size: int = 64
    ):
        self.checkpoint = checkpoint
        self.seed = seed
        self.batch_size = batch_size
        # Typed Any: transformers' Auto* factories return a union that varies with the
        # installed extras, so pinning a class here would be a fiction.
        self._head: Any = None
        self._tokenizer: Any = None
        self._encoder: Any = None

    def _load(self) -> None:
        if self._encoder is not None:
            return
        os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")
        from transformers import AutoModel, AutoTokenizer

        warnings.filterwarnings("ignore")
        self._tokenizer = AutoTokenizer.from_pretrained(self.checkpoint)
        self._encoder = AutoModel.from_pretrained(self.checkpoint)
        self._encoder.eval()

    def _embed(self, smiles: Sequence[str]):
        import numpy as np
        import torch

        self._load()
        vectors = []
        for start in range(0, len(smiles), self.batch_size):
            batch = list(smiles[start : start + self.batch_size])
            encoded = self._tokenizer(
                batch, return_tensors="pt", padding=True, truncation=True, max_length=512
            )
            with torch.no_grad():
                hidden = self._encoder(**encoded).last_hidden_state
            # Mean-pool over real tokens only; padding must not dilute the embedding.
            mask = encoded["attention_mask"].unsqueeze(-1).float()
            pooled = (hidden * mask).sum(dim=1) / mask.sum(dim=1).clamp(min=1e-9)
            vectors.append(pooled.numpy())
        return np.vstack(vectors)

    def fit(self, smiles: Sequence[str], y: Sequence[float]) -> None:
        from sklearn.linear_model import RidgeCV

        features = self._embed(smiles)
        self._head = RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0]).fit(features, list(y))

    def predict(self, smiles: Sequence[str]) -> list[float]:
        if self._head is None:
            raise RuntimeError("fit() before predict()")
        return [float(value) for value in self._head.predict(self._embed(smiles))]


# ---------------------------------------------------------------------------------------


def build(name: str, **params: object) -> Model:
    """Construct a model by name. Backends are imported lazily inside each branch."""
    if name not in MODEL_NAMES:
        raise ValueError(f"unknown model {name!r}; expected one of {MODEL_NAMES}")
    if name == "ecfp_svm":
        return EcfpSvm(**params)  # type: ignore[arg-type]
    if name == "chemprop":
        return Chemprop(**params)  # type: ignore[arg-type]
    return Foundation(**params)  # type: ignore[arg-type]


def available(name: str) -> tuple[bool, str]:
    """Whether a model's backend is importable, and why not if it is missing.

    Used so a run reports "chemprop was not installed" rather than silently comparing two
    models and calling it three.
    """
    try:
        if name == "ecfp_svm":
            import rdkit  # noqa: F401
            import sklearn  # noqa: F401
        elif name == "chemprop":
            import chemprop  # noqa: F401
            import lightning  # noqa: F401
        elif name == "foundation":
            import torch  # noqa: F401
            import transformers  # noqa: F401
        else:
            return False, f"unknown model {name!r}"
    except ImportError as exc:
        return False, str(exc)
    return True, ""
