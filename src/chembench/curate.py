"""ChEMBL retrieval and structure curation.

Follows the three-component pipeline of Bento et al. (*J Cheminform* 2020):

``check``
    Flag structures with serious errors before they are trusted.
``standardize``
    Apply one set of formatting conventions to every record.
``get_parent``
    Strip salts and solvents so that a compound is one compound.

Order matters and is the point. ``get_parent`` must run *before* any splitting, or the
same parent compound recorded as two different salts can be split across train and test.

RDKit is an optional extra and is imported lazily inside the functions that need it, so
that the splitting and evaluation code -- and CI -- run without a chemistry toolkit.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

#: ChEMBL stores tautomers as deposited; PubChem canonicalises them. Merging the two
#: without a single convention lets a model learn which database a record came from.
TAUTOMER_CONVENTION = "rdkit-v1"


@dataclass(frozen=True)
class Record:
    """One curated bioactivity measurement."""

    record_id: str
    smiles: str
    inchikey: str
    target_id: str
    activity: float
    activity_type: str


def fetch_target(target_id: str, out_path: Path, *, activity_type: str = "pChEMBL") -> Path:
    """Download one target's bioactivity records to ``out_path``, cached.

    Assay heterogeneity is the first thing to filter on: mixing IC50 measured under
    different conditions produces a target that no model can fit and that looks like a
    hard benchmark.
    """
    raise NotImplementedError("milestone 1: ChEMBL retrieval")


def check(smiles: Iterable[str]) -> dict[str, list[str]]:
    """Return, per SMILES, the list of validity problems found. Empty list means clean."""
    raise NotImplementedError("milestone 1: structure checker")


def standardize(smiles: str) -> str:
    """Normalise charges, isotopes and tautomers to a single stated convention."""
    raise NotImplementedError("milestone 1: standardizer")


def get_parent(smiles: str) -> str:
    """Strip salts and solvates, returning the parent structure."""
    raise NotImplementedError("milestone 1: parent extraction")


def murcko_scaffold(smiles: str, *, generic: bool = False) -> str:
    """Bemis-Murcko scaffold as canonical SMILES.

    ``generic=True`` erases atom types and bond orders, which merges scaffold groups and
    makes the resulting split materially harder. Whichever is used gets written into
    ``results/findings.json`` -- "scaffold split" alone is not a reproducible statement.
    """
    raise NotImplementedError("milestone 1: scaffold extraction")
