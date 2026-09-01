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

import json
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

#: ChEMBL stores tautomers as deposited; PubChem canonicalises them. Merging the two
#: without a single convention lets a model learn which database a record came from.
TAUTOMER_CONVENTION = "rdkit-v1"

#: The panel. Chosen for assay depth and for being standard in the activity-cliff
#: literature, so the numbers here can be argued with rather than taken on trust.
DEFAULT_TARGETS: tuple[tuple[str, str, str], ...] = (
    ("CHEMBL204", "Thrombin", "Ki"),
    ("CHEMBL214", "Serotonin 5-HT1a receptor", "Ki"),
    ("CHEMBL228", "Serotonin transporter", "Ki"),
    ("CHEMBL233", "mu opioid receptor", "Ki"),
    ("CHEMBL244", "Coagulation factor X", "Ki"),
)


@dataclass(frozen=True)
class Record:
    """One curated bioactivity measurement."""

    record_id: str
    smiles: str
    inchikey: str
    target_id: str
    activity: float
    activity_type: str


@dataclass
class CurationReport:
    """What curation kept, and what it refused.

    Refusals are counted by reason and written out with the data. A curated set whose
    losses are not stated cannot be argued with, and the losses here are large enough to
    change a conclusion.
    """

    target_id: str
    fetched: int = 0
    kept: int = 0
    rejected: dict[str, int] = None  # type: ignore[assignment]
    duplicates_collapsed: int = 0

    def __post_init__(self) -> None:
        if self.rejected is None:
            self.rejected = defaultdict(int)

    def reject(self, reason: str) -> None:
        self.rejected[reason] += 1

    def as_dict(self) -> dict:
        return {
            "target_id": self.target_id,
            "fetched": self.fetched,
            "kept": self.kept,
            "duplicates_collapsed": self.duplicates_collapsed,
            "rejected": dict(sorted(self.rejected.items())),
        }


# ---------------------------------------------------------------------------------------
# Retrieval
# ---------------------------------------------------------------------------------------


def _get(url: str, attempts: int = 5) -> dict:
    """One GET with exponential backoff.

    EBI throttles, and a run that dies halfway through a target leaves a cache file that
    looks complete. Retrying here is what makes the cache trustworthy.
    """
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            with urllib.request.urlopen(url, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"ChEMBL request failed after {attempts} attempts: {url}") from last


def fetch_target(
    target_id: str,
    out_path: Path,
    *,
    activity_type: str = "Ki",
    page_size: int = 1000,
) -> Path:
    """Download one target's bioactivity records to ``out_path``, cached.

    Assay heterogeneity is the first thing to filter on: mixing IC50 measured under
    different conditions produces a target that no model can fit and that looks like a
    hard benchmark. Only one ``standard_type`` is requested per target, and only records
    carrying a ``pchembl_value`` -- which ChEMBL computes only for an exact relation
    against a defined concentration -- are asked for.
    """
    if out_path.exists():
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)

    query = {
        "target_chembl_id": target_id,
        "standard_type": activity_type,
        "pchembl_value__isnull": "false",
        "limit": str(page_size),
    }
    url = f"{CHEMBL_API}/activity.json?{urllib.parse.urlencode(query)}"
    activities: list[dict] = []
    while url:
        payload = _get(url)
        activities.extend(payload["activities"])
        nxt = payload["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else ""

    out_path.write_text(json.dumps(activities, indent=1))
    return out_path


# ---------------------------------------------------------------------------------------
# Structure curation
# ---------------------------------------------------------------------------------------


def check(smiles: Iterable[str]) -> dict[str, list[str]]:
    """Return, per SMILES, the list of validity problems found. Empty list means clean."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    problems: dict[str, list[str]] = {}
    for text in smiles:
        found: list[str] = []
        mol = Chem.MolFromSmiles(text)
        if mol is None:
            found.append("unparseable")
        else:
            if mol.GetNumHeavyAtoms() == 0:
                found.append("no heavy atoms")
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    found.append("query or dummy atom")
                    break
            if "." in text and not found:
                found.append("multiple components")
        problems[text] = found
    return problems


def standardize(smiles: str) -> str:
    """Normalise charges, isotopes and tautomers to a single stated convention.

    The tautomer step is the one that matters for merging sources. It is also the slow
    one, which is why the enumerator is bounded -- an unbounded canonicaliser on a large
    flexible molecule can take minutes.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse {smiles!r}")

    mol = rdMolStandardize.Cleanup(mol)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)

    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetMaxTautomers(64)
    mol = enumerator.Canonicalize(mol)
    return Chem.MolToSmiles(mol)


def get_parent(smiles: str) -> str:
    """Strip salts and solvates, returning the parent structure.

    Runs before splitting, never after: two records of the same parent as different salts
    look like two molecules to a fingerprint and to a scaffold splitter, so they can land
    on opposite sides of a split and be scored as generalisation.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse {smiles!r}")
    parent = rdMolStandardize.FragmentParent(mol)
    return Chem.MolToSmiles(parent)


def inchikey(smiles: str) -> str:
    """InChIKey of a structure, used as the deduplication identity."""
    from rdkit import Chem, RDLogger

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse {smiles!r}")
    return Chem.MolToInchiKey(mol)


def murcko_scaffold(smiles: str, *, generic: bool = False) -> str:
    """Bemis-Murcko scaffold as canonical SMILES.

    ``generic=True`` erases atom types and bond orders, which merges scaffold groups and
    makes the resulting split materially harder. Whichever is used gets written into
    ``results/findings.json`` -- "scaffold split" alone is not a reproducible statement.
    """
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold

    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"cannot parse {smiles!r}")
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if generic:
        scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    return Chem.MolToSmiles(scaffold)


# ---------------------------------------------------------------------------------------
# The pipeline
# ---------------------------------------------------------------------------------------


def curate_target(
    target_id: str,
    raw_path: Path,
    *,
    activity_type: str = "Ki",
    max_records: int = 0,
) -> tuple[list[Record], CurationReport]:
    """Check, standardise, take parents, and collapse duplicates for one target.

    Duplicate measurements of the same parent structure are collapsed to their median
    pChEMBL rather than dropped or kept. Keeping them is leakage; dropping all but one is
    arbitrary; the median is the only choice that uses the replicates without letting a
    heavily-measured compound appear on both sides of a split.
    """
    activities = json.loads(raw_path.read_text())
    if max_records:
        activities = activities[:max_records]

    report = CurationReport(target_id=target_id, fetched=len(activities))
    by_structure: dict[str, list[tuple[str, float]]] = defaultdict(list)
    smiles_for_key: dict[str, str] = {}

    for activity in activities:
        smiles = activity.get("canonical_smiles")
        pchembl = activity.get("pchembl_value")
        if not smiles or pchembl is None:
            report.reject("missing structure or pChEMBL")
            continue
        # ChEMBL flags records its own curators distrust. Keeping them would put known-bad
        # measurements into a benchmark that exists to be trusted.
        if activity.get("data_validity_comment"):
            report.reject(f"data validity: {activity['data_validity_comment']}")
            continue
        if activity.get("standard_relation") != "=":
            report.reject("censored measurement (relation is not '=')")
            continue

        problems = check([smiles])[smiles]
        if any(
            problem in {"unparseable", "no heavy atoms", "query or dummy atom"}
            for problem in problems
        ):
            report.reject(problems[0])
            continue

        try:
            parent = get_parent(smiles)
            standardized = standardize(parent)
            key = inchikey(standardized)
        except Exception as exc:
            report.reject(f"curation failed: {type(exc).__name__}")
            continue

        smiles_for_key[key] = standardized
        by_structure[key].append((activity.get("molecule_chembl_id", key), float(pchembl)))

    records: list[Record] = []
    for key, measurements in by_structure.items():
        if len(measurements) > 1:
            report.duplicates_collapsed += len(measurements) - 1
        records.append(
            Record(
                record_id=measurements[0][0],
                smiles=smiles_for_key[key],
                inchikey=key,
                target_id=target_id,
                activity=median(value for _, value in measurements),
                activity_type=activity_type,
            )
        )

    report.kept = len(records)
    records.sort(key=lambda record: record.record_id)
    return records, report


def write_curated(records: list[Record], report: CurationReport, out_dir: Path) -> Path:
    """Write curated records and their curation report side by side."""
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{report.target_id}.json"
    path.write_text(
        json.dumps(
            {
                "report": report.as_dict(),
                "tautomer_convention": TAUTOMER_CONVENTION,
                "records": [asdict(record) for record in records],
            },
            indent=1,
        )
    )
    return path


def read_curated(path: Path) -> tuple[list[Record], dict]:
    """Read back what :func:`write_curated` produced."""
    payload = json.loads(path.read_text())
    records = [Record(**row) for row in payload["records"]]
    return records, payload["report"]


# ---------------------------------------------------------------------------------------
# Activity cliffs
# ---------------------------------------------------------------------------------------


def ecfp_matrix(smiles: list[str], *, n_bits: int = 2048, radius: int = 2):
    """ECFP fingerprints for a list of SMILES, as RDKit bit vectors."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    return [generator.GetFingerprint(Chem.MolFromSmiles(text)) for text in smiles]


def find_activity_cliffs(
    records: list[Record],
    *,
    similarity_threshold: float = 0.9,
    activity_threshold: float = 1.0,
) -> tuple[set[str], int]:
    """Identify compounds involved in an activity cliff.

    A cliff is a pair that is structurally near-identical (Tanimoto >= 0.9 on ECFP4) but
    differs by at least an order of magnitude in potency (>= 1 log unit of pChEMBL). Both
    members of such a pair are returned, along with the number of pairs found.

    Thresholds follow MoleculeACE. They are arguments rather than constants because the
    cliff definition changes the difficulty of the split, so a result is only meaningful
    alongside the thresholds that produced it -- which is why both are written into
    ``findings.json``.
    """
    from rdkit import DataStructs

    fingerprints = ecfp_matrix([record.smiles for record in records])
    members: set[str] = set()
    pairs = 0
    for i in range(len(records)):
        # BulkTanimoto over the upper triangle only: the relation is symmetric, and the
        # quadratic term is what dominates on a few thousand compounds.
        if i + 1 >= len(records):
            break
        similarities = DataStructs.BulkTanimotoSimilarity(
            fingerprints[i], fingerprints[i + 1 :]
        )
        for offset, similarity in enumerate(similarities):
            j = i + 1 + offset
            if similarity < similarity_threshold:
                continue
            if abs(records[i].activity - records[j].activity) >= activity_threshold:
                members.add(records[i].record_id)
                members.add(records[j].record_id)
                pairs += 1
    return members, pairs
