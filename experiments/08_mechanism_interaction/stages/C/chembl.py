"""ChEMBL `max_phase_for_ind` cross-check (PREREG §Measured variables; cross-check only).

Reads the ChEMBL 37 SQLite database. Per hypothesis: the maximum `max_phase_for_ind` over
`drug_indication` rows whose molecule is a hypothesis program or a child of one (ChEMBL
`molecule_hierarchy.parent_molecule_chembl_id`) and whose `efo_id` equals the indication id
(ChEMBL writes `EFO:0000270`; Open Targets writes `EFO_0000270`). If the dump is unavailable
the cross-check is reported as not run and no label changes (PREREG pre-stage amendments).

The database lives at `/chembl_37/chembl_37.db` on the Modal volume
`proteome-mr-claim-audit-inputs`, mounted read-only; its sha256 is pinned in
`/chembl_37/VERIFIED.json` beside it. `verified_chembl_db` checks that pin before the database
is opened, and raises if VERIFIED.json is missing, names no sha256 (or more than one) for the
database, or disagrees with the file.
"""
import hashlib
import json
import sqlite3
from pathlib import Path

from models import StageCError

DB_NAME = "chembl_37.db"
VERIFIED_NAME = "VERIFIED.json"
NAME_FIELDS = ("path", "name", "file", "filename")

QUERY = """
SELECT md.chembl_id AS molecule, pmd.chembl_id AS parent, di.efo_id, di.max_phase_for_ind
FROM drug_indication di
JOIN molecule_dictionary md ON md.molregno = di.molregno
LEFT JOIN molecule_hierarchy mh ON mh.molregno = di.molregno
LEFT JOIN molecule_dictionary pmd ON pmd.molregno = mh.parent_molregno
WHERE di.efo_id IS NOT NULL
"""


class ChemblNotVerified(StageCError):
    """The ChEMBL database is missing, or its sha256 is not the one VERIFIED.json pins."""


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pinned_sha256(verified: object, name: str = DB_NAME) -> set[str]:
    """Every sha256 that VERIFIED.json attaches to `name`: a key equal to `name` whose value is
    a hex string or an object with `sha256`, or an object with a path/name field ending in
    `name` and a `sha256` field, anywhere in the document."""
    found: set[str] = set()
    if isinstance(verified, dict):
        for k, v in verified.items():
            if k == name and isinstance(v, str):
                found.add(v.lower())
            if k == name and isinstance(v, dict) and isinstance(v.get("sha256"), str):
                found.add(v["sha256"].lower())
        named = any(isinstance(verified.get(f), str) and verified[f].endswith(name) for f in NAME_FIELDS)
        if named and isinstance(verified.get("sha256"), str):
            found.add(verified["sha256"].lower())
        for v in verified.values():
            found |= pinned_sha256(v, name)
    elif isinstance(verified, list):
        for v in verified:
            found |= pinned_sha256(v, name)
    return found


def verified_chembl_db(chembl_dir: Path) -> tuple[Path, str]:
    """(database path, sha256) after checking the database against VERIFIED.json."""
    verified, db = chembl_dir / VERIFIED_NAME, chembl_dir / DB_NAME
    if not verified.is_file():
        raise ChemblNotVerified(f"{verified} is missing; the ChEMBL database is not opened")
    if not db.is_file():
        raise ChemblNotVerified(f"{db} is missing")
    pins = pinned_sha256(json.loads(verified.read_text()))
    if len(pins) != 1:
        raise ChemblNotVerified(f"{verified} pins {len(pins)} sha256 values for {DB_NAME}: {sorted(pins)}")
    got = sha256_file(db)
    if got != next(iter(pins)):
        raise ChemblNotVerified(f"{db}: sha256 {got} differs from {next(iter(pins))} in {verified}")
    return db, got


def norm_ontology_id(x: str) -> str:
    return x.replace(":", "_")


def load_chembl_index(db: Path) -> dict[tuple[str, str], float]:
    """(program chembl id, indication id) -> max max_phase_for_ind. Opened read-only and
    immutable, so a read-only volume mount needs no lock or journal file."""
    out: dict[tuple[str, str], float] = {}
    with sqlite3.connect(f"file:{db}?mode=ro&immutable=1", uri=True) as con:
        for molecule, parent, efo, phase in con.execute(QUERY):
            if phase is None:
                continue
            key = (parent or molecule, norm_ontology_id(efo))
            out[key] = max(out.get(key, float("-inf")), float(phase))
    return out


def chembl_max_phase(index: dict[tuple[str, str], float], programs: tuple[str, ...],
                     indication_id: str) -> float | None:
    vals = [index[(p, indication_id)] for p in programs if (p, indication_id) in index]
    return max(vals) if vals else None
