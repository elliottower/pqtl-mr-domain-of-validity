import hashlib
import json
import sqlite3

import pytest

from chembl import ChemblNotVerified, load_chembl_index, pinned_sha256, verified_chembl_db


def make_db(path) -> str:
    with sqlite3.connect(path) as con:
        con.executescript("""
            CREATE TABLE molecule_dictionary (molregno INTEGER, chembl_id TEXT);
            CREATE TABLE molecule_hierarchy (molregno INTEGER, parent_molregno INTEGER);
            CREATE TABLE drug_indication (molregno INTEGER, efo_id TEXT, max_phase_for_ind REAL);
            INSERT INTO molecule_dictionary VALUES (1, 'P1'), (2, 'S1');
            INSERT INTO molecule_hierarchy VALUES (1, 1), (2, 1);
            INSERT INTO drug_indication VALUES (2, 'EFO:0000001', 3.0), (1, 'EFO:0000001', 2.0);
        """)
    return hashlib.sha256(path.read_bytes()).hexdigest()


@pytest.fixture
def chembl_dir(tmp_path):
    d = tmp_path / "chembl_37"
    d.mkdir()
    return d, make_db(d / "chembl_37.db")


@pytest.mark.parametrize("layout", [
    lambda sha: {"chembl_37.db": sha},
    lambda sha: {"chembl_37.db": {"sha256": sha, "bytes": 1}},
    lambda sha: {"files": [{"path": "/chembl_37/chembl_37.tar.gz", "sha256": "0" * 64},
                           {"path": "/chembl_37/chembl_37.db", "sha256": sha}]},
    lambda sha: {"verified": {"db": {"file": "chembl_37.db", "sha256": sha.upper()}}},
])
def test_verified_db_opens_when_the_pin_matches(chembl_dir, layout):
    d, sha = chembl_dir
    (d / "VERIFIED.json").write_text(json.dumps(layout(sha)))
    db, got = verified_chembl_db(d)
    assert (db, got) == (d / "chembl_37.db", sha)
    assert load_chembl_index(db) == {("P1", "EFO_0000001"): 3.0}


def test_refuses_without_verified_json(chembl_dir):
    d, _ = chembl_dir
    with pytest.raises(ChemblNotVerified, match="VERIFIED.json is missing"):
        verified_chembl_db(d)


def test_refuses_a_database_that_differs_from_the_pin(chembl_dir):
    d, sha = chembl_dir
    (d / "VERIFIED.json").write_text(json.dumps({"chembl_37.db": sha}))
    with sqlite3.connect(d / "chembl_37.db") as con:
        con.execute("INSERT INTO drug_indication VALUES (1, 'EFO:0000002', 4.0)")
    with pytest.raises(ChemblNotVerified, match="differs"):
        verified_chembl_db(d)


@pytest.mark.parametrize("doc", [{"chembl_37.tar.gz": "a" * 64}, {"chembl_37.db": "a" * 64, "x": {"chembl_37.db": "b" * 64}}])
def test_refuses_no_pin_or_two_pins(chembl_dir, doc):
    d, _ = chembl_dir
    (d / "VERIFIED.json").write_text(json.dumps(doc))
    with pytest.raises(ChemblNotVerified, match="pins"):
        verified_chembl_db(d)


def test_refuses_a_missing_database(tmp_path):
    (tmp_path / "VERIFIED.json").write_text(json.dumps({"chembl_37.db": "a" * 64}))
    with pytest.raises(ChemblNotVerified, match="missing"):
        verified_chembl_db(tmp_path)


def test_pin_for_another_file_is_not_taken():
    assert pinned_sha256({"files": [{"path": "chembl_37_sqlite.tar.gz", "sha256": "a" * 64}]}) == set()
