import hashlib
import json
import sqlite3
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ctgov import AACT_SNAPSHOT_DATE, STUDIES_COLUMNS, AactTrials
from hypotheses import read_hypotheses
from run_stage_c import build_rows, write_outcomes

# An AACT flat file carries more columns than stage C reads; two extra ones, one holding quotes.
AACT_HEADER = ("nct_id", "brief_title", *STUDIES_COLUMNS[1:], "enrollment")


def aact_row(nct: str, status: str, phase: str, start: str = "", end: str = "", why: str = "",
             submitted: str = "2010-01-01", updated: str = "2020-01-01", title: str = "A study") -> dict:
    """One `studies` row; `start` / `end` are `*_month_year` as ClinicalTrials.gov gives them."""
    return {"nct_id": nct, "brief_title": title, "overall_status": status, "why_stopped": why, "phase": phase,
            "start_month_year": start, "start_date": "", "start_date_type": "ACTUAL" if start else "",
            "completion_month_year": end, "completion_date": "", "completion_date_type": "ACTUAL" if end else "",
            "study_first_submitted_date": submitted, "last_update_submitted_date": updated, "enrollment": "100"}


def write_aact(root: Path, rows: list[dict], snapshot_date: str = AACT_SNAPSHOT_DATE, pin: bool = True) -> Path:
    """An AACT archive directory as download_aact leaves it: studies.txt and VERIFIED.json."""
    root.mkdir(parents=True, exist_ok=True)
    lines = ["|".join(AACT_HEADER)] + ["|".join(r[c] for c in AACT_HEADER) for r in rows]
    (root / "studies.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    sha = hashlib.sha256((root / "studies.txt").read_bytes()).hexdigest()
    if pin:
        (root / "VERIFIED.json").write_text(json.dumps({"snapshot_date": snapshot_date,
                                                        "files": {"studies.txt": {"sha256": sha}}}))
    return root


# Synthetic universe. Programs P1..P6, a salt S1 of P1, indications D1, D2.
HYPOTHESES = [
    # id, indication, programs, expected status_24
    ("h_adv", "D1", "P1", "advanced"),                 # salt S1 reached PHASE_3 on D1
    ("h_zero", "D2", "P1", "no_observed_advancement"),  # Phase II terminated for futility in 2018
    ("h_active", "D1", "P2", "active"),                # one Phase I recruiting
    ("h_business", "D1", "P3", "business_only"),       # only stop is a sponsor decision
    ("h_undated", "D1", "P4", "undated"),              # no linked ClinicalTrials.gov trial
    ("h_young", "D2", "P5", "active"),                 # last Phase II ended 2025-03: < 24 months
    ("h_nophase", "D2", "P6", "no_phase"),             # no clinical_indication row
]
CI = [
    ("P1", "D1", "PHASE_2", ["nct00000001"]),
    ("S1", "D1", "PHASE_3", ["nct00000002", "ttd-1"]),
    ("P1", "D2", "PHASE_2", ["nct00000003"]),
    ("P2", "D1", "PHASE_2", ["nct00000004", "nct00000005"]),
    ("P3", "D1", "PHASE_2", ["nct00000006"]),
    ("P4", "D1", "PHASE_2", ["dailymed-1"]),
    ("P5", "D2", "PHASE_2", ["nct00000007"]),
    ("P6", "D1", "PHASE_2", []),
]
AACT_ROWS = [
    aact_row("NCT00000001", "COMPLETED", "PHASE2", "2012-01", "2014-01"),
    aact_row("NCT00000002", "COMPLETED", "PHASE3", "2015-06-01", "2018-01", updated="2026-10-01"),
    aact_row("NCT00000003", "TERMINATED", "PHASE2", "2016-01", "2018-05", "Futility at interim",
             title='The "ALPHA" trial'),
    aact_row("NCT00000004", "Completed", "Phase 2", "January 2015", "January 2017"),   # legacy vocabulary
    aact_row("NCT00000005", "RECRUITING", "PHASE1", "2025-01"),
    aact_row("NCT00000006", "TERMINATED", "PHASE2", "2014-01", "2015-01", "Sponsor decision"),
    aact_row("NCT00000007", "COMPLETED", "PHASE2", "2023-01", "2025-03"),
    aact_row("NCT09999999", "RECRUITING", "PHASE2", "2024-01"),   # in the archive, linked to nothing
]


@pytest.fixture
def e2e(tmp_path: Path) -> dict:
    a = tmp_path / "stages" / "A" / "output"
    a.mkdir(parents=True)
    pd.DataFrame({"hypothesis_id": [h[0] for h in HYPOTHESES], "gene_symbol": "G",
                  "indication_id": [h[1] for h in HYPOTHESES], "direction": "decrease",
                  "mechanism_class": "blocking", "drug_program_ids": [h[2] for h in HYPOTHESES],
                  "eligible": True}).to_csv(a / "hypotheses.csv", index=False)
    # Unreadable for any user: a directory in the file's place, which no read can open. A mode-0 file
    # stops only an unprivileged user, and the Modal test container runs as root.
    b = tmp_path / "stages" / "B" / "output" / "evidence.csv"
    b.mkdir(parents=True)

    ot = tmp_path / "ot"
    ot.mkdir()
    pq.write_table(pa.table({"drugId": [r[0] for r in CI], "diseaseId": [r[1] for r in CI],
                             "maxClinicalStage": [r[2] for r in CI],
                             "clinicalReportIds": [r[3] for r in CI]}), ot / "clinical_indication.parquet")
    pq.write_table(pa.table({"id": ["P1", "S1", "P2"], "parentId": [None, "P1", None]}), ot / "drug_molecule.parquet")
    reports = sorted({r for row in CI for r in row[3]})
    pq.write_table(pa.table({"id": reports, "source": ["ClinicalTrials.gov" if r.startswith("nct") else "TTD"
                                                        for r in reports]}), ot / "clinical_report.parquet")

    db = tmp_path / "chembl.db"
    with sqlite3.connect(db) as con:
        con.executescript("""
            CREATE TABLE molecule_dictionary (molregno INTEGER, chembl_id TEXT);
            CREATE TABLE molecule_hierarchy (molregno INTEGER, parent_molregno INTEGER);
            CREATE TABLE drug_indication (molregno INTEGER, efo_id TEXT, max_phase_for_ind REAL);
            INSERT INTO molecule_dictionary VALUES (1, 'P1'), (2, 'S1'), (3, 'P2');
            INSERT INTO molecule_hierarchy VALUES (1, 1), (2, 1), (3, 3);
            INSERT INTO drug_indication VALUES (2, 'D1', 3.0), (1, 'D1', 2.0), (3, 'D1', 2.0), (1, 'D9', 4.0);
        """)

    aact = write_aact(tmp_path / "aact_20260930", AACT_ROWS)
    source = AactTrials(aact, tmp_path / "snap")
    hyps = read_hypotheses(a / "hypotheses.csv")
    built = build_rows(hyps, ot, ot / "clinical_report.parquet", source, db)
    out = tmp_path / "outcomes.csv"
    write_outcomes(built.rows, out)
    return {"rows": built.rows, "links": built.links, "diagnostics": built.post_freeze, "built": built, "csv": out,
            "source": source, "aact": aact, "expected": {h[0]: h[3] for h in HYPOTHESES}}
