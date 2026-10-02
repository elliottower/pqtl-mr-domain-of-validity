import hashlib
import json
from pathlib import Path

import pandas as pd

from models import OUTCOME_COLUMNS, Status

INTERFACES = Path(__file__).resolve().parents[2] / "INTERFACES.md"


def interfaces_stage_c_columns() -> list[str]:
    section = INTERFACES.read_text().split("## Stage C", 1)[1].split("## Stage D", 1)[0]
    cols = []
    for line in section.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if not line.startswith("| ") or cells[0] in ("column",) or set(cells[0]) <= {"-"}:
            continue
        cols += [c.strip() for c in cells[0].split(",")]
    return cols


def interfaces_statuses() -> set[str]:
    section = INTERFACES.read_text().split("## Stage C", 1)[1]
    line = next(x for x in section.splitlines() if x.startswith("| status_24 |"))
    return {s.strip() for s in line.split("|")[2].split("/")}


def test_outcome_columns_match_interfaces_exactly():
    assert OUTCOME_COLUMNS == interfaces_stage_c_columns()


def test_status_values_are_the_interfaces_values_including_no_phase():
    assert {s.value for s in Status} == interfaces_statuses()
    assert "no_phase" in interfaces_statuses()


def test_csv_header_and_statuses(e2e):
    df = pd.read_csv(e2e["csv"], dtype=str, keep_default_na=False)
    assert list(df.columns) == interfaces_stage_c_columns()
    assert dict(zip(df["hypothesis_id"], df["status_24"])) == e2e["expected"]


def test_csv_values(e2e):
    df = pd.read_csv(e2e["csv"], dtype=str, keep_default_na=False).set_index("hypothesis_id")
    assert df.loc["h_adv", "advanced_24"] == "1"
    assert df.loc["h_zero", "advanced_24"] == "0"
    assert all(df.loc[h, "advanced_24"] == "" for h in ("h_active", "h_business", "h_undated", "h_young", "h_nophase"))
    # h_young: last Phase II end 2025-03 is mature only under the 12-month window.
    assert (df.loc["h_young", "status_12"], df.loc["h_young", "advanced_12"]) == ("no_observed_advancement", "0")
    assert df.loc["h_young", "status_36"] == "active"
    assert (df.loc["h_zero", "stop_code"], df.loc["h_zero", "why_stopped_texts"],
            df.loc["h_zero", "efficacy_coded"]) == ("efficacy", "Futility at interim", "true")
    assert df.loc["h_business", "stop_code"] == "business"
    assert df.loc["h_young", "stop_code"] == "completed_no_successor"
    assert df.loc["h_adv", "stop_code"] == ""
    # Time to Phase III on h_adv: Phase II start 2012-01-01 to Phase III start 2015-06-01.
    assert (df.loc["h_adv", "phase2_start_date"], df.loc["h_adv", "phase3_start_date"],
            df.loc["h_adv", "phase3_event"]) == ("2012-01-01", "2015-06-01", "1")
    assert df.loc["h_adv", "time_to_phase3_days"] == str((pd.Timestamp("2015-06-01") - pd.Timestamp("2012-01-01")).days)
    assert df.loc["h_adv", "earliest_phase2_start"] == "2012-01-01"
    assert (df.loc["h_zero", "phase3_event"], df.loc["h_zero", "time_to_phase3_days"]) == (
        "0", str((pd.Timestamp("2026-09-30") - pd.Timestamp("2016-01-01")).days))
    assert df.loc["h_undated", "time_to_phase3_days"] == ""
    # h_adv: Phase III completed 2018-01, not approved -> Phase III success 0.
    assert (df.loc["h_adv", "phase3_success"], df.loc["h_adv", "approved"]) == ("0", "false")
    assert df.loc["h_zero", "phase3_success"] == ""
    # ChEMBL: S1 (child of P1) has 3 on D1, P1 has 2 on D1 and 4 on D9 -> 3 for h_adv; none for h_zero (D2).
    assert df.loc["h_adv", "chembl_max_phase_for_ind"] == "3"
    assert df.loc["h_active", "chembl_max_phase_for_ind"] == "2"
    assert df.loc["h_zero", "chembl_max_phase_for_ind"] == ""
    # Last Phase II completion or termination: the date the maturation rule reads.
    assert df.loc["h_adv", "last_phase2_end_date"] == "2014-01-01"
    assert df.loc["h_zero", "last_phase2_end_date"] == "2018-05-01"
    assert df.loc["h_young", "last_phase2_end_date"] == "2025-03-01"
    assert df.loc["h_business", "last_phase2_end_date"] == "2015-01-01"
    assert df.loc["h_active", "last_phase2_end_date"] == "2017-01-01"
    assert df.loc["h_undated", "last_phase2_end_date"] == ""


def test_links_follow_the_rule(e2e):
    links = {}
    for h, n in e2e["links"]:
        links.setdefault(h, []).append(n)
    assert links == {"h_adv": ["NCT00000001", "NCT00000002"], "h_zero": ["NCT00000003"],
                     "h_active": ["NCT00000004", "NCT00000005"], "h_business": ["NCT00000006"],
                     "h_young": ["NCT00000007"]}


def test_every_linked_record_is_in_the_snapshot_and_nothing_else(e2e):
    src = e2e["source"]
    used = pd.read_csv(src.rows, sep="\t", dtype=str, keep_default_na=False)
    assert set(used["nct_id"]) == {n for _, n in e2e["links"]}
    assert "NCT09999999" not in set(used["nct_id"])   # in the archive, linked to no hypothesis
    manifest = json.loads(src.manifest.read_text())
    assert manifest["source_sha256"] == hashlib.sha256((e2e["aact"] / "studies.txt").read_bytes()).hexdigest()
    assert (manifest["requested"], manifest["found"], manifest["not_found"]) == (7, 7, 0)


def test_post_freeze_update_diagnostic_counts_linked_trials_per_hypothesis(e2e):
    diag = {d.hypothesis_id: d for d in e2e["diagnostics"]}
    # NCT00000002 (h_adv) was last updated 2026-10-01, after the freeze; every other record before.
    assert (diag["h_adv"].n_linked_trials, diag["h_adv"].n_last_update_after_freeze) == (2, 1)
    assert diag["h_adv"].n_first_submitted_after_freeze == 0
    assert all(d.n_last_update_after_freeze == 0 for h, d in diag.items() if h != "h_adv")
    assert (diag["h_undated"].n_linked_trials, diag["h_nophase"].n_linked_trials) == (0, 0)
