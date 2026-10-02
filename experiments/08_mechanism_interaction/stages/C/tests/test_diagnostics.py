"""The three stage C diagnostics on hand-worked fixtures. Freeze 2026-09-30, so the maturation
thresholds (last Phase II end on or before) are 2025-09-30 (12 months), 2024-09-30 (24), 2023-09-30
(36) and 2022-09-30 (48)."""
from datetime import date

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from ctgov import AactTrials, parse_partial_date, parse_study
from diagnostics import (PHASE3_QUANTITIES, STATUS_QUANTITIES, differences, partial_date_boundary, phase23_as_phase3,
                         phase23_diagnostic)
from models import OUTCOME_FREEZE_DATE, Trial
from ot_phase import link_hypotheses, linkage_audit
from rules import code_hypothesis
from run_stage_c import build_rows
from tests.conftest import AACT_ROWS, CI, HYPOTHESES, aact_row, write_aact

F = OUTCOME_FREEZE_DATE
NOA, ACTIVE, UNDATED = "no_observed_advancement", "active", "undated"
_n = iter(range(20_000_000, 99_999_999))


def study(phase="PHASE2", status="COMPLETED", start="", end="", why="") -> dict:
    return aact_row(f"NCT{next(_n):08d}", status, phase, start, end, why)


def both_readings(rows: list[dict]) -> tuple[list[Trial], list[Trial]]:
    return [parse_study(r) for r in rows], [parse_study(r, last_day=True) for r in rows]


def boundary(rows: list[dict], stages=("PHASE_2",)) -> list[tuple[str, str, str]]:
    first, last = both_readings(rows)
    return [(d.quantity, d.registered, d.alternative) for d in partial_date_boundary("h", list(stages), first, last, F)]


def phase23(trials: list[Trial], stages=("PHASE_2_3",)) -> dict[str, tuple[str, str]]:
    return {d.quantity: (d.registered, d.alternative) for d in phase23_diagnostic("h", list(stages), trials, F)}


def trial(phases=("PHASE2",), status="COMPLETED", start=None, end=None) -> Trial:
    return Trial(nct_id=f"NCT{next(_n):08d}", phases=tuple(phases), overall_status=status, start_date=start,
                 completion_date=end)


# ---- partial dates on the last day ----------------------------------------------------------------

@pytest.mark.parametrize("raw,first,last", [
    ("2024", date(2024, 1, 1), date(2024, 12, 31)),
    ("2024-02", date(2024, 2, 1), date(2024, 2, 29)),
    ("2023-02", date(2023, 2, 1), date(2023, 2, 28)),
    ("2024-09", date(2024, 9, 1), date(2024, 9, 30)),
    ("2024-09-15", date(2024, 9, 15), date(2024, 9, 15)),
    ("February 2024", date(2024, 2, 1), date(2024, 2, 29)),
    ("December, 2019", date(2019, 12, 1), date(2019, 12, 31)),
    ("March 5, 2020", date(2020, 3, 5), date(2020, 3, 5)),
])
def test_a_partial_date_reads_as_the_first_day_and_in_the_diagnostic_as_the_last(raw, first, last):
    assert (parse_partial_date(raw), parse_partial_date(raw, last_day=True)) == (first, last)
    assert parse_partial_date("", last_day=True) is None


def test_the_last_day_reading_changes_partial_start_and_completion_dates_and_nothing_else():
    row = aact_row("NCT00000001", "COMPLETED", "PHASE2", "2015", "2017-06", submitted="2014-12-20", updated="2018-01-02")
    first, last = parse_study(row), parse_study(row, last_day=True)
    assert (first.start_date, first.completion_date) == (date(2015, 1, 1), date(2017, 6, 1))
    assert (last.start_date, last.completion_date) == (date(2015, 12, 31), date(2017, 6, 30))
    assert last.model_dump(exclude={"start_date", "completion_date"}) == first.model_dump(exclude={"start_date", "completion_date"})


@pytest.mark.parametrize("end,quantity,registered,alternative", [
    # 2025-01-01 is on or before 2025-09-30 (mature at 12 months); 2025-12-31 is not
    ("2025", "status_12", NOA, ACTIVE),
    # 2024-01-01 is on or before 2024-09-30; 2024-12-31 is not, and is still before 2025-09-30
    ("2024", "status_24", NOA, ACTIVE),
    ("2023", "status_36", NOA, ACTIVE),
    ("2022", "status_48", NOA, ACTIVE),
])
def test_a_year_only_phase2_end_moves_exactly_the_window_whose_threshold_falls_in_that_year(end, quantity, registered,
                                                                                            alternative):
    assert boundary([study(end=end)]) == [(quantity, registered, alternative)]


def test_no_month_level_date_can_move_a_window_because_every_threshold_is_the_last_day_of_september():
    for year in range(2020, 2027):
        for month in range(1, 13):
            assert boundary([study(end=f"{year}-{month:02d}")]) == [], (year, month)
    for raw in ("September 2024", "September 2025", "October 2024", "August 2023"):
        assert boundary([study(end=raw)]) == []


@pytest.mark.parametrize("rows", [
    [study(end="2024-09-30")], [study(end="2024-10-01")],                    # full dates are not partial
    [study(end="2021")],                                                     # mature in every window either way
    [study(end="2026")],             # 2026-01-01 ended and immature; 2026-12-31 is after the freeze: active either way
    [study(end="2024"), study(end="2025-03-01")],                            # a later full-date end decides
    [study(end="2024"), study(phase="PHASE1", status="RECRUITING")],         # active whatever the date
    [study(phase="PHASE3", end="2024")],                                     # not a Phase II end: undated either way
    [],
])
def test_hypotheses_whose_status_does_not_depend_on_the_reading_give_no_row(rows):
    assert boundary(rows) == []


def test_an_advanced_hypothesis_gives_no_row_whatever_its_dates():
    assert boundary([study(end="2024")], stages=("PHASE_2", "PHASE_3")) == []


def test_a_year_only_end_in_the_freeze_year_differs_where_the_registered_reading_was_mature_elsewhere():
    # two Phase II ends: 2021 (mature everywhere) and year-only 2026. Registered: 2026-01-01 is the last end, immature
    # in every window, so `active`. Last day: 2026-12-31 is after the freeze, so that trial is ongoing: `active` too.
    assert boundary([study(end="2021"), study(end="2026")]) == []


# ---- Phase II/III registrations read as Phase III ---------------------------------------------------

def test_phase23_as_phase3_rewrites_only_registrations_listing_both_phases():
    both = trial(phases=("PHASE2", "PHASE3"))
    assert phase23_as_phase3(both).phases == ("PHASE3",)
    assert phase23_as_phase3(both).model_dump(exclude={"phases"}) == both.model_dump(exclude={"phases"})
    for phases in (("PHASE2",), ("PHASE3",), ("PHASE1", "PHASE2"), ("PHASE1",), ()):
        t = trial(phases=phases)
        assert phase23_as_phase3(t) is t


def test_a_phase23_start_after_a_phase2_start_becomes_the_phase3_event():
    trials = [trial(start=date(2014, 1, 1), end=date(2016, 1, 1)),
              trial(phases=("PHASE2", "PHASE3"), start=date(2017, 6, 1), end=date(2019, 6, 1))]
    # registered: both are Phase II, no Phase III start, censored at the freeze. Alternative: event on 2017-06-01,
    # 1096 days (2014, 2015, 2016) + 151 days (January to May 2017) after 2014-01-01. Every window is mature either way.
    assert phase23(trials) == {"phase3_start_date": ("", "2017-06-01"),
                               "time_to_phase3_days": (str((F - date(2014, 1, 1)).days), "1247"),
                               "phase3_event": ("0", "1")}
    assert (date(2017, 6, 1) - date(2014, 1, 1)).days == 1096 + 151 == 1247


def test_a_recent_phase23_end_no_longer_delays_maturation():
    trials = [trial(end=date(2016, 1, 1)), trial(phases=("PHASE2", "PHASE3"), end=date(2025, 6, 1))]
    # registered: last Phase II end 2025-06-01, mature at 12 months only. Alternative: last Phase II end 2016-01-01.
    assert phase23(trials) == {"status_24": (ACTIVE, NOA), "status_36": (ACTIVE, NOA), "status_48": (ACTIVE, NOA)}


def test_a_hypothesis_whose_only_phase2_trial_is_a_phase23_registration_loses_its_phase2_dates():
    trials = [trial(phases=("PHASE2", "PHASE3"), start=date(2015, 1, 1), end=date(2018, 1, 1))]
    assert phase23(trials) == {
        **{q: (NOA, UNDATED) for q in STATUS_QUANTITIES},
        "phase2_start_date": ("2015-01-01", ""), "time_to_phase3_days": (str((F - date(2015, 1, 1)).days), ""),
        "phase3_event": ("0", "")}


@pytest.mark.parametrize("trials", [
    [],
    [trial(start=date(2014, 1, 1), end=date(2016, 1, 1))],
    [trial(start=date(2014, 1, 1), end=date(2016, 1, 1)), trial(phases=("PHASE3",), start=date(2017, 1, 1))],
    [trial(phases=("PHASE1", "PHASE2"), start=date(2014, 1, 1), end=date(2016, 1, 1))],
])
def test_a_hypothesis_without_a_phase23_registration_gives_no_row(trials):
    assert phase23(trials, stages=("PHASE_2",)) == {}


def test_the_diagnostic_leaves_the_open_targets_stage_as_registered():
    # PHASE_2_3 is not advanced under the registered rule, and the diagnostic does not change that
    trials = [trial(phases=("PHASE2", "PHASE3"), start=date(2015, 1, 1), end=date(2018, 1, 1))]
    reg = code_hypothesis(["PHASE_2_3"], trials, F)
    alt = code_hypothesis(["PHASE_2_3"], [phase23_as_phase3(t) for t in trials], F)
    assert set(reg.advanced.values()) == {0} and set(alt.advanced.values()) == {None}
    assert "advanced" not in {s.value for s in alt.statuses.values()}


def test_differences_lists_only_the_quantities_asked_for_and_only_where_they_differ():
    reg = code_hypothesis(["PHASE_2"], [trial(start=date(2014, 1, 1), end=date(2024, 1, 1))], F)
    alt = code_hypothesis(["PHASE_2"], [trial(start=date(2015, 1, 1), end=date(2024, 12, 31))], F)
    rows = differences("h", reg, alt, STATUS_QUANTITIES)
    assert [(r.hypothesis_id, r.quantity, r.registered, r.alternative) for r in rows] == [("h", "status_24", NOA, ACTIVE)]
    assert [r.quantity for r in differences("h", reg, alt, PHASE3_QUANTITIES)] == ["phase2_start_date", "time_to_phase3_days"]
    assert differences("h", reg, reg, STATUS_QUANTITIES + PHASE3_QUANTITIES) == []


# ---- linkage audit ----------------------------------------------------------------------------------

def counts(rows) -> dict[str, int]:
    return {r.step: r.count for r in rows}


def test_linkage_audit_counts_each_step_on_a_hand_worked_table():
    hyps = pd.DataFrame([("h1", "D", ("P1", "P2")), ("h2", "D", ("P2",)), ("h3", "E", ("P3",)), ("h4", "E", ("P9",))],
                        columns=["hypothesis_id", "indication_id", "programs"])
    ci = pd.DataFrame([("P1", "D", "PHASE_2", ["a", "b"]),          # h1
                       ("P2", "D", "PHASE_2", ["b", "c", "x"]),     # h1, h2; x is not a ClinicalTrials.gov report
                       ("SALT", "D", "PHASE_3", ["d"]),             # a salt of P2: h1, h2
                       ("P3", "E", "PHASE_2", ["e"]),               # h3
                       ("P1", "E", "PHASE_2", ["z"])],              # no hypothesis has P1 on E
                      columns=["drugId", "diseaseId", "maxClinicalStage", "clinicalReportIds"])
    parent = {"SALT": "P2"}
    report_nct = {"a": "NCT00000001", "b": "NCT00000002", "c": "NCT00000003", "d": "NCT00000004", "e": "NCT00000005",
                  "z": "NCT00000006"}
    in_snapshot = {"NCT00000001", "NCT00000002", "NCT00000003", "NCT00000006"}
    rows = linkage_audit(hyps, ci, parent, report_nct, in_snapshot)
    assert counts(rows) == {
        "hypotheses": 4,
        "hypotheses_with_clinical_indication_rows": 3,
        "hypotheses_without_clinical_indication_rows": 1,           # h4
        "clinical_indication_rows": 4,                              # the P1-on-E row belongs to no hypothesis
        "report_ids": 6,                                            # a, b, c, x, d, e
        "report_ids_not_clinicaltrials_gov": 1,                     # x
        "clinicaltrials_gov_ids": 5,
        "clinicaltrials_gov_ids_in_aact_snapshot": 3,
        "clinicaltrials_gov_ids_not_in_aact_snapshot": 2,           # NCT00000004, NCT00000005
        "hypotheses_with_zero_linked_trials": 1,                    # h4
        "hypotheses_with_zero_trials_in_aact_snapshot": 2,          # h3 (its one trial is not in the snapshot), h4
    }
    assert [r.step for r in rows] == list(counts(rows)) and all(r.description for r in rows)
    linked = {h.hypothesis_id: h.nct_ids for h in link_hypotheses(hyps, ci, parent, report_nct)}
    assert len(set().union(*linked.values())) == counts(rows)["clinicaltrials_gov_ids"]
    assert sum(1 for n in linked.values() if not n) == counts(rows)["hypotheses_with_zero_linked_trials"]


def test_linkage_audit_of_the_end_to_end_universe(e2e):
    # 7 hypotheses; h_nophase has no clinical_indication row; h_undated's only report is a DailyMed one
    assert counts(e2e["built"].linkage_audit) == {
        "hypotheses": 7, "hypotheses_with_clinical_indication_rows": 6, "hypotheses_without_clinical_indication_rows": 1,
        "clinical_indication_rows": 7, "report_ids": 9, "report_ids_not_clinicaltrials_gov": 2,
        "clinicaltrials_gov_ids": 7, "clinicaltrials_gov_ids_in_aact_snapshot": 7,
        "clinicaltrials_gov_ids_not_in_aact_snapshot": 0, "hypotheses_with_zero_linked_trials": 2,
        "hypotheses_with_zero_trials_in_aact_snapshot": 2}
    assert len(HYPOTHESES) == 7 and sum(len(r[3]) for r in CI) == 9 and len(AACT_ROWS) == 8
    # every Phase II end in that universe is a month-level or full date, and no registration is Phase II/III
    assert e2e["built"].partial_date_boundary == [] and e2e["built"].phase23 == []


def test_the_diagnostics_reach_the_run_from_the_archived_rows_and_change_no_outcome_row(tmp_path):
    a = tmp_path / "stages" / "A" / "output"
    a.mkdir(parents=True)
    hyps = pd.DataFrame({"hypothesis_id": ["h1"], "indication_id": ["D1"], "programs": [("P1",)]})
    ot = tmp_path / "ot"
    ot.mkdir()
    pq.write_table(pa.table({"drugId": ["P1"], "diseaseId": ["D1"], "maxClinicalStage": ["PHASE_2"],
                             "clinicalReportIds": [["nct00000001", "nct00000002", "nct00000003"]]}),
                   ot / "clinical_indication.parquet")
    pq.write_table(pa.table({"id": ["P1"], "parentId": pa.array([None], pa.string())}), ot / "drug_molecule.parquet")
    pq.write_table(pa.table({"id": ["nct00000001", "nct00000002", "nct00000003"], "source": ["ClinicalTrials.gov"] * 3}),
                   ot / "clinical_report.parquet")
    aact = write_aact(tmp_path / "aact_20260930", [
        aact_row("NCT00000001", "COMPLETED", "PHASE2", "2012", "2024"),
        aact_row("NCT00000002", "COMPLETED", "PHASE2/PHASE3", "2015-03", "2017"),
    ])
    built = build_rows(hyps, ot, ot / "clinical_report.parquet", AactTrials(aact, tmp_path / "snap"), None)
    (row,) = built.rows
    # registered: Phase II ends 2024-01-01 and 2017-01-01, origin 2012-01-01, no Phase III start
    assert (row.status_12.value, row.status_24.value, row.status_36.value, row.status_48.value) == (NOA, NOA, ACTIVE, ACTIVE)
    assert (row.last_phase2_end_date, row.phase2_start_date, row.phase3_event) == (date(2024, 1, 1), date(2012, 1, 1), 0)
    # last day: the last Phase II end is 2024-12-31, after the 24-month threshold
    assert [(d.hypothesis_id, d.quantity, d.registered, d.alternative) for d in built.partial_date_boundary] == [
        ("h1", "status_24", NOA, ACTIVE)]
    # Phase II/III as Phase III: the event is the start on 2015-03-01, 366 + 365 + 365 + 31 + 28 days after 2012-01-01
    assert {d.quantity: (d.registered, d.alternative) for d in built.phase23} == {
        "phase3_start_date": ("", "2015-03-01"),
        "time_to_phase3_days": (str((F - date(2012, 1, 1)).days), "1155"), "phase3_event": ("0", "1")}
    assert counts(built.linkage_audit)["clinicaltrials_gov_ids_not_in_aact_snapshot"] == 1
    assert counts(built.linkage_audit)["hypotheses_with_zero_trials_in_aact_snapshot"] == 0
