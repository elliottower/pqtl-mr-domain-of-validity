import pandas as pd

from ot_phase import link_hypotheses, nct_of_report, program_of


def hyps(*rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["hypothesis_id", "indication_id", "programs"])


def ci(*rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["drugId", "diseaseId", "maxClinicalStage", "clinicalReportIds"])


def test_nct_of_report_only_ctgov_source():
    assert nct_of_report("nct01234567", "ClinicalTrials.gov") == "NCT01234567"
    assert nct_of_report("nct01234567", "TTD") is None
    assert nct_of_report("0a1b2c", "ClinicalTrials.gov") is None
    assert nct_of_report("nct0123456", "ClinicalTrials.gov") is None


def test_program_of_uses_parent_where_present():
    parent = {"SALT": "PARENT", "PARENT": None}
    assert program_of("SALT", parent) == "PARENT"
    assert program_of("PARENT", parent) == "PARENT"
    assert program_of("UNLISTED", parent) == "UNLISTED"


def test_link_takes_child_rows_exact_indication_and_ctgov_reports_only():
    parent = {"SALT": "P1"}
    table = ci(("P1", "MONDO_1", "PHASE_2", ["nct00000001", "ttd1"]),
               ("SALT", "MONDO_1", "PHASE_3", ["nct00000002"]),
               ("P1", "MONDO_1_CHILD", "APPROVAL", ["nct00000003"]),
               ("P2", "MONDO_1", "APPROVAL", ["nct00000004"]),
               ("P1", "MONDO_9", "PHASE_2", None))
    report_nct = {"nct00000001": "NCT00000001", "nct00000002": "NCT00000002",
                  "nct00000003": "NCT00000003", "nct00000004": "NCT00000004"}
    (h,) = link_hypotheses(hyps(("h1", "MONDO_1", ("P1",))), table, parent, report_nct)
    assert h.max_stages == ("PHASE_2", "PHASE_3")
    assert h.nct_ids == ("NCT00000001", "NCT00000002")


def test_link_unions_programs_and_dedupes_trials():
    table = ci(("P1", "D", "PHASE_2", ["a"]), ("P2", "D", "PHASE_2", ["a", "b"]))
    (h,) = link_hypotheses(hyps(("h", "D", ("P1", "P2"))), table, {}, {"a": "NCT00000001", "b": "NCT00000002"})
    assert h.nct_ids == ("NCT00000001", "NCT00000002")


def test_link_with_no_rows_gives_no_phase():
    (h,) = link_hypotheses(hyps(("h", "D", ("P1",))), ci(("P1", "E", "PHASE_2", [])), {}, {})
    assert (h.max_stages, h.nct_ids) == ((), ())
