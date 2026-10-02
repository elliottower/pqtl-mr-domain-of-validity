"""Open Targets 26.09 phase and the trial-linking rule.

Phase (PREREG §Measured variables): `clinical_indication.maxClinicalStage` for the
hypothesis's drugs and the indication. A clinical_indication `drugId` belongs to a program
through `drug_molecule.parentId` (PREREG §Study design, Drug programs), so a salt's row counts
for its parent's program.

TRIAL-LINKING RULE (the plan does not fix one; see README.md): the ClinicalTrials.gov trials of
a hypothesis are the NCT registrations among the `clinicalReportIds` of the same
clinical_indication rows that give its phase, that is, rows whose `diseaseId` equals the
indication id exactly and whose drug's program is one of the hypothesis's programs. A report
is an NCT registration when its `clinical_report.source` is "ClinicalTrials.gov" and its `id` is
an NCT id (in 26.09 these are the 230,989 AACT rows, ids written `nct` + 8 digits). Reports from
other sources (TTD, DailyMed, EMA, FDA, ...) are not trials here.
"""
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from models import HypothesisDrugs, LinkageAuditRow

NCT_RE = re.compile(r"NCT\d{8}", re.IGNORECASE)
CTGOV_SOURCE = "ClinicalTrials.gov"


def program_of(drug_id: str, parent: dict[str, str | None]) -> str:
    p = parent.get(drug_id)
    return p if isinstance(p, str) and p else drug_id


def nct_of_report(report_id: str, source: str | None) -> str | None:
    """The NCT id of a ClinicalTrials.gov clinical_report row, else None."""
    if source != CTGOV_SOURCE or not NCT_RE.fullmatch(report_id or ""):
        return None
    return report_id.upper()


def load_ot(ot_dir: Path, clinical_report: Path) -> tuple[pd.DataFrame, dict[str, str | None], dict[str, str]]:
    """clinical_indication rows, drugId -> parentId, report id -> NCT id."""
    ci = pq.read_table(ot_dir / "clinical_indication.parquet",
                       columns=["drugId", "diseaseId", "maxClinicalStage", "clinicalReportIds"]).to_pandas()
    mol = pq.read_table(ot_dir / "drug_molecule.parquet", columns=["id", "parentId"]).to_pandas()
    rep = pq.read_table(clinical_report, columns=["id", "source"]).to_pandas()
    report_nct = {}
    for rid, src in zip(rep["id"], rep["source"]):
        n = nct_of_report(rid, src)
        if n is not None:
            report_nct[rid] = n
    return ci, dict(zip(mol["id"], mol["parentId"])), report_nct


def rows_by_program_indication(ci: pd.DataFrame, parent: dict[str, str | None]
                               ) -> dict[tuple[str, str], list[tuple[int, str, list[str]]]]:
    """(program, diseaseId) -> its clinical_indication rows as (row position, maxClinicalStage,
    clinicalReportIds)."""
    by_key: dict[tuple[str, str], list[tuple[int, str, list[str]]]] = defaultdict(list)
    for i, (drug, disease, stage, reports) in enumerate(zip(ci["drugId"], ci["diseaseId"], ci["maxClinicalStage"],
                                                            ci["clinicalReportIds"])):
        rids = list(reports) if reports is not None else []
        by_key[(program_of(drug, parent), disease)].append((i, stage, rids))
    return by_key


def link_hypotheses(hyps: pd.DataFrame, ci: pd.DataFrame, parent: dict[str, str | None],
                    report_nct: dict[str, str]) -> list[HypothesisDrugs]:
    """Per hypothesis: every maxClinicalStage and every linked NCT id (trial-linking rule)."""
    by_key = rows_by_program_indication(ci, parent)
    out = []
    for hid, ind, programs in zip(hyps["hypothesis_id"], hyps["indication_id"], hyps["programs"]):
        rows = [r for p in programs for r in by_key.get((p, ind), [])]
        stages = tuple(sorted({s for _, s, _ in rows if isinstance(s, str)}))
        ncts = tuple(sorted({report_nct[r] for _, _, rids in rows for r in rids if r in report_nct}))
        out.append(HypothesisDrugs(hypothesis_id=hid, indication_id=ind, drug_program_ids=tuple(programs),
                                   max_stages=stages, nct_ids=ncts))
    return out


def linkage_audit(hyps: pd.DataFrame, ci: pd.DataFrame, parent: dict[str, str | None], report_nct: dict[str, str],
                  in_snapshot: set[str]) -> list[LinkageAuditRow]:
    """Counts at each step of the trial-linking rule, over all hypotheses: clinical_indication rows,
    their report ids, the ClinicalTrials.gov registrations among them, those present in the AACT
    snapshot (`in_snapshot`), and the hypotheses left with no trial. Descriptive; no label reads it."""
    by_key = rows_by_program_indication(ci, parent)
    ci_rows: set[int] = set()
    reports: set[str] = set()
    ncts: set[str] = set()
    with_rows = no_trial = no_trial_in_snapshot = 0
    for ind, programs in zip(hyps["indication_id"], hyps["programs"]):
        rows = [r for p in programs for r in by_key.get((p, ind), [])]
        own_reports = {r for _, _, rids in rows for r in rids}
        own_ncts = {report_nct[r] for r in own_reports if r in report_nct}
        ci_rows |= {i for i, _, _ in rows}
        reports |= own_reports
        ncts |= own_ncts
        with_rows += bool(rows)
        no_trial += not own_ncts
        no_trial_in_snapshot += not (own_ncts & in_snapshot)
    n = len(hyps)
    steps = [
        ("hypotheses", n, "rows of hypotheses.csv"),
        ("hypotheses_with_clinical_indication_rows", with_rows,
         "hypotheses with at least one Open Targets clinical_indication row for one of their programs and their "
         "indication id"),
        ("hypotheses_without_clinical_indication_rows", n - with_rows, "hypotheses with none (status no_phase)"),
        ("clinical_indication_rows", len(ci_rows), "distinct clinical_indication rows of those hypotheses"),
        ("report_ids", len(reports), "distinct clinicalReportIds on those rows"),
        ("report_ids_not_clinicaltrials_gov", len({r for r in reports if r not in report_nct}),
         "report ids that are not ClinicalTrials.gov registrations in clinical_report (another source, or absent)"),
        ("clinicaltrials_gov_ids", len(ncts), "distinct NCT ids among the report ids"),
        ("clinicaltrials_gov_ids_in_aact_snapshot", len(ncts & in_snapshot), "NCT ids with a row in the AACT studies table"),
        ("clinicaltrials_gov_ids_not_in_aact_snapshot", len(ncts - in_snapshot), "NCT ids with no row there (unmatched)"),
        ("hypotheses_with_zero_linked_trials", no_trial, "hypotheses with no NCT id among their report ids"),
        ("hypotheses_with_zero_trials_in_aact_snapshot", no_trial_in_snapshot,
         "hypotheses with no NCT id that has a row in the AACT studies table (includes the previous row)"),
    ]
    return [LinkageAuditRow(step=s, count=int(c), description=d) for s, c, d in steps]
