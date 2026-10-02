"""Stage C: outcome coding (PREREG §Measured variables; INTERFACES.md §Stage C).

Reads only: A/output/hypotheses.csv (hypothesis_id, indication_id, drug_program_ids), Open
Targets 26.09 clinical_indication, drug_molecule and clinical_report, the ChEMBL 37 SQLite
database (cross-check, opened only after its sha256 matches VERIFIED.json), and ClinicalTrials.gov
as of the outcome-freeze date: `studies.txt` of the AACT daily flat-file snapshot of 2026-09-30,
read only after its sha256 matches the VERIFIED.json the archive check wrote (ctgov.py). Never reads
stages/B or any evidence value.

Inputs are resolved under `--inputs-root`, a directory laid out like
experiments/08_mechanism_interaction: the Open Targets tables at
feasibility/v2_all_indications/inputs/ot_26.09/ and inputs/ot_26.09/clinical_report.parquet.
Real runs use the Modal volume `pqtl-v8-inputs` (modal_stage_c.py mounts it); the local tree
is the default only for development.

Writes, in order, before printing anything:
  <out>/trials_snapshot/studies_linked.tsv, not_found.jsonl, manifest.json  (the archived rows used)
  <out>/trial_links.csv          hypothesis_id, nct_id  (the trial-linking rule, for audit)
  <out>/outcomes.csv             INTERFACES.md columns
  <out>/post_freeze_updates.csv  per hypothesis, linked trials first submitted / last updated after
                                 the freeze (reported diagnostic; no label reads it)
  <out>/linkage_audit.csv        counts at each trial-linking step (diagnostic; ot_phase.linkage_audit)
  <out>/partial_date_boundary.csv  maturation statuses that would differ with partial dates on the
                                 last day of the stated month or year (diagnostic; diagnostics.py)
  <out>/phase23_diagnostic.csv   statuses and time-to-Phase-III fields that would differ if Phase
                                 II/III registrations counted as Phase III (diagnostic; diagnostics.py)
  <out>/run_info.json           run token, repository commit, code digest, ChEMBL cross-check run or
                                 not, freeze date, AACT source
  <out>/INPUTS.tsv, MANIFEST.tsv shared format (stages/run_guard/v8_manifest.py): inputs with
                                 repo- or volume-relative paths; one row per output file

Refuses to run unless PREREG.md passes `prereg check` and its log holds the entry
`RUN_START stage=C token=<token>` for the `--run-token` given and, before it, `SEAL stage=A
manifest_sha256=<sha256>` and the same for stage B (shared guard, stages/run_guard). The A seal
must equal the sha256 of the MANIFEST.tsv beside hypotheses.csv, and hypotheses.csv must have the
sha256 and row count that manifest lists. Stage B's seal is required in the log only; none of its
files is opened. The run records one repository commit (v8_run_guard.run_commit: the value baked
into the Modal image, or HEAD of a clean working tree) in run_info.json and in the `stage_code`
row of INPUTS.tsv.

Usage (from this directory):
  uv run python run_stage_c.py --run-token <token> [--inputs-root <dir>] [--aact-dir <dir>] \\
      [--chembl-dir <dir holding chembl_37.db and VERIFIED.json>]
"""
import argparse
import csv
import hashlib
import json
from collections.abc import Callable, Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd
from v8_manifest import (MANIFEST_NAME, RUN_INFO_NAME, InputRecord, code_record, code_sha256, guard_files,
                         portable_path, relative_files, verify_listed, write_manifest)
from v8_run_guard import require_run, run_commit

from chembl import chembl_max_phase, load_chembl_index, verified_chembl_db
from ctgov import AACT_DIR_NAME, AACT_SNAPSHOT_DATE, AactTrials, parse_study, verified_studies
from diagnostics import partial_date_boundary, phase23_diagnostic
from hypotheses import assert_not_stage_b, read_hypotheses
from models import (DIFFERENCE_COLUMNS, FREEZE_DIAGNOSTIC_COLUMNS, LINKAGE_AUDIT_COLUMNS, OUTCOME_COLUMNS,
                    OUTCOME_FREEZE_DATE, DifferenceRow, FreezeDiagnosticRow, LinkageAuditRow, OutcomeRow, StageCError,
                    StopCode)
from ot_phase import link_hypotheses, linkage_audit, load_ot
from rules import code_hypothesis, post_freeze_counts

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
REPO = EXP.parents[1]
OUT = HERE / "output"
DEFAULT_HYPOTHESES = HERE.parent / "A" / "output" / "hypotheses.csv"
OT_DIR_REL = Path("feasibility") / "v2_all_indications" / "inputs" / "ot_26.09"
CLINICAL_REPORT_REL = Path("inputs") / "ot_26.09" / "clinical_report.parquet"
PINS = EXP / "feasibility" / "v2_all_indications" / "coverage_all_indications.json"
PREREG = EXP / "PREREG.md"
DEFAULT_AACT = EXP / "inputs" / AACT_DIR_NAME


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_pins(ot_dir: Path, pins: Path) -> dict[str, tuple[Path, str]]:
    """sha256 of the Open Targets tables against coverage_all_indications.json; name -> (path, sha256)."""
    pinned = json.loads(pins.read_text())["inputs"]
    out = {}
    for name in ("clinical_indication", "drug_molecule"):
        p = ot_dir / f"{name}.parquet"
        got = sha256(p)
        if got != pinned[f"ot_{name}"]["sha256"]:
            raise StageCError(f"{p}: sha256 {got} differs from the pin in {pins.name}")
        out[f"ot_{name}"] = (p, got)
    return out


def fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, Enum):
        return v.value
    if isinstance(v, float):
        return f"{v:g}"
    return str(v)


def stage_code_sha256() -> str:
    """Digest of the code a stage C run executes: the modules of this directory and the shared
    guard. In the Modal image the guard modules are copied beside these; they are hashed under
    `run_guard/`, as in the repository, so the digest is the same in both places."""
    own = sorted(p for p in HERE.glob("*.py") if not p.name.startswith("v8_"))
    return code_sha256({**relative_files(HERE, own), **guard_files()})


@dataclass(frozen=True)
class Built:
    """What `build_rows` computes: the outcome rows, and the diagnostics no label reads."""
    rows: list[OutcomeRow]
    links: list[tuple[str, str]]
    post_freeze: list[FreezeDiagnosticRow]
    linkage_audit: list[LinkageAuditRow]
    partial_date_boundary: list[DifferenceRow]
    phase23: list[DifferenceRow]


def build_rows(hyps: pd.DataFrame, ot_dir: Path, clinical_report: Path, trials_source: AactTrials,
               chembl_db: Path | None) -> Built:
    ci, parent, report_nct = load_ot(ot_dir, clinical_report)
    linked = link_hypotheses(hyps, ci, parent, report_nct)
    trials_by_nct = trials_source.load(sorted({n for h in linked for n in h.nct_ids}))
    last_day_by_nct = {n: parse_study(r, last_day=True) for n, r in trials_source.raw.items()}
    audit = linkage_audit(hyps, ci, parent, report_nct, {n for n, t in trials_by_nct.items() if t is not None})
    chembl_index = load_chembl_index(chembl_db) if chembl_db is not None else None
    rows, links, diagnostics, boundary, phase23 = [], [], [], [], []
    for h in linked:
        trials = [t for t in (trials_by_nct[n] for n in h.nct_ids) if t is not None]
        o = code_hypothesis(list(h.max_stages), trials, OUTCOME_FREEZE_DATE)
        boundary += partial_date_boundary(h.hypothesis_id, list(h.max_stages), trials,
                                          [last_day_by_nct[t.nct_id] for t in trials], OUTCOME_FREEZE_DATE)
        phase23 += phase23_diagnostic(h.hypothesis_id, list(h.max_stages), trials, OUTCOME_FREEZE_DATE)
        links += [(h.hypothesis_id, n) for n in h.nct_ids]
        first_after, update_after = post_freeze_counts(trials, OUTCOME_FREEZE_DATE)
        diagnostics.append(FreezeDiagnosticRow(
            hypothesis_id=h.hypothesis_id, n_linked_trials=len(h.nct_ids),
            n_not_in_archive=len(h.nct_ids) - len(trials), n_first_submitted_after_freeze=first_after,
            n_last_update_after_freeze=update_after))
        rows.append(OutcomeRow(
            hypothesis_id=h.hypothesis_id,
            status_24=o.statuses[24], advanced_24=o.advanced[24],
            status_12=o.statuses[12], status_36=o.statuses[36], status_48=o.statuses[48],
            advanced_12=o.advanced[12], advanced_36=o.advanced[36], advanced_48=o.advanced[48],
            stop_code=o.stop_code, why_stopped_texts=o.why_stopped_texts,
            efficacy_coded=o.stop_code == StopCode.efficacy,
            phase2_start_date=o.phase2_start, phase3_start_date=o.phase3_start,
            time_to_phase3_days=o.time_to_phase3_days, phase3_event=o.phase3_event,
            phase3_success=o.phase3_success, approved=o.approved,
            chembl_max_phase_for_ind=(chembl_max_phase(chembl_index, h.drug_program_ids, h.indication_id)
                                      if chembl_index is not None else None),
            earliest_phase2_start=o.phase2_start, last_phase2_end_date=o.last_phase2_end))
    return Built(rows=rows, links=links, post_freeze=diagnostics, linkage_audit=audit, partial_date_boundary=boundary,
                 phase23=phase23)


def write_rows(rows: list, columns: list[str], path: Path) -> None:
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(columns)
        for r in rows:
            w.writerow([fmt(getattr(r, c)) for c in columns])


def write_outcomes(rows: list[OutcomeRow], path: Path) -> None:
    write_rows(rows, OUTCOME_COLUMNS, path)


def run(run_token: str, hypotheses: Path, inputs_root: Path, aact_dir: Path, chembl_dir: Path | None, out_dir: Path,
        prereg: Path = PREREG, pins: Path = PINS, commit: Callable[[], None] = lambda: None,
        roots: Sequence[tuple[str, Path]] = (("", REPO),), *, repo_commit: str) -> Path:
    """`roots` make INPUTS.tsv paths repo- or volume-relative (v8_manifest.portable_path).
    `repo_commit` is v8_run_guard.run_commit of the caller: the commit the run is recorded under."""
    ot_dir, clinical_report = inputs_root / OT_DIR_REL, inputs_root / CLINICAL_REPORT_REL
    for p in (hypotheses, ot_dir, clinical_report, aact_dir) + ((chembl_dir,) if chembl_dir else ()):
        assert_not_stage_b(p)
    require_run(prereg, "C", run_token, manifests={"A": hypotheses.parent / MANIFEST_NAME})
    hypotheses_sha = verify_listed(hypotheses.parent, [hypotheses.name])[hypotheses.name]
    studies, studies_sha = verified_studies(aact_dir)
    inputs = {"hypotheses": (hypotheses, hypotheses_sha), **check_pins(ot_dir, pins),
              "ot_clinical_report": (clinical_report, sha256(clinical_report)),
              "aact_studies": (studies, studies_sha)}
    chembl_db = None
    if chembl_dir is not None:
        chembl_db, chembl_sha = verified_chembl_db(chembl_dir)
        inputs["chembl_37"] = (chembl_db, chembl_sha)
    script = stage_code_sha256()
    records = [InputRecord(name=k, path=portable_path(p, roots), sha256=s) for k, (p, s) in inputs.items()]
    records.append(code_record("C", script, repo_commit))

    out_dir.mkdir(parents=True, exist_ok=True)
    trials_source = AactTrials(aact_dir, out_dir / "trials_snapshot")
    hyps = read_hypotheses(hypotheses)
    built = build_rows(hyps, ot_dir, clinical_report, trials_source, chembl_db)
    rows, links = built.rows, built.links
    commit()

    links_path = out_dir / "trial_links.csv"
    pd.DataFrame(links, columns=["hypothesis_id", "nct_id"]).to_csv(links_path, index=False)
    outcomes_path = out_dir / "outcomes.csv"
    write_outcomes(rows, outcomes_path)
    diagnostics_path = out_dir / "post_freeze_updates.csv"
    write_rows(built.post_freeze, FREEZE_DIAGNOSTIC_COLUMNS, diagnostics_path)
    audit_path, boundary_path, phase23_path = (out_dir / n for n in ("linkage_audit.csv", "partial_date_boundary.csv",
                                                                     "phase23_diagnostic.csv"))
    write_rows(built.linkage_audit, LINKAGE_AUDIT_COLUMNS, audit_path)
    write_rows(built.partial_date_boundary, DIFFERENCE_COLUMNS, boundary_path)
    write_rows(built.phase23, DIFFERENCE_COLUMNS, phase23_path)
    run_info = out_dir / RUN_INFO_NAME
    run_info.write_text(json.dumps({
        "stage": "C", "run_token": run_token, "repo_commit": repo_commit, "code_sha256": script,
        "chembl_cross_check": "run" if chembl_db else "not run (no ChEMBL directory given)",
        "outcome_freeze_date": OUTCOME_FREEZE_DATE.isoformat(),
        "clinicaltrials_gov_source": f"AACT flat files {AACT_SNAPSHOT_DATE}, studies.txt sha256 {studies_sha}"},
        indent=1))
    files = [outcomes_path, links_path, diagnostics_path, audit_path, boundary_path, phase23_path, run_info,
             trials_source.rows, trials_source.missing, trials_source.manifest]
    manifest = write_manifest(out_dir, files, script, records)
    commit()
    print(f"wrote {outcomes_path} ({len(rows)} hypotheses), {links_path}, {diagnostics_path}, {manifest}")
    return manifest


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-token", required=True, help="token of the PREREG.md entry 'RUN_START stage=C token=<token>'")
    ap.add_argument("--hypotheses", type=Path, default=DEFAULT_HYPOTHESES)
    ap.add_argument("--inputs-root", type=Path, default=EXP)
    ap.add_argument("--aact-dir", type=Path, default=DEFAULT_AACT)
    ap.add_argument("--chembl-dir", type=Path, default=None)
    ap.add_argument("--out-dir", type=Path, default=OUT)
    ap.add_argument("--prereg", type=Path, default=PREREG)
    a = ap.parse_args()
    run(a.run_token, a.hypotheses, a.inputs_root, a.aact_dir, a.chembl_dir, a.out_dir, prereg=a.prereg,
        repo_commit=run_commit(REPO))


if __name__ == "__main__":
    main()
