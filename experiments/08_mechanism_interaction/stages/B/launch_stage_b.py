"""Stage B launcher: plan, collect, analyze, status, assemble. Every step but `status` is part of the
stage B run and needs `prereg log "SEAL stage=A manifest_sha256=<sha256>"`, then
`prereg log "RUN_START stage=B token=<token>"`, a commit, and `modal deploy stages/B/modal_stage_b.py`
after it (the image bakes PREREG.md and the sealed stage A files).

    B="uv run --project stages/B --with modal==1.4.3 python stages/B/launch_stage_b.py"
    # 1. units, on Modal, from the sealed stage A files and the pinned instrument tables
    $B plan --run-token <token>          # writes B/output/unit_plan.json and inputs/stage_b/units.jsonl; commit the plan
    # 2. collect: one call per whole file (deCODE files, UKB-PPP tars and rsID maps, GWAS Catalog
    #    files without an index), each downloaded once to the stage B volume
    $B collect --run-token <token> [--source decode ...]
    #    re-collect GWAS Catalog accessions recorded absent under an earlier rule: their `absent` record
    #    is moved to /stage_b/superseded/collect/ and the task runs again (stage_b.collect.supersede_absent)
    $B collect --run-token <token> --accessions-file <file>   # or --accessions GCST... GCST...
    $B status                            # until every file has a record
    #    pre-analysis validation of the GWAS Catalog files read without a harmonised copy (stage_b/validate.py):
    #    aggregate parser diagnostics, and a file failing a fixed rule made `unreadable`
    $B validate --run-token <token>      # writes inputs/stage_b/validation/outcome_validation_<utc>.json
    #    a file made unreadable by a reader defect, after the reader is fixed: its collected record moved back
    #    (refused unless the bytes match the recorded sha256; nothing deleted), then `validate` again
    $B restore-validated --run-token <token> --accessions GCST... [--accessions-file <file>]
    #    the classification of the last report (files by uncertainty mode, failed, error) against an expected one;
    #    accessions and counts only, written to inputs/stage_b/validation/classification_check_<utc>.json
    $B check-validation --expected <json> [--report inputs/stage_b/validation/outcome_validation_<utc>.json]
    # 3. analyze: one call per instrument unit; refused while a file has no collect record, or while the
    #    validation report is missing or does not cover the current GWAS Catalog collect records
    $B spawn --run-token <token>
    $B status                            # until every unit is done
    # 4. after `modal volume get pqtl-v8-stage-b /stage_b/units inputs/stage_b/units` and
    #    `modal volume get pqtl-v8-stage-b /stage_b/collect inputs/stage_b/collect` (the gitignored
    #    inputs directory: the unit directories hold regional extracts, which are never committed)
    $B assemble --run-token <token> --units-dir inputs/stage_b/units --collect-dir inputs/stage_b/collect \
        --hypotheses stages/A/output/hypotheses.csv --st29 <Eldjarn 2023 MOESM3.xlsx>
    # 5. after `prereg log "SEAL stage=B manifest_sha256=<sha256>"`, a commit and a new deploy (the
    #    image bakes the log): delete the collected whole files; the records and extracts stay
    $B purge

`plan` reads no local table: the Modal function reads the stage A files baked into its image and
the five tables pinned in modal_inputs_manifest.json on the `pqtl-v8-inputs` volume, and returns
the plan (stage_b/plan.py). unit_plan.json records the sha256 of every file it read. units.jsonl
holds identities only (file names, sizes, ETags; no link, no token); it lives in the gitignored
inputs directory because B/output/ holds only what MANIFEST.tsv lists and unit_plan.json. A plan
already on disk that differs is moved to inputs/stage_b/superseded/ first.

The deCODE folder link is a credential: its token is in the Modal secret `pqtl-stage-b`
(DECODE_FOLDER_TOKEN) and nowhere else. `collect` can be given `--source decode` to fetch the
deCODE files first, while a link is valid; a link re-issued later changes no record, no unit and
no fingerprint.

`plan`, `collect`, `spawn` and `assemble` refuse to start unless `git status --porcelain` is empty
for stages/ and PREREG.md (v8_run_guard.clean_commit; commit unit_plan.json after `plan`, and
launch_log.jsonl after `collect` and `spawn`) and PREREG.md authorizes the run for the `--run-token`
given (stage_b.launch.authorize, the shared guard in stages/run_guard); each Modal call makes the
same check on its baked copies. `collect` and `spawn` refuse units planned from files other than
the sealed ones. `assemble` accepts a unit directory only under the fingerprint of this run
(stage_b/checkpoint.py), recomputed from the unit record, the pins, the code and the plan held here,
the tool versions the unit recorded and the collect digest recomputed from the collect records
given (so a unit computed from other bytes than the recorded ones is refused), only when every unit
ran under the same tool versions, and only when every planned file has a collect record; it writes
the commit to run_info.json and INPUTS.tsv and the collect records to collected_files.tsv. It
writes B/output/ only: evidence.csv, regional_manifest.tsv, collected_files.tsv, run_info.json,
INPUTS.tsv and MANIFEST.tsv, beside the unit_plan.json `plan` left there. No regional extract and
no whole file is written to B/output/; `--units-dir` and `--collect-dir` are refused inside the
repository unless under the gitignored inputs directory (stage_b.launch.private_copy).

`status` writes inputs/stage_b/status/status_<utc>.json and prints its summary: per source, files
collected, absent, pending and failed (with the error class), and units done, pending and failed.
It asks Modal for the state of the deployed app and of every call the launcher spawned, so work
that stopped because the app is gone is reported as stranded, not as pending.
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import modal

from v8_manifest import MANIFEST_NAME, sha256_file, verify_output_dir
from v8_run_guard import PLAN_SHA256, clean_commit

from stage_b.assemble import (build_evidence, collect_unit_dir, load_collect_records, load_st29, regional_rows,
                              write_outputs)
from stage_b.checkpoint import collect_digest, common_tools, package_sha256, source_pins, unit_fingerprint, unit_tools
from stage_b.launch import (authorize, check_plan_inputs, planned_tasks, private_copy, sealed_stage_b, spawn_collect,
                            spawn_units)
from stage_b.schemas import CollectRecord, InputContractError, InstrumentUnit
from stage_b.status import error_of, status_report
from stage_b.units import load_hypotheses
from stage_b.validate import check_classification

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
REPO = EXP.parents[1]
PREREG = EXP / "PREREG.md"
A_OUTPUT = HERE.parent / "A" / "output"
INPUTS_MANIFEST = EXP / "modal_inputs_manifest.json"
WORK = EXP / "inputs" / "stage_b"
UNITS = WORK / "units.jsonl"
PLAN_JSON = HERE / "output" / "unit_plan.json"
LOG = HERE / "launch_log.jsonl"
APP = "pqtl-v8-stage-b"
VOLUME = "pqtl-v8-stage-b"
STAGE_ROOT = "stage_b"     # the stage B root on the volume (modal_stage_b.ROOT without the mount point)
FUNCTIONS = ("plan_remote", "collect_file", "validate_outcomes", "validate_outcome_file", "restore_validated_records",
             "run_unit", "stage_state", "purge_raw_files")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def log_entry(entry: dict) -> None:
    with LOG.open("a") as fh:
        fh.write(json.dumps(entry) + "\n")


def cmd_plan(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    authorize(PREREG, a.run_token, A_OUTPUT)
    out = modal.Function.from_name(APP, "plan_remote").remote(a.run_token)
    units_text, plan = out["units_jsonl"], json.dumps(out["unit_plan"], indent=1, sort_keys=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    for path, text in ((UNITS, units_text), (PLAN_JSON, plan)):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and path.read_text() != text:
            old = WORK / "superseded" / f"{path.stem}_{stamp}{path.suffix}"
            old.parent.mkdir(parents=True, exist_ok=True)
            path.rename(old)
        path.write_text(text)
    if sha256_file(UNITS) != out["unit_plan"]["units_sha256"]:
        raise InputContractError(f"{UNITS} is not the units file the plan recorded")
    log_entry({"what": "plan", "repo_commit": commit, "units_sha256": out["unit_plan"]["units_sha256"], "utc": utc_now()})
    print(f"wrote {UNITS} ({out['unit_plan']['units']} units: {out['unit_plan']['units_by_source']}) and {PLAN_JSON}; "
          f"{len(out['unit_plan']['unresolved'])} hypotheses without a regional file; "
          f"files to collect: {out['unit_plan']['collect_tasks_by_source']}")


def accession_list(a: argparse.Namespace) -> list[str]:
    """--accessions, and the non-empty lines of --accessions-file (first tab-separated field)."""
    listed = list(a.accessions or [])
    if a.accessions_file is not None:
        listed += [line.split("\t")[0].strip() for line in a.accessions_file.read_text().splitlines() if line.strip()]
    return sorted(set(listed))


def cmd_collect(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    fn = modal.Function.from_name(APP, "collect_file")
    accessions = accession_list(a)
    supersede = bool(accessions)          # only an `absent` GWAS Catalog record is ever superseded
    calls = spawn_collect(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT,
                          lambda task, token: fn.spawn(task, token, supersede).object_id, a.source or (), accessions)
    log_entry({"what": "collect", "repo_commit": commit, "calls": calls, "supersede_absent": supersede,
               "accessions": accessions, "utc": utc_now()})
    print(f"{len(calls)} collect calls spawned at {utc_now()}" + (f" (supersede absent records of {len(accessions)} "
                                                                  "accessions)" if supersede else ""))


# The pre-analysis validation on Modal. Files already validated are not read again, so an interrupted
# call is resumed by running the command again. The guard is the first statement (tests/test_launch_guard.py).
def cmd_validate(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    authorize(PREREG, a.run_token, A_OUTPUT)
    call = modal.Function.from_name(APP, "validate_outcomes").spawn(a.run_token)
    log_entry({"what": "validate", "repo_commit": commit, "call_id": call.object_id, "utc": utc_now()})
    report = call.get()
    out = WORK / "validation" / f"outcome_validation_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print(f"{report['files_validated']} files validated: {report['passed']} passed, failed {report['failed']}, "
          f"errors {report['errors']}; unreadable at collect {report['unreadable_at_collect']}; full report: {out}")


# Undo the validation's `unreadable` for accessions whose failure was a reader defect (stage_b/validate.py
# `restorable`): refused, nothing moved, unless every file's bytes still match its recorded sha256. Run
# `validate` afterwards. The guard is the first statement (tests/test_launch_guard.py).
def cmd_restore_validated(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    authorize(PREREG, a.run_token, A_OUTPUT)
    accessions = accession_list(a)
    if not accessions:
        raise InputContractError("restore-validated needs --accessions or --accessions-file")
    out = modal.Function.from_name(APP, "restore_validated_records").remote(accessions, a.run_token)
    log_entry({"what": "restore_validated", "repo_commit": commit, "accessions": accessions, **out, "utc": utc_now()})
    print(f"{len(out['restored'])} records restored to collected: {[r['key'] for r in out['restored']]}; "
          "run `validate` again before `spawn`")


def cmd_check_validation(a: argparse.Namespace) -> None:
    """The classification of a validation report written by `validate` (the newest unless --report)
    against the expected one (--expected: class -> accessions, stage_b.validate.CLASSES); refuses on
    any difference. Reads local files only."""
    reports = sorted((WORK / "validation").glob("outcome_validation_*.json"))
    report = a.report or (reports[-1] if reports else None)
    if report is None:
        raise InputContractError(f"no validation report under {WORK / 'validation'}; run `validate` first")
    got = check_classification(json.loads(report.read_text()), json.loads(a.expected.read_text()))
    out = WORK / "validation" / f"classification_check_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
    out.write_text(json.dumps({"report": str(report), "report_sha256": sha256_file(report), "expected": str(a.expected),
                               "expected_sha256": sha256_file(a.expected), **got, "utc": utc_now()}, indent=1))
    print(f"counts {got['counts']}; expected {got['expected_counts']}; written to {out}")
    if not got["agrees"]:
        raise InputContractError(f"the validation report differs from the expected classification: {got['mismatches']}")


def cmd_spawn(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    state = modal.Function.from_name(APP, "stage_state").remote()
    fn = modal.Function.from_name(APP, "run_unit")
    calls = spawn_units(PREREG, a.run_token, UNITS, PLAN_JSON, A_OUTPUT, lambda line, token: fn.spawn(line, token).object_id,
                        recorded={(r["source"], r["key"]) for r in state["records"]},
                        validation=(state["outcome_validation"],
                                    [CollectRecord.model_validate(r) for r in state["gwas_catalog_records"]]))
    log_entry({"what": "spawn", "repo_commit": commit, "calls": calls, "utc": utc_now()})
    print(f"{len(calls)} units spawned at {utc_now()}")


def app_state() -> dict:
    """Whether the deployed app answers, and each function's backlog and running containers."""
    functions = {}
    for name in FUNCTIONS:
        try:
            stats = modal.Function.from_name(APP, name).get_current_stats()
            functions[name] = {"backlog": stats.backlog, "runners": stats.num_total_runners}
        except modal.exception.NotFoundError:
            return {"deployed": False, "functions": functions}
    return {"deployed": True, "functions": functions}


def volume_state_without_the_app() -> dict:
    """stage_b.status.volume_state read through the volume API, for when the app is not deployed and
    its `stage_state` function cannot be called. One listing, then one read per small JSON record."""
    vol = modal.Volume.from_name(VOLUME)
    files = {e.path: e for e in vol.listdir(f"/{STAGE_ROOT}", recursive=True) if e.type == modal.volume.FileEntryType.FILE}

    def read(path: str) -> dict:
        return json.loads(b"".join(vol.read_file(path)))

    def under(prefix: str, suffix: str) -> list[str]:
        return sorted(p for p in files if p.startswith(f"{STAGE_ROOT}/{prefix}/") and p.endswith(suffix))

    records = [{k: read(p)[k] for k in ("source", "key", "status")} for p in under("collect", ".json")]
    partials = [{"source": Path(p).parent.name, "name": Path(p).name, "bytes": files[p].size} for p in under("raw", ".part")]
    in_units = [Path(p).relative_to(f"{STAGE_ROOT}/units").parts for p in under("units", "")]
    errors = {kind: {m["name"]: m for m in (read(p) for p in under(f"errors/{kind}", ".json"))} for kind in ("collect", "units")}
    return {"records": records, "partials": partials, "errors": errors,
            "units_done": sorted({parts[0] for parts in in_units if parts[1:] == ("result.json",)}),
            "units_started": sorted({parts[0] for parts in in_units})}


def call_state(call_id: str) -> dict:
    """What Modal says of one spawned call: done, failed (with the error class) or pending."""
    try:
        modal.FunctionCall.from_id(call_id).get(timeout=0)
    except TimeoutError:
        return {"call_id": call_id, "state": "pending"}
    except Exception as err:   # whatever the call raised, remotely or when its result was fetched
        return {"call_id": call_id, "state": "failed", **error_of(err)}
    return {"call_id": call_id, "state": "done"}


def latest_calls() -> dict[str, dict[str, str]]:
    """The last call id the launch log holds for each file (`collect`) and unit (`units`)."""
    out: dict[str, dict[str, str]] = {"collect": {}, "units": {}}
    for line in LOG.read_text().splitlines() if LOG.exists() else []:
        entry = json.loads(line)
        for c in entry.get("calls", []):
            if entry["what"] == "collect":
                out["collect"][f"{c['source']}/{c['key']}"] = c["call_id"]
            elif entry["what"] == "spawn":
                out["units"][c["unit_key"]] = c["call_id"]
    return out


def cmd_status(a: argparse.Namespace) -> None:
    lines = UNITS.read_text().splitlines()
    units = [InstrumentUnit.model_validate_json(line) for line in lines]
    app = app_state()
    volume = modal.Function.from_name(APP, "stage_state").remote() if app["deployed"] else volume_state_without_the_app()
    finished_files = {f"{r['source']}/{r['key']}" for r in volume["records"]}
    finished = {"collect": finished_files, "units": set(volume["units_done"])}
    calls = {kind: {name: call_state(cid) for name, cid in ids.items() if name not in finished[kind]}
             for kind, ids in latest_calls().items()}
    report = {"utc": utc_now(), **status_report(planned_tasks(lines), units, volume, calls, app)}
    out = WORK / "status" / f"status_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(report, indent=1))
    print(json.dumps({k: report[k] for k in ("app", "files_by_source", "units_by_source", "failures_by_error_class")}, indent=1))
    print(f"partial downloads: {len(report['partial_downloads'])}; full report: {out}")
    if not app["deployed"]:
        print(f"APP {APP} IS NOT DEPLOYED: nothing that is not finished can make progress (stranded above).")


def cmd_assemble(a: argparse.Namespace) -> None:
    commit = clean_commit(REPO)
    sealed = authorize(PREREG, a.run_token, A_OUTPUT)
    check_plan_inputs(PLAN_JSON, UNITS, sealed)
    if sha256_file(a.hypotheses) != sealed["hypotheses.csv"]:
        raise InputContractError(f"{a.hypotheses} is not the hypotheses.csv stage A sealed")
    hyps = load_hypotheses(a.hypotheses)
    plan = json.loads(PLAN_JSON.read_text())
    lines = UNITS.read_text().splitlines()
    units = [InstrumentUnit.model_validate_json(line) for line in lines]
    units_dir = private_copy(a.units_dir, REPO, EXP / "inputs")
    collected = load_collect_records(private_copy(a.collect_dir, REPO, EXP / "inputs"))
    records = {(r.source, r.key): r for r in collected}
    missing = {(t.source, t.key) for t in planned_tasks(lines)} - set(records)
    if missing:
        raise InputContractError(f"{len(missing)} planned files have no collect record in {a.collect_dir} "
                                 f"(first: {sorted(missing)[0]}); stage B is not complete")
    pins, code = source_pins(INPUTS_MANIFEST), package_sha256(HERE / "stage_b")
    results, metas, tools = {}, {}, {}
    for u in units:
        tools[u.unit_key] = unit_tools(units_dir / u.unit_key)
        digest = collect_digest(u, lambda source, key: records[(source, key)])
        results[u.unit_key], metas[u.unit_key] = collect_unit_dir(
            u, units_dir / u.unit_key, unit_fingerprint(u, pins, code, PLAN_SHA256, tools[u.unit_key], digest), digest)
    unit_ids = {k: tuple(v) for k, v in plan["unit_ids"].items()}
    rows = build_evidence(hyps, plan["hypothesis_unit"], unit_ids, results, load_st29(a.st29),
                          plan["hypothesis_source_units"])
    scripts = sorted((HERE / "stage_b").glob("*.py")) + [HERE / "stage_b" / "coloc_run.R", Path(__file__).resolve()]
    inputs = {"hypotheses": a.hypotheses, "st29_workbook": a.st29, "unit_plan": PLAN_JSON, "units": UNITS}
    manifest = write_outputs(HERE / "output", rows, regional_rows(metas), scripts, inputs, [("", REPO)],
                             script_root=HERE, run_token=a.run_token, repo_commit=commit, tools=common_tools(tools),
                             collected=collected)
    print(f"wrote {HERE / 'output' / 'evidence.csv'} ({len(rows)} rows) and {manifest}")


def cmd_purge(a: argparse.Namespace) -> None:
    """Delete the whole files on the volume, once stage B's output here is the sealed, intact one."""
    commit = clean_commit(REPO)
    seal = sealed_stage_b(PREREG, HERE / "output" / MANIFEST_NAME)
    verify_output_dir(HERE / "output", ["evidence.csv", "regional_manifest.tsv", "collected_files.tsv"], ["unit_plan.json"])
    out = modal.Function.from_name(APP, "purge_raw_files").remote(UNITS.read_text())
    log_entry({"what": "purge", "repo_commit": commit, "stage_b_seal": seal, **out, "utc": utc_now()})
    print(f"deleted {out['files_deleted']} whole files ({out['bytes_deleted']} bytes); {out['records_kept']} records kept")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    token_help = "token of the PREREG.md entry 'RUN_START stage=B token=<token>'"
    for name in ("plan", "collect", "validate", "restore-validated", "spawn", "assemble"):
        sub.add_parser(name).add_argument("--run-token", required=True, help=token_help)
    sub.choices["collect"].add_argument("--source", nargs="*", help="collect only these sources (e.g. decode)")
    sub.choices["collect"].add_argument("--accessions", nargs="*",
                                        help="collect only these GWAS Catalog accessions again, superseding an `absent` record")
    sub.choices["restore-validated"].add_argument("--accessions", nargs="*",
                                                  help="GWAS Catalog accessions made unreadable by a reader defect")
    for name in ("collect", "restore-validated"):
        sub.choices[name].add_argument("--accessions-file", type=Path,
                                       help="a file of accessions, one per line (first tab-separated field), as --accessions")
    for name in ("hypotheses", "units-dir", "collect-dir", "st29"):
        sub.choices["assemble"].add_argument(f"--{name}", type=Path, required=True)
    check = sub.add_parser("check-validation")
    check.add_argument("--expected", type=Path, required=True, help="JSON: class -> accessions (stage_b.validate.CLASSES)")
    check.add_argument("--report", type=Path, help="a report `validate` wrote; the newest when not given")
    sub.add_parser("status")
    sub.add_parser("purge")
    a = ap.parse_args()
    {"plan": cmd_plan, "collect": cmd_collect, "validate": cmd_validate, "restore-validated": cmd_restore_validated,
     "check-validation": cmd_check_validation, "spawn": cmd_spawn, "status": cmd_status, "assemble": cmd_assemble,
     "purge": cmd_purge}[a.cmd](a)


if __name__ == "__main__":
    main()
