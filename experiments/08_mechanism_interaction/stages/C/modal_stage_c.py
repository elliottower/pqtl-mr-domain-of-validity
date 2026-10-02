"""Thin Modal wrapper for stage C (logic in run_stage_c.py, aact_download.py and the modules beside it).

Volumes:
  /inputs  pqtl-v8-inputs, read-only for the run: /inputs/08_mechanism_interaction mirrors the
           experiment directory (Open Targets 26.09 tables under feasibility/v2_all_indications/
           inputs/ot_26.09 and inputs/ot_26.09/clinical_report.parquet; scripts/upload_inputs_to_modal.py),
           and /inputs/aact_20260930/ holds the AACT daily flat-file snapshot of 2026-09-30
           (20260930_export_ctgov.zip, uploaded), its extracted studies.txt and VERIFIED.json,
           written by `verify_uploaded_aact` (or `download_aact`), the only writable mounts.
  /chembl  proteome-mr-claim-audit-inputs, read-only: /chembl/chembl_37/chembl_37.db, opened only
           after its sha256 matches /chembl/chembl_37/VERIFIED.json.
  /vol     pqtl-v8-stage-c: output/ (outcomes.csv, trial_links.csv, post_freeze_updates.csv,
           trials_snapshot/, MANIFEST.tsv).

Baked into the image: this directory's modules, the shared run guard (stages/run_guard/
v8_run_guard.py, v8_manifest.py, v8_test_report.py), PREREG.md, coverage_all_indications.json (the Open Targets pins),
A/output/hypotheses.csv and A/output/MANIFEST.tsv (the sealed stage A manifest the guard checks
hypotheses.csv against). No stage B file is baked.

1. Once. The daily snapshot of 2026-09-30 was downloaded by the author through a signed-in AACT
   session and uploaded to pqtl-v8-inputs:/aact_20260930/; AACT deletes daily snapshots at month
   end, so that copy is the record. Verify it against the sha256 recorded at download:
    uv run --project experiments/08_mechanism_interaction/stages/C --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/C/modal_stage_c.py::verify_uploaded_aact \\
        --name 20260930_export_ctgov.zip --expected-sha256 <sha256>
   (`download_aact --url <link>` is the URL route to the same directory and the same check.)
2. Real run only after `SEAL stage=A manifest_sha256=<sha256>`, `SEAL stage=B manifest_sha256=<sha256>`
   and then `RUN_START stage=C token=<token>` are logged below the line of PREREG.md, and the stage
   code and PREREG.md are committed: the launcher refuses unless `git status --porcelain` is empty
   for stages/ and PREREG.md (v8_run_guard.clean_commit), and the image carries that commit as
   V8_REPO_COMMIT, which the run writes to run_info.json and INPUTS.tsv.
   The guard runs locally before the call and again in the container on the baked PREREG.md, so
   the image must be built after the entry is logged:
    uv run --project experiments/08_mechanism_interaction/stages/C --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/C/modal_stage_c.py --run-token <token>
    modal volume get pqtl-v8-stage-c /output experiments/08_mechanism_interaction/stages/C/output

Wiring probe (no study data; allowed before the run is logged), in the run image with the volumes of
`run_remote`, so a start of it is a container start of the real run's configuration. It imports
every module the run imports and reports the baked files (PREREG.md, the Open Targets pins, stage
A's hypotheses.csv against its sealed manifest), the frozen-plan check of the shared guard, and the
mounted volumes: top-level listings, the Open Targets tables by size and column names, the AACT
snapshot by VERIFIED.json, size and the header line of studies.txt, and the ChEMBL database by
VERIFIED.json, size and its 16-byte file header. It hashes nothing on a volume, reads no trial row
and never calls `require_run`. The report is printed and written to /vol/probe/probe.json on
pqtl-v8-stage-c, the previous one moved to /vol/probe/superseded/ first:
    uv run --project experiments/08_mechanism_interaction/stages/C --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/C/modal_stage_c.py::probe

Synthetic test suite (no study data), in the run image plus pytest, with tests/ and stages/run_guard
baked beside the modules at /root/exp/stages/. The summary is printed and written to
/vol/tests/modal_tests.json on pqtl-v8-stage-c, the previous report moved to /vol/tests/superseded/
first:
    uv run --project experiments/08_mechanism_interaction/stages/C --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/C/modal_stage_c.py::tests
"""
import json
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

import modal
import pandas
import pyarrow
import pyarrow.parquet
import pydantic
import requests
import v8_manifest
import v8_run_guard
import v8_test_report
from v8_manifest import MANIFEST_NAME, sha256_file, verify_listed
from v8_run_guard import (COMMIT_ENV, PLAN_SHA256, RUN_START_RE, baked_commit, check_plan, clean_commit, logged_seals,
                          parse_log, require_run, run_commit)
from v8_test_report import run_pytest, write_report

import aact_download
import chembl
import ctgov
import diagnostics
import hypotheses
import models
import ot_phase
import rules
import run_stage_c
from aact_download import UPLOADED_SOURCE, download_archive, verify_archive
from chembl import DB_NAME, pinned_sha256
from ctgov import AACT_DIR_NAME, AACT_SNAPSHOT_DATE, STUDIES_COLUMNS, STUDIES_FILE, VERIFIED_NAME
from run_stage_c import CLINICAL_REPORT_REL, OT_DIR_REL, run, stage_code_sha256

HERE = Path(__file__).resolve().parent
A_OUTPUT = HERE.parent / "A" / "output"
A_FILES = ("hypotheses.csv", "MANIFEST.tsv")
ROOT = Path("/root")
INPUTS = Path("/inputs/08_mechanism_interaction")
AACT = Path("/inputs") / AACT_DIR_NAME
CHEMBL = Path("/chembl/chembl_37")
OUT = Path("/vol/output")
PREREG_IMAGE = ROOT / "exp" / "PREREG.md"
PY_PACKAGES = ["pandas==2.2.3", "pyarrow==21.0.0", "pydantic==2.11.7", "requests==2.32.3", "prereg==0.4.2",
               "provenance-core==0.4.2"]
TEST_PACKAGES = ["pytest==8.4.2"]
TEST_IGNORE = [".venv", "**/__pycache__", "**/.pytest_cache", "output", "work", "launch_log.jsonl"]
STAGES_IMAGE = ROOT / "exp" / "stages"
A_IMAGE = ROOT / "stages" / "A" / "output"
PINS_IMAGE = ROOT / "exp" / "coverage_all_indications.json"
PROBE_REPORT = Path("/vol/probe/probe.json")
RUN_MODULES = (aact_download, chembl, ctgov, diagnostics, hypotheses, models, ot_phase, rules, run_stage_c, v8_manifest,
               v8_run_guard, v8_test_report, pandas, pyarrow, pydantic, requests)
OT_TABLES = ("clinical_indication", "drug_molecule")
SQLITE_MAGIC = b"SQLite format 3\x00"

image = modal.Image.debian_slim(python_version="3.12").pip_install(*PY_PACKAGES).env({"PYTHONPATH": "/root/exp/stages/C"})
test_image = image
if modal.is_local():
    image = (image.add_local_dir(str(HERE), "/root/exp/stages/C", copy=True,
                                 ignore=["output", "tests", "__pycache__", ".venv", ".pytest_cache"])
             .add_local_file(str(HERE.parent / "run_guard" / "v8_run_guard.py"), "/root/exp/stages/C/v8_run_guard.py", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_manifest.py"), "/root/exp/stages/C/v8_manifest.py", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_test_report.py"), "/root/exp/stages/C/v8_test_report.py",
                             copy=True)
             .add_local_file(str(HERE.parents[1] / "PREREG.md"), str(PREREG_IMAGE), copy=True)
             .add_local_file(str(HERE.parents[1] / "feasibility" / "v2_all_indications" / "coverage_all_indications.json"),
                             "/root/exp/coverage_all_indications.json", copy=True))
    for _name in A_FILES:
        if (A_OUTPUT / _name).is_file():
            image = image.add_local_file(str(A_OUTPUT / _name), f"/root/stages/A/output/{_name}", copy=True)
    image = image.env({COMMIT_ENV: baked_commit(HERE.parents[3])})
    test_image = (image.pip_install(*TEST_PACKAGES)
                  .env({"PYTHONPATH": f"/root/exp/stages/C:{STAGES_IMAGE / 'run_guard'}"})
                  .add_local_dir(str(HERE / "tests"), str(STAGES_IMAGE / "C" / "tests"), copy=True, ignore=["**/__pycache__"])
                  .add_local_file(str(HERE / "pyproject.toml"), str(STAGES_IMAGE / "C" / "pyproject.toml"), copy=True)
                  .add_local_dir(str(HERE.parent / "run_guard"), str(STAGES_IMAGE / "run_guard"), copy=True,
                                 ignore=TEST_IGNORE)
                  .add_local_file(str(HERE.parent / "INTERFACES.md"), str(STAGES_IMAGE / "INTERFACES.md"), copy=True))

app = modal.App("pqtl-v8-stage-c", image=image)
inputs_vol = modal.Volume.from_name("pqtl-v8-inputs")
chembl_vol = modal.Volume.from_name("proteome-mr-claim-audit-inputs").with_mount_options(read_only=True)
vol = modal.Volume.from_name("pqtl-v8-stage-c", create_if_missing=True)
# What `run_remote` mounts. `probe` takes the same image and volumes.
RUN_VOLUMES = {"/inputs": inputs_vol.with_mount_options(read_only=True), "/chembl": chembl_vol, "/vol": vol}


@app.function(cpu=2, memory=4096, timeout=24 * 3600, retries=3, volumes={"/inputs": inputs_vol})
def download_aact(url: str, expected_sha256: str = "") -> dict:
    """AACT archive -> /inputs/aact_20260930/ with VERIFIED.json; resumes a partial download."""
    inputs_vol.reload()
    return download_archive(url, AACT, commit=inputs_vol.commit, expected_sha256=expected_sha256 or None)


@app.function(cpu=2, memory=4096, timeout=24 * 3600, retries=1, volumes={"/inputs": inputs_vol})
def verify_uploaded_aact(name: str, expected_sha256: str) -> dict:
    """/inputs/aact_20260930/<name>, uploaded -> studies.txt and VERIFIED.json beside it; the zip is never deleted."""
    inputs_vol.reload()
    return verify_archive(AACT / name, AACT, UPLOADED_SOURCE, expected_sha256, commit=inputs_vol.commit)


@app.function(cpu=2, memory=8192, timeout=24 * 3600, retries=2, volumes=RUN_VOLUMES)
def run_remote(run_token: str) -> str:
    vol.reload()
    manifest = run(run_token, A_IMAGE / "hypotheses.csv", INPUTS, AACT, CHEMBL, OUT,
                   prereg=PREREG_IMAGE, pins=PINS_IMAGE, commit=vol.commit,
                   roots=[("experiments/08_mechanism_interaction/stages", ROOT / "stages"), ("pqtl-v8-inputs:", Path("/inputs")),
                          ("proteome-mr-claim-audit-inputs:", Path("/chembl"))], repo_commit=run_commit(None))
    vol.commit()
    return str(manifest)


def _attempt(check: Callable[[], object]) -> dict:
    """One probe check: its value, or the exception it raised, so one failure does not hide the rest."""
    try:
        return {"ok": True, "value": check()}
    except Exception as err:  # a probe reports every failure, whatever its type
        return {"ok": False, "error": f"{type(err).__name__}: {err}"[:2000]}


def _file(path: Path, sha256: bool = False) -> dict:
    if not path.is_file():
        return {"path": str(path), "exists": False}
    return {"path": str(path), "exists": True, "bytes": path.stat().st_size,
            **({"sha256": sha256_file(path)} if sha256 else {})}


def _top_level(path: Path) -> dict:
    if not path.is_dir():
        return {"path": str(path), "exists": False}
    return {"path": str(path), "exists": True, "entries": sorted(p.name + ("/" if p.is_dir() else "") for p in path.iterdir())}


def _first_line(path: Path) -> str:
    with path.open(encoding="utf-8") as f:
        return f.readline().rstrip("\r\n")


def _magic(path: Path) -> bool:
    with path.open("rb") as f:
        return f.read(len(SQLITE_MAGIC)) == SQLITE_MAGIC


def _log_state() -> dict:
    """The baked PREREG.md under the guard's frozen-plan check (not its authorization), the seals
    its log holds, and the stages with a RUN_START entry."""
    text = PREREG_IMAGE.read_text(encoding="utf-8")
    check_plan(text)
    entries = parse_log(text)
    starts = sorted({m.group(1) for m in (RUN_START_RE.fullmatch(e.event) for e in entries) if m is not None})
    return {"plan_sha256": PLAN_SHA256, "log_entries": len(entries), "seals": logged_seals(entries, len(entries) + 1),
            "run_start_stages": starts}


@app.function(cpu=2, memory=8192, timeout=1800, volumes=RUN_VOLUMES)
def probe() -> dict:
    """WIRING ONLY, in the container configuration of `run_remote`: imports, baked files, mounted
    volumes, package versions. Reads no trial row, hashes no volume file, never calls `require_run`."""
    vol.reload()
    ot_dir = INPUTS / OT_DIR_REL
    tables = {**{name: ot_dir / f"{name}.parquet" for name in OT_TABLES}, "clinical_report": INPUTS / CLINICAL_REPORT_REL}
    aact_verified = _attempt(lambda: json.loads((AACT / VERIFIED_NAME).read_text()))
    chembl_verified = _attempt(lambda: json.loads((CHEMBL / VERIFIED_NAME).read_text()))
    aact = aact_verified.get("value", {})
    studies, archive = _file(AACT / STUDIES_FILE), _file(AACT / aact.get("archive", {}).get("name", "archive name unknown"))
    header = _attempt(lambda: _first_line(AACT / STUDIES_FILE).split("|"))
    report = {
        "stage": "C", "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), COMMIT_ENV: _attempt(lambda: run_commit(None)),
        "imports": {m.__name__: {"file": str(m.__file__), "version": str(getattr(m, "__version__", ""))} for m in RUN_MODULES},
        "baked": {"prereg": _file(PREREG_IMAGE, sha256=True), "pins": _file(PINS_IMAGE, sha256=True),
                  "pinned_tables": _attempt(lambda: sorted(json.loads(PINS_IMAGE.read_text())["inputs"])),
                  "stage_a": {n: _file(A_IMAGE / n, sha256=True) for n in A_FILES},
                  "stage_a_hypotheses_match_manifest": _attempt(lambda: verify_listed(A_IMAGE, ["hypotheses.csv"])),
                  "stage_code_sha256": _attempt(stage_code_sha256)},
        "prereg_log": _attempt(_log_state),
        "volumes": {
            "/inputs": _top_level(Path("/inputs")), "ot_dir": _top_level(ot_dir),
            "ot_tables": {name: {**_file(path), "columns": _attempt(lambda path=path: pyarrow.parquet.read_schema(path).names)}
                          for name, path in tables.items()},
            "aact_dir": _top_level(AACT), "aact_verified": aact_verified, "studies": studies, "studies_header": header,
            "aact_archive": archive,
            "/chembl": _top_level(Path("/chembl")), "chembl_dir": _top_level(CHEMBL), "chembl_verified": chembl_verified,
            "chembl_db": _file(CHEMBL / DB_NAME), "chembl_db_is_sqlite": _attempt(lambda: _magic(CHEMBL / DB_NAME)),
            "chembl_pins": _attempt(lambda: sorted(pinned_sha256(chembl_verified["value"]))),
            "/vol": _top_level(Path("/vol")), "output": _top_level(OUT)},
    }
    log = report["prereg_log"].get("value", {})
    seal = report["baked"]["stage_a"][MANIFEST_NAME].get("sha256")
    pinned = aact.get("files", {}).get(STUDIES_FILE, {})
    vols = report["volumes"]
    problems = [f"baked file missing: {r['path']}" for r in (report["baked"]["prereg"], report["baked"]["pins"],
                                                              *report["baked"]["stage_a"].values()) if not r["exists"]]
    problems += [f"{k}: {report['baked'][k]['error']}" for k in ("pinned_tables", "stage_a_hypotheses_match_manifest",
                                                                 "stage_code_sha256") if not report["baked"][k]["ok"]]
    if report["baked"]["pinned_tables"]["ok"]:
        problems += [f"{PINS_IMAGE.name} pins no ot_{name}" for name in OT_TABLES
                     if f"ot_{name}" not in report["baked"]["pinned_tables"]["value"]]
    if not report["prereg_log"]["ok"]:
        problems.append(f"frozen-plan check: {report['prereg_log']['error']}")
    elif log["seals"].get("A") != seal:
        problems.append(f"baked stage A MANIFEST.tsv sha256 {seal} is not the logged seal {log['seals'].get('A')}")
    problems += [f"Open Targets table {name}: {rec}" for name, rec in vols["ot_tables"].items()
                 if not rec["exists"] or not rec["columns"]["ok"]]
    if not aact_verified["ok"] or aact.get("snapshot_date") != AACT_SNAPSHOT_DATE:
        problems.append(f"{AACT / VERIFIED_NAME}: {aact_verified.get('error', aact.get('snapshot_date'))}")
    if not studies["exists"] or studies["bytes"] != pinned.get("bytes"):
        problems.append(f"{STUDIES_FILE}: on the volume {studies.get('bytes')} bytes, VERIFIED.json records {pinned.get('bytes')}")
    if not header["ok"] or not set(STUDIES_COLUMNS) <= set(header["value"]) or header["value"][0] != "nct_id":
        problems.append(f"{STUDIES_FILE} header: {header}")
    if not archive["exists"] or archive["bytes"] != aact.get("archive", {}).get("bytes"):
        problems.append(f"AACT archive: {archive}, VERIFIED.json records {aact.get('archive')}")
    if not vols["chembl_db"]["exists"] or vols["chembl_db_is_sqlite"].get("value") is not True:
        problems.append(f"ChEMBL database: {vols['chembl_db']}, SQLite header {vols['chembl_db_is_sqlite']}")
    if not vols["chembl_pins"]["ok"] or len(vols["chembl_pins"]["value"]) != 1:
        problems.append(f"{CHEMBL / VERIFIED_NAME} pins: {vols['chembl_pins']}")
    report["sealed_before_stage_c_may_start"] = {s: s in log.get("seals", {}) for s in ("A", "B")}
    report["problems"], report["passed"] = problems, not problems
    write_report(PROBE_REPORT, report)
    vol.commit()
    print(json.dumps(report, indent=1))
    print(f"stage C probe: {'passed' if report['passed'] else 'FAILED'}; {len(problems)} problem(s)")
    return report


@app.function(image=test_image, cpu=2, memory=8192, timeout=3600, volumes={"/vol": vol})
def tests() -> dict:
    """SYNTHETIC: the stage C pytest suite on the modules baked at /root/exp/stages/C."""
    report = run_pytest(STAGES_IMAGE / "C", [STAGES_IMAGE / "C", STAGES_IMAGE / "run_guard"])
    vol.reload()
    write_report(Path("/vol/tests/modal_tests.json"), report)
    vol.commit()
    print(f"stage C tests: {report['summary']} (exit {report['returncode']})")
    if report["returncode"] != 0:
        print(report["stdout_tail"], report["stderr_tail"])
    return report


@app.local_entrypoint()
def main(run_token: str):
    clean_commit(HERE.parents[3])
    require_run(HERE.parents[1] / "PREREG.md", "C", run_token, manifests={"A": A_OUTPUT / "MANIFEST.tsv"})
    print(run_remote.remote(run_token))
