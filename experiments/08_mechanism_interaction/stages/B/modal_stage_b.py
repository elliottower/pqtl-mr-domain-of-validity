"""Thin Modal wrapper for stage B (logic in stage_b/). Nothing here decides anything.

Log `SEAL stage=A manifest_sha256=<sha256>` and then `RUN_START stage=B token=<token>` below the
line of PREREG.md first (`prereg log`): the image bakes PREREG.md and stage A's MANIFEST.tsv,
hypotheses.csv and outcome_trait_coding.tsv, and `plan_remote`, `collect_file` and `run_unit` refuse
unless both entries are in the baked log, the seal is the baked manifest's sha256 and the two files
match it (stage_b.launch.authorize over the shared guard, stages/run_guard/v8_run_guard.py). Commit
the stage code and PREREG.md before deploying: the image carries the commit of a clean tree as
V8_REPO_COMMIT and those functions refuse in an image built from an unclean one. Then deploy, and
start work with the launcher, whose calls survive the client disconnecting:
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal deploy experiments/08_mechanism_interaction/stages/B/modal_stage_b.py
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        python experiments/08_mechanism_interaction/stages/B/launch_stage_b.py plan --run-token <token>
    ... launch_stage_b.py collect --run-token <token>;  status;  spawn --run-token <token>;  status;  assemble ...

Three phases (launch_stage_b.py):
- `plan_remote` reads the baked stage A files and the instrument tables pinned on `pqtl-v8-inputs`
  (stage_b/plan.py) and returns units.jsonl and unit_plan.json. A unit names its files by
  identity (file name, size, ETag); no link or token is in a unit or a fingerprint.
- `collect_file`, one call per whole file (stage_b/collect.py): deCODE per-SeqId files, UKB-PPP
  tars and rsID maps, GWAS Catalog files without a tabix index (for an accession without a
  harmonised file, its GWAS-SSF file or reviewed author format, stage_b/outcome_files.py). Each is downloaded once to
  /vol/stage_b/raw/<source>/, resumed by HTTP Range, the volume committed every 256 MB, and gets a
  record under /vol/stage_b/collect/<source>/ (size, sha256, ETag, Last-Modified, the address
  without query or token, UTC). A file the source does not hold gets a record of `absent`.
- `run_unit`, one call per instrument unit (stage_b/pipeline.py). Whole files are read only from
  the volume, after they match their record; a missing record raises. Regional queries (OpenGWAS,
  tabix on FinnGen, indexed GWAS Catalog files and the eQTL Catalogue, bcftools on 1000 Genomes,
  Ensembl, GTEx) stay remote. Each step of a unit (positions, pQTL region, LD, each outcome region,
  each colocalization, each VEP lookup, splicing, S16) is written to /vol/stage_b/units/<unit_key>/
  and the volume committed before the next step; a restarted call resumes at the first missing step.

Errors (stage_b/remote.py): only a definitive absence (HTTP 404/410, a file the source's listing
does not name, a GWAS Catalog study directory with no summary-statistics file) is recorded as
unavailable, and a GWAS Catalog file the source holds but no fixed rule of
stage_b/outcome_files.py reads is recorded `unreadable` with its reason; both take the registered
consequence of an unavailable file. A refused
credential (401/403), a rate limit, a server error, a timeout, a dropped connection and a truncated
or corrupt stream raise: the call fails, leaves a note under /vol/stage_b/errors/, and the file or
unit stays unfinished until a later call finishes it. `launch_stage_b.py status` lists them. A source
that keeps failing is never turned into "unavailable": Modal's `retries=3` only runs the same call
again, each run raising the same way, and nothing counts failures or writes a record after them.

Secret `pqtl-stage-b`: SYNAPSE_PAT (UKB-PPP), OPENGWAS, DECODE_FOLDER_TOKEN (the token of the deCODE
folder link `https://download.decode.is/folder/<token>`; the link expires and is requested again
through deCODE's form, and a new token only has to replace the old one in the secret; it is needed
by the deCODE collect calls only). Optional: DECODE_SMP_FOLDER_TOKEN (the folder of the
SMP-normalized release, S16), DECODE_FILE_URL_TEMPLATE and DECODE_SMP_FILE_URL_TEMPLATE (the address
of one file of a folder link, with `{token}` and `{key}`, where it differs from
stage_b.fetch.DECODE_FILE_URL).

The whole files under /vol/stage_b/raw/ are working copies, held between collection and the stage B
seal: `purge_raw_files` deletes them after stage B is sealed, and the records and the regional
extracts of the unit directories stay on the volume. The extracts hold rows of the downloaded files
and stay private: they are never written to B/output/ or committed; only their hashes
(regional_manifest.tsv) and the derived summary results are.

A unit checkpoint is resumed only in a directory bound to this run's fingerprint
(stage_b/checkpoint.py: the unit record, the pinned deCODE annotation files as hashed in the
container, the stage_b/ code digest, the frozen plan hash, the versions of Python, the Python
packages, R, coloc, susieR, jsonlite, bcftools/htslib and tabix read in the container at the start
of the call, and the collect digest: the sha256 over the collect records of the unit's whole files,
each with the size and sha256 of the bytes on the volume); any other directory raises. Pinned: R 4.4.3 (rocker/r-ver, dated CRAN snapshot for
transitive packages), coloc 5.2.3, susieR 0.12.35, jsonlite 2.0.0, htslib/bcftools 1.21, and the
Python wheels below.

Pre-analysis validation (stage_b/validate.py), after collect and before spawn: `validate_outcomes`
reads every collected GWAS Catalog whole file (remote-indexed files are queried by region and not
validated) through the analysis reader, one
`validate_outcome_file` call per file, and writes /vol/stage_b/validation/outcome_validation.json
(aggregate parser diagnostics only; the previous report moved to validation/superseded/ first). A
file that fails a rule fixed in stage_b/validate.py is made `unreadable`; a file that passes has
its record bound to the uncertainty mode the validation chose for it (native_se, or_ci_derived_se,
pvalue_coloc):
    ... launch_stage_b.py validate --run-token <token>
A file made `unreadable` for its standard error alone has its `collected` record moved back
(`restore_validated_records`; refused unless its bytes match the recorded sha256 and its result
failed for that reason only; nothing deleted), and is validated again:
    ... launch_stage_b.py restore-validated --run-token <token> --accessions-file <file> ;  validate --run-token <token>

Wiring probe and dry run (no study data; allowed before the run is logged). Both run in the run
image with the volumes and the secret of the real-run functions, so a start of either is a
container start of the real run's configuration:
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/modal_stage_b.py::probe
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/modal_stage_b.py::dry_run
`probe` imports every module the run imports and reports the baked files (PREREG.md, the inputs
manifest, stage A's sealed files against the logged seal), the frozen-plan check of the shared
guard, the mounted volumes (top level; the deCODE annotation files and the tables `plan` reads by
size against their pins and by header only), the tool versions against the pins, which names the
secret sets, whether Synapse and OpenGWAS accept their credentials, whether Synapse's file-handle
service answers as stage_b/synapse_source.py expects (on an rsID map; the link is not followed),
and whether tabix reads a remote header. It never calls `authorize`, makes no request to deCODE and
opens no regional file. The report is printed and written to /vol/probe/probe.json on
pqtl-v8-stage-b, the previous one moved to /vol/probe/superseded/ first. `dry_run` serves the
synthetic sources of synthetic_unit.py from a server inside the container (fake_remote.py) and runs
the real collect and analyze code against it under /vol/dry_run/<utc>/, with faults injected:
broken connections, an expired folder token, a rate limit, a re-issued folder link between collect
and analyze, an expired OpenGWAS token during analyze. Its report is /vol/dry_run/<utc>/report.json.

Synthetic test suite (no study data), in the run image plus pytest, with this directory and
stages/run_guard baked at /root/exp/stages/; the R coloc tests and the tabix and bcftools tests
skipped on the laptop run here. The summary is printed and written to /vol/tests/modal_tests.json
on pqtl-v8-stage-b, the previous report moved to /vol/tests/superseded/ first. `tests` belongs to
its own app, pqtl-v8-stage-b-tests, because the suite needs no credentials:
    uv run --project experiments/08_mechanism_interaction/stages/B --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/B/modal_stage_b.py::tests
"""
import base64
import gzip
import importlib
import json
import os
import subprocess
import tempfile
import time
from collections import Counter
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

import modal

from v8_manifest import MANIFEST_NAME, sha256_file, verify_listed
from v8_run_guard import (COMMIT_ENV, PLAN_SHA256, RUN_START_RE, baked_commit, check_plan, logged_seals, parse_log,
                          run_commit)
from v8_test_report import run_pytest, write_report

R_VERSION = "4.4.3"
HTSLIB_VERSION = "1.21"
R_PACKAGES = {"jsonlite": "2.0.0", "susieR": "0.12.35", "coloc": "5.2.3"}
PY_PACKAGES = ["numpy==2.1.3", "pandas==2.2.3", "pydantic==2.13.5", "requests==2.32.3", "openpyxl==3.1.5",
               "synapseclient==4.12.0", "prereg==0.4.2", "provenance-core==0.4.2", "PyYAML==6.0.3"]
HERE = Path(__file__).resolve().parent
PREREG_IMAGE = Path("/root/exp/PREREG.md")
INPUTS_MANIFEST_IMAGE = Path("/root/exp/modal_inputs_manifest.json")
A_IMAGE = Path("/root/stages/A/output")
PACKAGE_IMAGE = Path("/root/stageb/stage_b")
A_BAKED = ("MANIFEST.tsv", "hypotheses.csv", "outcome_trait_coding.tsv")
RUN_MODULES = ("numpy", "pandas", "pydantic", "requests", "openpyxl", "synapseclient", "yaml", "prereg.log", "prereg.plan",
               "v8_run_guard", "v8_manifest", "v8_test_report", "synthetic_unit", "fake_remote", "stage_b.assemble",
               "stage_b.checkpoint", "stage_b.collect", "stage_b.coloc_backend", "stage_b.evidence", "stage_b.fetch",
               "stage_b.harmonize", "stage_b.launch", "stage_b.ld", "stage_b.parsers", "stage_b.pipeline", "stage_b.plan",
               "stage_b.outcome_files", "stage_b.author_formats",
               "stage_b.remote", "stage_b.schemas", "stage_b.sentinels", "stage_b.status", "stage_b.synapse_source",
               "stage_b.units", "stage_b.validate")
SECRET_NAMES = ("SYNAPSE_PAT", "OPENGWAS")
COLLECT_SECRET_NAMES = ("DECODE_FOLDER_TOKEN",)     # set just before collect; its absence blocks `collect --source decode` only
OPTIONAL_SECRET_NAMES = ("DECODE_SMP_FOLDER_TOKEN", "DECODE_FILE_URL_TEMPLATE", "DECODE_SMP_FILE_URL_TEMPLATE")
# Containers at once. Collect: parallel whole-file downloads from one source. Analyze: every unit
# container queries OpenGWAS, Ensembl and the eQTL Catalogue, whose rate limit is per client.
COLLECT_CONTAINERS = 16
ANALYZE_CONTAINERS = 24
TEST_PACKAGES = ["pytest==9.1.1"]
TEST_IGNORE = [".venv", "**/__pycache__", "**/.pytest_cache", "output", "work", "launch_log.jsonl"]
STAGES_IMAGE = Path("/root/exp/stages")

_r_install = "; ".join(
    f"remotes::install_version('{p}', version='{v}', upgrade='never')" for p, v in R_PACKAGES.items())
image = (
    modal.Image.from_registry(f"rocker/r-ver:{R_VERSION}", add_python="3.12")
    .apt_install("wget", "bzip2", "make", "gcc", "libcurl4-openssl-dev", "libbz2-dev", "liblzma-dev", "zlib1g-dev",
                 "libssl-dev", "libncurses-dev")
    .run_commands(
        f"wget -q https://github.com/samtools/bcftools/releases/download/{HTSLIB_VERSION}/bcftools-{HTSLIB_VERSION}.tar.bz2"
        f" && tar xjf bcftools-{HTSLIB_VERSION}.tar.bz2"
        f" && cd bcftools-{HTSLIB_VERSION}/htslib-{HTSLIB_VERSION} && ./configure --enable-libcurl && make -j4 && make install"
        f" && cd .. && ./configure --with-htslib=system && make -j4 && make install && ldconfig",
        "Rscript -e \"install.packages('remotes')\"",
        f"Rscript -e \"{_r_install}\"",
        "Rscript -e \"stopifnot(packageVersion('coloc') == '5.2.3', packageVersion('susieR') == '0.12.35')\"",
    )
    .pip_install(*PY_PACKAGES)
    .env({"PYTHONPATH": "/root/stageb"})
)
test_image = image
if modal.is_local():
    EXP = HERE.parents[1]   # local only: in the container this file sits directly under /root
    image = (image.add_local_dir(str(HERE / "stage_b"), "/root/stageb/stage_b", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_run_guard.py"), "/root/stageb/v8_run_guard.py", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_manifest.py"), "/root/stageb/v8_manifest.py", copy=True)
             # this file imports v8_test_report at module level, so the run container needs it too
             .add_local_file(str(HERE.parent / "run_guard" / "v8_test_report.py"), "/root/stageb/v8_test_report.py",
                             copy=True)
             .add_local_file(str(HERE / "synthetic_unit.py"), "/root/stageb/synthetic_unit.py", copy=True)
             .add_local_file(str(HERE / "fake_remote.py"), "/root/stageb/fake_remote.py", copy=True)
             .add_local_file(str(EXP / "PREREG.md"), str(PREREG_IMAGE), copy=True)
             .add_local_file(str(EXP / "modal_inputs_manifest.json"), str(INPUTS_MANIFEST_IMAGE), copy=True))
    for _name in A_BAKED:
        if (HERE.parent / "A" / "output" / _name).is_file():
            image = image.add_local_file(str(HERE.parent / "A" / "output" / _name), str(A_IMAGE / _name), copy=True)
    image = image.env({COMMIT_ENV: baked_commit(EXP.parents[1])})
    test_image = (image.pip_install(*TEST_PACKAGES)
                  .env({"PYTHONPATH": f"/root/stageb:{STAGES_IMAGE / 'run_guard'}"})
                  .add_local_dir(str(HERE), str(STAGES_IMAGE / "B"), copy=True, ignore=TEST_IGNORE)
                  .add_local_dir(str(HERE.parent / "run_guard"), str(STAGES_IMAGE / "run_guard"), copy=True,
                                 ignore=TEST_IGNORE)
                  .add_local_file(str(HERE.parent / "INTERFACES.md"), str(STAGES_IMAGE / "INTERFACES.md"), copy=True))

app = modal.App("pqtl-v8-stage-b", image=image)
test_app = modal.App("pqtl-v8-stage-b-tests")
vol = modal.Volume.from_name("pqtl-v8-stage-b", create_if_missing=True)
inputs_vol = modal.Volume.from_name("pqtl-v8-inputs").with_mount_options(read_only=True)
secret = modal.Secret.from_name("pqtl-stage-b")
ROOT = Path("/vol/stage_b")
PROBE_REPORT = Path("/vol/probe/probe.json")
DRY_RUN_ROOT = Path("/vol/dry_run")
# What the real-run functions run in. `probe` and `dry_run` take the same image, volumes and secret.
RUN_CONTAINER = {"cpu": 2, "memory": 8192, "volumes": {"/vol": vol, "/inputs": inputs_vol}, "secrets": [secret]}
EXP_INPUTS = Path("/inputs/08_mechanism_interaction")
DECODE_INPUTS = EXP_INPUTS / "inputs" / "decode"

with image.imports():
    import numpy as np
    import requests

    from fake_remote import FakeRemote
    from stage_b.checkpoint import (package_sha256, source_pins, tool_versions, unit_fingerprint, unit_tools,
                                    verified_source_pins, volume_collect_digest)
    from stage_b.collect import collect_one, purge_raw, read_record, supersede_absent
    from stage_b.coloc_backend import ColocDataset, ColocTask, RscriptColoc
    from stage_b.fetch import (DECODE_FILE_URL, ENSEMBL_REST, EQTLCAT_PATHS, GWASCAT_FTP, KG_PANEL, KG_VCF, OPENGWAS_API,
                               SYNAPSE_RSID_MAPS, SYNAPSE_UKBPPP_EUR, Endpoints, RemoteSources, VolumeFetcher, https_path)
    from stage_b.launch import A_FILES, authorize, sealed_stage_b
    from stage_b.pipeline import DirStore, process_unit, reset_unavailable
    from stage_b.plan import PLAN_TABLES, make_plan, pinned_tables
    from stage_b.schemas import CollectTask, InstrumentUnit, StageBError
    from stage_b.sentinels import DECODE_ST02, INTERVAL_ST4, UKBPPP_ST9, sheet_header
    from stage_b.status import marked, volume_state
    from stage_b.synapse_source import SynapseSource
    from stage_b.validate import (REPORT_NAME, VALIDATION_DIR, apply_validation, gwas_catalog_records, in_scope,
                                  restorable, restore_validated, validate_file, validation_report)
    from synthetic_unit import build_world, scenario


def endpoints() -> "Endpoints":
    """The registered sources; the address of one file under a deCODE folder link from the
    environment where the secret sets it."""
    return Endpoints(decode_file=os.environ.get("DECODE_FILE_URL_TEMPLATE") or DECODE_FILE_URL,
                     decode_smp_file=os.environ.get("DECODE_SMP_FILE_URL_TEMPLATE") or DECODE_FILE_URL)


def checkpointed_unit(unit: "InstrumentUnit", fetcher: "VolumeFetcher", pins: dict[str, str], units_root: Path) -> dict:
    """One unit under its fingerprint, every step committed to the volume: what `run_unit` does
    after its guard, and what `dry_run` does with the synthetic units. The collect digest is taken
    from the records under the root the fetcher reads whole files from; a missing record raises."""
    tools = tool_versions()
    collected = volume_collect_digest(fetcher.root, unit)
    fingerprint = unit_fingerprint(unit, pins, package_sha256(PACKAGE_IMAGE), PLAN_SHA256, tools, collected)
    store = DirStore(units_root / unit.unit_key, commit=vol.commit)
    result = process_unit(unit, fetcher, RscriptColoc(), store, fingerprint, tools, collected)
    vol.commit()
    return result


@app.function(**{**RUN_CONTAINER, "memory": 16384}, timeout=2 * 3600)
def plan_remote(run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    units_text, plan = make_plan(A_IMAGE / "hypotheses.csv", A_IMAGE / "outcome_trait_coding.tsv",
                                 pinned_tables(INPUTS_MANIFEST_IMAGE, EXP_INPUTS))
    vol.reload()
    out = ROOT / "plan" / plan["units_sha256"]          # one directory per distinct plan; none is overwritten
    out.mkdir(parents=True, exist_ok=True)
    (out / "units.jsonl").write_text(units_text)
    (out / "unit_plan.json").write_text(json.dumps(plan, indent=1, sort_keys=True))
    vol.commit()
    return {"unit_plan": plan, "units_jsonl": units_text}


# One collect task. With `supersede_absent_record`, an `absent` GWAS Catalog record is first moved to
# /stage_b/superseded/collect/ (stage_b.collect.supersede_absent; any other record is kept), so the task
# runs again under the current rules. The guard stays the first statement (tests/test_launch_guard.py).
@app.function(**RUN_CONTAINER, timeout=24 * 3600, retries=3, max_containers=COLLECT_CONTAINERS)
def collect_file(task_json: str, run_token: str, supersede_absent_record: bool = False) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    task = CollectTask.model_validate_json(task_json)
    superseded = supersede_absent(ROOT, task, vol.commit) if supersede_absent_record else None
    sources = RemoteSources(ROOT / "cache", os.environ.get("DECODE_FOLDER_TOKEN", ""),
                            lambda: SynapseSource(os.environ["SYNAPSE_PAT"]), endpoints(),
                            decode_smp_token=os.environ.get("DECODE_SMP_FOLDER_TOKEN", ""))
    record = marked(ROOT, "collect", f"{task.source}__{task.key}", f"{task.source}/{task.key}",
                    lambda: collect_one(task, sources, ROOT, vol.commit), vol.commit)
    return {"source": task.source, "key": task.key, "status": record.status, "bytes": record.bytes,
            "superseded": None if superseded is None else str(superseded.relative_to(ROOT))}


# Pre-analysis validation (stage_b/validate.py): every collected GWAS Catalog whole file (harmonised,
# GWAS-SSF or reviewed author format) is read whole by the analysis reader, one call per file, checkpointed every 1M rows; a file that fails
# a fixed rule is made `unreadable`. Aggregate parser diagnostics only. Run after collect, before spawn.
@app.function(**{**RUN_CONTAINER, "memory": 16384}, timeout=24 * 3600, retries=3, max_containers=COLLECT_CONTAINERS)
def validate_outcome_file(key: str, run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    record = read_record(ROOT, "gwas_catalog", key)
    result = validate_file(ROOT, record, vol.commit)
    apply_validation(ROOT, record, result, vol.commit)
    return result


@app.function(**RUN_CONTAINER, timeout=24 * 3600)
def validate_outcomes(run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    keys = [r.key for r in gwas_catalog_records(ROOT) if in_scope(r)]
    results = {}
    for key, got in zip(keys, validate_outcome_file.map(keys, kwargs={"run_token": run_token}, return_exceptions=True),
                        strict=True):
        results[key] = got if isinstance(got, dict) else f"{type(got).__name__}: {got}"[:500]
    vol.reload()
    report = {**validation_report(gwas_catalog_records(ROOT), results),
              "modal_image_id": os.environ.get("MODAL_IMAGE_ID", ""), COMMIT_ENV: os.environ.get(COMMIT_ENV)}
    write_report(ROOT / VALIDATION_DIR / REPORT_NAME, report)
    vol.commit()
    return report


# Undo the validation's `unreadable` for files whose failure was a reader defect (stage_b/validate.py
# `restorable`, `restore_validated`). Every key is checked, its bytes against the recorded sha256, before
# anything is moved; one refusal moves nothing. Then run `validate` again. The guard stays first.
@app.function(**RUN_CONTAINER, timeout=3600)
def restore_validated_records(keys: list[str], run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    plans = [restorable(ROOT, key) for key in keys]
    restored = [restore_validated(ROOT, plan, vol.commit) for plan in plans]
    vol.commit()
    return {"restored": restored}


@app.function(**RUN_CONTAINER, timeout=24 * 3600, retries=3, max_containers=ANALYZE_CONTAINERS)
def run_unit(unit_json: str, run_token: str) -> dict:
    authorize(PREREG_IMAGE, run_token, A_IMAGE)
    run_commit(None)
    vol.reload()
    unit = InstrumentUnit.model_validate_json(unit_json)
    pins = verified_source_pins(INPUTS_MANIFEST_IMAGE, EXP_INPUTS)
    fetcher = VolumeFetcher(ROOT, os.environ["OPENGWAS"], DECODE_INPUTS / "assocvariants.annotated.txt.gz",
                            DECODE_INPUTS / "assocvariants.excluded.txt.gz", endpoints())
    result = marked(ROOT, "units", unit.unit_key, unit.unit_key,
                    lambda: checkpointed_unit(unit, fetcher, pins, ROOT / "units"), vol.commit)
    return {"unit_key": unit.unit_key, "pqtl_available": result["pqtl_available"], "outcomes": len(result["outcomes"])}


@app.function(timeout=600, volumes={"/vol": vol})
def reset_unit(unit_key: str) -> list[str]:
    vol.reload()
    return reset_unavailable(DirStore(ROOT / "units" / unit_key, commit=vol.commit))


@app.function(timeout=3600, volumes={"/vol": vol})
def purge_raw_files(units_jsonl: str) -> dict:
    """Delete the collected whole files once stage B is sealed (stage_b.collect.purge_raw): refuses
    unless the baked PREREG.md log holds stage B's SEAL, every unit has its result and every planned
    file its record. The records, the unit directories and their regional extracts stay."""
    sealed_stage_b(PREREG_IMAGE)
    vol.reload()
    units = [InstrumentUnit.model_validate_json(line) for line in units_jsonl.splitlines()]
    out = purge_raw(ROOT, units, vol.commit)
    return {"files_deleted": len(out["deleted"]), "bytes_deleted": out["bytes_deleted"], "records_kept": out["records_kept"]}


@app.function(timeout=900, volumes={"/vol": vol})
def stage_state() -> dict:
    """What the volume holds, for `launch_stage_b.py status` and `spawn`: collect records, partial
    downloads, finished units and error notes (stage_b.status.volume_state). No regional file is opened."""
    vol.reload()
    return volume_state(ROOT)


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


def _gz_header(path: Path) -> list[str]:
    with gzip.open(path, "rt") as fh:
        return fh.readline().rstrip("\n").split("\t")


def _http_status(url: str, headers: dict[str, str] | None = None) -> int:
    """The status line of a GET; the body is not read."""
    with requests.get(url, headers=headers or {}, stream=True, timeout=60) as r:
        return r.status_code


def _jwt_expiry(token: str) -> str:
    """The `exp` claim of a JSON Web Token as a UTC time; no other claim is returned."""
    payload = token.split(".")[1]
    exp = json.loads(base64.urlsafe_b64decode(payload + "=" * (-len(payload) % 4)))["exp"]
    return datetime.fromtimestamp(exp, timezone.utc).isoformat(timespec="seconds")


def _log_state() -> dict:
    """The baked PREREG.md under the guard's frozen-plan check (not its authorization), the seals
    its log holds, and the stages with a RUN_START entry."""
    text = PREREG_IMAGE.read_text(encoding="utf-8")
    check_plan(text)
    entries = parse_log(text)
    starts = sorted({m.group(1) for m in (RUN_START_RE.fullmatch(e.event) for e in entries) if m is not None})
    return {"plan_sha256": PLAN_SHA256, "log_entries": len(entries), "seals": logged_seals(entries, len(entries) + 1),
            "run_start_stages": starts}


def _remote_header(url: str) -> dict:
    """`tabix -H` on a remote indexed file: the header lines only, read through htslib's libcurl."""
    with tempfile.TemporaryDirectory() as tmp:
        res = subprocess.run(["tabix", "-H", url], cwd=tmp, capture_output=True, text=True, timeout=600)
    lines = res.stdout.splitlines()
    return {"returncode": res.returncode, "header_lines": len(lines), "first": lines[0][:80] if lines else "",
            "stderr": res.stderr[-300:]}


def _credentials() -> dict:
    """Whether Synapse and OpenGWAS accept their credentials, and whether Synapse's file-handle
    service returns what stage_b/synapse_source.py reads, asked for one rsID map (a position-to-rsID
    table). The link it returns is not followed and not reported. No protein file is listed."""
    synapse = SynapseSource(os.environ["SYNAPSE_PAT"])
    permissions = synapse.syn.restGET(f"/entity/{SYNAPSE_UKBPPP_EUR}/permissions")
    maps = synapse.children(SYNAPSE_RSID_MAPS)
    handle = synapse.file(maps[0]["id"])
    bearer = {"Authorization": f"Bearer {os.environ['OPENGWAS']}"}
    return {"synapse_login": True,
            "synapse_ukbppp_folder": {k: permissions.get(k) for k in ("canView", "canDownload")},
            "synapse_rsid_maps_listed": len(maps),
            "synapse_file_handle": {"name": handle["name"], "size": handle["size"], "md5_declared": bool(handle["md5"]),
                                    "link_returned": handle["url"].startswith("https://")},
            "opengwas_user_status": _http_status(f"{OPENGWAS_API}/user", bearer),
            "opengwas_token_expires_utc": _jwt_expiry(os.environ["OPENGWAS"])}


def _eqtl_catalogue_access() -> dict:
    """How the eQTL Catalogue's GTEx files are addressed in its path table, and whether the first
    of them and its index answer over HTTPS (HEAD requests: no file is read)."""
    with requests.get(EQTLCAT_PATHS, timeout=120) as r:
        rows = [line.split("\t") for line in r.text.splitlines()]
    cols = {name: i for i, name in enumerate(rows[0])}
    paths = sorted(row[cols["ftp_path"]] for row in rows[1:] if row[cols["study_label"]] == "GTEx")
    first = https_path(paths[0])
    head = {name: requests.head(url, allow_redirects=True, timeout=60).status_code
            for name, url in (("file", first), ("index", first + ".tbi"))}
    return {"gtex_datasets": len(paths), "schemes_in_path_table": dict(Counter(p.split("://")[0] for p in paths)),
            "read_as": first.split("://")[0], "all_on_ebi_ftp_host": all(https_path(p).startswith("https://") for p in paths),
            "head_status_of_first": head}


def _plan_tables() -> dict:
    """The tables `plan` reads, on the inputs volume: size against the pin, and the header row of
    the sheet each workbook is read from (the listing and the OpenGWAS file by size only)."""
    manifest = {f["path"]: f for f in json.loads(INPUTS_MANIFEST_IMAGE.read_text())["files"]}
    sheets = {"ukbppp_st9": UKBPPP_ST9, "decode_st02": DECODE_ST02, "interval_st4": INTERVAL_ST4}
    out = {}
    for name, rel in PLAN_TABLES.items():
        rec = {**_file(EXP_INPUTS / rel), "pinned": rel in manifest, "pinned_bytes": manifest.get(rel, {}).get("bytes")}
        if name in sheets and rec["exists"]:
            sheet, row, columns = sheets[name]
            rec["sheet"], rec["columns_read"] = sheet, list(columns)
            rec["header"] = _attempt(lambda p=EXP_INPUTS / rel, s=sheet, r=row: sheet_header(p, s, r))
        out[name] = rec
    return out


@app.function(**RUN_CONTAINER, timeout=1800)
def probe() -> dict:
    """WIRING ONLY, in the container configuration of the real-run functions: imports, baked files,
    mounted volumes, tool versions, credentials. Reads no regional statistics, makes no request to
    deCODE and never calls `authorize`."""
    vol.reload()
    pinned_bytes = {f["path"]: f["bytes"] for f in json.loads(INPUTS_MANIFEST_IMAGE.read_text())["files"]}
    decode = {}
    for name, columns in (("assocvariants.annotated.txt.gz", ["Name", "effectAllele", "otherAllele", "effectAlleleFreq"]),
                          ("assocvariants.excluded.txt.gz", ["Name"])):
        rec = _file(DECODE_INPUTS / name)
        header = _attempt(lambda name=name: _gz_header(DECODE_INPUTS / name))
        decode[name] = {**rec, "pinned_bytes": pinned_bytes.get(f"inputs/decode/{name}"), "header": header,
                        "columns_the_parser_needs": columns}
    tools = _attempt(tool_versions)
    pinned_tools = {**{f"python:{p.split('==')[0]}": p.split("==")[1] for p in PY_PACKAGES}, "R": R_VERSION,
                    **{f"R:{p}": v for p, v in R_PACKAGES.items()}, "bcftools": HTSLIB_VERSION,
                    "bcftools:htslib": HTSLIB_VERSION, "tabix:htslib": HTSLIB_VERSION}
    report = {
        "stage": "B", "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "modal_image_id": os.environ.get("MODAL_IMAGE_ID", ""), COMMIT_ENV: os.environ.get(COMMIT_ENV),
        "imports": {m: _attempt(lambda m=m: str(importlib.import_module(m).__file__)) for m in RUN_MODULES},
        "baked": {"prereg": _file(PREREG_IMAGE, sha256=True), "inputs_manifest": _file(INPUTS_MANIFEST_IMAGE, sha256=True),
                  "coloc_run.R": _file(PACKAGE_IMAGE / "coloc_run.R", sha256=True),
                  "stage_a": {n: _file(A_IMAGE / n, sha256=True) for n in A_BAKED},
                  "stage_a_files_match_manifest": _attempt(lambda: verify_listed(A_IMAGE, A_FILES)),
                  "source_pins": _attempt(lambda: source_pins(INPUTS_MANIFEST_IMAGE)),
                  "package_sha256": _attempt(lambda: package_sha256(PACKAGE_IMAGE))},
        "prereg_log": _attempt(_log_state),
        "volumes": {"/inputs": _top_level(Path("/inputs")), "experiment": _top_level(EXP_INPUTS),
                    "decode_dir": _top_level(DECODE_INPUTS), "decode": decode, "plan_tables": _attempt(_plan_tables),
                    "/vol": _top_level(Path("/vol")), "stage_root": _top_level(ROOT)},
        "tools": tools, "pinned_tools": pinned_tools,
        "secret": {name: bool(os.environ.get(name)) for name in SECRET_NAMES},
        "secret_for_collect": {name: bool(os.environ.get(name)) for name in COLLECT_SECRET_NAMES},
        "secret_optional": {name: bool(os.environ.get(name)) for name in OPTIONAL_SECRET_NAMES},
        "decode_file_url_template": endpoints().decode_file.replace("{token}", "<token>"),
        "credentials": _attempt(_credentials),
        "eqtl_catalogue_access": _attempt(_eqtl_catalogue_access),
        "network": {"ensembl_ping": {b: _attempt(lambda u=u: _http_status(f"{u}/info/ping?content-type=application/json"))
                                     for b, u in ENSEMBL_REST.items()},
                    "opengwas_status": _attempt(lambda: _http_status(f"{OPENGWAS_API}/status")),
                    "gwas_catalog_ftp": _attempt(lambda: _http_status(f"{GWASCAT_FTP}/")),
                    "eqtl_catalogue_paths": _attempt(lambda: _http_status(EQTLCAT_PATHS)),
                    "1000g_panel": _attempt(lambda: _http_status(KG_PANEL)),
                    "1000g_vcf_header_chr22": _attempt(lambda: _remote_header(KG_VCF.format(chrom="22")))},
    }
    log = report["prereg_log"].get("value", {})
    seal = report["baked"]["stage_a"][MANIFEST_NAME].get("sha256")
    creds = report["credentials"].get("value", {})
    vcf = report["network"]["1000g_vcf_header_chr22"].get("value", {})
    problems = [f"import {m}: {r['error']}" for m, r in report["imports"].items() if not r["ok"]]
    problems += [f"baked file missing: {r['path']}" for r in (report["baked"]["prereg"], report["baked"]["inputs_manifest"],
                                                              report["baked"]["coloc_run.R"],
                                                              *report["baked"]["stage_a"].values()) if not r["exists"]]
    problems += [f"{k}: {report['baked'][k]['error']}" for k in ("stage_a_files_match_manifest", "source_pins", "package_sha256")
                 if not report["baked"][k]["ok"]]
    if not report["prereg_log"]["ok"]:
        problems.append(f"frozen-plan check: {report['prereg_log']['error']}")
    elif log["seals"].get("A") != seal:
        problems.append(f"baked stage A MANIFEST.tsv sha256 {seal} is not the logged seal {log['seals'].get('A')}")
    for name, rec in decode.items():
        if not rec["exists"] or rec["bytes"] != rec["pinned_bytes"]:
            problems.append(f"{name}: on the volume {rec.get('bytes')} bytes, pinned {rec['pinned_bytes']}")
        elif not rec["header"]["ok"] or not set(rec["columns_the_parser_needs"]) <= set(rec["header"]["value"]):
            problems.append(f"{name}: header {rec['header']} lacks {rec['columns_the_parser_needs']}")
    tables = report["volumes"]["plan_tables"]
    if not tables["ok"]:
        problems.append(f"plan tables: {tables['error']}")
    for name, rec in tables.get("value", {}).items():
        if not rec["exists"] or not rec["pinned"] or rec["bytes"] != rec["pinned_bytes"]:
            problems.append(f"plan table {name}: on the volume {rec.get('bytes')} bytes, pinned {rec['pinned_bytes']}")
        elif "header" in rec and (not rec["header"]["ok"] or not set(rec["columns_read"]) <= set(rec["header"]["value"])):
            problems.append(f"plan table {name}: sheet {rec['sheet']} header lacks the columns read: {rec['header']}")
    if not tools["ok"]:
        problems.append(f"tool versions: {tools['error']}")
    else:
        problems += [f"{k}: installed {tools['value'].get(k)}, pinned {v}" for k, v in pinned_tools.items()
                     if tools["value"].get(k) != v]
    problems += [f"secret pqtl-stage-b does not set {name}" for name, present in report["secret"].items() if not present]
    if not report["credentials"]["ok"]:
        problems.append(f"credentials: {report['credentials']['error']}")
    elif creds["opengwas_user_status"] != 200:
        problems.append(f"OpenGWAS /user answered {creds['opengwas_user_status']} to the token")
    eqtl = report["eqtl_catalogue_access"]
    if not eqtl["ok"]:
        problems.append(f"eQTL Catalogue access: {eqtl['error']}")
    elif not eqtl["value"]["all_on_ebi_ftp_host"] or set(eqtl["value"]["head_status_of_first"].values()) != {200}:
        problems.append(f"eQTL Catalogue files do not answer over HTTPS as expected: {eqtl['value']}")
    report["warnings"] = [f"secret pqtl-stage-b does not set {name} yet: this blocks `collect --source decode` only"
                          for name, present in report["secret_for_collect"].items() if not present]
    if vcf.get("returncode") != 0 or not vcf.get("header_lines"):
        problems.append(f"tabix could not read a remote header: {report['network']['1000g_vcf_header_chr22']}")
    report["problems"], report["passed"] = problems, not problems
    write_report(PROBE_REPORT, report)
    vol.commit()
    print(json.dumps(report, indent=1))
    print(f"stage B probe: {'passed' if report['passed'] else 'FAILED'}; {len(problems)} problem(s); "
          f"{len(report['warnings'])} warning(s) {report['warnings']}")
    return report


@app.function(**RUN_CONTAINER, timeout=3600)
def dry_run() -> dict:
    """SYNTHETIC: collect and analyze the units of synthetic_unit.py against a server inside this
    container, through `collect_one` and `checkpointed_unit` with the R backend, under injected
    faults (synthetic_unit.scenario). The pins are read from the baked manifest without hashing the
    mounted files, and the annotation files are synthetic, so no input file is opened."""
    vol.reload()
    started = time.monotonic()
    root = DRY_RUN_ROOT / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    pins = source_pins(INPUTS_MANIFEST_IMAGE)
    with tempfile.TemporaryDirectory() as tmp, FakeRemote() as remote:
        report = scenario(build_world(remote, Path(tmp)), root, lambda unit, fetcher, units_root:
                          checkpointed_unit(unit, fetcher, pins, units_root), vol.commit)
    decode_dir = root / "units" / next(k for k in report["units"] if k.startswith("decode__"))
    sessions = {k: u["coloc_session"].get("coloc") for k, u in report["units"].items()}
    report["checks"]["r_session_is_the_pinned_coloc"] = set(sessions.values()) == {R_PACKAGES["coloc"]}
    report.update({"synthetic": True, "stage": "B", "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                   "modal_image_id": os.environ.get("MODAL_IMAGE_ID", ""), "root": str(root),
                   "tools": unit_tools(decode_dir), "passed": all(report["checks"].values()),
                   "seconds": round(time.monotonic() - started, 1)})
    (root / "report.json").write_text(json.dumps(report, indent=1))
    vol.commit()
    print(json.dumps({k: report[k] for k in ("checks", "events", "status", "passed", "seconds", "root")}, indent=1))
    return report


@test_app.function(image=test_image, cpu=2, memory=8192, timeout=3600, volumes={"/vol": vol})
def tests() -> dict:
    """SYNTHETIC: the stage B pytest suite on the copy of this directory baked at /root/exp/stages/B."""
    report = run_pytest(STAGES_IMAGE / "B", [STAGES_IMAGE / "B", STAGES_IMAGE / "B" / "tests", STAGES_IMAGE / "run_guard"])
    vol.reload()
    write_report(Path("/vol/tests/modal_tests.json"), report)
    vol.commit()
    print(f"stage B tests: {report['summary']} (exit {report['returncode']})")
    if report["returncode"] != 0:
        print(report["stdout_tail"], report["stderr_tail"])
    return report


@app.function(cpu=2, memory=4096, timeout=1800)
def smoke() -> dict:
    """coloc.abf and coloc.susie on a tiny synthetic region with a known answer, plus tool
    versions. Region: 60 variants, AR(1) LD r = 0.6^|i-j|; both traits' z-scores are the LD
    column of variant 30 scaled to z = 10 (one shared causal variant), so PP.H4 must dominate
    and the lead variant must be v30. A second pair puts the outcome signal on v10 (distinct
    causal variants), so PP.H3 must dominate."""
    m, causal, other = 60, 30, 10
    idx = np.arange(m)
    ld = 0.6 ** np.abs(idx[:, None] - idx[None, :])
    snp = [f"v{i}" for i in range(m)]
    se = np.full(m, 0.02)

    def ds(zc: np.ndarray, kind: str) -> ColocDataset:
        extra = {"sdY": 1.0} if kind == "quant" else {"s": 0.3}
        return ColocDataset(snp=snp, beta=list(zc * se), varbeta=list(se ** 2), N=20000.0, type=kind,
                            MAF=[0.3] * m, **extra)

    z_same, z_other = 10 * ld[causal], 10 * ld[other]
    tasks = [ColocTask(id="shared", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"), d2=ds(z_same, "cc")),
             ColocTask(id="distinct", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"), d2=ds(z_other, "cc")),
             ColocTask(id="susie", method="susie", p1=1e-4, p2=1e-4, p12=5e-6, d1=ds(z_same, "quant"),
                       d2=ds(z_same, "cc"), LD=ld.tolist())]
    backend = RscriptColoc()
    res = backend.run(tasks)
    shared, distinct = res["shared"], res["distinct"]
    if not (shared.pp[4] > 0.9 and shared.lead_variant() == f"v{causal}" and distinct.pp[3] > 0.9):
        raise StageBError(f"coloc smoke failed: shared {shared.pp}, lead {shared.lead_variant()}, distinct {distinct.pp}")
    tools = {t: subprocess.run([t, "--version"], capture_output=True, text=True).stdout.splitlines()[0]
             for t in ("tabix", "bcftools")}
    return {"shared_pp": shared.pp, "lead": shared.lead_variant(), "distinct_pp": distinct.pp,
            "susie": res["susie"].model_dump(), "session": backend.session, "tools": tools}
