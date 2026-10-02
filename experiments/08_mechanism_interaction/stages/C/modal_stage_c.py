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

Synthetic test suite (no study data), in the run image plus pytest, with tests/ and stages/run_guard
baked beside the modules at /root/exp/stages/. The summary is printed and written to
/vol/tests/modal_tests.json on pqtl-v8-stage-c, the previous report moved to /vol/tests/superseded/
first:
    uv run --project experiments/08_mechanism_interaction/stages/C --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/C/modal_stage_c.py::tests
"""
from pathlib import Path

import modal

from aact_download import UPLOADED_SOURCE, download_archive, verify_archive
from ctgov import AACT_DIR_NAME
from run_stage_c import run
from v8_run_guard import COMMIT_ENV, baked_commit, clean_commit, require_run, run_commit
from v8_test_report import run_pytest, write_report

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


@app.function(cpu=2, memory=8192, timeout=24 * 3600, retries=2,
              volumes={"/inputs": inputs_vol.with_mount_options(read_only=True), "/chembl": chembl_vol, "/vol": vol})
def run_remote(run_token: str) -> str:
    vol.reload()
    manifest = run(run_token, ROOT / "stages" / "A" / "output" / "hypotheses.csv", INPUTS, AACT, CHEMBL, OUT,
                   prereg=PREREG_IMAGE, pins=ROOT / "exp" / "coverage_all_indications.json", commit=vol.commit,
                   roots=[("experiments/08_mechanism_interaction/stages", ROOT / "stages"), ("pqtl-v8-inputs:", Path("/inputs")),
                          ("proteome-mr-claim-audit-inputs:", Path("/chembl"))], repo_commit=run_commit(None))
    vol.commit()
    return str(manifest)


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
