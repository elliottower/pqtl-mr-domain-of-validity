"""Thin Modal wrapper for stage A (logic in stage_a/ and run_stage_a.py).

Volumes:
  /inputs  pqtl-v8-inputs, read-only: /inputs/08_mechanism_interaction mirrors the experiment
           directory (inputs/{ukbppp,karim2026,...} and feasibility/v2_all_indications/inputs/);
           see scripts/upload_inputs_to_modal.py.
  /vol     pqtl-v8-stage-a: output/ (the stage A tables and MANIFEST.tsv).

Baked into the image: stage_a/, run_stage_a.py, the shared run guard (stages/run_guard/
v8_run_guard.py), PREREG.md at /root/exp/PREREG.md, and the repository files stage A reads that
are not on the volume (deCODE supplementary tables, the UKB-PPP Olink map, the EpiGraphDB gene
list, frozen_candidates_v34.csv, classification_v5.csv, classification_v5_1.csv), at their
repository-relative paths under /root/exp and /root/repo.

Real run only after the OSF registration is approved, `RUN_START stage=A token=<token>` is logged
below the line of PREREG.md (`prereg log`), stage_a.flags.PUBLICATION_DATES is set by a
logged amendment, and the stage code and PREREG.md are committed: the launcher refuses unless
`git status --porcelain` is empty for stages/ and PREREG.md (v8_run_guard.clean_commit), and the
image carries that commit as V8_REPO_COMMIT, which the run writes to run_info.json and INPUTS.tsv.
/vol/output must be empty or hold an unfinished run of the same token. The guard runs twice:
locally before the call, and inside the container on the
PREREG.md baked into the image, so the image must be built after the RUN_START entry is logged:
    uv run --project experiments/08_mechanism_interaction/stages/A --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/A/modal_stage_a.py --run-token <token>
    modal volume get pqtl-v8-stage-a /output experiments/08_mechanism_interaction/stages/A/output

Synthetic test suite (no study data; allowed before the run is logged), in the run image plus
pytest, with this directory and stages/run_guard baked at /root/exp/stages/. The summary is
printed and written to /vol/tests/modal_tests.json on pqtl-v8-stage-a, the previous report moved to
/vol/tests/superseded/ first:
    uv run --project experiments/08_mechanism_interaction/stages/A --with modal==1.4.3 \\
        modal run experiments/08_mechanism_interaction/stages/A/modal_stage_a.py::tests
"""
from pathlib import Path

import modal

from run_stage_a import main as run_stage_a_main
from v8_run_guard import COMMIT_ENV, baked_commit, clean_commit, require_run

HERE = Path(__file__).resolve().parent
INPUTS = Path("/inputs/08_mechanism_interaction")
EXP = Path("/root/exp")
REPO = Path("/root/repo")
OUT = Path("/vol/output")
PY_PACKAGES = ["numpy==2.2.6", "pandas==2.2.3", "pyarrow==21.0.0", "openpyxl==3.1.5", "pydantic==2.11.7",
               "prereg==0.4.2", "provenance-core==0.4.2"]
GUARD = HERE.parent / "run_guard" / "v8_run_guard.py"
PREREG_IMAGE = EXP / "PREREG.md"
TEST_PACKAGES = ["pytest==8.4.1"]
TEST_IGNORE = [".venv", "**/__pycache__", "**/.pytest_cache", "output", "work", "launch_log.jsonl"]
STAGES_IMAGE = EXP / "stages"
EXP_FILES = ("feasibility/inputs/ferkingstad2021_MOESM4_ESM.xlsx",)
REPO_FILES = ("planning/perplexity_v8_design/F_data_readmes/olink_protein_map_3k_v1.tsv",
              "zenodo_export/data/epigraphdb/v34_mr_catalog.csv", "data/frozen_candidates_v34.csv",
              "results/v5/classification_v5.csv", "results/v5_1/classification_v5_1.csv")

image = modal.Image.debian_slim(python_version="3.12").pip_install(*PY_PACKAGES).env({"PYTHONPATH": "/root/stagea"})
test_image = image
if modal.is_local():
    EXP_LOCAL = HERE.parents[1]      # local only: in the container this file sits directly under /root
    REPO_LOCAL = EXP_LOCAL.parents[1]
    image = (image.add_local_dir(str(HERE / "stage_a"), "/root/stagea/stage_a", copy=True, ignore=["__pycache__"])
             .add_local_file(str(HERE / "run_stage_a.py"), "/root/stagea/run_stage_a.py", copy=True)
             .add_local_file(str(GUARD), "/root/stagea/v8_run_guard.py", copy=True)
             .add_local_file(str(GUARD.with_name("v8_manifest.py")), "/root/stagea/v8_manifest.py", copy=True)
             # the wrapper imports v8_test_report at module level, so the run container needs it too
             .add_local_file(str(GUARD.with_name("v8_test_report.py")), "/root/stagea/v8_test_report.py", copy=True)
             .add_local_file(str(EXP_LOCAL / "PREREG.md"), str(PREREG_IMAGE), copy=True))
    for rel in EXP_FILES:
        image = image.add_local_file(str(EXP_LOCAL / rel), str(EXP / rel), copy=True)
    for rel in REPO_FILES:
        image = image.add_local_file(str(REPO_LOCAL / rel), str(REPO / rel), copy=True)
    image = image.env({COMMIT_ENV: baked_commit(REPO_LOCAL)})
    test_image = (image.pip_install(*TEST_PACKAGES)
                  .env({"PYTHONPATH": f"/root/stagea:{STAGES_IMAGE / 'run_guard'}"})
                  .add_local_dir(str(HERE), str(STAGES_IMAGE / "A"), copy=True, ignore=TEST_IGNORE)
                  .add_local_dir(str(HERE.parent / "run_guard"), str(STAGES_IMAGE / "run_guard"), copy=True,
                                 ignore=TEST_IGNORE)
                  .add_local_file(str(HERE.parent / "INTERFACES.md"), str(STAGES_IMAGE / "INTERFACES.md"), copy=True))

with test_image.imports():
    from v8_test_report import run_pytest, write_report

app = modal.App("pqtl-v8-stage-a", image=image)
inputs_vol = modal.Volume.from_name("pqtl-v8-inputs").with_mount_options(read_only=True)
vol = modal.Volume.from_name("pqtl-v8-stage-a", create_if_missing=True)


@app.function(cpu=4, memory=32768, timeout=6 * 3600, volumes={"/inputs": inputs_vol, "/vol": vol})
def run_remote(run_token: str) -> str:
    """run_stage_a.main applies the run guard to the PREREG.md baked into the image."""
    vol.reload()
    manifest = run_stage_a_main(["--inputs-root", str(INPUTS), "--exp-dir", str(EXP), "--repo-root", str(REPO),
                                 "--out-dir", str(OUT), "--run-token", run_token, "--prereg", str(PREREG_IMAGE)])
    vol.commit()
    return str(manifest)


@app.function(image=test_image, cpu=2, memory=8192, timeout=3600, volumes={"/vol": vol})
def tests() -> dict:
    """SYNTHETIC: the stage A pytest suite on the copy of this directory baked at /root/exp/stages/A."""
    report = run_pytest(STAGES_IMAGE / "A", [STAGES_IMAGE / "A", STAGES_IMAGE / "run_guard"])
    vol.reload()
    write_report(Path("/vol/tests/modal_tests.json"), report)
    vol.commit()
    print(f"stage A tests: {report['summary']} (exit {report['returncode']})")
    if report["returncode"] != 0:
        print(report["stdout_tail"], report["stderr_tail"])
    return report


@app.local_entrypoint()
def main(run_token: str):
    clean_commit(HERE.parents[3])
    require_run(HERE.parents[1] / "PREREG.md", "A", run_token)
    print(run_remote.remote(run_token))
