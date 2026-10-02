"""Thin Modal wrapper for stage D (logic in stage_d/ and run_stage_d.py).

Image: the calibration image of power/power_v9/modal_calibrate_v9.py (PyMC 5.28.5, PyTensor
2.38.2, ArviZ 0.23.4, NumPy 2.5.3, SciPy 1.18.1, pandas 3.0.6, statsmodels 0.15.0, OpenBLAS)
plus matplotlib, pydantic and pytest.

BLAS: PyTensor links against OpenBLAS (Debian bookworm libopenblas-dev 0.3.21+ds-4, pthread
variant) through PYTENSOR_FLAGS=blas__ldflags=-lopenblas, so its C BLAS ops (CGemv, Dot22)
are used instead of the degraded fallback. Modal sets OPENBLAS_NUM_THREADS to the function's CPU
request at container start, so the image does not set it. Checked by pytensor_blas_check.check()
inside the image.

Synthetic only (allowed before stages A-C exist). This file imports the shared guard, so it is run
in the stage D project (or with `--with-editable experiments/08_mechanism_interaction/stages/run_guard`):
    uv run --project experiments/08_mechanism_interaction/stages/D --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/D/modal_stage_d.py --what smoke
    uv run --project experiments/08_mechanism_interaction/stages/D --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/D/modal_stage_d.py --what tests
    uv run --project experiments/08_mechanism_interaction/stages/D --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/D/modal_stage_d.py --what blas
    uv run --project experiments/08_mechanism_interaction/stages/D --with modal==1.4.3 \
        modal run experiments/08_mechanism_interaction/stages/D/modal_stage_d.py --what env

Real run (only after `SEAL stage=A manifest_sha256=<sha256>`, the same for B and C, and then
`RUN_START stage=D token=<token>` are logged below the line of PREREG.md): the stage outputs under
stages/{A,B,C}/output and PREREG.md are baked into the image, so the image is built after those
entries are logged. `prepare_remote`, `frequentist`, `fit_one` and `assemble_remote` each take
the run token and apply the shared guard (stages/run_guard/v8_run_guard.py) to the baked
PREREG.md and baked MANIFEST.tsv files before reading anything; the expected manifest hashes are
the logged seals. They read and write the volume `pqtl-v8-stage-d`, committing after every
checkpoint. Each computes the run fingerprint (run_stage_d.run_fingerprint: plan hash, logged
seals, script digest, package versions in the image, run token, environment digest);
`prepare_remote` refuses a non-empty /vol/work written under another fingerprint and the others
refuse unless /vol/work, its plan and every checkpoint they resume carry this run's. The
environment digest (stage_d/environment.py) is computed in the container from MODAL_IMAGE_ID, the
OS release, the architecture, the installed libopenblas* packages, the BLAS PyTensor and NumPy
link and every installed distribution, so a run resumed in a rebuilt image that differs in any of
them is refused. `--what env` writes the components and the digest of the current image.

Commit identity: the stage code and PREREG.md are committed before the image is built. The image
carries the commit of a clean tree as V8_REPO_COMMIT (v8_run_guard.baked_commit; "" when
`git status --porcelain` is not empty for stages/ and PREREG.md), and each of the four real-run
functions calls `run_commit(None)`, which refuses in an image built from an unclean tree.
`assemble_remote` writes the commit to results.json and INPUTS.tsv.
"""
import json
import os
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import modal

from v8_run_guard import COMMIT_ENV, baked_commit

OPENBLAS_APT = ("libopenblas-dev=0.3.21+ds-4", "libopenblas-pthread-dev=0.3.21+ds-4", "libopenblas0-pthread=0.3.21+ds-4")
BLAS_ENV = {"PYTENSOR_FLAGS": "blas__ldflags=-lopenblas"}
image = modal.Image.debian_slim(python_version="3.12").apt_install(*OPENBLAS_APT).pip_install(
    "pymc==5.28.5", "pytensor==2.38.2", "arviz==0.23.4", "numpy==2.5.3", "scipy==1.18.1",
    "pandas==3.0.6", "statsmodels==0.15.0", "tqdm==4.70.1", "matplotlib==3.11.2", "pydantic==2.13.5",
    "pytest==9.1.1", "prereg==0.4.2", "provenance-core==0.4.2",
).env({"PYTHONPATH": "/root/stage_d_root:/root/power", "V8_PREREG": "/root/exp/PREREG.md", **BLAS_ENV})
PREREG_IMAGE = Path("/root/exp/PREREG.md")
if modal.is_local():
    HERE = Path(__file__).resolve().parent
    POWER = HERE.parents[1] / "power" / "power_v9"
    image = (image.add_local_dir(str(HERE / "stage_d"), "/root/stage_d_root/stage_d", copy=True,
                                 ignore=["__pycache__"])
             .add_local_dir(str(HERE / "tests"), "/root/stage_d_root/tests", copy=True, ignore=["__pycache__"])
             .add_local_file(str(HERE / "run_stage_d.py"), "/root/stage_d_root/run_stage_d.py", copy=True)
             .add_local_file(str(HERE / "pyproject.toml"), "/root/stage_d_root/pyproject.toml", copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_run_guard.py"), "/root/stage_d_root/v8_run_guard.py",
                             copy=True)
             .add_local_file(str(HERE.parent / "run_guard" / "v8_manifest.py"), "/root/stage_d_root/v8_manifest.py",
                             copy=True)
             .add_local_file(str(HERE.parents[1] / "PREREG.md"), str(PREREG_IMAGE), copy=True)
             .add_local_file(str(POWER / "power_v9_logic.py"), "/root/power/power_v9_logic.py", copy=True)
             .add_local_file(str(POWER / "power_v9_fast.py"), "/root/power/power_v9_fast.py", copy=True)
             .add_local_file(str(POWER / "power_v9_calibrate.py"), "/root/power/power_v9_calibrate.py", copy=True)
             .add_local_file(str(POWER / "pytensor_blas_check.py"), "/root/power/pytensor_blas_check.py", copy=True))
    for stage in ("A", "B", "C"):
        if (HERE.parent / stage / "output" / "MANIFEST.tsv").is_file():
            image = image.add_local_dir(str(HERE.parent / stage / "output"), f"/root/stages/{stage}/output", copy=True)
    image = image.env({COMMIT_ENV: baked_commit(HERE.parents[3])})

with image.imports():
    from pytensor_blas_check import check as blas_check_in_image
    from run_stage_d import package_versions, run_fingerprint, run_fits
    from stage_d.bayes import fit
    from stage_d.constants import SamplerSettings
    from stage_d.environment import collect as collect_environment
    from stage_d.fingerprint import canonical_sha256
    from stage_d.fitting import run_fit
    from stage_d.guard import logged_seals
    from stage_d.pipeline import assemble, load_plan, load_planned_design, prepare, run_frequentist
    from stage_d.synthetic import make_tables, synthetic_fingerprint, write_stage_outputs
    from v8_run_guard import run_commit

app = modal.App("pqtl-v8-stage-d", image=image)
vol = modal.Volume.from_name("pqtl-v8-stage-d", create_if_missing=True)
WORK = Path("/vol/work")
OUT = Path("/vol/output")
STAGES = Path("/root/stages")
SMOKE_SAMPLER = {"chains": 2, "draws": 300, "tune": 300, "target_accept": 0.95}


@app.function(cpu=4, memory=8192, timeout=3600)
def smoke() -> dict:
    """SYNTHETIC: fit the H1 model on a tiny simulated dataset through the full attempt sequence."""
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        expected = write_stage_outputs(root / "stages", *make_tables(seed=7, n_genes=40, hyps_per_gene=6))
        fp = synthetic_fingerprint(expected)
        plan = prepare(root / "stages", expected, root / "work", fp, samplers=(SamplerSettings(**SMOKE_SAMPLER),) * 2)
        spec = next(s for s in plan.fits if s.fit_id == "h1__S1__normal15")
        design, design_sha = load_planned_design(root / "work", plan, spec.design_id)
        final = run_fit(spec, design, root / "work" / "fits", fit, fingerprint=fp.digest, design_sha256=design_sha)
        att = final["attempts"][-1]
        components = collect_environment()      # after the fit, so PyTensor has compiled and loaded its BLAS
        return {"synthetic": True, "n": design.n, "columns": list(design.columns), "outcome": final["outcome"],
                "attempts": [a["attempt"] for a in final["attempts"]], "diagnostics": att["diagnostics"],
                "focal": att["summary"]["focal"], "marginal_rdi": att["summary"]["marginal"]["risk_difference_interaction"],
                "seconds": att["seconds"], "versions": package_versions(),
                "environment_sha256": canonical_sha256(components), "environment": components}


@app.function(cpu=4, memory=16384, timeout=3 * 3600)
def tests() -> str:
    """SYNTHETIC: the full pytest suite, including the PyMC tests skipped locally."""
    res = subprocess.run(["python", "-m", "pytest", "-q", "-p", "no:cacheprovider", "tests"], cwd="/root/stage_d_root",
                         capture_output=True, text=True)
    return res.stdout[-6000:] + res.stderr[-2000:]


@app.function(cpu=2, memory=4096, timeout=900)
def blas() -> dict:
    """PyTensor's BLAS link in this image (no data)."""
    return blas_check_in_image()


@app.function(cpu=2, memory=4096, timeout=900)
def runtime_environment() -> dict:
    """The environment components of the run fingerprint in this image and their digest (no data),
    with the names of the MODAL_* variables the container sets."""
    components = collect_environment()
    return {"environment_sha256": canonical_sha256(components), "environment": components,
            "modal_variables": sorted(k for k in os.environ if k.startswith("MODAL_"))}


@app.function(cpu=2, memory=8192, timeout=3600, volumes={"/vol": vol})
def prepare_remote(run_token: str) -> dict:
    expected = logged_seals(PREREG_IMAGE, run_token, STAGES)
    run_commit(None)
    vol.reload()
    plan = prepare(STAGES, expected, WORK, run_fingerprint(expected, run_token))
    vol.commit()
    return {"designs": list(plan.designs), "fits": [s.fit_id for s in plan.fits]}


@app.function(cpu=2, memory=8192, timeout=24 * 3600, retries=2, volumes={"/vol": vol})
def frequentist(design_id: str, run_token: str) -> str:
    expected = logged_seals(PREREG_IMAGE, run_token, STAGES)
    run_commit(None)
    vol.reload()
    run_frequentist(WORK, run_fingerprint(expected, run_token), [design_id], on_commit=vol.commit)
    vol.commit()
    return design_id


@app.function(cpu=4, memory=16384, timeout=24 * 3600, retries=2, volumes={"/vol": vol})
def fit_one(fit_id: str, run_token: str) -> dict:
    expected = logged_seals(PREREG_IMAGE, run_token, STAGES)
    run_commit(None)
    vol.reload()
    out = run_fits(WORK, run_fingerprint(expected, run_token), [fit_id], on_commit=vol.commit)
    vol.commit()
    return out


@app.function(cpu=2, memory=16384, timeout=4 * 3600, volumes={"/vol": vol})
def assemble_remote(run_token: str) -> dict:
    expected = logged_seals(PREREG_IMAGE, run_token, STAGES)
    commit = run_commit(None)
    vol.reload()
    res = assemble(STAGES, expected, WORK, OUT, fp=run_fingerprint(expected, run_token),
                   runner_path=Path("/root/stage_d_root/run_stage_d.py"), package_versions=package_versions(),
                   repo_commit=commit)
    vol.commit()
    return res["decisions"]


@app.function(cpu=1, memory=2048, timeout=600, volumes={"/vol": vol})
def plan_status(run_token: str) -> dict:
    expected = logged_seals(PREREG_IMAGE, run_token, STAGES)
    vol.reload()
    plan = load_plan(WORK, run_fingerprint(expected, run_token))
    return {s.fit_id: (WORK / "fits" / s.fit_id / "final.json").exists() for s in plan.fits}


@app.local_entrypoint()
def main(what: str = "smoke"):
    """SYNTHETIC checks only: `what` is `smoke`, `tests`, `blas` or `env`; the result is written to
    work/, after any earlier result of the same name is moved to work/superseded/."""
    result = {"smoke": smoke.remote, "blas": blas.remote, "env": runtime_environment.remote,
              "tests": lambda: {"pytest": tests.remote()}}[what]()
    out = Path(__file__).resolve().parent / "work" / f"modal_{what}.json"
    out.parent.mkdir(exist_ok=True)
    if out.exists():
        written = datetime.fromtimestamp(out.stat().st_mtime, timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        old = out.parent / "superseded" / f"{out.stem}_{written}.json"
        old.parent.mkdir(exist_ok=True)
        if old.exists():
            raise FileExistsError(f"{old} already exists; not overwriting an archived result")
        out.rename(old)
    out.write_text(json.dumps(result, indent=2, default=str))
    print(out)
