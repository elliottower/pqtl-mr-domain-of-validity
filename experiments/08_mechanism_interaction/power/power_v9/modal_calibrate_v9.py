"""Thin Modal wrapper for the exact-model check (logic in power_v9_calibrate.py).

    modal deploy experiments/08_mechanism_interaction/power/power_v9/modal_calibrate_v9.py
    uv run --with modal==1.4.3 python -c "import modal; modal.Function.from_name('pqtl-v8-calibrate-v9','smoke').remote()"
    uv run --with modal==1.4.3 python experiments/08_mechanism_interaction/power/power_v9/launch_calibrate_v9.py

One repetition per call; each writes its own JSON to the volume and commits, and a call whose
file already exists returns without refitting.

BLAS: PyTensor links against OpenBLAS (Debian bookworm libopenblas-dev 0.3.21+ds-4, pthread
variant) through PYTENSOR_FLAGS=blas__ldflags=-lopenblas, so its C BLAS ops (CGemv, Dot22)
are used instead of the degraded fallback. Modal sets OPENBLAS_NUM_THREADS to the function's CPU
request at container start, so the image does not set it. Checked by pytensor_blas_check.check()
inside the image.

    uv run --with modal==1.4.3 modal run experiments/08_mechanism_interaction/power/power_v9/modal_calibrate_v9.py::blas
"""
from pathlib import Path

import modal

OPENBLAS_APT = ("libopenblas-dev=0.3.21+ds-4", "libopenblas-pthread-dev=0.3.21+ds-4", "libopenblas0-pthread=0.3.21+ds-4")
BLAS_ENV = {"PYTENSOR_FLAGS": "blas__ldflags=-lopenblas"}
image = modal.Image.debian_slim(python_version="3.12").apt_install(*OPENBLAS_APT).pip_install(
    "pymc==5.28.5", "pytensor==2.38.2", "arviz==0.23.4", "numpy==2.5.3", "scipy==1.18.1",
    "pandas==3.0.6", "statsmodels==0.15.0", "tqdm==4.70.1",
).env({"PYTHONPATH": "/root/power", **BLAS_ENV})
if modal.is_local():
    HERE = Path(__file__).resolve().parent
    FEAS = HERE.parents[1] / "feasibility" / "v4_round3"
    image = (
        image.add_local_file(str(HERE / "power_v9_logic.py"), "/root/power/power_v9_logic.py", copy=True)
        .add_local_file(str(HERE / "power_v9_fast.py"), "/root/power/power_v9_fast.py", copy=True)
        .add_local_file(str(HERE / "power_v9_calibrate.py"), "/root/power/power_v9_calibrate.py", copy=True)
        .add_local_file(str(HERE / "pytensor_blas_check.py"), "/root/power/pytensor_blas_check.py", copy=True)
        .add_local_file(str(FEAS / "h1_membership.csv"), "/root/data/h1_membership.csv", copy=True)
        .add_local_file(str(FEAS / "s1_membership.csv"), "/root/data/s1_membership.csv", copy=True)
    )

app = modal.App("pqtl-v8-calibrate-v9", image=image)
vol = modal.Volume.from_name("pqtl-v8-power", create_if_missing=True)
OUT = Path("/vol/calibrate_v9")


@app.function(cpu=4, memory=8192, timeout=4 * 3600, retries=2, volumes={"/vol": vol})
def run_rep(label: str, rep: int):
    from power_v9_calibrate import run_one

    vol.reload()
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{label}_{rep:03d}.json"
    if out.exists():
        return {"label": label, "rep": rep, "skipped": True}
    rec = run_one(label, rep, Path("/root/data"), out)
    vol.commit()
    return {"label": label, "rep": rep, "exact_reject": rec["exact_reject"],
            "fast_z": rec["fast_z"], "divergences": rec["exact"]["divergences"]}


@app.function(cpu=2, memory=4096, timeout=900)
def blas() -> dict:
    """PyTensor's BLAS link in this image (no data)."""
    from pytensor_blas_check import check

    return check()


@app.function(cpu=4, memory=8192, timeout=3600)
def smoke():
    import time

    import numpy as np

    from power_v9_calibrate import calibration_scenarios, fit_exact
    from power_v9_logic import load_structure, simulate

    st = load_structure(Path("/root/data/h1_membership.csv"))
    y, s = simulate(np.random.default_rng(1), st, calibration_scenarios()["or3"], "h1")
    t = time.time()
    r = fit_exact(y, s, st, draws=100, tune=100, chains=2, seed=1)
    return {**r, "seconds": round(time.time() - t)}
