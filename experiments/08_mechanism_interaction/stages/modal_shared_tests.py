"""Synthetic test suites that belong to no single stage image (no study data is read).

    run_guard    the shared run guard and manifest suite (stages/run_guard/tests), in an image
                 holding only its pinned dependencies (run_guard/pyproject.toml)
    cross_stage  stages/D/tests/test_cross_stage_seal.py, which runs the real stage A, B and C
                 writers through `uv run --project` and then stage D's guard, so it needs all four
                 stage projects at once: a combined image with uv and each project synced from its
                 own uv.lock

    uv run --with modal==1.4.3 modal run experiments/08_mechanism_interaction/stages/modal_shared_tests.py --what run_guard
    uv run --with modal==1.4.3 modal run experiments/08_mechanism_interaction/stages/modal_shared_tests.py --what cross_stage

The report is printed and written to run_guard/work/modal_tests.json or D/work/modal_cross_stage.json;
a report already there is first moved to superseded/ beside it.
"""
import json
from datetime import datetime, timezone
from pathlib import Path

import modal

UV_VERSION = "0.11.3"
GUARD_PACKAGES = ["prereg==0.4.2", "provenance-core==0.4.2", "pytest==8.4.2"]
EXP = Path("/root/exp")
STAGES = EXP / "stages"
IGNORE = [".venv", "**/__pycache__", "**/.pytest_cache", "output", "work", "launch_log.jsonl"]
PROJECTS = ("A", "B", "C", "D")

guard_image = modal.Image.debian_slim(python_version="3.12").pip_install(*GUARD_PACKAGES).env(
    {"PYTHONPATH": str(STAGES / "run_guard")})
cross_image = modal.Image.debian_slim(python_version="3.12").pip_install(f"uv=={UV_VERSION}", *GUARD_PACKAGES).env(
    {"PYTHONPATH": str(STAGES / "run_guard"), "UV_LINK_MODE": "copy"})
if modal.is_local():
    HERE = Path(__file__).resolve().parent
    guard_image = (guard_image.add_local_dir(str(HERE / "run_guard"), str(STAGES / "run_guard"), copy=True, ignore=IGNORE)
                   .add_local_file(str(HERE.parent / "PREREG.md"), str(EXP / "PREREG.md"), copy=True))
    cross_image = (cross_image.add_local_dir(str(HERE / "run_guard"), str(STAGES / "run_guard"), copy=True, ignore=IGNORE)
                   .add_local_file(str(HERE.parent / "PREREG.md"), str(EXP / "PREREG.md"), copy=True)
                   .add_local_file(str(HERE / "INTERFACES.md"), str(STAGES / "INTERFACES.md"), copy=True))
    for _stage in PROJECTS:
        cross_image = cross_image.add_local_dir(str(HERE / _stage), str(STAGES / _stage), copy=True, ignore=IGNORE)
    cross_image = cross_image.run_commands(*(f"uv sync --frozen --project {STAGES / s}" for s in PROJECTS))

with guard_image.imports():
    from v8_test_report import run_pytest

with cross_image.imports():
    from v8_test_report import run_pytest  # noqa: F811  (the same helper, in the other image)

app = modal.App("pqtl-v8-shared-tests")


@app.function(image=guard_image, cpu=2, memory=4096, timeout=1800)
def run_guard() -> dict:
    """SYNTHETIC: the run guard and manifest suite."""
    return run_pytest(STAGES / "run_guard", [STAGES / "run_guard"])


@app.function(image=cross_image, cpu=4, memory=16384, timeout=3600)
def cross_stage() -> dict:
    """SYNTHETIC: the cross-stage seal test, in stage D's own uv environment."""
    return run_pytest(STAGES / "D", [STAGES / "D"], args=("tests/test_cross_stage_seal.py",),
                      command=("uv", "run", "--frozen", "--project", str(STAGES / "D"), "python", "-m", "pytest"))


@app.local_entrypoint()
def main(what: str):
    here = Path(__file__).resolve().parent
    fn, out = {"run_guard": (run_guard, here / "run_guard" / "work" / "modal_tests.json"),
               "cross_stage": (cross_stage, here / "D" / "work" / "modal_cross_stage.json")}[what]
    report = fn.remote()
    out.parent.mkdir(exist_ok=True)
    if out.exists():
        written = datetime.fromtimestamp(out.stat().st_mtime, timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        old = out.parent / "superseded" / f"{out.stem}_{written}.json"
        old.parent.mkdir(exist_ok=True)
        if old.exists():
            raise FileExistsError(f"{old} already exists; not overwriting an archived report")
        out.rename(old)
    out.write_text(json.dumps(report, indent=2))
    print(f"{what}: {report['summary']} (exit {report['returncode']}) -> {out}")
