"""Stage D entry point. Runs on Linux in the image of modal_stage_d.py (PyMC 5.28.5, PyTensor
2.38.2); PyMC does not install on the local Intel Mac.

Refuses to run unless PREREG.md passes `prereg check`, its log holds the entry
`RUN_START stage=D token=<token>` for the `--run-token` given and, before it, the entries
`SEAL stage=A manifest_sha256=<sha256>` (and B, C) whose values are the sha256 of the
MANIFEST.tsv files under --stages-root (shared guard, stages/run_guard/v8_run_guard.py). The
logged seals, not a caller-supplied value, are the expected hashes stage_d.guard then checks every
listed output file against.

All work under --work is bound to the run fingerprint (stage_d/fingerprint.py): the frozen plan
hash, the three logged seals, the digest of this file, the stage_d modules, the shared guard
modules and power_v9_fast.py (each by relative path, length and bytes), the installed package
versions, the run token, and the digest of the runtime environment (stage_d/environment.py: Modal
image id, OS release, architecture, OpenBLAS apt packages, linked BLAS, every installed
distribution), so work written in one image is refused in a rebuilt image that differs. The run also records one repository commit
(v8_run_guard.run_commit: the value baked into the Modal image, or HEAD of a clean working tree) in
results.json and INPUTS.tsv; the commit is not part of the fingerprint. `prepare` refuses a non-empty work directory
written under another fingerprint; the later phases refuse unless the directory, the plan and
every checkpoint they resume carry this run's.

Phases (each resumable; `all` runs them in order):
  prepare      guard, join, sets S1-S22, gates, designs, fit list        -> work/plan.json
  frequentist  two-way cluster-robust check, size correction, bootstrap  -> work/freq/
  fit          Bayesian fits with the rerun / reduced-model sequence    -> work/fits/
  assemble     decisions, results.json, tables, figures, MANIFEST.tsv    -> output/

    PYTHONPATH=.:../../power/power_v9 python run_stage_d.py --stages-root .. --run-token <token> --phase all
"""
import argparse
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import power_v9_fast

from stage_d.bayes import fit
from stage_d.constants import PLAN_SHA256
from stage_d.environment import collect as collect_environment
from stage_d.fingerprint import RunFingerprint
from stage_d.fitting import run_fit
from stage_d.guard import logged_seals
from stage_d.pipeline import assemble, load_plan, load_planned_design, prepare, run_frequentist, script_sha256
from v8_run_guard import run_commit

HERE = Path(__file__).resolve().parent
PREREG = HERE.parent.parent / "PREREG.md"
REPO = PREREG.parent.parent.parent
PACKAGES = ("pymc", "pytensor", "arviz", "numpy", "scipy", "pandas", "statsmodels", "matplotlib", "pydantic")


def run_fits(work_dir: Path, fp: RunFingerprint, fit_ids: list[str] | None = None, on_commit=None) -> dict:
    plan = load_plan(work_dir, fp)
    specs = [s for s in plan.fits if fit_ids is None or s.fit_id in fit_ids]
    out = {}
    for spec in specs:
        design, design_sha = load_planned_design(work_dir, plan, spec.design_id)
        final = run_fit(spec, design, work_dir / "fits", fit, on_commit, fingerprint=fp.digest, design_sha256=design_sha)
        out[spec.fit_id] = final["outcome"]
        print(f"{datetime.now(timezone.utc).isoformat()} {spec.fit_id}: {final['outcome']}")
    return out


def package_versions() -> dict[str, str]:
    return {p: version(p) for p in PACKAGES}


def run_fingerprint(seals: dict[str, str], run_token: str) -> RunFingerprint:
    """The fingerprint of this run: frozen plan hash, logged seals, the digest of this file, the
    stage_d modules and power_v9_fast.py, the Python and package versions installed here, the
    run token, and the runtime environment collected here (stage_d.environment.collect)."""
    digest = script_sha256(HERE / "stage_d", Path(__file__).resolve(),
                           extra={"power_v9/power_v9_fast.py": Path(power_v9_fast.__file__).resolve()})
    return RunFingerprint(plan_sha256=PLAN_SHA256, seals=seals, script_sha256=digest,
                          pins={"python": platform.python_version(), **package_versions()}, run_token=run_token,
                          environment=collect_environment())


def main(argv: list[str] | None = None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stages-root", required=True, type=Path)
    ap.add_argument("--run-token", required=True, help="token of the PREREG.md entry 'RUN_START stage=D token=<token>'")
    ap.add_argument("--prereg", type=Path, default=PREREG)
    ap.add_argument("--work", type=Path, default=HERE / "work" / "run",
                    help="work directory of this run; bound to the run fingerprint")
    ap.add_argument("--out", type=Path, default=HERE / "output")
    ap.add_argument("--phase", choices=("prepare", "frequentist", "fit", "assemble", "all"), default="all")
    ap.add_argument("--fit-id", nargs="*", default=None)
    args = ap.parse_args(argv)
    expected = logged_seals(args.prereg, args.run_token, args.stages_root)
    commit = run_commit(REPO)
    fp = run_fingerprint(expected, args.run_token)
    if args.phase in ("prepare", "all"):
        plan = prepare(args.stages_root, expected, args.work, fp)
        print(f"{datetime.now(timezone.utc).isoformat()} prepared {len(plan.designs)} designs, {len(plan.fits)} fits")
    if args.phase in ("frequentist", "all"):
        run_frequentist(args.work, fp)
    if args.phase in ("fit", "all"):
        run_fits(args.work, fp, args.fit_id)
    if args.phase in ("assemble", "all"):
        res = assemble(args.stages_root, expected, args.work, args.out, fp=fp, runner_path=Path(__file__).resolve(),
                       package_versions=package_versions(), repo_commit=commit)
        print(json.dumps(res["decisions"], indent=2))


if __name__ == "__main__":
    main()
