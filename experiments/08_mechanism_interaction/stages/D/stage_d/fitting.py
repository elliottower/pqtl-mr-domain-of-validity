"""Runs one fit specification through the attempt sequence (fitspec.ATTEMPTS), writing one JSON
per attempt and resuming from the attempts already on disk. The sampler is passed in as
`fit_fn` (stage_d.bayes.fit on Linux), so this module does not import PyMC.

A NUTS run cannot be resumed mid-chain, so the checkpoint unit is one attempt. Every attempt file
and final.json carry the run fingerprint and the sha256 of the design file fitted
(stage_d/fingerprint.py); an attempt on disk under another fingerprint or another design raises
StaleWork and is not resumed."""
import json
import time
from pathlib import Path

import numpy as np

from stage_d.constants import SEED
from stage_d.designs import Design
from stage_d.fingerprint import check_record
from stage_d.fitspec import FitSpec, fit_outcome, next_attempt
from stage_d.frequentist import design_seed
from stage_d.gates import sampler_passed
from stage_d.posterior import PRIORS, summarize_fit


def _write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=2, default=float))
    tmp.replace(path)


def load_attempts(fit_dir: Path, fingerprint: str, design_sha256: str) -> list[dict]:
    out = []
    for name in ("full", "full_rerun", "reduced", "reduced_rerun"):
        p = fit_dir / f"{name}.json"
        if not p.exists():
            break
        rec = json.loads(p.read_text())
        check_record(rec, str(p), fingerprint, design_sha256)
        out.append(rec)
    return out


def run_fit(spec: FitSpec, design: Design, fits_dir: Path, fit_fn, on_commit=None, *, fingerprint: str,
            design_sha256: str) -> dict:
    """`fingerprint` is the run fingerprint digest, `design_sha256` the sha256 of the design file
    `design` was loaded from."""
    fit_dir = fits_dir / spec.fit_id
    done = load_attempts(fit_dir, fingerprint, design_sha256)
    prior = PRIORS[spec.prior]
    while (step := next_attempt(done)) is not None:
        name, which, keep_slope = step
        sampler = spec.sampler(which)
        seed = design_seed(spec.fit_id) + len(done)
        t0 = time.time()
        draws, diag, re = fit_fn(design, prior, sampler, keep_slope, seed)
        rng = np.random.default_rng([SEED, design_seed(spec.fit_id), len(done)])
        rec = {"attempt": name, "fit_id": spec.fit_id, "fingerprint": fingerprint, "design_sha256": design_sha256,
               "sampler": sampler.model_dump(), "keep_slope": keep_slope,
               "random_effects": re, "seed": seed, "seconds": round(time.time() - t0, 1),
               "diagnostics": diag.model_dump(), "passed": sampler_passed(diag),
               "summary": summarize_fit(draws, design, prior, re, rng)}
        _write_json(fit_dir / f"{name}.json", rec)
        if on_commit:
            on_commit()
        done.append(rec)
    final = {"spec": spec.model_dump(), "fingerprint": fingerprint, "design_sha256": design_sha256,
             "outcome": fit_outcome(done), "attempts": done}
    _write_json(fit_dir / "final.json", final)
    if on_commit:
        on_commit()
    return final
