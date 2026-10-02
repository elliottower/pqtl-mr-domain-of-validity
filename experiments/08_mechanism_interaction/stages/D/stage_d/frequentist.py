"""Frequentist check beside every Bayesian model (PREREG §Statistical models, Frequentist check).

Logistic regression with the model's fixed effects and two-way cluster-robust standard errors
(clusters: connected components of the gene-drug-program graph, and indications), reusing
power_v9_fast (`_fit_logit`, `_meat`, the Cameron-Gelbach-Miller combination of fit_two_way).
Its one-sided z is compared with a size-corrected critical value: the 95th percentile of z in a
matched no-effect simulation on the realized structure and realized evidence regressor, with the
data-generating model of power_v9_logic.simulate (gene intercept, gene slope on the evidence
regressor, indication intercept, multiple-membership program intercept; base advancement 0.30)
over the power_v9 null grid (random-effect SD 0.3/0.7/1.0 x slope SD 0.3/0.7 x main-effect OR
0.8/1.0/1.25; H2 kinds OR 1.0 only). The reported critical value is the largest over the grid;
the reference cell (0.7, 0.7, 1.0) is reported beside it. The simulation reads no outcome.

A bootstrap over whole gene-drug-program components (10,000 resamples) gives an interval for
the focal quantity. Both loops checkpoint to JSONL every chunk and resume from disk.
"""
import hashlib
import json
from math import erf, log, sqrt
from pathlib import Path

import numpy as np
from power_v9_fast import Z_ONE_SIDED, _fit_logit, _meat

from stage_d.constants import (BOOTSTRAP_CHUNK, MAX_ABS_LOGIT_COEF, N_BOOTSTRAP, N_NULL_REPS, NULL_BASE_ADVANCE,
                               NULL_CHUNK, NULL_MAIN_ORS, NULL_RE_SDS, NULL_REFERENCE_CELL, NULL_SLOPE_SDS, SEED)
from stage_d.designs import Design


class FitFailure(ArithmeticError):
    pass


def _cdf(z: float) -> float:
    return 0.5 * (1 + erf(z / sqrt(2)))


def design_seed(design_id: str) -> int:
    return int(hashlib.sha256(design_id.encode()).hexdigest()[:8], 16)


def two_way_covariance(x, y, g1, g2):
    """power_v9_fast.fit_two_way, returning the full covariance instead of its diagonal."""
    beta, p = _fit_logit(x, y)
    bread = np.linalg.inv(x.T @ (x * (p * (1 - p))[:, None]))
    s = x * (y - p)[:, None]
    g12 = g1.astype(np.int64) * (g2.max() + 1) + g2
    v = bread @ (_meat(s, g1) + _meat(s, g2) - _meat(s, g12)) @ bread
    evals, evecs = np.linalg.eigh(v)
    return beta, evecs @ np.diag(np.clip(evals, 0, None)) @ evecs.T


def _with_intercept(design: Design) -> tuple[np.ndarray, np.ndarray]:
    return np.column_stack([np.ones(design.n), design.X]), np.concatenate([[0.0], design.focal])


def focal_z(design: Design, y: np.ndarray) -> dict:
    """Estimate, SE and signed z of the focal contrast; raises FitFailure as power_v9_fast does
    (constant outcome, singular fit, non-finite or |estimate| > 15)."""
    if y.min() == y.max():
        raise FitFailure("constant outcome")
    x, c = _with_intercept(design)
    try:
        beta, v = two_way_covariance(x, y, design.component, design.indication)
    except np.linalg.LinAlgError as err:
        raise FitFailure(str(err)) from err
    est = float(c @ beta)
    se = float(np.sqrt(c @ v @ c))
    if not np.isfinite(est) or not np.isfinite(se) or se <= 0 or abs(est) > MAX_ABS_LOGIT_COEF:
        raise FitFailure(f"estimate {est}, se {se}")
    return {"estimate": est, "se": se, "z": design.focal_sign * est / se,
            "coefficients": dict(zip(("intercept",) + design.columns, map(float, beta))),
            "se_all": dict(zip(("intercept",) + design.columns, map(float, np.sqrt(np.diag(v)))))}


def main_regressor(design: Design) -> np.ndarray:
    return design.X[:, 0] + design.X[:, 1] if design.kind == "h3" else design.X[:, design.ev_col]


def simulate_null(rng, design: Design, re_sd: float, slope_sd: float, main_log_or: float,
                  base: float = NULL_BASE_ADVANCE) -> np.ndarray:
    """Outcome under no focal effect on the realized structure; draws in the order of
    power_v9_logic.simulate after S (u, w, v, program effects, outcome)."""
    u = rng.normal(0, re_sd, design.n_gene)[design.gene]
    w = rng.normal(0, slope_sd, design.n_gene)[design.gene]
    v = rng.normal(0, re_sd, design.n_indication)[design.indication]
    z = design.W @ rng.normal(0, re_sd, design.W.shape[1])
    lo = log(base / (1 - base)) + u + w * design.slope + v + z + main_regressor(design) * main_log_or
    return (rng.random(design.n) < 1 / (1 + np.exp(-lo))).astype(float)


def null_cells(kind: str) -> list[tuple[float, float, float]]:
    mains = (1.0,) if kind in ("h2", "h2_E") else NULL_MAIN_ORS
    return [(r, s, m) for r in NULL_RE_SDS for s in NULL_SLOPE_SDS for m in mains]


def _resume(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.open()] if path.exists() else []


def _append(path: Path, recs: list[dict], on_commit) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write("\n".join(json.dumps(r) for r in recs) + "\n")
    if on_commit:
        on_commit()


def run_null_cell(design: Design, cell_index: int, cell, out_dir: Path, reps: int = N_NULL_REPS, on_commit=None) -> list[dict]:
    path = out_dir / f"{design.design_id}__null_{cell_index:02d}.jsonl"
    done = _resume(path)
    buf = []
    for rep in range(len(done), reps):
        rng = np.random.default_rng([SEED, design_seed(design.design_id), cell_index, rep])
        y = simulate_null(rng, design, cell[0], cell[1], log(cell[2]))
        try:
            z = focal_z(design, y)["z"]
        except FitFailure:
            z = None
        buf.append({"rep": rep, "z": z})
        if len(buf) == NULL_CHUNK or rep == reps - 1:
            _append(path, buf, on_commit)
            buf = []
    return _resume(path)


def size_correction(design: Design, out_dir: Path, reps: int = N_NULL_REPS, on_commit=None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    cells = []
    for i, cell in enumerate(null_cells(design.kind)):
        recs = run_null_cell(design, i, cell, out_dir, reps, on_commit)
        zs = np.array([r["z"] for r in recs if r["z"] is not None])
        cells.append({"re_sd": cell[0], "slope_sd": cell[1], "main_or": cell[2], "reps": len(recs),
                      "failed": len(recs) - len(zs),
                      "size_nominal": float(np.sum(zs >= Z_ONE_SIDED) / len(recs)),
                      "crit_z": float(np.quantile(zs, 0.95)) if len(zs) else float("nan"), "_z": zs})
    crits = [c["crit_z"] for c in cells]
    ref = next((c for c in cells if (c["re_sd"], c["slope_sd"], c["main_or"]) == NULL_REFERENCE_CELL), None)
    return {"cells": cells, "crit_z": float(np.nanmax(crits)) if not np.all(np.isnan(crits)) else float("nan"),
            "crit_z_reference_cell": ref["crit_z"] if ref else None}


def component_bootstrap(design: Design, out_path: Path, n_boot: int = N_BOOTSTRAP, on_commit=None) -> dict:
    x, c = _with_intercept(design)
    comps = np.unique(design.component)
    members = [np.flatnonzero(design.component == k) for k in comps]
    done = _resume(out_path)
    buf = []
    for rep in range(len(done), n_boot):
        rng = np.random.default_rng([SEED, design_seed(design.design_id), 999, rep])
        idx = np.concatenate([members[j] for j in rng.integers(0, len(comps), len(comps))])
        yb = design.y[idx]
        est = None
        if yb.min() != yb.max():
            try:
                beta, _ = _fit_logit(x[idx], yb)
                e = float(c @ beta)
                est = e if np.isfinite(e) and abs(e) <= MAX_ABS_LOGIT_COEF else None
            except np.linalg.LinAlgError:
                est = None
        buf.append({"rep": rep, "estimate": est})
        if len(buf) == BOOTSTRAP_CHUNK or rep == n_boot - 1:
            _append(out_path, buf, on_commit)
            buf = []
    vals = np.array([r["estimate"] for r in _resume(out_path) if r["estimate"] is not None])
    q = (lambda p: float(np.quantile(vals, p))) if len(vals) else (lambda p: float("nan"))
    return {"resamples": n_boot, "failed": n_boot - len(vals), "components": len(comps),
            "q025": q(0.025), "q05": q(0.05), "median": q(0.5), "q95": q(0.95), "q975": q(0.975)}


def frequentist_check(design: Design, work_dir: Path, reps: int = N_NULL_REPS, n_boot: int = N_BOOTSTRAP,
                      on_commit=None) -> dict:
    sc = size_correction(design, work_dir / "null", reps, on_commit)
    null_z = [c.pop("_z") for c in sc["cells"]]
    try:
        obs = focal_z(design, design.y)
    except FitFailure as err:
        return {"fitted": False, "reason": str(err), "size_correction": sc}
    z = obs["z"]
    p_corr = max(float(np.sum(zs >= z) / c["reps"]) for zs, c in zip(null_z, sc["cells"]))
    boot = component_bootstrap(design, work_dir / "bootstrap" / f"{design.design_id}.jsonl", n_boot, on_commit)
    return {"fitted": True, **obs, "p_one_sided_unadjusted": 1 - _cdf(z),
            "reject_nominal": z >= Z_ONE_SIDED, "reject_size_corrected": z >= sc["crit_z"],
            "p_one_sided_size_corrected": p_corr, "size_correction": sc, "bootstrap": boot,
            "clusters": {"components": int(len(np.unique(design.component))),
                         "indications": int(design.n_indication)}}
