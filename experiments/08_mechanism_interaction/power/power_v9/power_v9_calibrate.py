"""Exact-model check for the V8 power simulation (pure logic; no Modal dependency).

For selected scenarios, each repetition simulates one dataset (power_v9_logic.simulate) and
analyses it twice: with the fast two-way cluster-robust logistic model (power_v9_fast) and
with the registered H1 model fitted by NUTS:

    advance ~ S + A + S:A + (1 + S | gene, uncorrelated) + (1 | indication)
              + multiple-membership drug intercept (weights 1/k over the hypothesis's programs)

Priors as registered: Normal(0, 1.5) on fixed coefficients, Normal(0, 2.5) on the intercept,
Half-Normal(1) on random-effect SDs; non-centered random effects; target acceptance 0.95.
Covariates (platform, N_eff, oncology) are omitted because the simulation does not generate
them. The output per repetition is both decisions, so agreement is measured dataset by dataset.
"""
import json
from dataclasses import asdict, replace
from math import log
from pathlib import Path

import numpy as np
import pymc as pm

from power_v9_fast import Z_ONE_SIDED, component_clusters, fit_two_way
from power_v9_logic import load_structure, scenario_grid, simulate

SEED = 20260930
LOG_SESOI = log(1.5)


def calibration_scenarios():
    base = scenario_grid()[0]
    common = dict(re_sd=0.7, slope_sd=0.7, p_support=0.06, aligned_support_mult=1.0,
                  or_blocking=1.0, base_advance=0.30)
    return {
        "null": replace(base, scenario_id=9001, or_ratio=1.0, **common),
        "or2": replace(base, scenario_id=9002, or_ratio=2.0, **common),
        "or3": replace(base, scenario_id=9003, or_ratio=3.0, **common),
        # round-4 representative cells
        "or1p5": replace(base, scenario_id=9004, or_ratio=1.5, **common),
        "or3_ps03": replace(base, scenario_id=9005, or_ratio=3.0, **{**common, "p_support": 0.03}),
        "or3_ps10": replace(base, scenario_id=9006, or_ratio=3.0, **{**common, "p_support": 0.10}),
        "or3_lowhet": replace(base, scenario_id=9007, or_ratio=3.0, **{**common, "re_sd": 0.3, "slope_sd": 0.3}),
        "or3_highhet": replace(base, scenario_id=9008, or_ratio=3.0, **{**common, "re_sd": 1.0, "slope_sd": 0.7}),
        "null_ps03": replace(base, scenario_id=9009, or_ratio=1.0, **{**common, "p_support": 0.03}),
        "null_highhet": replace(base, scenario_id=9010, or_ratio=1.0, **{**common, "re_sd": 1.0, "slope_sd": 0.7}),
    }


def fit_exact(y, s, st, draws=1000, tune=1000, chains=4, seed=0, mode="h1"):
    w = st.drug_matrix.toarray()
    with pm.Model():
        b0 = pm.Normal("b0", 0, 2.5)
        k = 1 if mode == "h2" else 3
        b = pm.Normal("b", 0, 1.5, shape=k)
        other = st.cns if mode == "h4" else st.aligned
        sd = pm.HalfNormal("sd", 1.0, shape=4)
        zg = pm.Normal("zg", 0, 1, shape=st.n_gene)
        zw = pm.Normal("zw", 0, 1, shape=st.n_gene)
        zi = pm.Normal("zi", 0, 1, shape=st.n_indication)
        zd = pm.Normal("zd", 0, 1, shape=w.shape[1])
        fixed = b[0] * s if mode == "h2" else b[0] * s + b[1] * other + b[2] * s * other
        eta = (b0 + fixed
               + sd[0] * zg[st.gene] + sd[1] * zw[st.gene] * s
               + sd[2] * zi[st.indication] + sd[3] * pm.math.dot(w, zd))
        pm.Bernoulli("y", logit_p=eta, observed=y)
        idata = pm.sample(draws=draws, tune=tune, chains=chains, cores=chains,
                          target_accept=0.95, random_seed=seed, progressbar=False)
    post = idata.posterior["b"].values[..., -1].ravel()
    if mode == "h4":
        post = -post   # H4 predicts a weaker association in CNS indications
    diverging = int(idata.sample_stats["diverging"].values.sum())
    return {"pr_gt0": float(np.mean(post > 0)), "pr_gt_sesoi": float(np.mean(post > LOG_SESOI)),
            "median": float(np.median(post)), "divergences": diverging}


NULLS = {"h2_null": ("h2", dict(or_blocking=1.0, or_ratio=1.0)),
         "h4_null": ("h4", dict(or_blocking=1.25, cns_ratio=1.0))}


def run_one(label, rep, membership_dir: Path, out_path: Path):
    if label in NULLS:
        mode, kw = NULLS[label]
        sc = replace(calibration_scenarios()["null"], scenario_id=9100 + list(NULLS).index(label), **kw)
        fname = "s1_membership.csv"
    else:
        mode, sc, fname = "h1", calibration_scenarios()[label], "h1_membership.csv"
    st = load_structure(membership_dir / fname)
    comp = component_clusters(membership_dir / fname)
    rng = np.random.default_rng([SEED, {"h1": 1, "h2": 2, "h4": 4}[mode], sc.scenario_id, rep])
    y, s = simulate(rng, st, sc, mode)
    one = np.ones_like(s)
    other = st.cns if mode == "h4" else st.aligned
    x = np.column_stack([one, s]) if mode == "h2" else np.column_stack([one, s, other, s * other])
    beta, se = fit_two_way(x, y, comp, st.indication)
    z = float(beta[-1] / se[-1]) * (-1 if mode == "h4" else 1)
    exact = fit_exact(y, s, st, seed=rep, mode=mode)
    rec = {"label": label, "mode": mode, "scenario": asdict(sc), "rep": rep,
           "fast_z": z, "fast_reject_nominal": z >= Z_ONE_SIDED,
           "exact": exact, "exact_reject": exact["pr_gt0"] >= 0.95,
           "n_supportive": int(s.sum())}
    out_path.write_text(json.dumps(rec, indent=2))
    return rec
