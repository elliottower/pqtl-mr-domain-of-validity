# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy==2.1.3",
#     "statsmodels==0.14.4",
#     "pandas==2.2.3",
#     "scipy==1.14.1",
# ]
# ///
"""Design power for experiment 08 (V8). Uses no outcome, MR estimate or evidence state.

Inputs are counts from the outcome-blind feasibility file and assumed effect sizes. The
data-generating model has a gene random intercept; the analysis model is a logistic
regression with gene-clustered standard errors, a stand-in for the registered Bayesian
hierarchical model that is cheap enough to simulate.

    H1 (interaction): advance ~ S + A + S:A, one-sided test that the S:A coefficient > 0
    H2 (main effect): advance ~ S, one-sided test that the S coefficient > 0

S = supportive evidence state (1/0). A = abundance-aligned (1) vs function-blocking (0).

Usage:
    uv run experiments/08_mechanism_interaction/power/power_v8.py
"""
import itertools
import json
import math
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

HERE = Path(__file__).resolve().parent
FEASIBILITY = HERE.parent / "feasibility" / "coverage_counts.json"
OUT_JSON = HERE / "power_v8.json"
CHECKPOINT = HERE / "power_v8_checkpoint.jsonl"

SEED = 20260929
REPS = 1000
ALPHA_ONE_SIDED = 0.05
PAIRS_PER_GENE = 2.5
GENE_RE_SD = 0.7
SHARE_ALIGNED = 0.40

# Grid. n_h1 is the number of hypotheses entering the H1 model (aligned + blocking).
N_H1 = [200, 275, 350]
P_SUPPORT = [0.05, 0.08, 0.12]
BASE_ADVANCE = [0.20, 0.30]
OR_ALIGNED = [2.0, 3.0, 5.0]          # supportive-vs-not OR in the aligned stratum
H1_ELIGIBLE_FRACTION = 0.70           # share of H2 hypotheses that are aligned or blocking


def simulate(rng, n, p_support, base, or_aligned, or_blocking):
    n_genes = max(1, int(round(n / PAIRS_PER_GENE)))
    gene = rng.integers(0, n_genes, size=n)
    u = rng.normal(0.0, GENE_RE_SD, size=n_genes)[gene]
    a = (rng.random(n) < SHARE_ALIGNED).astype(float)
    s = (rng.random(n) < p_support).astype(float)
    b0 = np.log(base / (1 - base))
    lo = b0 + u + s * np.where(a == 1, np.log(or_aligned), np.log(or_blocking))
    y = (rng.random(n) < 1 / (1 + np.exp(-lo))).astype(float)
    return pd.DataFrame({"y": y, "s": s, "a": a, "gene": gene})


def one_sided_p(df, terms, coef):
    x = sm.add_constant(df[terms], has_constant="add")
    if df["s"].sum() == 0 or df["y"].nunique() < 2:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fit = sm.GLM(df["y"], x, family=sm.families.Binomial()).fit(
                cov_type="cluster", cov_kwds={"groups": df["gene"]})
        except Exception:
            return None
    b, se = fit.params[coef], fit.bse[coef]
    if not np.isfinite(b) or not np.isfinite(se) or se <= 0 or abs(b) > 15:
        return None
    return 0.5 * math.erfc((b / se) / math.sqrt(2))


def scenario(rng, n_h1, p_support, base, or_aligned):
    h1_hits = h1_fail = h2_hits = h2_fail = 0
    n_h2 = int(round(n_h1 / H1_ELIGIBLE_FRACTION))
    for _ in range(REPS):
        d1 = simulate(rng, n_h1, p_support, base, or_aligned, 1.0)
        d1["sa"] = d1["s"] * d1["a"]
        p1 = one_sided_p(d1, ["s", "a", "sa"], "sa")
        if p1 is None:
            h1_fail += 1
        elif p1 < ALPHA_ONE_SIDED:
            h1_hits += 1
        # H2 population: aligned + blocking + other, other carrying the average effect
        d2 = simulate(rng, n_h2, p_support, base, np.sqrt(or_aligned), np.sqrt(or_aligned))
        p2 = one_sided_p(d2, ["s"], "s")
        if p2 is None:
            h2_fail += 1
        elif p2 < ALPHA_ONE_SIDED:
            h2_hits += 1
    return {
        "n_h1": n_h1, "n_h2": n_h2, "p_support": p_support, "base_advance": base,
        "or_aligned": or_aligned, "or_blocking": 1.0,
        "expected_supportive_aligned": n_h1 * SHARE_ALIGNED * p_support,
        "power_h1_interaction": h1_hits / REPS, "h1_fits_failed": h1_fail,
        "h2_true_or": float(np.sqrt(or_aligned)),
        "power_h2_main": h2_hits / REPS, "h2_fits_failed": h2_fail,
    }


def main():
    started = datetime.now(timezone.utc).isoformat()
    done = {}
    if CHECKPOINT.exists():
        for line in CHECKPOINT.read_text().splitlines():
            r = json.loads(line)
            done[(r["n_h1"], r["p_support"], r["base_advance"], r["or_aligned"])] = r
    grid = list(itertools.product(N_H1, P_SUPPORT, BASE_ADVANCE, OR_ALIGNED))
    for i, (n, ps, base, ora) in enumerate(grid):
        if (n, ps, base, ora) in done:
            continue
        rng = np.random.default_rng([SEED, i])
        r = scenario(rng, n, ps, base, ora)
        with CHECKPOINT.open("a") as fh:
            fh.write(json.dumps(r) + "\n")
        done[(n, ps, base, ora)] = r
        print(f"{datetime.now(timezone.utc).isoformat()} {i + 1}/{len(grid)} n={n} ps={ps} "
              f"base={base} OR={ora}: H1 {r['power_h1_interaction']:.2f} H2 {r['power_h2_main']:.2f}")
    feas = json.loads(FEASIBILITY.read_text())
    p3 = feas["counts"]["phase3"]["by_source"]["pqtl_union"]
    p2 = feas["counts"]["phase2_only"]["by_source"]["pqtl_union"]
    result = {
        "experiment": "08_mechanism_interaction design power (V8)",
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED, "reps_per_scenario": REPS, "alpha_one_sided": ALPHA_ONE_SIDED,
        "assumptions": {
            "pairs_per_gene": PAIRS_PER_GENE, "gene_random_intercept_sd_logit": GENE_RE_SD,
            "share_aligned_among_h1": SHARE_ALIGNED, "or_blocking": 1.0,
            "h1_eligible_fraction_of_h2": H1_ELIGIBLE_FRACTION,
            "h2_true_or": "sqrt(or_aligned) for every hypothesis",
            "analysis_model": "GLM binomial, gene-clustered SE, one-sided Wald",
            "failed_fit_counts_as": "no rejection",
        },
        "feasibility_inputs": {
            "file": str(FEASIBILITY),
            "heldout_phase3_pqtl_pairs_with_outcome_gwas": p3["pairs_with_outcome_gwas_new"],
            "heldout_phase2only_pqtl_pairs_with_outcome_gwas": p2["pairs_with_outcome_gwas_new"],
            "outcome_information_used": ("none per pair; the Phase III vs Phase II-only split "
                                         "of the universe is a marginal count only"),
        },
        "scenarios": [done[k] for k in grid],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
