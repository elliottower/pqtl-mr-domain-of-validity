# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "numpy==2.1.3",
#     "statsmodels==0.14.4",
#     "pandas==2.2.3",
#     "scipy==1.14.1",
# ]
# ///
"""Design power for experiment 08 (V8) on the realized strict H1 structure (S1). Uses no outcome, MR estimate
or evidence state. Copy of power_v8.py with: gene-level mechanism class, the observed
pairs-per-gene distribution of the held-out set, and a grid sized for the expanded universe.

Inputs are counts from the outcome-blind feasibility file and assumed effect sizes. The
data-generating model has a gene random intercept; the analysis model is a logistic
regression with gene-clustered standard errors, a stand-in for the registered Bayesian
hierarchical model that is cheap enough to simulate.

    H1 (interaction): advance ~ S + A + S:A, one-sided test that the S:A coefficient > 0
    H2 (main effect): advance ~ S, one-sided test that the S coefficient > 0

S = supportive evidence state (1/0). A = abundance-aligned (1) vs function-blocking (0).

Usage:
    uv run experiments/08_mechanism_interaction/power/power_v8b.py
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
FEASIBILITY = HERE.parent / "feasibility" / "v2_all_indications" / "coverage_all_indications.json"
HELDOUT_PAIRS = HERE.parent / "feasibility" / "v3_strict" / "h1_gene_class_sizes.csv"
OUT_JSON = HERE / "power_v8c.json"
CHECKPOINT = HERE / "power_v8c_checkpoint.jsonl"

SEED = 202609302
REPS = 1000
ALPHA_ONE_SIDED = 0.05
GENE_RE_SD = 0.7

# Grid. n_h1 is the number of hypotheses entering the H1 model (aligned + blocking).
N_H1 = [0]                            # unused: H1 size is the realized S1 structure
P_SUPPORT = [0.03, 0.06, 0.10]
BASE_ADVANCE = [0.25, 0.40]
OR_ALIGNED = [1.5, 2.0, 3.0]          # supportive-vs-not OR in the aligned stratum
H1_ELIGIBLE_FRACTION = 0.70           # share of H2 hypotheses that are aligned or blocking

# Observed pairs per gene in the held-out set (gene column only; no phase, outcome or evidence).
SIZES = pd.read_csv(HELDOUT_PAIRS)
ALIGNED_SIZES = SIZES.loc[SIZES["class"] == "aligned", "n_hypotheses"].to_numpy()
BLOCKING_SIZES = SIZES.loc[SIZES["class"] == "blocking", "n_hypotheses"].to_numpy()


def simulate(rng, n, p_support, base, or_aligned, or_blocking):
    """One dataset on the realized S1 H1 structure; n is ignored."""
    sizes = np.concatenate([ALIGNED_SIZES, BLOCKING_SIZES])
    a_gene = np.concatenate([np.ones(len(ALIGNED_SIZES)), np.zeros(len(BLOCKING_SIZES))])
    gene = np.repeat(np.arange(len(sizes)), sizes)
    a = a_gene[gene]
    u = rng.normal(0.0, GENE_RE_SD, size=len(sizes))[gene]
    s = (rng.random(len(gene)) < p_support).astype(float)
    b0 = np.log(base / (1 - base))
    lo = b0 + u + s * np.where(a == 1, np.log(or_aligned), np.log(or_blocking))
    y = (rng.random(len(gene)) < 1 / (1 + np.exp(-lo))).astype(float)
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
        "expected_supportive_aligned": float(ALIGNED_SIZES.sum() * p_support),
        "h1_hypotheses": int(ALIGNED_SIZES.sum() + BLOCKING_SIZES.sum()),
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
    result = {
        "experiment": "08_mechanism_interaction design power (V8)",
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED, "reps_per_scenario": REPS, "alpha_one_sided": ALPHA_ONE_SIDED,
        "assumptions": {
            "pairs_per_gene": "realized S1 H1 gene sizes by class (feasibility/v3_strict/h1_gene_class_sizes.csv)",
            "mechanism_class_level": "gene", "gene_random_intercept_sd_logit": GENE_RE_SD,
            "aligned_genes": int(len(ALIGNED_SIZES)), "aligned_hypotheses": int(ALIGNED_SIZES.sum()),
            "blocking_genes": int(len(BLOCKING_SIZES)), "blocking_hypotheses": int(BLOCKING_SIZES.sum()),
            "or_blocking": 1.0, "h2_note": "H2 simulated on the same H1 structure (smaller than the full S1; conservative)",
            "h1_eligible_fraction_of_h2": H1_ELIGIBLE_FRACTION,
            "h2_true_or": "sqrt(or_aligned) for every hypothesis",
            "analysis_model": "GLM binomial, gene-clustered SE, one-sided Wald",
            "failed_fit_counts_as": "no rejection",
        },
        "feasibility_inputs": {
            "file": str(HELDOUT_PAIRS),
            "source": "strict S1 recount, per-gene hypothesis counts by class; no gene names, phase, outcome or evidence",
            "outcome_information_used": "none; gene column only",
        },
        "scenarios": [done[k] for k in grid],
    }
    OUT_JSON.write_text(json.dumps(result, indent=2))
    print(f"Wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
