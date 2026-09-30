"""Design power for experiment 08 (V8) under the round-3 dependence structure.

Pure logic; no Modal dependency. Uses no outcome, evidence state, phase or disease
association: only the opaque membership structure written by the v4 feasibility count
(gene, indication and drug-program membership, mechanism class, CNS and secretion flags).

Data-generating model per hypothesis i:

    logit P(advance_i) = b0 + u_gene + w_gene * S_i + v_indication + mean_{d in drugs(i)} z_d
                         + S_i * (log OR_blocking + A_i * log(OR ratio))

with S_i ~ Bernoulli(p_support, scaled by class), A_i the mechanism class, u, w, v, z normal.

Analysis model per simulated dataset: statsmodels BinomialBayesMixedGLM fitted by variational
Bayes, with fixed effects S, A, S:A and random effects for gene intercept, gene slope on S,
indication intercept and a multiple-membership drug term (weights 1/k over a hypothesis's k
drug programs). Decision rules match the registration: Pr(beta > 0) >= 0.95 for the focal
coefficient, with Pr(beta > log 1.5) recorded. The variational posterior is a stand-in for
the registered NUTS fit.
"""
from dataclasses import dataclass
from math import erf, log, sqrt

import numpy as np
import pandas as pd
import scipy.sparse as sp
from statsmodels.genmod.bayes_mixed_glm import BinomialBayesMixedGLM

FE_PRIOR_SD = 1.5
PROB_THRESHOLD = 0.95
LOG_SESOI = log(1.5)


@dataclass(frozen=True)
class Scenario:
    scenario_id: int
    re_sd: float            # SD of gene, indication and drug intercepts
    slope_sd: float         # SD of the gene random slope on S
    p_support: float        # supportive share among blocking hypotheses
    aligned_support_mult: float  # supportive share in aligned = p_support * this
    or_blocking: float      # supportive-evidence OR in the blocking stratum
    or_ratio: float         # aligned OR / blocking OR (H1 interaction)
    base_advance: float
    cns_ratio: float        # H4: supportive OR in CNS / non-CNS (< 1 means weaker in CNS)


@dataclass(frozen=True)
class Structure:
    gene: np.ndarray
    indication: np.ndarray
    aligned: np.ndarray
    cns: np.ndarray
    drug_matrix: sp.csr_matrix   # rows = hypotheses, cols = drug programs, weights 1/k
    n_gene: int
    n_indication: int


def load_structure(path) -> Structure:
    df = pd.read_csv(path)
    gene_codes, gene = np.unique(df["gene_id"].to_numpy(), return_inverse=True)
    ind_codes, ind = np.unique(df["indication_id"].to_numpy(), return_inverse=True)
    programs = [str(x).split(";") for x in df["drug_program_ids"].fillna("")]
    all_programs = sorted({p for ps in programs for p in ps if p})
    col = {p: j for j, p in enumerate(all_programs)}
    rows, cols, vals = [], [], []
    for i, ps in enumerate(programs):
        ps = [p for p in ps if p]
        for p in ps:
            rows.append(i)
            cols.append(col[p])
            vals.append(1.0 / len(ps))
    drug = sp.csr_matrix((vals, (rows, cols)), shape=(len(df), len(all_programs)))
    aligned = (df["class"] == "aligned").to_numpy().astype(float) if "class" in df else np.zeros(len(df))
    cns = df["cns"].to_numpy().astype(float)
    return Structure(gene, ind, aligned, cns, drug, len(gene_codes), len(ind_codes))


def simulate(rng, st: Structure, sc: Scenario, mode: str) -> tuple[np.ndarray, np.ndarray]:
    n = len(st.gene)
    p_s = np.where(st.aligned == 1, sc.p_support * sc.aligned_support_mult, sc.p_support)
    s = (rng.random(n) < np.clip(p_s, 0, 1)).astype(float)
    u = rng.normal(0, sc.re_sd, st.n_gene)[st.gene]
    w = rng.normal(0, sc.slope_sd, st.n_gene)[st.gene]
    v = rng.normal(0, sc.re_sd, st.n_indication)[st.indication]
    z = st.drug_matrix @ rng.normal(0, sc.re_sd, st.drug_matrix.shape[1])
    b0 = log(sc.base_advance / (1 - sc.base_advance))
    effect = np.log(sc.or_blocking) + st.aligned * np.log(sc.or_ratio)
    if mode == "h4":
        effect = np.log(sc.or_blocking) + st.cns * np.log(sc.cns_ratio)
    lo = b0 + u + w * s + v + z + s * effect
    y = (rng.random(n) < 1 / (1 + np.exp(-lo))).astype(float)
    return y, s


def _vc_design(st: Structure, s: np.ndarray) -> tuple[sp.csr_matrix, np.ndarray]:
    n = len(st.gene)
    rows = np.arange(n)
    g_int = sp.csr_matrix((np.ones(n), (rows, st.gene)), shape=(n, st.n_gene))
    g_slope = sp.csr_matrix((s, (rows, st.gene)), shape=(n, st.n_gene))
    ind = sp.csr_matrix((np.ones(n), (rows, st.indication)), shape=(n, st.n_indication))
    vc = sp.hstack([g_int, g_slope, ind, st.drug_matrix]).tocsr()
    ident = np.concatenate([
        np.zeros(st.n_gene, int), np.ones(st.n_gene, int),
        np.full(st.n_indication, 2), np.full(st.drug_matrix.shape[1], 3),
    ])
    return vc, ident


def _norm_cdf(x: float) -> float:
    return 0.5 * (1 + erf(x / sqrt(2)))


def fit_focal(y, s, st: Structure, mode: str) -> dict | None:
    if s.sum() == 0 or y.min() == y.max():
        return None
    if mode == "h1":
        exog = np.column_stack([np.ones_like(s), s, st.aligned, s * st.aligned])
    elif mode == "h4":
        exog = np.column_stack([np.ones_like(s), s, st.cns, s * st.cns])
    else:
        exog = np.column_stack([np.ones_like(s), s])
    vc, ident = _vc_design(st, s)
    model = BinomialBayesMixedGLM(y, exog, vc, ident, fe_p=FE_PRIOR_SD)
    try:
        res = model.fit_vb()
    except (np.linalg.LinAlgError, ValueError, FloatingPointError):
        return None
    mean, sd = float(res.fe_mean[-1]), float(res.fe_sd[-1])
    if not np.isfinite(mean) or not np.isfinite(sd) or sd <= 0:
        return None
    sign = -1.0 if mode == "h4" else 1.0   # H4 predicts a weaker association in CNS
    return {
        "mean": mean, "sd": sd,
        "pr_gt0": _norm_cdf(sign * mean / sd),
        "pr_gt_sesoi": _norm_cdf((sign * mean - LOG_SESOI) / sd),
    }


def run_rep(rng, st: Structure, sc: Scenario, mode: str) -> dict:
    y, s = simulate(rng, st, sc, mode)
    r = fit_focal(y, s, st, mode)
    base = {"n_supportive": int(s.sum()),
            "n_supportive_aligned": int((s * st.aligned).sum()),
            "n_supportive_genes_aligned": int(len(np.unique(st.gene[(s == 1) & (st.aligned == 1)]))),
            "n_supportive_genes_blocking": int(len(np.unique(st.gene[(s == 1) & (st.aligned == 0)])))}
    if r is None:
        return {**base, "failed": True, "reject": False, "reject_sesoi": False}
    return {**base, "failed": False, **r,
            "reject": r["pr_gt0"] >= PROB_THRESHOLD,
            "reject_sesoi": r["pr_gt_sesoi"] >= PROB_THRESHOLD}


def scenario_grid() -> list[Scenario]:
    grid, k = [], 0
    for re_sd in (0.3, 0.7, 1.0):
        for slope_sd in (0.3, 0.7):
            for p_support in (0.03, 0.06, 0.10):
                for mult in (1.0, 1.5):
                    for or_blocking in (0.8, 1.0, 1.25):
                        for or_ratio in (1.5, 2.0, 3.0):
                            grid.append(Scenario(k, re_sd, slope_sd, p_support, mult,
                                                 or_blocking, or_ratio, 0.30, 1.0 / or_ratio))
                            k += 1
    return grid


def summarize(records: list[dict]) -> dict:
    n = len(records)
    rej = sum(r["reject"] for r in records)
    p = rej / n if n else float("nan")
    mc = 1.96 * sqrt(p * (1 - p) / n) if n else float("nan")
    return {"reps": n, "power": p, "mc_95": [max(0.0, p - mc), min(1.0, p + mc)],
            "power_sesoi": sum(r["reject_sesoi"] for r in records) / n if n else float("nan"),
            "failed_fits": sum(r["failed"] for r in records),
            "median_supportive_genes_aligned": float(np.median([r["n_supportive_genes_aligned"] for r in records])) if n else float("nan")}
