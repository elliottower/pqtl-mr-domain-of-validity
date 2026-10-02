"""Everything computed from posterior draws, in numpy (no PyMC): coefficient summaries, the
focal quantity and its registered probabilities, stratum odds ratios, marginal standardized
probabilities, and prior- and posterior-predictive checks (PREREG §Inference criteria,
Magnitude and Marginal quantities; §Statistical models, Checks).
"""
from dataclasses import dataclass, field
from math import sqrt
from typing import Literal

import numpy as np
from pydantic import BaseModel, ConfigDict

from stage_d.constants import (LOG_SESOI, PRIOR_SD_FIXED, PRIOR_SD_INTERCEPT, PRIOR_SD_RANDOM,
                               PRIOR_SENSITIVITY_NORMAL_SD, PRIOR_SENSITIVITY_T_NU, PRIOR_SENSITIVITY_T_SCALE)
from stage_d.designs import Design
from stage_d.gates import prior_dominated

RE_NAMES = ("gene", "slope", "indication", "program")


def random_effects(design: Design, keep_slope: bool) -> list[str]:
    """Random-effect SDs in the model: all four, minus the gene intercept for the within-gene
    diagnostic (gene fixed effect) and minus the gene slope in the reduced model."""
    return [n for n in RE_NAMES if not (n == "gene" and design.gene_fixed) and not (n == "slope" and not keep_slope)]


class PriorSpec(BaseModel):
    model_config = ConfigDict(frozen=True)
    name: str
    family: Literal["normal", "student_t"]
    scale: float
    nu: float = 0.0

    @property
    def sd(self) -> float:
        if self.family == "normal":
            return self.scale
        return self.scale * sqrt(self.nu / (self.nu - 2.0))


PRIORS = {
    "normal15": PriorSpec(name="normal15", family="normal", scale=PRIOR_SD_FIXED),
    "normal1": PriorSpec(name="normal1", family="normal", scale=PRIOR_SENSITIVITY_NORMAL_SD),
    "student_t3": PriorSpec(name="student_t3", family="student_t", scale=PRIOR_SENSITIVITY_T_SCALE,
                            nu=PRIOR_SENSITIVITY_T_NU),
}
PRIMARY_PRIOR = "normal15"


@dataclass
class Draws:
    """Flattened posterior draws. `sd` holds the random-effect SDs present in the fitted model
    (subset of RE_NAMES); `z` holds thinned non-centered level effects for predictive checks."""
    b0: np.ndarray
    b: np.ndarray
    sd: dict[str, np.ndarray]
    z: dict[str, np.ndarray] = field(default_factory=dict)
    z_index: np.ndarray | None = None


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))


def summarize(x: np.ndarray, sign: float = 1.0) -> dict:
    """Median, 90% interval, and the registered probabilities in the predicted direction."""
    s = sign * x
    return {"median": float(np.median(x)), "mean": float(np.mean(x)), "sd": float(np.std(x, ddof=1)),
            "q05": float(np.quantile(x, 0.05)), "q95": float(np.quantile(x, 0.95)),
            "pr_predicted": float(np.mean(s > 0)), "pr_opposite": float(np.mean(s < 0)),
            "pr_beyond_sesoi": float(np.mean(s > LOG_SESOI)),
            "mass_within_sesoi": float(np.mean(np.abs(x) < LOG_SESOI))}


def describe(x: np.ndarray) -> dict:
    """Median, mean, SD and 90% interval, for quantities with no predicted direction."""
    return {"median": float(np.median(x)), "mean": float(np.mean(x)), "sd": float(np.std(x, ddof=1)),
            "q05": float(np.quantile(x, 0.05)), "q95": float(np.quantile(x, 0.95))}


def directional(x: np.ndarray, sign: float = 1.0) -> dict:
    """describe() plus Pr(> 0) in the predicted direction; for probability-scale quantities, where
    the log(1.5) odds-ratio threshold does not apply."""
    return {**describe(x), "pr_predicted": float(np.mean(sign * x > 0)), "pr_opposite": float(np.mean(sign * x < 0))}


def coefficient_table(draws: Draws, design: Design) -> dict:
    out = {"intercept": summarize(draws.b0)}
    for j, name in enumerate(design.columns):
        out[name] = summarize(draws.b[:, j])
    for name, v in draws.sd.items():
        out[f"sd_{name}"] = summarize(v)
    return out


def focal_draws(draws: Draws, design: Design) -> np.ndarray:
    return draws.b @ design.focal


def stratum_odds_ratios(draws: Draws, design: Design) -> dict | None:
    if design.inter_col < 0 or design.kind.endswith("_E"):
        return None
    base = draws.b[:, design.ev_col]
    other = base + draws.b[:, design.inter_col]
    labels = ("neuro_psych", "other_indications") if design.kind == "h4" else ("aligned", "blocking")
    out = {}
    for label, log_or in zip(labels, (other, base)):
        o = np.exp(log_or)
        out[label] = {"or_median": float(np.median(o)), "or_q05": float(np.quantile(o, 0.05)),
                      "or_q95": float(np.quantile(o, 0.95)), "log_or": summarize(log_or)}
    return out


def _counterfactual_parts(draws: Draws, design: Design):
    roles = [design.ev_col, design.stratum_col, design.inter_col]
    rest = [j for j in range(design.X.shape[1]) if j not in roles]
    return rest, draws.b[:, design.ev_col], draws.b[:, design.stratum_col], draws.b[:, design.inter_col]


def standardized_probabilities(draws: Draws, design: Design, rng: np.random.Generator, chunk: int = 200) -> dict | None:
    """Standardized advancement probability with and without supportive evidence in each stratum,
    integrating over new random-effect levels: per posterior draw, every observed row gets a new
    gene intercept, gene slope, indication intercept and k new program effects averaged with
    weight 1/k (drawn as their exact sum, Normal(0, sd_program / sqrt(k))). The same new-level
    draws serve all four (S, stratum) cells, and predictions average over the observed covariates."""
    if design.inter_col < 0 or design.kind.endswith("_E") or design.gene_fixed:
        return None
    rest, b_ev, b_st, b_int = _counterfactual_parts(draws, design)
    n_draws, n = len(draws.b0), design.n
    p = {(s, a): np.empty(n_draws) for s in (0, 1) for a in (0, 1)}
    sqrt_k = np.sqrt(design.k)
    for lo in range(0, n_draws, chunk):
        sl = slice(lo, min(lo + chunk, n_draws))
        m = sl.stop - sl.start
        base = draws.b0[sl][None, :] + design.X[:, rest] @ draws.b[sl][:, rest].T
        noise = np.zeros((n, m))
        slope = np.zeros((n, m))
        if "gene" in draws.sd:
            noise += rng.standard_normal((n, m)) * draws.sd["gene"][sl][None, :]
        if "slope" in draws.sd:
            slope = rng.standard_normal((n, m)) * draws.sd["slope"][sl][None, :]
        noise += rng.standard_normal((n, m)) * draws.sd["indication"][sl][None, :]
        noise += rng.standard_normal((n, m)) * draws.sd["program"][sl][None, :] / sqrt_k[:, None]
        for s in (0, 1):
            for a in (0, 1):
                eta = base + noise + s * slope + (s * b_ev[sl] + a * b_st[sl] + s * a * b_int[sl])[None, :]
                p[(s, a)][sl] = _sigmoid(eta).mean(axis=0)
    labels = {1: "neuro_psych", 0: "other_indications"} if design.kind == "h4" else {1: "aligned", 0: "blocking"}
    rd = {a: p[(1, a)] - p[(0, a)] for a in (0, 1)}
    rdi = rd[1] - rd[0]
    return {
        "probabilities": {labels[a]: {"supportive": describe(p[(1, a)]), "not_supportive": describe(p[(0, a)])}
                          for a in (0, 1)},
        "risk_difference": {labels[a]: directional(rd[a]) for a in (0, 1)},
        "risk_difference_interaction": directional(rdi, sign=design.focal_sign),
        "definition": f"RD({labels[1]}) - RD({labels[0]})",
    }


def _linear_predictor_in_sample(draws: Draws, design: Design, idx: np.ndarray) -> np.ndarray:
    eta = draws.b0[idx][None, :] + design.X @ draws.b[idx].T
    zi = draws.z
    if "gene" in draws.sd:
        eta += (zi["gene"] * draws.sd["gene"][idx][:, None])[:, design.gene].T
    if "slope" in draws.sd:
        eta += (zi["slope"] * draws.sd["slope"][idx][:, None])[:, design.gene].T * design.slope[:, None]
    eta += (zi["indication"] * draws.sd["indication"][idx][:, None])[:, design.indication].T
    eta += design.W @ (zi["program"] * draws.sd["program"][idx][:, None]).T
    return eta


def _group_rates(y_rep: np.ndarray, groups: np.ndarray) -> dict[str, np.ndarray]:
    return {g: y_rep[groups == g].mean(axis=0) for g in np.unique(groups)}


def posterior_predictive(draws: Draws, design: Design, rng: np.random.Generator) -> dict:
    """Observed vs replicated advancement rate per mechanism stratum x evidence state, using the
    in-sample random effects of the thinned draws in `draws.z` (index `draws.z_index`)."""
    idx = draws.z_index
    eta = _linear_predictor_in_sample(draws, design, idx)
    y_rep = (rng.random(eta.shape) < _sigmoid(eta)).astype(float)
    out = {}
    for g, rep in _group_rates(y_rep, design.group).items():
        obs = float(design.y[design.group == g].mean())
        out[g] = {"n": int((design.group == g).sum()), "observed": obs, "replicated_mean": float(rep.mean()),
                  "replicated_q05": float(np.quantile(rep, 0.05)), "replicated_q95": float(np.quantile(rep, 0.95)),
                  "ppp_ge_observed": float(np.mean(rep >= obs))}
    return out


def sample_prior(prior: PriorSpec, n_fixed: int, re_names, n_sims: int, rng: np.random.Generator) -> Draws:
    b = (rng.normal(0.0, prior.scale, (n_sims, n_fixed)) if prior.family == "normal"
         else prior.scale * rng.standard_t(prior.nu, (n_sims, n_fixed)))
    return Draws(b0=rng.normal(0.0, PRIOR_SD_INTERCEPT, n_sims), b=b,
                 sd={name: np.abs(rng.normal(0.0, PRIOR_SD_RANDOM, n_sims)) for name in re_names})


def prior_predictive(design: Design, prior: PriorSpec, re_names, rng: np.random.Generator, n_sims: int = 1000) -> dict:
    """Distribution of stratum x evidence-state advancement rates under the priors, on the
    observed design structure."""
    d = sample_prior(prior, design.X.shape[1], re_names, n_sims, rng)
    sizes = {"gene": design.n_gene, "slope": design.n_gene, "indication": design.n_indication,
             "program": design.W.shape[1]}
    d.z = {name: rng.standard_normal((n_sims, sizes[name])) for name in RE_NAMES}
    d.z_index = np.arange(n_sims)
    for name in RE_NAMES:
        if name not in re_names:
            d.sd[name] = np.zeros(n_sims)
    eta = _linear_predictor_in_sample(d, design, d.z_index)
    y_rep = (rng.random(eta.shape) < _sigmoid(eta)).astype(float)
    return {g: {"n": int((design.group == g).sum()), "q05": float(np.quantile(r, 0.05)),
                "median": float(np.median(r)), "q95": float(np.quantile(r, 0.95))}
            for g, r in _group_rates(y_rep, design.group).items()}


def focal_label(design: Design) -> str:
    return " - ".join(n for n, c in zip(design.columns, design.focal) if c != 0) if design.kind == "h3" \
        else design.columns[int(np.flatnonzero(design.focal)[0])]


def summarize_fit(draws: Draws, design: Design, prior: PriorSpec, re_names, rng: np.random.Generator) -> dict:
    """Everything reported for one fit. The focal prior SD is the prior SD of the focal contrast
    (prior SD x ||contrast||), which is the prior SD of the coefficient itself for H1 and H4."""
    f = focal_draws(draws, design)
    post_sd = float(np.std(f, ddof=1))
    focal_prior_sd = prior.sd * float(np.sqrt(np.sum(design.focal ** 2)))
    return {
        "focal": {"quantity": focal_label(design), "predicted_sign": design.focal_sign,
                  **summarize(f, design.focal_sign)},
        "focal_posterior_sd": post_sd, "focal_prior_sd": focal_prior_sd,
        "prior_dominated": prior_dominated(post_sd, focal_prior_sd),
        "coefficients": coefficient_table(draws, design),
        "stratum_odds_ratios": stratum_odds_ratios(draws, design),
        "marginal": standardized_probabilities(draws, design, rng),
        "posterior_predictive": posterior_predictive(draws, design, rng),
        "prior_predictive": prior_predictive(design, prior, re_names, rng),
    }
