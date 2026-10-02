"""The registered Bayesian hierarchical logistic model in PyMC (Linux only; PyMC 5.28.5,
PyTensor 2.38.2, as in the calibration image).

power/power_v9/power_v9_calibrate.py `fit_exact` is the reference implementation; this is the
same model with the design's fixed effects (focal terms plus platform, z(log10 N_eff) and
oncology), the prior family as an argument, and the gene slope and gene intercept optional
(the reduced model drops the slope; the within-gene diagnostic replaces the gene intercept by
fixed effects). Priors: Normal(0, 2.5) intercept; Normal(0, 1.5) fixed effects (sensitivity:
Normal(0, 1), Student-t(3, 0, 1.5)); Half-Normal(1) random-effect SDs; non-centered effects.
"""
import arviz as az
import numpy as np
import pymc as pm

from stage_d.constants import PRIOR_SD_INTERCEPT, PRIOR_SD_RANDOM, SamplerSettings
from stage_d.designs import Design
from stage_d.gates import Diagnostics
from stage_d.posterior import Draws, PriorSpec, random_effects

N_PREDICTIVE_DRAWS = 1000


def build_model(design: Design, prior: PriorSpec, keep_slope: bool) -> tuple[pm.Model, list[str]]:
    re = random_effects(design, keep_slope)
    w = design.W.toarray()
    with pm.Model() as model:
        b0 = pm.Normal("b0", 0, PRIOR_SD_INTERCEPT)
        if prior.family == "normal":
            b = pm.Normal("b", 0, prior.scale, shape=design.X.shape[1])
        else:
            b = pm.StudentT("b", nu=prior.nu, mu=0, sigma=prior.scale, shape=design.X.shape[1])
        sd = pm.HalfNormal("sd", PRIOR_SD_RANDOM, shape=len(re))
        eta = b0 + pm.math.dot(design.X, b)
        for j, name in enumerate(re):
            if name == "gene":
                eta = eta + sd[j] * pm.Normal("z_gene", 0, 1, shape=design.n_gene)[design.gene]
            elif name == "slope":
                eta = eta + sd[j] * pm.Normal("z_slope", 0, 1, shape=design.n_gene)[design.gene] * design.slope
            elif name == "indication":
                eta = eta + sd[j] * pm.Normal("z_indication", 0, 1, shape=design.n_indication)[design.indication]
            else:
                eta = eta + sd[j] * pm.math.dot(w, pm.Normal("z_program", 0, 1, shape=w.shape[1]))
        pm.Bernoulli("y", logit_p=eta, observed=design.y)
    return model, re


def diagnostics(idata) -> Diagnostics:
    rhat = az.rhat(idata)
    ess = az.ess(idata, method="bulk")
    return Diagnostics(rhat_max=float(max(np.nanmax(v.values) for v in rhat.data_vars.values())),
                       ess_bulk_min=float(min(np.nanmin(v.values) for v in ess.data_vars.values())),
                       divergences=int(idata.sample_stats["diverging"].values.sum()))


def fit(design: Design, prior: PriorSpec, sampler: SamplerSettings, keep_slope: bool, seed: int) -> tuple[Draws, Diagnostics, list[str]]:
    model, re = build_model(design, prior, keep_slope)
    with model:
        idata = pm.sample(draws=sampler.draws, tune=sampler.tune, chains=sampler.chains, cores=sampler.chains,
                          target_accept=sampler.target_accept, random_seed=seed, progressbar=False)
    post = idata.posterior
    n = sampler.chains * sampler.draws
    sd = post["sd"].values.reshape(n, len(re))
    thin = np.unique(np.linspace(0, n - 1, min(N_PREDICTIVE_DRAWS, n)).astype(int))
    z = {name: post[f"z_{name}"].values.reshape(n, -1)[thin] for name in re}
    draws = Draws(b0=post["b0"].values.reshape(n), b=post["b"].values.reshape(n, design.X.shape[1]),
                  sd={name: sd[:, j] for j, name in enumerate(re)}, z=z, z_index=thin)
    return draws, diagnostics(idata), re
