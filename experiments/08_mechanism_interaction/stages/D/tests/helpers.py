"""SYNTHETIC helpers shared by tests: small frames, and a numpy stand-in for the PyMC sampler
with the signature of stage_d.bayes.fit."""
import numpy as np
import pandas as pd

from stage_d.gates import Diagnostics
from stage_d.posterior import Draws, random_effects


def small_frame(rng, n=200, n_genes=20, n_ind=10, k_max=3, programs_per_gene=4):
    genes = rng.integers(0, n_genes, n)
    rows = []
    for i, g in enumerate(genes):
        k = int(rng.integers(1, k_max + 1))
        progs = sorted(set(f"p{g}_{j}" for j in rng.integers(0, programs_per_gene, k)))
        rows.append({"hypothesis_id": f"h{i}", "gene_ensembl": f"g{g}", "indication_id": f"i{rng.integers(0, n_ind)}",
                     "programs": progs, "cls": "aligned" if g % 2 == 0 else "blocking"})
    f = pd.DataFrame(rows)
    f["S"] = (rng.random(n) < 0.3).astype(float)
    f["state"] = np.where(f["S"] == 1, "supportive", np.where(rng.random(n) < 0.2, "contradictory", "inconclusive"))
    f["E"] = np.where(f["S"] == 1, rng.uniform(0.8, 1, n), rng.uniform(-1, 0.8, n))
    f["C"] = (rng.random(n) < 0.3).astype(float)
    f["platform_somascan"] = (rng.random(n) < 0.4).astype(float)
    f["z_log10_neff"] = rng.normal(0, 1, n)
    f["oncology_f"] = (rng.random(n) < 0.2).astype(float)
    f["y"] = (rng.random(n) < 0.35).astype(float)
    return f


def fake_draws(design, rng, n_draws=400, focal_mean=0.0, focal_sd=0.3, sd_values=None, keep_slope=True, n_z=100):
    re = random_effects(design, keep_slope)
    b = rng.normal(0.0, 0.3, (n_draws, design.X.shape[1]))
    nz = np.flatnonzero(design.focal)
    b[:, nz[0]] = rng.normal(focal_mean, focal_sd, n_draws)
    if len(nz) > 1:
        b[:, nz[1]] = 0.0
    sizes = {"gene": design.n_gene, "slope": design.n_gene, "indication": design.n_indication,
             "program": design.W.shape[1]}
    idx = np.arange(min(n_z, n_draws))
    sd = {name: np.full(n_draws, (sd_values or {}).get(name, 0.5)) for name in re}
    return Draws(b0=rng.normal(-0.5, 0.1, n_draws), b=b, sd=sd,
                 z={name: rng.standard_normal((len(idx), sizes[name])) for name in re}, z_index=idx), re


def stub_fit_factory(fail_full=False, focal_mean=0.0, focal_sd=0.3, calls=None):
    """Stand-in for stage_d.bayes.fit: numpy draws, diagnostics that pass unless `fail_full`
    and the gene slope is kept (so the sequence reaches the reduced model)."""
    def fit_fn(design, prior, sampler, keep_slope, seed):
        if calls is not None:
            calls.append((design.design_id, sampler.target_accept, keep_slope))
        rng = np.random.default_rng(seed)
        draws, re = fake_draws(design, rng, focal_mean=focal_mean, focal_sd=focal_sd, keep_slope=keep_slope)
        bad = fail_full and keep_slope
        return draws, Diagnostics(rhat_max=1.05 if bad else 1.001, ess_bulk_min=900.0, divergences=0), re
    return fit_fn
