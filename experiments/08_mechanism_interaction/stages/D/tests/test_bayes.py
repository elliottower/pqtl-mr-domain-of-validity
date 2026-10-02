"""PyMC tests; they run in the Modal image (modal_stage_d.py::tests) and skip where PyMC is absent."""
import numpy as np
import pytest

pm = pytest.importorskip("pymc")

from helpers import small_frame  # noqa: E402
from power_v9_calibrate import fit_exact  # noqa: E402
from power_v9_logic import Structure  # noqa: E402

from stage_d.bayes import build_model, fit  # noqa: E402
from stage_d.constants import SamplerSettings  # noqa: E402
from stage_d.designs import build_design  # noqa: E402
from stage_d.fitspec import FitSpec  # noqa: E402
from stage_d.fitting import run_fit  # noqa: E402
from stage_d.posterior import PRIORS  # noqa: E402

SMALL = SamplerSettings(chains=2, draws=500, tune=500, target_accept=0.95)


def _frame(seed, n=400, constant_covariates=False):
    rng = np.random.default_rng(seed)
    f = small_frame(rng, n=n, n_genes=40, n_ind=20)
    if constant_covariates:
        f["platform_somascan"] = 0.0
        f["z_log10_neff"] = 0.0
        f["oncology_f"] = 0.0
    lo = -0.8 + f["S"] * (0.1 + 1.0 * (f["cls"] == "aligned"))
    f["y"] = (rng.random(n) < 1 / (1 + np.exp(-lo))).astype(float)
    return f


def test_model_variables_by_variant():
    d = build_design("h1__T", "h1", _frame(0))
    names = lambda m: {v.name for v in m.free_RVs}  # noqa: E731
    full, re = build_model(d, PRIORS["normal15"], keep_slope=True)
    assert names(full) == {"b0", "b", "sd", "z_gene", "z_slope", "z_indication", "z_program"}
    assert re == ["gene", "slope", "indication", "program"]
    reduced, re = build_model(d, PRIORS["normal15"], keep_slope=False)
    assert "z_slope" not in names(reduced) and re == ["gene", "indication", "program"]
    f = _frame(0)
    odd = f.index[f["gene_ensembl"].isin(["g1", "g3"])]
    f.loc[odd[::2], "cls"] = "aligned"   # two genes now carry both classes
    within = build_design("within__T", "within", f)
    assert within.n_gene == 2 and within.gene_fixed
    m, re = build_model(within, PRIORS["student_t3"], keep_slope=True)
    assert "z_gene" not in names(m) and re == ["slope", "indication", "program"]


def test_matches_reference_fit_exact_without_covariates():
    d = build_design("h1__R", "h1", _frame(1, constant_covariates=True))
    assert d.columns == ("S", "A", "SxA")
    st = Structure(gene=d.gene, indication=d.indication, aligned=d.X[:, 1], cns=np.zeros(d.n), drug_matrix=d.W,
                   n_gene=d.n_gene, n_indication=d.n_indication)
    ref = fit_exact(d.y, d.X[:, 0], st, draws=1000, tune=1000, chains=4, seed=3, mode="h1")
    draws, diag, _ = fit(d, PRIORS["normal15"], SamplerSettings(chains=4, draws=1000, tune=1000, target_accept=0.95),
                         True, seed=4)
    focal = draws.b[:, d.inter_col]
    assert np.median(focal) == pytest.approx(ref["median"], abs=0.15)
    assert np.mean(focal > 0) == pytest.approx(ref["pr_gt0"], abs=0.05)


def test_fit_through_attempt_sequence(tmp_path):
    d = build_design("h1__F", "h1", _frame(2))
    spec = FitSpec(fit_id="h1__F__normal15", design_id="h1__F", kind="h1", set_id="T", prior="normal15", role="H1",
                   primary=SMALL, rerun=SMALL)
    final = run_fit(spec, d, tmp_path, fit, fingerprint="f" * 64, design_sha256="d" * 64)
    att = final["attempts"][0]
    assert att["diagnostics"]["rhat_max"] > 0.99
    s = att["summary"]
    assert set(s["coefficients"]) >= {"intercept", "S", "A", "SxA", "sd_gene", "sd_slope", "sd_indication", "sd_program"}
    assert 0 <= s["focal"]["pr_predicted"] <= 1
    assert s["marginal"]["risk_difference_interaction"]["median"] > 0
    assert (tmp_path / spec.fit_id / "full.json").exists()
