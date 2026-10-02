from math import log

import numpy as np
import pytest
from helpers import fake_draws, small_frame
from scipy.stats import norm

from stage_d.designs import build_design
from stage_d.posterior import (PRIORS, Draws, posterior_predictive, prior_predictive, standardized_probabilities,
                               summarize, summarize_fit)

LOG15 = log(1.5)


def _sigmoid(x):
    return 1 / (1 + np.exp(-x))


def _gauss_hermite_mean_sigmoid(mu, sigma, n=80):
    x, w = np.polynomial.hermite_e.hermegauss(n)
    return float(np.sum(w * _sigmoid(mu + sigma * x)) / np.sqrt(2 * np.pi))


@pytest.fixture
def design():
    return build_design("h1__T", "h1", small_frame(np.random.default_rng(0)))


def _const_draws(design, b0, b, sds, n):
    return Draws(b0=np.full(n, b0), b=np.tile(np.asarray(b, float), (n, 1)),
                 sd={k: np.full(n, v) for k, v in sds.items()})


def test_summarize_registered_probabilities_match_normal_theory():
    x = np.random.default_rng(1).normal(0.5, 0.3, 400_000)
    s = summarize(x)
    assert s["pr_predicted"] == pytest.approx(norm.cdf(0.5 / 0.3), abs=3e-3)
    assert s["pr_beyond_sesoi"] == pytest.approx(norm.cdf((0.5 - LOG15) / 0.3), abs=3e-3)
    assert s["mass_within_sesoi"] == pytest.approx(norm.cdf((LOG15 - 0.5) / 0.3) - norm.cdf((-LOG15 - 0.5) / 0.3), abs=3e-3)
    flipped = summarize(x, sign=-1.0)
    assert flipped["pr_predicted"] == pytest.approx(s["pr_opposite"])
    assert flipped["q05"] == s["q05"]


def test_standardized_probabilities_without_random_effects_are_exact(design):
    b = np.linspace(-0.4, 0.6, design.X.shape[1])
    d = _const_draws(design, -0.7, b, {"gene": 0.0, "slope": 0.0, "indication": 0.0, "program": 0.0}, 3)
    m = standardized_probabilities(d, design, np.random.default_rng(2))
    rest = [j for j in range(design.X.shape[1]) if j not in (design.ev_col, design.stratum_col, design.inter_col)]
    lin = -0.7 + design.X[:, rest] @ b[rest]
    for a, label in ((1, "aligned"), (0, "blocking")):
        for s, key in ((1, "supportive"), (0, "not_supportive")):
            expect = _sigmoid(lin + s * b[design.ev_col] + a * b[design.stratum_col] + s * a * b[design.inter_col]).mean()
            assert m["probabilities"][label][key]["median"] == pytest.approx(expect, rel=1e-12)
    rd = {lab: m["probabilities"][lab]["supportive"]["median"] - m["probabilities"][lab]["not_supportive"]["median"]
          for lab in ("aligned", "blocking")}
    assert m["risk_difference_interaction"]["median"] == pytest.approx(rd["aligned"] - rd["blocking"], rel=1e-12)


@pytest.mark.parametrize("k", [1, 4])
def test_new_program_effects_average_k_draws(k):
    rng = np.random.default_rng(3)
    f = small_frame(rng, n=40)
    f["programs"] = [[f"x{i}_{j}" for j in range(k)] for i in range(len(f))]
    d = build_design("h1__K", "h1", f)
    assert np.all(d.k == k)
    b = np.zeros(d.X.shape[1])
    draws = _const_draws(d, 1.0, b, {"gene": 0.0, "slope": 0.0, "indication": 0.0, "program": 2.0}, 20_000)
    m = standardized_probabilities(draws, d, np.random.default_rng(4))
    got = m["probabilities"]["blocking"]["not_supportive"]["mean"]
    assert got == pytest.approx(_gauss_hermite_mean_sigmoid(1.0, 2.0 / np.sqrt(k)), abs=2e-3)


def test_all_random_effect_sds_enter_new_level_integration():
    rng = np.random.default_rng(5)
    d = build_design("h1__R", "h1", small_frame(rng, n=40))
    b = np.zeros(d.X.shape[1])
    sds = {"gene": 0.6, "slope": 0.8, "indication": 0.5, "program": 0.9}
    m = standardized_probabilities(_const_draws(d, 0.8, b, sds, 20_000), d, np.random.default_rng(6))
    k = d.k
    for s, key in ((0, "not_supportive"), (1, "supportive")):
        total = np.sqrt(sds["gene"] ** 2 + s * sds["slope"] ** 2 + sds["indication"] ** 2 + sds["program"] ** 2 / k)
        expect = np.mean([_gauss_hermite_mean_sigmoid(0.8, t) for t in total])
        assert m["probabilities"]["aligned"][key]["mean"] == pytest.approx(expect, abs=2e-3)


def test_h4_marginal_labels_and_sign():
    d = build_design("h4__T", "h4", small_frame(np.random.default_rng(7)))
    draws, _ = fake_draws(d, np.random.default_rng(8), focal_mean=-0.8)
    m = standardized_probabilities(draws, d, np.random.default_rng(9))
    assert set(m["probabilities"]) == {"neuro_psych", "other_indications"}
    assert m["risk_difference_interaction"]["pr_predicted"] > 0.5


def test_posterior_predictive_matches_model_rates(design):
    n = 3000
    b = np.full(design.X.shape[1], 0.2)
    d = _const_draws(design, -0.3, b, {"gene": 0.0, "slope": 0.0, "indication": 0.0, "program": 0.0}, n)
    d.z = {"gene": np.zeros((n, design.n_gene)), "slope": np.zeros((n, design.n_gene)),
           "indication": np.zeros((n, design.n_indication)), "program": np.zeros((n, design.W.shape[1]))}
    d.z_index = np.arange(n)
    ppc = posterior_predictive(d, design, np.random.default_rng(10))
    p = _sigmoid(-0.3 + design.X @ b)
    for g, v in ppc.items():
        m = design.group == g
        assert v["observed"] == pytest.approx(design.y[m].mean())
        assert v["replicated_mean"] == pytest.approx(p[m].mean(), abs=4 * np.sqrt(0.25 / (n * m.sum())) + 1e-9)


def test_prior_predictive_is_symmetric_about_one_half(design):
    pp = prior_predictive(design, PRIORS["normal15"], ["gene", "slope", "indication", "program"],
                          np.random.default_rng(11), n_sims=4000)
    for v in pp.values():
        assert v["q05"] < v["median"] < v["q95"]
        assert v["median"] == pytest.approx(0.5, abs=0.06)


def test_prior_dominance_flag_from_focal_posterior_sd(design):
    for sd, dominated in ((1.4, True), (0.4, False)):
        draws, re = fake_draws(design, np.random.default_rng(12), n_draws=4000, focal_sd=sd)
        s = summarize_fit(draws, design, PRIORS["normal15"], re, np.random.default_rng(13))
        assert s["prior_dominated"] is dominated
        assert s["focal_prior_sd"] == pytest.approx(1.5)
        assert s["focal"]["quantity"] == "SxA"


def test_h3_focal_is_contrast_with_contrast_prior_sd():
    d = build_design("h3__T", "h3", small_frame(np.random.default_rng(14)))
    draws, re = fake_draws(d, np.random.default_rng(15))
    s = summarize_fit(draws, d, PRIORS["normal15"], re, np.random.default_rng(16))
    assert s["focal"]["quantity"] == "supportive - contradictory"
    assert s["focal_prior_sd"] == pytest.approx(1.5 * np.sqrt(2))
    assert s["marginal"] is None
    assert PRIORS["student_t3"].sd == pytest.approx(1.5 * np.sqrt(3))
