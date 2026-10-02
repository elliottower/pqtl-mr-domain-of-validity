import numpy as np
import pytest
from statsmodels.duration.hazard_regression import PHReg

from stage_d.cox import cox_two_way, fit_cox


@pytest.fixture
def survival():
    rng = np.random.default_rng(0)
    n = 600
    X = np.column_stack([rng.random(n) < 0.3, rng.normal(0, 1, n), rng.random(n) < 0.5]).astype(float)
    hazard = np.exp(X @ np.array([0.7, -0.4, 0.2]))
    t_event = rng.exponential(1 / hazard)
    t_cens = rng.exponential(1.2, n)
    time = np.round(np.minimum(t_event, t_cens), 1) + 0.1   # rounding creates ties
    event = (t_event <= t_cens).astype(float)
    groups = rng.integers(0, 40, n)
    return time, event, X, groups


def test_cox_breslow_coefficients_match_statsmodels(survival):
    time, event, X, _ = survival
    assert len(np.unique(time)) < len(time)
    beta, info = fit_cox(time, event, X)
    ref = PHReg(time, X, status=event, ties="breslow").fit()
    np.testing.assert_allclose(beta, ref.params, rtol=1e-6)
    np.testing.assert_allclose(np.linalg.inv(info), ref.cov_params(), rtol=1e-5)


def test_two_way_reduces_to_statsmodels_one_way_robust(survival):
    time, event, X, groups = survival
    unique = np.arange(len(time))
    ours = cox_two_way(time, event, X, groups, unique, ["a", "b", "c"])
    ref = PHReg(time, X, status=event, ties="breslow").fit(groups=groups)
    np.testing.assert_allclose([ours[k]["se"] for k in "abc"], ref.bse, rtol=1e-5)
