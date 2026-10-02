"""Cox proportional hazards for time to Phase III with two-way cluster-robust standard errors
(PREREG §Other planned analysis, Further secondary analyses): Breslow ties, Newton-Raphson on the
partial likelihood, Lin-Wei score residuals, and the Cameron-Gelbach-Miller combination
V = I^-1 (M_components + M_indications - M_intersection) I^-1 with negative eigenvalues clipped,
as in power_v9_fast.fit_two_way.
"""
import numpy as np
from power_v9_fast import _meat


class CoxFailure(ArithmeticError):
    pass


def _risk_sums(time, xb, X):
    """For each row, S0, S1, S2 over its risk set {j: t_j >= t_i} (Breslow)."""
    order = np.argsort(-time, kind="mergesort")
    t = time[order]
    r = np.exp(xb[order])
    s0 = np.cumsum(r)
    s1 = np.cumsum(r[:, None] * X[order], axis=0)
    s2 = np.cumsum(r[:, None, None] * X[order][:, :, None] * X[order][:, None, :], axis=0)
    # rows tied in time share the risk set ending at the last tied position in descending order
    last = np.searchsorted(-t, -t, side="right") - 1
    inv = np.empty_like(order)
    inv[order] = np.arange(len(order))
    pos = last[inv]
    return s0[pos], s1[pos], s2[pos]


def fit_cox(time, event, X, iters: int = 50, tol: float = 1e-9):
    beta = np.zeros(X.shape[1])
    ev = event.astype(bool)
    for _ in range(iters):
        s0, s1, s2 = _risk_sums(time, X @ beta, X)
        xbar = s1 / s0[:, None]
        grad = (X[ev] - xbar[ev]).sum(axis=0)
        info = (s2[ev] / s0[ev][:, None, None] - xbar[ev][:, :, None] * xbar[ev][:, None, :]).sum(axis=0)
        try:
            step = np.linalg.solve(info, grad)
        except np.linalg.LinAlgError as err:
            raise CoxFailure(str(err)) from err
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    s0, s1, s2 = _risk_sums(time, X @ beta, X)
    xbar = s1 / s0[:, None]
    info = (s2[ev] / s0[ev][:, None, None] - xbar[ev][:, :, None] * xbar[ev][:, None, :]).sum(axis=0)
    if not np.all(np.isfinite(beta)):
        raise CoxFailure("non-finite coefficients")
    return beta, info


def score_residuals(time, event, X, beta) -> np.ndarray:
    """Lin-Wei score residuals: delta_i (x_i - xbar(t_i)) - sum_{k event, t_k <= t_i} r_i / S0(t_k) (x_i - xbar(t_k))."""
    ev = event.astype(bool)
    r = np.exp(X @ beta)
    s0, s1, _ = _risk_sums(time, X @ beta, X)
    xbar = s1 / s0[:, None]
    et, inv_s0, xbar_s0 = time[ev], 1.0 / s0[ev], xbar[ev] / s0[ev][:, None]
    order = np.argsort(et, kind="mergesort")
    et, c0, c1 = et[order], np.cumsum(inv_s0[order]), np.cumsum(xbar_s0[order], axis=0)
    upto = np.searchsorted(et, time, side="right") - 1
    a0 = np.where(upto >= 0, c0[np.clip(upto, 0, None)], 0.0)
    a1 = np.where(upto[:, None] >= 0, c1[np.clip(upto, 0, None)], 0.0)
    return ev[:, None] * (X - xbar) - r[:, None] * (X * a0[:, None] - a1)


def cox_two_way(time, event, X, g1, g2, names) -> dict:
    beta, info = fit_cox(time, event, X)
    bread = np.linalg.inv(info)
    s = score_residuals(time, event, X, beta)
    g12 = g1.astype(np.int64) * (g2.max() + 1) + g2
    v = bread @ (_meat(s, g1) + _meat(s, g2) - _meat(s, g12)) @ bread
    evals, evecs = np.linalg.eigh(v)
    v = evecs @ np.diag(np.clip(evals, 0, None)) @ evecs.T
    se = np.sqrt(np.diag(v))
    return {name: {"log_hr": float(b), "se": float(e), "hr": float(np.exp(b)),
                   "hr_ci95": [float(np.exp(b - 1.959964 * e)), float(np.exp(b + 1.959964 * e))]}
            for name, b, e in zip(names, beta, se)}
