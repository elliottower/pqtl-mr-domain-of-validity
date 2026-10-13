"""Q6: surrogate decision calibration, not power (PREREG §Statistical models, Q6), on V8's planning
simulation (`power_v9_logic.simulate` and the variational surrogate `fit_focal`, mode h1) over the
structure of the realized S1 H1 set. Runs in the V8 power image (numpy 2.1.3, statsmodels 0.14.4)."""
import json
from pathlib import Path

import numpy as np

from followup.common import SEED, wilson
from power_v9_logic import FE_PRIOR_SD, Scenario, Structure, fit_focal, simulate

K_GRID = (14, 20, 30, 40, 60, 80, 120, 160, 240)
RATIOS = (1.0, 1.5, 2.0, 3.0)
FIRST_PASS = 1000
SECOND_PASS = 1000
DECISION = 0.95
TARGET = 0.8
GATE_HYPOTHESES = 20
GATE_GENES = 10
PRIOR_DOMINANCE = 0.8
CHECKPOINT_EVERY = 50


def scenario(k: int, ratio: float, n_aligned: int, n_blocking: int) -> Scenario:
    """Expected k supportive hypotheses in each stratum; V8 reference cell otherwise."""
    return Scenario(scenario_id=-1, re_sd=0.7, slope_sd=0.7, p_support=k / n_blocking,
                    aligned_support_mult=n_blocking / n_aligned, or_blocking=1.0, or_ratio=ratio,
                    base_advance=0.30, cns_ratio=1.0)


def dataset(st: Structure, sc: Scenario, rng: np.random.Generator) -> dict:
    y, s = simulate(rng, st, sc, "h1")
    al = st.aligned == 1
    sup = s == 1
    counts = {"aligned": int((sup & al).sum()), "blocking": int((sup & ~al).sum()),
              "genes_aligned": int(len(np.unique(st.gene[sup & al]))),
              "genes_blocking": int(len(np.unique(st.gene[sup & ~al])))}
    gate = (min(counts["aligned"], counts["blocking"]) >= GATE_HYPOTHESES
            and min(counts["genes_aligned"], counts["genes_blocking"]) >= GATE_GENES)
    r = fit_focal(y, s, st, "h1")
    if r is None:
        return {**counts, "gate": gate, "failed": True, "prior_dominated": None, "decision": False}
    dominated = r["sd"] > PRIOR_DOMINANCE * FE_PRIOR_SD
    return {**counts, "gate": gate, "failed": False, "sd": r["sd"], "pr_gt0": r["pr_gt0"],
            "prior_dominated": dominated, "decision": bool(gate and not dominated and r["pr_gt0"] >= DECISION)}


def cell_rng(k_idx: int, r_idx: int, pass_no: int, rep: int) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(k_idx, r_idx, pass_no, rep)))


def run_cell(st: Structure, k_idx: int, r_idx: int, pass_no: int, n: int, path: Path, on_commit=None) -> None:
    n_al, n_bl = int((st.aligned == 1).sum()), int((st.aligned == 0).sum())
    sc = scenario(K_GRID[k_idx], RATIOS[r_idx], n_al, n_bl)
    done = len(path.read_text().splitlines()) if path.exists() else 0
    buf = []
    for rep in range(done, n):
        buf.append({"rep": rep, **dataset(st, sc, cell_rng(k_idx, r_idx, pass_no, rep))})
        if len(buf) == CHECKPOINT_EVERY or rep == n - 1:
            path.parent.mkdir(parents=True, exist_ok=True)
            with path.open("a") as fh:
                fh.write("".join(json.dumps(r) + "\n" for r in buf))
            buf = []
            if on_commit:
                on_commit()


def summarize_cell(records: list[dict]) -> dict:
    n = len(records)

    def prob(key):
        k = sum(1 for r in records if r[key] is True)
        return {"p": k / n, "wilson95": wilson(k, n)}
    ok = [r for r in records if not r["failed"]]
    return {"datasets": n, "gate": prob("gate"), "prior_dominated": prob("prior_dominated"), "decision": prob("decision"),
            "failed_share": (n - len(ok)) / n,
            "median_interval_width_90": float(np.median([3.29 * r["sd"] for r in ok])) if ok else None,
            "median_supportive": {s: float(np.median([r[s] for r in records])) for s in ("aligned", "blocking")}}


def second_pass_cells(first: dict[tuple[int, int], float]) -> list[tuple[int, int]]:
    """first: (k_idx, r_idx) -> decision probability. For each ratio, k* = smallest grid k reaching 0.8;
    k* and its two lower and two upper grid neighbours (truncated); if none, the two largest k."""
    cells = []
    for r_idx in range(len(RATIOS)):
        hit = [k for k in range(len(K_GRID)) if first[(k, r_idx)] >= TARGET]
        ks = range(max(0, hit[0] - 2), min(len(K_GRID), hit[0] + 3)) if hit else (len(K_GRID) - 2, len(K_GRID) - 1)
        cells += [(k, r_idx) for k in ks]
    return cells
