import json

import numpy as np
import pandas as pd
import pytest

from followup.q6 import K_GRID, RATIOS, run_cell, scenario, second_pass_cells, summarize_cell
from power_v9_logic import load_structure


def test_second_pass_takes_k_star_and_two_neighbours_each_side_truncated_at_the_grid_ends():
    first = {(k, r): 0.0 for k in range(len(K_GRID)) for r in range(len(RATIOS))}
    for k in range(4, len(K_GRID)):
        first[(k, 0)] = 0.9                     # k* = index 4
    first[(0, 1)] = 0.85                        # k* = index 0 -> 0, 1, 2
    cells = second_pass_cells(first)
    assert [k for k, r in cells if r == 0] == [2, 3, 4, 5, 6]
    assert [k for k, r in cells if r == 1] == [0, 1, 2]
    assert [k for k, r in cells if r == 2] == [len(K_GRID) - 2, len(K_GRID) - 1]


def test_scenario_gives_k_expected_supportive_in_each_stratum():
    sc = scenario(40, 2.0, n_aligned=500, n_blocking=2500)
    assert sc.p_support * 2500 == pytest.approx(40) and sc.p_support * sc.aligned_support_mult * 500 == pytest.approx(40)


def tiny_structure(tmp_path):
    rng = np.random.default_rng(5)
    n = 400
    pd.DataFrame({"gene_id": [f"g{i % 60}" for i in range(n)], "indication_id": [f"d{i % 25}" for i in range(n)],
                  "drug_program_ids": [f"p{i % 40}" for i in range(n)],
                  "class": np.where(rng.random(n) < 0.3, "aligned", "blocking"), "cns": 0}).to_csv(tmp_path / "s.csv", index=False)
    return load_structure(tmp_path / "s.csv")


def test_cells_resume_to_identical_datasets_and_summaries_count_failures(tmp_path):
    st = tiny_structure(tmp_path)
    run_cell(st, 3, 2, 1, 4, tmp_path / "a.jsonl")
    run_cell(st, 3, 2, 1, 2, tmp_path / "b.jsonl")
    run_cell(st, 3, 2, 1, 4, tmp_path / "b.jsonl")
    a = [json.loads(x) for x in (tmp_path / "a.jsonl").read_text().splitlines()]
    b = [json.loads(x) for x in (tmp_path / "b.jsonl").read_text().splitlines()]
    simulated = ("rep", "aligned", "blocking", "genes_aligned", "genes_blocking", "gate", "failed")
    assert [{k: r[k] for k in simulated} for r in a] == [{k: r[k] for k in simulated} for r in b]
    for ra, rb in zip(a, b):           # the variational fit is reproducible to floating-point noise only
        if not ra["failed"]:
            assert (ra["sd"], ra["pr_gt0"]) == pytest.approx((rb["sd"], rb["pr_gt0"]), abs=1e-5)
    s = summarize_cell(a)
    assert s["datasets"] == 4 and 0 <= s["decision"]["p"] <= 1
    assert all(r["decision"] is False or (r["gate"] and not r["prior_dominated"]) for r in a)
