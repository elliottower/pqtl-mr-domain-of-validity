from dataclasses import replace
from math import log

import numpy as np
import pandas as pd
import pytest
from helpers import small_frame
from power_v9_fast import component_clusters, fit_two_way
from power_v9_logic import Scenario, load_structure, simulate

from stage_d.designs import Design, build_design, component_labels
from stage_d.frequentist import (component_bootstrap, focal_z, null_cells, run_null_cell, simulate_null,
                                 size_correction, two_way_covariance)


@pytest.fixture
def membership_csv(tmp_path):
    """SYNTHETIC membership file in the power_v9 format."""
    rng = np.random.default_rng(0)
    n = 400
    genes = rng.integers(0, 40, n)
    progs = [";".join(sorted({f"d{g}_{j}" for j in rng.integers(0, 3, rng.integers(1, 3))} |
                             ({f"shared{rng.integers(0, 4)}"} if rng.random() < 0.1 else set()))) for g in genes]
    df = pd.DataFrame({"hypothesis_id": range(n), "gene_id": genes, "indication_id": rng.integers(0, 30, n),
                       "class": np.where(genes % 3 == 0, "aligned", "blocking"), "cns": (rng.random(n) < 0.2).astype(int),
                       "secreted": 1, "drug_program_ids": progs})
    path = tmp_path / "membership.csv"
    df.to_csv(path, index=False)
    return path


def _design_from_structure(st, s, kind="h1"):
    X = np.column_stack([s, st.aligned, s * st.aligned])
    return Design(design_id="sim", kind=kind, y=np.zeros(len(s)), X=X, columns=("S", "A", "SxA"),
                  focal=np.array([0.0, 0.0, 1.0]), focal_sign=1.0, ev_col=0, stratum_col=1, inter_col=2, slope=s,
                  gene=st.gene, n_gene=st.n_gene, indication=st.indication, n_indication=st.n_indication,
                  W=st.drug_matrix, k=np.ones(len(s)), component=np.zeros(len(s), int), group=np.array(["g"] * len(s)),
                  hypothesis_ids=np.arange(len(s)).astype(str), gene_fixed=False)


def test_component_labels_equal_power_v9_clusters(membership_csv):
    df = pd.read_csv(membership_csv)
    ours = component_labels(list(df["gene_id"]), [str(p).split(";") for p in df["drug_program_ids"]])
    np.testing.assert_array_equal(ours, component_clusters(membership_csv))


def test_simulate_null_reproduces_power_v9_data_generating_model(membership_csv):
    st = load_structure(membership_csv)
    sc = Scenario(scenario_id=0, re_sd=0.7, slope_sd=0.3, p_support=0.1, aligned_support_mult=1.5,
                  or_blocking=1.25, or_ratio=1.0, base_advance=0.3, cns_ratio=1.0)
    for seed in range(5):
        y_ref, s_ref = simulate(np.random.default_rng(seed), st, sc, "h1")
        rng = np.random.default_rng(seed)
        p_s = np.where(st.aligned == 1, sc.p_support * sc.aligned_support_mult, sc.p_support)
        s = (rng.random(len(st.gene)) < p_s).astype(float)
        y = simulate_null(rng, _design_from_structure(st, s), sc.re_sd, sc.slope_sd, log(sc.or_blocking), 0.3)
        np.testing.assert_array_equal(s, s_ref)
        np.testing.assert_array_equal(y, y_ref)


def test_two_way_covariance_diagonal_equals_power_v9_fit_two_way():
    f = small_frame(np.random.default_rng(1), n=300)
    d = build_design("h1__F", "h1", f)
    x = np.column_stack([np.ones(d.n), d.X])
    beta, se = fit_two_way(x, d.y, d.component, d.indication)
    beta2, v = two_way_covariance(x, d.y, d.component, d.indication)
    np.testing.assert_allclose(beta2, beta, rtol=1e-12)
    np.testing.assert_allclose(np.sqrt(np.diag(v)), se, rtol=1e-10)
    z = focal_z(d, d.y)
    assert z["estimate"] == pytest.approx(beta[1 + d.inter_col])
    assert z["z"] == pytest.approx(beta[1 + d.inter_col] / se[1 + d.inter_col])


def test_h3_and_h4_focal_statistics():
    f = small_frame(np.random.default_rng(2), n=300)
    d3 = build_design("h3__F", "h3", f)
    x = np.column_stack([np.ones(d3.n), d3.X])
    beta, v = two_way_covariance(x, d3.y, d3.component, d3.indication)
    z = focal_z(d3, d3.y)
    assert z["estimate"] == pytest.approx(beta[1] - beta[2])
    assert z["se"] == pytest.approx(np.sqrt(v[1, 1] + v[2, 2] - 2 * v[1, 2]))
    d4 = build_design("h4__F", "h4", f)
    z4 = focal_z(d4, d4.y)
    assert z4["z"] == pytest.approx(-z4["estimate"] / z4["se"])


def test_null_grid_matches_power_v9():
    assert len(null_cells("h1")) == 18 and len(null_cells("h2")) == 6
    assert {c[2] for c in null_cells("h2")} == {1.0}


def test_null_cell_checkpoint_resume_is_exact(tmp_path):
    d = build_design("h1__N", "h1", small_frame(np.random.default_rng(3), n=250))
    full = run_null_cell(d, 4, (0.7, 0.7, 1.25), tmp_path / "a", reps=25)
    run_null_cell(d, 4, (0.7, 0.7, 1.25), tmp_path / "b", reps=10)
    resumed = run_null_cell(d, 4, (0.7, 0.7, 1.25), tmp_path / "b", reps=25)
    assert resumed == full and len(full) == 25


def test_bootstrap_checkpoint_resume_is_exact(tmp_path):
    d = build_design("h1__B", "h1", small_frame(np.random.default_rng(4), n=250))
    a = component_bootstrap(d, tmp_path / "a.jsonl", n_boot=40)
    component_bootstrap(d, tmp_path / "b.jsonl", n_boot=15)
    b = component_bootstrap(d, tmp_path / "b.jsonl", n_boot=40)
    assert a == b
    assert a["q05"] <= a["median"] <= a["q95"]


def test_size_corrected_critical_value_controls_null_rejections(tmp_path):
    d = build_design("h1__C", "h1", small_frame(np.random.default_rng(5), n=300))
    sc = size_correction(replace(d, kind="h2"), tmp_path, reps=200)
    for cell in sc["cells"]:
        zs = cell["_z"]
        assert np.mean(zs >= cell["crit_z"]) == pytest.approx(0.05, abs=0.011)
    assert sc["crit_z"] == max(c["crit_z"] for c in sc["cells"])
    assert sc["crit_z_reference_cell"] == next(c["crit_z"] for c in sc["cells"] if (c["re_sd"], c["slope_sd"]) == (0.7, 0.7))
