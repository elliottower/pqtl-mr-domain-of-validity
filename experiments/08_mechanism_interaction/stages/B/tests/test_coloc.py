"""The Python side of the coloc interface everywhere; the R side where Rscript and coloc exist.

The R tests compare coloc.abf against the closed form it implements (Wakefield approximate Bayes
factors, prior SD 0.15 x sdY for a quantitative trait and 0.2 for case-control, combined over
H0-H4 with p1, p2, p12), so they pass only if the JSON bridge hands R the right numbers and
reads back the right ones. Point R at a library with STAGE_B_R_LIBS if coloc is not in the
default one.

Every R process started here runs under conftest.py's profile: `options(warnPartialMatchDollar =
TRUE)`, and a partial-match warning is an error.
"""
import json
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd
import pytest

from stage_b.coloc_backend import (R_SCRIPT, AbfResult, ColocDataset, ColocTask, RscriptColoc, SusieResult, outcome_dataset,
                                   parse_response, pqtl_dataset, qtl_dataset)
from stage_b.schemas import ColocBackendError

R_LIBS = os.environ.get("STAGE_B_R_LIBS")


def _r_has_coloc() -> bool:
    if shutil.which("Rscript") is None:
        return False
    env = {**os.environ, **({"R_LIBS_USER": R_LIBS} if R_LIBS else {})}
    return subprocess.run(["Rscript", "-e", "library(coloc)"], env=env, capture_output=True).returncode == 0


needs_r = pytest.mark.skipif(not _r_has_coloc(), reason="Rscript with coloc not available")
needs_rscript = pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript not available")


def _logsum(x: np.ndarray) -> float:
    m = float(np.max(x))
    return m + float(np.log(np.sum(np.exp(x - m))))


def closed_form_abf(z1, v1, w1, z2, v2, w2, p1, p2, p12):
    l1 = 0.5 * (np.log(1 - w1 / (w1 + v1)) + (w1 / (w1 + v1)) * z1 ** 2)
    l2 = 0.5 * (np.log(1 - w2 / (w2 + v2)) + (w2 / (w2 + v2)) * z2 ** 2)
    s1, s2, s12 = _logsum(l1), _logsum(l2), _logsum(l1 + l2)
    a, b = s1 + s2, s12
    h3 = np.log(p1) + np.log(p2) + max(a, b) + np.log(np.exp(a - max(a, b)) - np.exp(b - max(a, b)))
    lh = np.array([0.0, np.log(p1) + s1, np.log(p2) + s2, h3, np.log(p12) + s12])
    pp = np.exp(lh - _logsum(lh))
    snp = np.exp((l1 + l2) - s12)
    return pp, snp


def region(m: int, rho: float):
    idx = np.arange(m)
    return rho ** np.abs(idx[:, None] - idx[None, :])


def datasets(z1: np.ndarray, z2: np.ndarray, se1: float, se2: float):
    m = len(z1)
    snp = [f"rs{i}" for i in range(m)]
    d1 = ColocDataset(snp=snp, beta=list(z1 * se1), varbeta=[se1 ** 2] * m, N=30000.0, type="quant",
                      MAF=[0.25] * m, sdY=1.0)
    d2 = ColocDataset(snp=snp, beta=list(z2 * se2), varbeta=[se2 ** 2] * m, N=50000.0, type="cc",
                      MAF=[0.25] * m, s=0.2)
    return d1, d2


def test_parse_response_reads_abf_and_susie_and_rejects_failures():
    payload = {"tasks": [{"id": "a", "ok": True, "method": "abf", "pp": [0.1, 0.1, 0.1, 0.1, 0.6], "nsnps": 2,
                          "snp": ["rs1", "rs2"], "snp_pp_h4": [0.3, 0.7]},
                         {"id": "s", "ok": True, "method": "susie", "n_cs1": 2, "n_cs2": 1, "n_pairs": 2,
                          "max_pp_h4": None}]}
    out = parse_response(payload, {"a", "s"})
    assert isinstance(out["a"], AbfResult) and out["a"].lead_variant() == "rs2"
    assert isinstance(out["s"], SusieResult) and out["s"].max_pp_h4 is None
    with pytest.raises(ColocBackendError):
        parse_response({"tasks": [{"id": "a", "ok": False, "method": "abf", "error": "boom"}]}, {"a"})
    with pytest.raises(ColocBackendError):
        parse_response(payload, {"a", "s", "missing"})


def test_datasets_from_a_harmonized_table():
    h = pd.DataFrame({"rsid": ["rs1", "rs2"], "beta_p": [0.1, -0.2], "se_p": [0.01, 0.02], "eaf_p": [0.2, 0.9],
                      "n_p": [1000.0, 1200.0], "beta_o": [0.3, 0.1], "se_o": [0.05, 0.05], "eaf_o": [0.2, 0.9]})
    d1 = pqtl_dataset(h)
    d2 = outcome_dataset(h, 250, 750)
    assert d1.varbeta == pytest.approx([1e-4, 4e-4]) and d1.N == 1100.0 and d1.sdY == 1.0
    assert d1.MAF == pytest.approx([0.2, 0.1])
    assert d2.type == "cc" and d2.s == pytest.approx(0.25) and d2.N == 1000.0
    with pytest.raises(ColocBackendError):
        outcome_dataset(h, 0, 750)


@needs_r
def test_r_coloc_abf_matches_the_closed_form():
    m = 80
    ld = region(m, 0.7)
    z1 = 6.0 * ld[40] + 0.3 * np.cos(np.arange(m))
    z2 = 4.0 * ld[40] + 0.3 * np.sin(np.arange(m))
    se1, se2 = 0.015, 0.04
    d1, d2 = datasets(z1, z2, se1, se2)
    for p12 in (5e-6, 1e-6, 1e-5):
        res = RscriptColoc(r_libs=R_LIBS).run([ColocTask(id="t", method="abf", p1=1e-4, p2=1e-4, p12=p12, d1=d1, d2=d2)])["t"]
        pp, snp = closed_form_abf(z1, se1 ** 2, 0.15 ** 2, z2, se2 ** 2, 0.2 ** 2, 1e-4, 1e-4, p12)
        assert list(res.pp) == pytest.approx(list(pp), rel=1e-6, abs=1e-10)
        got = np.array([res.snp_pp_h4[f"rs{i}"] for i in range(m)])
        assert got == pytest.approx(snp, rel=1e-6, abs=1e-12)


@needs_r
def test_r_coloc_abf_estimates_sdy_for_a_quantitative_dataset_sent_without_one():
    m, n2, sd_y = 80, 600.0, 2.0
    ld = region(m, 0.7)
    maf = 0.1 + 0.3 * (np.arange(m) % 5) / 4
    v2 = sd_y ** 2 / (2 * n2 * maf * (1 - maf))          # the variance coloc's sdY estimate inverts exactly
    z1, z2, se1 = 6.0 * ld[40], 5.0 * ld[40], 0.015
    snp = [f"rs{i}" for i in range(m)]
    d1 = ColocDataset(snp=snp, beta=list(z1 * se1), varbeta=[se1 ** 2] * m, N=30000.0, type="quant", MAF=list(maf), sdY=1.0)
    d2 = ColocDataset(snp=snp, beta=list(z2 * np.sqrt(v2)), varbeta=list(v2), N=n2, type="quant", MAF=list(maf))
    res = RscriptColoc(r_libs=R_LIBS).run([ColocTask(id="t", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=d1, d2=d2)])["t"]
    pp, _ = closed_form_abf(z1, se1 ** 2, 0.15 ** 2, z2, v2, (0.15 * sd_y) ** 2, 1e-4, 1e-4, 5e-6)
    assert list(res.pp) == pytest.approx(list(pp), rel=1e-6, abs=1e-10)


@needs_r
def test_r_coloc_abf_known_answers_shared_and_distinct_signals():
    ld = region(60, 0.6)
    d1, d2 = datasets(10 * ld[30], 10 * ld[30], 0.02, 0.02)
    d1b, d2b = datasets(10 * ld[30], 10 * ld[10], 0.02, 0.02)
    out = RscriptColoc(r_libs=R_LIBS).run([
        ColocTask(id="shared", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=d1, d2=d2),
        ColocTask(id="distinct", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=d1b, d2=d2b)])
    assert out["shared"].pp[4] > 0.95 and out["shared"].lead_variant() == "rs30"
    assert out["distinct"].pp[3] > 0.95


@needs_r
def test_r_coloc_susie_finds_two_shared_signals():
    ld = region(60, 0.5)
    z = 9 * ld[12] + 8 * ld[47]
    d1, d2 = datasets(z, z, 0.02, 0.02)
    out = RscriptColoc(r_libs=R_LIBS).run([ColocTask(id="g", method="susie", p1=1e-4, p2=1e-4, p12=5e-6, d1=d1, d2=d2,
                                                     LD=ld.tolist())])["g"]
    assert out.n_cs1 >= 2 and out.n_cs2 >= 2
    assert out.max_pp_h4 is not None and out.max_pp_h4 > 0.8


# ---- exact field access in coloc_run.R ------------------------------------------------------------------

def test_the_r_script_reads_no_field_with_the_dollar_operator():
    code = "\n".join(line.split("#", 1)[0] for line in R_SCRIPT.read_text().splitlines())
    assert "$" not in code                                         # `$` matches names partially; `[[ ]]` does not
    for field in ("snp", "beta", "varbeta", "N", "type", "MAF", "sdY", "s", "d1", "d2", "p1", "p2", "p12", "id", "method",
                  "LD", "tasks"):
        assert f'[["{field}"]]' in code, field
    assert re.findall(r"p1 = t\[\[\"p1\"\]\], p2 = t\[\[\"p2\"\]\], p12 = t\[\[\"p12\"\]\]", code) and "coloc.abf(" in code
    assert "coloc.susie(" in code and "runsusie(" in code          # the two registered methods, priors passed through


@needs_rscript
def test_r_started_by_the_tests_treats_a_partial_dollar_match_as_an_error():
    assert "warnPartialMatchDollar = TRUE" in open(os.environ["R_PROFILE_USER"]).read()
    partial = subprocess.run(["Rscript", "-e", 'd <- list(snp = c("a", "b"), beta = 1); x <- d$s; cat("matched")'],
                             capture_output=True, text=True)
    assert partial.returncode != 0 and "partial match of 's' to 'snp'" in partial.stderr and "matched" not in partial.stdout
    exact = subprocess.run(["Rscript", "-e", 'd <- list(snp = c("a", "b"), beta = 1); cat(is.null(d[["s"]]))'],
                           capture_output=True, text=True)
    assert exact.returncode == 0 and exact.stdout == "TRUE"


def _request() -> dict:
    snp = [f"rs{i}" for i in range(60)]
    d = {"snp": snp, "beta": [0.01] * 60, "varbeta": [1e-4] * 60, "N": 30000.0, "type": "quant", "MAF": [0.25] * 60, "sdY": 1.0}
    return {"tasks": [{"id": "t", "method": "abf", "p1": 1e-4, "p2": 1e-4, "p12": 5e-6, "d1": dict(d), "d2": dict(d)}]}


def _drop(request: dict, where: tuple) -> dict:
    target = request
    for key in where[:-1]:
        target = target[key]
    del target[where[-1]]
    return request


@needs_r
@pytest.mark.parametrize("where,message", [
    (("tasks",), "the request lacks the required field(s): tasks"),
    (("tasks", 0, "p12"), "task 1 lacks the required field(s): p12"),
    (("tasks", 0, "id"), "task 1 lacks the required field(s): id"),
    (("tasks", 0, "method"), "task 1 lacks the required field(s): method"),
    (("tasks", 0, "d2"), "task 1 lacks the required field(s): d2"),
    (("tasks", 0, "d1", "snp"), "task 1, dataset d1 lacks the required field(s): snp"),
    (("tasks", 0, "d2", "varbeta"), "task 1, dataset d2 lacks the required field(s): varbeta"),
    (("tasks", 0, "d2", "N"), "task 1, dataset d2 lacks the required field(s): N"),
])
def test_r_stops_with_a_clear_error_when_the_request_lacks_a_required_field(tmp_path, where, message):
    req, resp = tmp_path / "request.json", tmp_path / "response.json"
    env = {**os.environ, **({"R_LIBS_USER": R_LIBS} if R_LIBS else {})}
    req.write_text(json.dumps(_request()))
    whole = subprocess.run(["Rscript", str(R_SCRIPT), str(req), str(resp)], env=env, capture_output=True, text=True)
    assert whole.returncode == 0 and json.loads(resp.read_text())["tasks"][0]["ok"] is True      # the request is otherwise valid
    resp.unlink()
    req.write_text(json.dumps(_drop(_request(), where)))
    res = subprocess.run(["Rscript", str(R_SCRIPT), str(req), str(resp)], env=env, capture_output=True, text=True)
    assert res.returncode != 0 and message in res.stderr and not resp.exists()


@needs_r
def test_r_requires_the_ld_matrix_of_a_susie_task_and_refuses_another_method(tmp_path):
    req, resp = tmp_path / "request.json", tmp_path / "response.json"
    env = {**os.environ, **({"R_LIBS_USER": R_LIBS} if R_LIBS else {})}
    for method, message in (("susie", "task 1 lacks the required field(s): LD"), ("abc", "task 1 has a method other than abf or susie")):
        request = _request()
        request["tasks"][0]["method"] = method
        req.write_text(json.dumps(request))
        res = subprocess.run(["Rscript", str(R_SCRIPT), str(req), str(resp)], env=env, capture_output=True, text=True)
        assert res.returncode != 0 and message in res.stderr and not resp.exists()


@needs_r
def test_r_reads_each_optional_field_exactly_on_the_qtl_side_under_the_strict_profile():
    m = 80
    ld = region(m, 0.7)
    h = pd.DataFrame({"rsid": [f"rs{i}" for i in range(m)], "beta_p": 6.0 * ld[40] * 0.015, "se_p": 0.015, "eaf_p": 0.25,
                      "n_p": 30000.0, "beta_o": 5.0 * ld[40] * 0.05, "se_o": 0.05, "eaf_o": 0.25, "n_o": 600.0})
    d2 = qtl_dataset(h)
    assert d2.s is None and d2.sdY is None and d2.MAF is not None      # the dataset on which `d$s` returned `snp`
    res = RscriptColoc(r_libs=R_LIBS).run([ColocTask(id="t", method="abf", p1=1e-4, p2=1e-4, p12=5e-6, d1=pqtl_dataset(h), d2=d2)])["t"]
    assert res.nsnps == m and sum(res.pp) == pytest.approx(1.0) and res.lead_variant() == "rs40"
