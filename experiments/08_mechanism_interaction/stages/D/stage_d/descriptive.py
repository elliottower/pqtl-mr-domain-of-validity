"""Descriptive tables 1-8 and 11-15 (PREREG §Other planned analysis). Tables 9 and 10 are built
from the fits in pipeline.py. Every function returns a list of records, which go into
results.json before any CSV is written.
"""
from datetime import date

import numpy as np
import pandas as pd

from stage_d.constants import KARIM_LAUNCHED_TARGETS, N_BOOTSTRAP, SEED, TABLE11_MIN_HYPOTHESES
from stage_d.designs import component_labels

OUTCOME_FREEZE_DATE = date(2026, 9, 30)   # PREREG header: frozen 2026-09-30 (§Starting and stopping rules)
CLASSES = ("aligned", "blocking", "other")
METRICS = ("rs", "or", "sensitivity", "specificity", "ppv", "npv", "balanced_accuracy")


def _records(df: pd.DataFrame) -> list[dict]:
    return df.reset_index(drop=True).astype(object).where(df.notna(), None).to_dict(orient="records")


def cells_2x2(S: np.ndarray, y: np.ndarray) -> np.ndarray:
    """[a, b, c, d] = [S1 advanced, S1 not, S0 advanced, S0 not]."""
    return np.array([np.sum((S == 1) & (y == 1)), np.sum((S == 1) & (y == 0)),
                     np.sum((S == 0) & (y == 1)), np.sum((S == 0) & (y == 0))], dtype=float)


def metrics_2x2(cells: np.ndarray) -> dict[str, np.ndarray]:
    """Vectorized over leading axes; cells[..., 0:4] = a, b, c, d. Zero denominators give NaN."""
    a, b, c, d = (cells[..., i] for i in range(4))
    with np.errstate(divide="ignore", invalid="ignore"):
        sens, spec = a / (a + c), d / (b + d)
        return {"rs": (a / (a + b)) / (c / (c + d)), "or": (a * d) / (b * c), "sensitivity": sens,
                "specificity": spec, "ppv": a / (a + b), "npv": d / (c + d), "balanced_accuracy": (sens + spec) / 2}


def stratum_2x2_table(frame: pd.DataFrame, n_boot: int = N_BOOTSTRAP, by: str = "cls", seed_tag: int = 8) -> list[dict]:
    """Per-stratum 2x2, RS, OR, sensitivity, specificity, PPV, NPV and balanced accuracy, with 90%
    intervals from a bootstrap over whole gene-drug-program components of `frame`."""
    frame = frame.loc[frame["y"].notna()]
    if frame.empty:
        return []
    comp = component_labels(list(frame["gene_ensembl"]), list(frame["programs"]))
    comps, cidx = np.unique(comp, return_inverse=True)
    strata = sorted(frame[by].astype(str).unique())
    per_comp = np.zeros((len(comps), len(strata), 4))
    for j, s in enumerate(strata):
        m = (frame[by].astype(str) == s).to_numpy()
        S, y = frame["S"].to_numpy()[m], frame["y"].to_numpy()[m]
        for i, cell in enumerate(((1, 1), (1, 0), (0, 1), (0, 0))):
            np.add.at(per_comp[:, j, i], cidx[m], ((S == cell[0]) & (y == cell[1])).astype(float))
    rng = np.random.default_rng([SEED, seed_tag])
    weights = rng.multinomial(len(comps), np.full(len(comps), 1 / len(comps)), size=n_boot)
    boot = np.einsum("bc,csk->bsk", weights, per_comp)
    point = per_comp.sum(axis=0)
    pm, bm = metrics_2x2(point), metrics_2x2(boot)
    out = []
    for j, s in enumerate(strata):
        rec = {"stratum": s, "a_supportive_advanced": int(point[j, 0]), "b_supportive_not": int(point[j, 1]),
               "c_not_supportive_advanced": int(point[j, 2]), "d_not_supportive_not": int(point[j, 3])}
        for k in METRICS:
            vals = bm[k][:, j]
            vals = vals[np.isfinite(vals)]
            rec[k] = float(pm[k][j])
            rec[f"{k}_q05"] = float(np.quantile(vals, 0.05)) if len(vals) else float("nan")
            rec[f"{k}_q95"] = float(np.quantile(vals, 0.95)) if len(vals) else float("nan")
        rec["bootstrap_finite"] = {k: int(np.isfinite(bm[k][:, j]).sum()) for k in ("rs", "or")}
        out.append(rec)
    return out


def table1_funnel(funnel_a: pd.DataFrame, s1_pre: pd.DataFrame) -> list[dict]:
    rows = _records(funnel_a.assign(stage="A"))
    n = len(s1_pre)
    rows.append({"stage": "D", "step": "S1 held-out hypotheses", "remaining": n, "excluded": 0, "reason": ""})
    for status in ("active", "business_only", "undated", "no_phase"):
        k = int((s1_pre["status_24"] == status).sum())
        n -= k
        rows.append({"stage": "D", "step": f"primary outcome: {status}", "remaining": n, "excluded": k,
                     "reason": f"status_24 = {status}"})
    other = int(((s1_pre["y"].notna()) & (s1_pre["cls"] == "other")).sum())
    rows.append({"stage": "D", "step": "H1 only: class other", "remaining": n - other, "excluded": other,
                 "reason": "mechanism class other"})
    return rows


def table3_by_class(s1_pre: pd.DataFrame) -> list[dict]:
    out = []
    for cls in CLASSES:
        f = s1_pre.loc[s1_pre["cls"] == cls]
        sup = f.loc[f["S"] == 1]
        per_gene = f.groupby("gene_ensembl").size()
        per_ind = f.groupby("indication_id").size()
        out.append({"class": cls, "n": len(f), "genes": int(f["gene_ensembl"].nunique()),
                    "indications": int(f["indication_id"].nunique()),
                    "somascan_share": float(f["platform_somascan"].mean()) if len(f) else None,
                    "median_neff": float(f["outcome_neff"].median()) if len(f) else None,
                    "oncology_share": float(f["oncology_f"].mean()) if len(f) else None,
                    "hypotheses_per_gene_median": float(per_gene.median()) if len(f) else None,
                    "hypotheses_per_gene_max": int(per_gene.max()) if len(f) else None,
                    "hypotheses_per_indication_median": float(per_ind.median()) if len(f) else None,
                    "hypotheses_per_indication_max": int(per_ind.max()) if len(f) else None,
                    "supportive_hypotheses": len(sup), "supportive_genes": int(sup["gene_ensembl"].nunique()),
                    "supportive_indications": int(sup["indication_id"].nunique()),
                    "supportive_programs": len({p for ps in sup["programs"] for p in ps})})
    return out


def table4_other_reasons(crosstab_a: pd.DataFrame, s1_pre: pd.DataFrame) -> dict:
    other = s1_pre.loc[s1_pre["cls"] == "other", ["hypothesis_id", "gene_symbol", "indication_id", "class_reason"]]
    return {"crosstab": _records(crosstab_a), "other_hypotheses": _records(other),
            "other_reason_counts": _records(other["class_reason"].value_counts().rename_axis("class_reason")
                                            .reset_index(name="n"))}


def _quantiles(x: pd.Series) -> dict:
    x = x.dropna().astype(float)
    if x.empty:
        return {"n": 0}
    return {"n": int(len(x)), **{f"q{int(q * 100):02d}": float(x.quantile(q)) for q in (0, .05, .25, .5, .75, .95, 1)}}


def table5_evidence(s1_pre: pd.DataFrame) -> dict:
    by_class = pd.crosstab(s1_pre["cls"], s1_pre["state"]).rename_axis(index="class").reset_index()
    by_source = pd.crosstab(s1_pre["instrument_source"], s1_pre["state"]).reset_index()
    not_run = s1_pre.loc[~s1_pre["coloc_run"], "not_run_reason"].fillna("").value_counts()
    dists = {cls: {col: _quantiles(f[col]) for col in ("pp_h4", "pp_h3", "n_shared", "frac_pqtl_retained",
                                                       "frac_outcome_retained")}
             for cls, f in s1_pre.groupby("cls")}
    return {"state_by_class": _records(by_class), "state_by_instrument_source": _records(by_source),
            "not_run_by_reason": {str(k): int(v) for k, v in not_run.items()}, "distributions_by_class": dists}


MISSINGNESS_COLUMNS = ("pp_h4", "protein_altering", "platform_concordant", "splicing_candidate", "low_coverage",
                       "pre_pqtl_publication", "outcome_neff", "platform")


def table6_missingness(s1_pre: pd.DataFrame) -> list[dict]:
    out = []
    for (cls, state), f in s1_pre.groupby(["cls", "state"]):
        rec = {"class": cls, "state": state, "n": len(f)}
        for col in MISSINGNESS_COLUMNS:
            rec[f"{col}_missing"] = int(f[col].isna().sum())
        out.append(rec)
    return out


def outcome_category(frame: pd.DataFrame) -> pd.Series:
    status = frame["status_24"].astype(str)
    stop = frame["stop_code"].astype(object).where(frame["stop_code"].notna(), "unclassified")
    cat = status.map({"advanced": "advanced", "active": "active", "business_only": "business", "undated": "undated",
                      "no_phase": "no_phase"})
    return cat.where(cat.notna(), stop.astype(str))


def table7_outcomes(s1_pre: pd.DataFrame) -> dict:
    cat = outcome_category(s1_pre)
    decomposition = pd.crosstab(s1_pre["cls"], cat).rename_axis(index="class").reset_index()
    if "last_phase2_end_date" in s1_pre and s1_pre["last_phase2_end_date"].notna().any():
        months = (pd.Timestamp(OUTCOME_FREEZE_DATE) - s1_pre["last_phase2_end_date"]).dt.days / 30.4375
        follow = {cls: _quantiles(months[s1_pre["cls"] == cls]) for cls in CLASSES}
    else:
        follow = {"not_available": "stage C writes no last Phase II completion or termination date"}
    y = s1_pre["y"]
    chembl3 = (s1_pre["chembl_max_phase_for_ind"] >= 3).where(s1_pre["chembl_max_phase_for_ind"].notna())
    both = y.notna() & chembl3.notna()
    agree = {"n_compared": int(both.sum()),
             "agree": int((y[both] == chembl3[both].astype(float)).sum()),
             "chembl_missing": int(s1_pre["chembl_max_phase_for_ind"].isna().sum())}
    approved = pd.crosstab(s1_pre["cls"], s1_pre["approved"]).rename_axis(index="class").reset_index()
    approved.columns = [str(c) for c in approved.columns]
    return {"decomposition_by_class": _records(decomposition), "followup_months_last_phase2": follow,
            "chembl_agreement_phase3": agree, "approved_by_class": _records(approved)}


def table11_per_indication(s1: pd.DataFrame) -> list[dict]:
    f = s1.loc[s1["y"].notna()].copy()
    sizes = f.groupby("indication_id")["hypothesis_id"].transform("size")
    f["row_indication"] = f["indication_id"].where(sizes >= TABLE11_MIN_HYPOTHESES, f"pooled (<{TABLE11_MIN_HYPOTHESES})")
    out = []
    for ind, g in f.groupby("row_indication"):
        rec = {"indication": ind, "n": len(g), "indications_in_row": int(g["indication_id"].nunique()),
               "supportive": int(g["S"].sum()), "advanced": int(g["y"].sum())}
        for cls in ("aligned", "blocking"):
            h = g.loc[g["cls"] == cls]
            rec[f"{cls}_n"] = len(h)
            rec[f"{cls}_or"] = float(metrics_2x2(cells_2x2(h["S"].to_numpy(), h["y"].to_numpy()))["or"])
        out.append(rec)
    return out


def cross_source_agreement(s1_pre: pd.DataFrame) -> tuple[dict, list[dict], list[dict]]:
    """Evidence-state agreement where a hypothesis has a state from both a UKB-PPP and a deCODE
    instrument (stage B `evidence_state_ukbppp`, `evidence_state_decode`): (summary, UKB-PPP x
    deCODE state crosstab, agreement by mechanism class)."""
    cols = ("evidence_state_ukbppp", "evidence_state_decode")
    if any(c not in s1_pre or s1_pre[c].isna().all() for c in cols):
        return {"not_available": "stage B wrote no hypothesis with an evidence state from both UKB-PPP and deCODE"}, [], []
    both = s1_pre.loc[s1_pre[cols[0]].notna() & s1_pre[cols[1]].notna()]
    agree = both[cols[0]] == both[cols[1]]
    xt = pd.crosstab(both[cols[0]], both[cols[1]]).rename_axis(index="ukbppp_state").reset_index()
    xt.columns = [c if c == "ukbppp_state" else f"decode_{c}" for c in xt.columns]
    by_class = [{"class": cls, "n": int((both["cls"] == cls).sum()), "agree": int(agree[both["cls"] == cls].sum())}
                for cls in CLASSES]
    summary = {"n_both_sources": len(both), "genes": int(both["gene_ensembl"].nunique()), "agree": int(agree.sum()),
               "agreement_fraction": float(agree.mean()) if len(both) else None,
               "supportive_either": int(((both[cols[0]] == "supportive") | (both[cols[1]] == "supportive")).sum()),
               "supportive_both": int(((both[cols[0]] == "supportive") & (both[cols[1]] == "supportive")).sum())}
    return summary, _records(xt), by_class


def table12_platform(s1_pre: pd.DataFrame) -> dict:
    ct = pd.crosstab(s1_pre["cls"], s1_pre["platform_concordant"].astype(object).fillna("missing"))
    summary, crosstab, by_class = cross_source_agreement(s1_pre)
    return {"concordance_by_class": _records(ct.rename_axis(index="class").reset_index()),
            "cross_source_state_agreement": summary, "cross_source_state_crosstab": crosstab,
            "cross_source_state_by_class": by_class}


def table13_neuro(s1: pd.DataFrame, n_boot: int = N_BOOTSTRAP) -> list[dict]:
    f = s1.copy()
    f["neuro_group"] = np.where(f["neurological"], "neurological",
                                np.where(f["psych_only"], "psychiatric_only", "neither"))
    f["group_stratum"] = f["neuro_group"] + "|" + f["cls"].astype(str)
    return stratum_2x2_table(f, n_boot=n_boot, by="group_stratum", seed_tag=13)


def table15_karim(df: pd.DataFrame) -> list[dict]:
    cols = ["hypothesis_id", "gene_symbol", "indication_id", "indication_name", "mechanism_class", "evidence_state",
            "status_24", "advanced_24", "heldout", "in_s1", "karim_launched"]
    primary = df["variant"].fillna("primary") == "primary"   # S8/S9-only rows are not universe hypotheses
    return _records(df.loc[primary & df["gene_symbol"].isin(KARIM_LAUNCHED_TARGETS), cols]
                    .sort_values(["gene_symbol", "indication_id"]))


def deviation_frozen_ci_rule(s1_pre: pd.DataFrame, sensitivity: pd.DataFrame | None, n_boot: int = N_BOOTSTRAP) -> dict:
    """The S1 hypotheses by the frozen CI-versus-p verdict of their outcome file (stage B
    `outcome_file_frozen_ci_p_pass`: pass, fail, or no check), the evidence states of the failing ones
    by class, and, where the sensitivity set S_frozen_ci_rule was formed, its per-stratum 2x2 table
    (as table 8)."""
    col = "outcome_file_frozen_ci_p_pass"
    if col not in s1_pre or s1_pre[col].isna().all():
        return {"not_available": "stage B wrote no frozen CI-versus-p verdict for any S1 hypothesis's outcome file"}
    verdict = s1_pre[col].astype(object).map({True: "pass", False: "fail"}).fillna("no_ci_p_check")
    by_verdict = pd.crosstab(s1_pre["cls"], verdict).rename_axis(index="class").reset_index()
    by_verdict.columns = [str(c) for c in by_verdict.columns]
    failing = s1_pre.loc[verdict == "fail"]
    states = pd.crosstab(failing["cls"], failing["state"]).rename_axis(index="class").reset_index()
    states.columns = [str(c) for c in states.columns]
    return {"outcome_file_verdict_by_class": _records(by_verdict),
            "failing_outcome_file_state_by_class": _records(states),
            "stratum_2x2": [] if sensitivity is None else stratum_2x2_table(sensitivity, n_boot=n_boot, seed_tag=16)}

