"""Orchestration of the follow-up on the work volume (no Modal here; modal_followup.py wraps it).

Layout under `root` (the work volume):
  work/frame.pkl                  the V8 join with every follow-up index (prepare)
  work/designs/<id>.npz           designs of the Bayesian fits (prepare; Q4 designs in prepare_q4)
  work/fits/<fit_id>/             stage_d.fitting attempts and final.json
  work/q1/<model>.jsonl           Q1 bootstrap replicates
  work/q6/<pass>/<k>_<r>.jsonl    Q6 datasets
  results/                        everything reported, each JSON with the declaration
"""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from followup.common import FROZEN_PLAN_SHA256, FollowupError, sha256_file, write_json
from followup.data import (assay_strength, b_and_c, check_unit_integrity, coloc_class, decode_strength, midrank_quintile,
                           s05_state, ukbppp_strength, unit_id)
from followup.q1 import MIN_C_RUNS, MODELS, class_share_figure, describe, make_spec, population, standardize
from followup.q1 import fit as q1_fit
from followup.q2_q3 import ba_bootstrap, q2_table, s17_phase3_frame
from followup.q4_q5 import augment, cell_table, check_evidence, closure, g_columns, pp_bin_table, support
from followup.q7 import registered_terms, survival_table
from stage_d.designs import build_design, save_design
from stage_d.fitspec import FitSpec
from stage_d.gates import h1_gate
from stage_d.pipeline import load_joined
from stage_d.posterior import PRIMARY_PRIOR, PRIORS
from stage_d.sets import _make

FITS = (
    *(FitSpec(fit_id=f"q3_s17_p3s__{p}", design_id="q3_s17_p3s", kind="h1", set_id="S17_phase3_success", prior=p,
              role="q3") for p in PRIORS),
    FitSpec(fit_id="q4_h2_G", design_id="q4_h2_G", kind="h2", set_id="S1", prior=PRIMARY_PRIOR, role="q4"),
    FitSpec(fit_id="q4_h1_G", design_id="q4_h1_G", kind="h1", set_id="S1", prior=PRIMARY_PRIOR, role="q4"),
    FitSpec(fit_id="q4_h2_Gpre", design_id="q4_h2_Gpre", kind="h2", set_id="S1_Gpre_known", prior=PRIMARY_PRIOR, role="q4"),
    FitSpec(fit_id="q4_h1_Gpre", design_id="q4_h1_Gpre", kind="h1", set_id="S1_Gpre_known", prior=PRIMARY_PRIOR, role="q4"),
    FitSpec(fit_id="q5_h1_s05", design_id="q5_h1_s05", kind="h1", set_id="S1_s05", prior=PRIMARY_PRIOR, role="q5"),
    FitSpec(fit_id="q5_h2_s05", design_id="q5_h2_s05", kind="h2", set_id="S1_s05", prior=PRIMARY_PRIOR, role="q5"),
)


def code_digest(package: Path) -> str:
    h = hashlib.sha256()
    for p in sorted(package.glob("*.py")) + sorted(package.glob("*.R")):
        h.update(p.name.encode() + b"\0" + p.read_bytes() + b"\0")
    return h.hexdigest()


def run_fingerprint(token: str, package: Path) -> str:
    return hashlib.sha256(f"{FROZEN_PLAN_SHA256}|{token}|{code_digest(package)}".encode()).hexdigest()


def save_frame(df: pd.DataFrame, root: Path) -> None:
    p = root / "work" / "frame.pkl"
    p.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(p)


def load_frame(root: Path) -> pd.DataFrame:
    return pd.read_pickle(root / "work" / "frame.pkl")


def _save(d, root: Path) -> str:
    path = root / "work" / "designs" / f"{d.design_id}.npz"
    save_design(d, path)
    return sha256_file(path)


def prepare(stages_root: Path, seals: dict[str, str], sun_xlsx: Path, ferk_xlsx: Path, root: Path) -> dict:
    """The join, every index, Q1 descriptive, Q2, the Q3 and Q5 designs, Q7's table, Q6's structure."""
    _, df, info, sets = load_joined(stages_root, seals)
    check_unit_integrity(df)
    df["unit_id"] = unit_id(df)
    df["coloc_class"] = coloc_class(df)
    df["B"], df["C"] = b_and_c(df)
    df["S_05"] = s05_state(df)
    strength, strength_report = assay_strength(df, {"ukbppp": ukbppp_strength(sun_xlsx), "decode": decode_strength(ferk_xlsx)})
    df["assay_cis_signal_strength"] = strength
    neff = df["outcome_neff"].astype(float).fillna(info.imputation.neff_median)
    df["log10_neff"] = np.log10(neff)
    df["platform"] = df["platform"].fillna(info.imputation.platform_mode)
    df["neff_quintile"] = midrank_quintile(df["outcome_neff"])
    df["strength_quintile"] = midrank_quintile(df["assay_cis_signal_strength"])
    save_frame(df, root)
    res = root / "results"

    write_json(res / "q1_descriptive.json", {"strength_matching": strength_report,
                                             "not_run_by_reason": df.loc[~df["coloc_run"].astype(bool), "not_run_reason"]
                                             .value_counts().to_dict(), **describe(df)})
    class_share_figure(df, res / "q1_class_share_by_neff_quintile")
    q2 = q2_table(df)
    q2.to_csv(res / "q2_pp_h4_ge_0.5.csv", index=False)

    designs = {}
    s17 = s17_phase3_frame(sets["S17"].frame)
    designs["q3_s17_p3s"] = _save(build_design("q3_s17_p3s", "h1", s17, "y"), root)
    base = (df["in_s1"] & df["heldout"]).to_numpy()
    s05 = _make(df, base, "S1_s05", "S1 with S_0.5 in place of S", state=df["S_05"].to_numpy())
    designs["q5_h1_s05"] = _save(build_design("q5_h1_s05", "h1", s05.frame), root)
    designs["q5_h2_s05"] = _save(build_design("q5_h2_s05", "h2", s05.frame), root)
    write_json(res / "q5_counts.json", {
        "s05_state_counts_S1": s05.frame["state"].value_counts().to_dict(),
        "s08_state_counts_S1": sets["S1"].frame["state"].value_counts().to_dict(),
        "gate_h1_s05": h1_gate(s05.frame).model_dump(), "gate_h1_s08": h1_gate(sets["S1"].frame).model_dump(),
        "pp_h4_bins_by_direction_S1": pp_bin_table(sets["S1"].frame)})
    write_json(res / "q3_population.json", {"set": "S17 with Phase III success", "pre_outcome": int(len(s17)),
                                            "analysed_h1": int(s17.loc[s17["cls"].isin(["aligned", "blocking"]), "y"]
                                                               .notna().sum()),
                                            "gate_h1": h1_gate(s17).model_dump()})
    write_json(res / "q3_balanced_accuracy.json", ba_bootstrap(s17))

    surv, counts = survival_table(sets["S1"].frame)
    (root / "work" / "q7").mkdir(parents=True, exist_ok=True)
    surv.to_csv(root / "work" / "q7" / "survival.csv", index=False)
    for model in ("H1", "H2"):
        sub = surv.loc[surv["cls"].isin(["aligned", "blocking"])] if model == "H1" else surv
        (root / "work" / "q7" / f"terms_{model}.txt").write_text("\n".join(registered_terms(sub, model)) + "\n")
    write_json(res / "q7_event_counts.json", counts)

    h1rows = sets["S1"].frame
    h1rows = h1rows.loc[h1rows["cls"].isin(["aligned", "blocking"]) & h1rows["y"].notna()]
    pd.DataFrame({"gene_id": h1rows["gene_ensembl"], "indication_id": h1rows["indication_id"],
                  "drug_program_ids": h1rows["drug_program_ids"], "class": h1rows["cls"],
                  "cns": h1rows["neuro_psych"].astype(int)}).to_csv(root / "work" / "q6_structure.csv", index=False)
    write_json(root / "work" / "designs.json", designs)
    return {"designs": designs, "q1_models": [m for m in MODELS], "q7": counts}


def prepare_q4(evidence: pd.DataFrame, disease: pd.DataFrame, root: Path) -> dict:
    df = load_frame(root)
    ev_report = check_evidence(evidence)
    g = support(df, evidence, closure(disease)).set_index("hypothesis_id")
    frame = df.set_index("hypothesis_id").join(g).reset_index()
    save_frame(frame, root)
    s1 = _make(frame, (frame["in_s1"] & frame["heldout"]).to_numpy(), "S1", "S1").frame
    designs = json.loads((root / "work" / "designs.json").read_text())
    for col in ("G", "G_pre"):
        sub = s1 if col == "G" else s1.loc[s1["G_pre"].notna()]
        suffix = "G" if col == "G" else "Gpre"
        gid = sub.set_index("hypothesis_id")[col].astype(float)
        h2 = build_design(f"q4_h2_{suffix}", "h2", sub)
        designs[f"q4_h2_{suffix}"] = _save(augment(h2, f"q4_h2_{suffix}", g_columns(h2, gid, False)), root)
        h1 = build_design(f"q4_h1_{suffix}", "h1", sub)
        designs[f"q4_h1_{suffix}"] = _save(augment(h1, f"q4_h1_{suffix}", g_columns(h1, gid, True)), root)
    write_json(root / "work" / "designs.json", designs)
    s1g = s1.assign(G_pre=s1["G_pre"].astype("float"))
    write_json(root / "results" / "q4_tables.json", {
        "evidence_check": ev_report,
        "G_counts_all": frame["G"].value_counts().to_dict(), "G_pre_counts_all": frame["G_pre"].value_counts(dropna=False).to_dict(),
        "S1_G_pre_missing": int(s1["G_pre"].isna().sum()),
        "cells_G": cell_table(s1, "G"), "cells_G_pre": cell_table(s1g, "G_pre"),
        "S_by_G_S1": pd.crosstab(s1["S"], s1["G"]).to_dict()})
    return designs


def q1_point(root: Path, model: str) -> tuple:
    df = load_frame(root)
    response, with_strength = MODELS[model]
    pop = population(df, model)
    if response == "C" and len(pop) < MIN_C_RUNS:
        return pop, None, None
    spec = make_spec(pop, response, with_strength)
    res = q1_fit(pop, spec)
    if res is None:
        raise FollowupError(f"{model}: the fit on the analysis population failed")
    return pop, spec, standardize(res, pop, spec)
