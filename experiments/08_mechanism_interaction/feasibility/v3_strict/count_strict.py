# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
#     "pandas==2.2.3",
#     "pyarrow==21.0.0",
#     "requests==2.32.3",
#     "openpyxl==3.1.5",
# ]
# ///
"""Outcome-blind strict count for experiment 08: the registered rules of PREREG.md applied
to the inputs of feasibility/v2_all_indications.

Sets counted:
- S1 (primary): indication rule; Phase II+ drug; cis-pQTL instrument; outcome GWAS selected by
  the joint instrument/outcome rule is tier 1 with overlap "no"; related-indication rule;
  held out from the 161 V5/V5.1 keys.
- S4: S1 plus hypotheses flagged overlap_unknown (selected candidate tier 1, overlap unknown).
- H1-eligible: S1 hypotheses classed abundance-aligned or function-blocking.

Blinding:
- clinical_indication is read with drugId, diseaseId, maxClinicalStage; maxClinicalStage is
  converted at once to a boolean Phase II+ flag and dropped. No phase value is kept, written
  or printed.
- drug_mechanism_of_action: mechanismOfAction, actionType, chemblIds, targets, targetType.
- drug_molecule: id, drugType, parentId.
- classification_v5*.csv: gene, disease only. frozen_candidates_v34.csv: disease,
  ot_disease_id only (the Open Targets disease name used by count_coverage.py).
- feasibility/coverage_counts.json is never opened. The pilot label -> ID map is rebuilt
  offline from frozen_candidates_v34.csv names by exact name match in the OT 26.09 disease
  index, with the multiple-sclerosis override of count_coverage.py.
- No MR, colocalization, disease-association, outcome or approval value is read.

Candidate mapping, tiers, qualification and overlap flags are the functions of
count_all_indications.py, loaded from that file (PREREG: "mapped to ontology IDs exactly as
in count_all_indications.py").

Usage:
    cd ~/Documents/GitHub/pqtl-mr-domain-of-validity
    uv run experiments/08_mechanism_interaction/feasibility/v3_strict/count_strict.py
"""
import importlib.util
import io
import json
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq
import requests

HERE = Path(__file__).resolve().parent
V2_DIR = HERE.parent / "v2_all_indications"
V2_SCRIPT = V2_DIR / "count_all_indications.py"
V2_JSON = V2_DIR / "coverage_all_indications.json"
OUT_JSON = HERE / "coverage_strict.json"

_spec = importlib.util.spec_from_file_location("count_all_indications", V2_SCRIPT)
v2 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(v2)

REPO = v2.REPO
FROZEN_CANDIDATES = REPO / "data" / "frozen_candidates_v34.csv"
OT_DIR = v2.OT_DIR

PHASE2_PLUS = {"PHASE_2", "PHASE_2_3", "PHASE_3", "PREAPPROVAL", "APPROVAL"}
ALLOWED_PREFIXES = {"MONDO", "EFO", "Orphanet"}
NON_DISEASE_TAS = {"GO_0008150", "EFO_0000651", "EFO_0001444", "EFO_0002571", "MONDO_0005583"}
MS_OVERRIDE = {"multiple sclerosis": "EFO_0803536"}  # count_coverage.py OT_ID_OVERRIDES

HPA_VERSION = "25.1"
HPA_VERSION_BASIS = ("proteinatlas.org/about/releases lists version 25.1 (release date "
                     "2026.05.25) as current; the zip entry is dated 2026-06-03 and the "
                     "live download's Last-Modified and Content-Length are checked below")

# Mechanism class table (PREREG §Measured variables), first matching row wins.
BIO_TYPES = {"Antibody", "Protein", "Enzyme"}
BIO_ACTIONS = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR"}
OLIGO_ACTIONS = {"ANTISENSE INHIBITOR", "RNAI INHIBITOR"}
SM_ACTIONS = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "INVERSE AGONIST"}
NEUTRALIZING_TOKENS = ("inhibitor", "antagonist", "blocker", "neutraliz", "neutralis", "sequestr")

# Direction tiers (PREREG_direction_concordance.md §3.3).
TIER1_OVERRIDES = {"INVERSE AGONIST": "BLOCKING", "NEGATIVE ALLOSTERIC MODULATOR": "BLOCKING",
                   "POSITIVE ALLOSTERIC MODULATOR": "ACTIVATING", "PARTIAL AGONIST": "ACTIVATING"}
BLOCKING_SUBSTR = ("INHIBITOR", "ANTAGONIST", "BLOCKER", "DEGRADER", "DISRUPT", "SUPPRESSOR")
ACTIVATING_SUBSTR = ("AGONIST", "ACTIVATOR", "OPENER", "STABILISER", "STABILIZER", "RELEASING AGENT")

SOURCE_ORDER = ["ukbppp", "decode", "interval"]
OVERLAP_COL = {"ukbppp": "includes_ukb", "decode": "includes_iceland", "interval": "includes_interval"}
OVERLAP_RANK = {"no_per_cohort_list": 0, "unknown": 1, "yes": 2}
CAND_SOURCE_RANK = {"gwas_catalog": 0, "finngen": 1, "opengwas": 2}


# ---- rules -------------------------------------------------------------------------------

def token_direction(token: str) -> str:
    t = " ".join(str(token).upper().split())
    if t in TIER1_OVERRIDES:
        return TIER1_OVERRIDES[t]
    if any(s in t for s in BLOCKING_SUBSTR):
        return "BLOCKING"
    if any(s in t for s in ACTIVATING_SUBSTR):
        return "ACTIVATING"
    return "AMBIGUOUS"


def is_neutralizing(action: str, moa_text: str) -> bool:
    return action == "BINDING AGENT" and any(k in str(moa_text).lower() for k in NEUTRALIZING_TOKENS)


def effective_action(action: str, moa_text: str) -> str:
    return "INHIBITOR" if is_neutralizing(action, moa_text) else action


def mechanism_class(mol_type: str, action: str, blood: bool) -> tuple[str, str]:
    """Return (class, matched table row)."""
    if mol_type in BIO_TYPES and action in BIO_ACTIONS and blood:
        return "aligned", "row1_biologic_inhibitor_blood"
    if mol_type == "Oligonucleotide" and action in OLIGO_ACTIONS and blood:
        return "aligned", "row2_oligo_blood"
    if action == "DEGRADER" and blood:
        return "aligned", "row3_degrader_blood"
    if mol_type == "Small molecule" and action in SM_ACTIONS:
        return "blocking", "row4_small_molecule_blocker"
    return "other", "row5_anything_else"


def direction_of(action_eff: str) -> str:
    d = token_direction(action_eff) if action_eff else "AMBIGUOUS"
    return {"BLOCKING": "decrease", "ACTIVATING": "increase"}.get(d, "ambiguous")


def indication_rule(dis: pd.DataFrame) -> dict[str, str]:
    """diseaseId -> 'pass' or the first failing reason."""
    out = {}
    for did, ont, tas in zip(dis["id"], dis["ontology"], dis["therapeuticAreas"]):
        tas = set(tas) if tas is not None else set()
        if did.split("_")[0] not in ALLOWED_PREFIXES:
            out[did] = "prefix_not_mondo_efo_orphanet"
        elif bool(ont["isTherapeuticArea"]):
            out[did] = "therapeutic_area_root"
        elif not (tas - NON_DISEASE_TAS):
            out[did] = "no_disease_therapeutic_area"
        else:
            out[did] = "pass"
    return out


def select_outcome(sources: list[str], pools: dict, cands: pd.DataFrame) -> dict:
    """Joint instrument / outcome-GWAS selection for one gene-indication pair."""
    per_source = []
    for s in sources:
        col = OVERLAP_COL[s]
        ranked = []
        for tier, idx in ((1, pools["t1"]), (2, pools["t2"]), (3, pools["t3"])):
            for i in idx:
                ov = cands.at[i, col]
                if s == "decode" and ov == "yes":
                    continue  # a deCODE instrument is never paired with a deCODE outcome GWAS
                ranked.append((tier, OVERLAP_RANK[ov], -float(cands.at[i, "neff"]),
                               CAND_SOURCE_RANK[cands.at[i, "source"]], str(cands.at[i, "accession"]), i))
        if ranked:
            best = min(ranked)
            per_source.append({"source": s, "tier": best[0],
                               "overlap": {0: "no", 1: "unknown", 2: "yes"}[best[1]], "idx": best[5]})
        else:
            per_source.append({"source": s, "tier": None, "overlap": None, "idx": None})
    for sel in per_source:
        if sel["tier"] == 1 and sel["overlap"] == "no":
            return {**sel, "rule": "first_source_tier1_no"}
    first = per_source[0] if per_source else {"source": None, "tier": None, "overlap": None, "idx": None}
    return {**first, "rule": "first_source_fallback",
            "tier1_unknown_available_other_source": any(
                p["tier"] == 1 and p["overlap"] == "unknown" for p in per_source[1:])}


def related_indication_rule(hyp: pd.DataFrame, descendants: dict, drug_col: str) -> pd.Series:
    """True = retained. Within (gene, direction, class), drop X when every Phase II+ drug of X
    is a Phase II+ drug of some retained descendant of X; most specific terms first."""
    keep = pd.Series(True, index=hyp.index)
    for _, g in hyp.groupby(["gene", "direction", "class"]):
        if len(g) < 2:
            continue
        order = sorted(g.index, key=lambda i: len(descendants.get(g.at[i, "diseaseId"], ())))
        retained: list[int] = []
        for i in order:
            desc = descendants.get(g.at[i, "diseaseId"], set())
            covered = set()
            for j in retained:
                if g.at[j, "diseaseId"] in desc:
                    covered |= g.at[j, drug_col]
            if g.at[i, drug_col] and g.at[i, drug_col] <= covered:
                keep[i] = False
            else:
                retained.append(i)
    return keep


def shares_program(hyp: pd.DataFrame) -> pd.Series:
    flag = pd.Series(False, index=hyp.index)
    for _, g in hyp.groupby("gene"):
        drug_count = Counter(d for s in g["drugs"] for d in s)
        for i in g.index:
            flag[i] = any(drug_count[d] > 1 for d in g.at[i, "drugs"])
    return flag


def block(df: pd.DataFrame) -> dict:
    pairs = df[["gene", "diseaseId"]].drop_duplicates()
    return {"hypotheses": int(len(df)), "pairs": int(len(pairs)), "genes": int(df["gene"].nunique()),
            "indications": int(df["diseaseId"].nunique())}


def class_block(df: pd.DataFrame) -> dict:
    return {c: {"hypotheses": int((df["class"] == c).sum()),
                "genes": int(df.loc[df["class"] == c, "gene"].nunique()),
                "pairs": int(len(df.loc[df["class"] == c, ["gene", "diseaseId"]].drop_duplicates())),
                "indications": int(df.loc[df["class"] == c, "diseaseId"].nunique()),
                "direction": dict(Counter(df.loc[df["class"] == c, "direction"]))}
            for c in ("aligned", "blocking", "other")}


def per_gene(df: pd.DataFrame) -> dict:
    return v2.distribution(df[["gene", "diseaseId"]].drop_duplicates().groupby("gene").size())


# ---- main --------------------------------------------------------------------------------

def main() -> None:
    run_ts = datetime.now(timezone.utc).isoformat()
    ens2sym, sym2ens, uni2sym = v2.load_symbol_maps()

    # Universe: Phase II+ drug x target rows. Phase is reduced to a boolean and dropped.
    ci = pq.read_table(OT_DIR / "clinical_indication.parquet",
                       columns=["drugId", "diseaseId", "maxClinicalStage"]).to_pandas()
    n_ci = len(ci)
    ci = ci[ci["maxClinicalStage"].isin(PHASE2_PLUS)][["drugId", "diseaseId"]].drop_duplicates()
    moa = pq.read_table(OT_DIR / "drug_mechanism_of_action.parquet",
                        columns=["mechanismOfAction", "actionType", "chemblIds", "targets",
                                 "targetType"]).to_pandas()
    moa = moa.explode("chemblIds").explode("targets").dropna(subset=["chemblIds", "targets"])
    moa = moa.rename(columns={"chemblIds": "moa_drug", "targets": "ensembl"})
    mol = pq.read_table(OT_DIR / "drug_molecule.parquet", columns=["id", "drugType", "parentId"]).to_pandas()
    parent = dict(zip(mol["id"], mol["parentId"]))
    dtype = dict(zip(mol["id"], mol["drugType"]))
    own = set(moa["moa_drug"])
    ci["moa_drug"] = [d if d in own else (parent.get(d) if parent.get(d) in own else None)
                      for d in ci["drugId"]]
    dt = ci.dropna(subset=["moa_drug"]).merge(moa, on="moa_drug", how="inner")
    dt["gene"] = dt["ensembl"].map(ens2sym)
    dt = dt.dropna(subset=["gene"]).copy()
    dt["drug_type"] = dt["drugId"].map(dtype).fillna(dt["moa_drug"].map(dtype)).fillna("Unknown")
    dt["drug_parent"] = [parent.get(d) if isinstance(parent.get(d), str) else d for d in dt["drugId"]]

    # Localization (HPA secreted to blood, by the target's Ensembl ID).
    hpa = pd.read_csv(io.BytesIO(zipfile.ZipFile(v2.HPA_FILE).read("proteinatlas.tsv")), sep="\t",
                      usecols=["Gene", "Ensembl", "Secretome location"], dtype=str)
    hpa_blood = set(hpa.loc[hpa["Secretome location"].fillna("") == "Secreted to blood", "Ensembl"])
    dt["blood"] = dt["ensembl"].isin(hpa_blood)

    # Mechanism class and direction per drug-target row.
    dt["neutralizing_ba"] = [is_neutralizing(a, m) for a, m in zip(dt["actionType"], dt["mechanismOfAction"])]
    dt["action_eff"] = [effective_action(a, m) for a, m in zip(dt["actionType"], dt["mechanismOfAction"])]
    cr = [mechanism_class(t, a, b) for t, a, b in zip(dt["drug_type"], dt["action_eff"], dt["blood"])]
    dt["class"] = [c for c, _ in cr]
    dt["class_row"] = [r for _, r in cr]
    dt["direction"] = dt["action_eff"].map(direction_of)
    multi = (dt.groupby(["drugId", "gene", "diseaseId"])[["class", "direction"]]
             .apply(lambda g: len(set(zip(g["class"], g["direction"])))))
    drug_gene_multi = int((multi > 1).sum())

    # Hypotheses: gene + indication + direction + class.
    hyp = (dt.groupby(["gene", "diseaseId", "direction", "class"])
           .agg(ensembl=("ensembl", "first"), drugs=("drugId", lambda s: frozenset(s)),
                drug_parents=("drug_parent", lambda s: frozenset(s)),
                any_single_protein=("targetType", lambda s: bool((s == "single protein").any())))
           .reset_index())

    # Indication rule.
    dis = pq.read_table(OT_DIR / "disease.parquet",
                        columns=["id", "name", "exactSynonyms", "parents", "children", "ancestors",
                                 "descendants", "obsoleteTerms", "therapeuticAreas", "ontology"]).to_pandas()
    ind_status = indication_rule(dis)
    hyp["ind_rule"] = hyp["diseaseId"].map(ind_status).fillna("not_in_ot_disease_index")
    descendants = {i: set(d) if d is not None else set() for i, d in zip(dis["id"], dis["descendants"])}

    # Instruments (count_all_indications.py flags).
    inst = v2.load_instruments(uni2sym)
    for s in ("ukbppp", "decode", "interval_st4", "interval_epigraphdb"):
        hyp[s] = hyp["gene"].isin(inst[s]["sym"]) | hyp["ensembl"].isin(inst[s].get("ens", set()))
    hyp["interval"] = hyp["interval_st4"] | hyp["interval_epigraphdb"]
    hyp["pqtl"] = hyp["ukbppp"] | hyp["decode"] | hyp["interval"]

    # Outcome-GWAS candidates (count_all_indications.py functions).
    labels = v2.label_index(dis)
    gcat, gcat_raw = v2.gwas_catalog_candidates()
    fgn, fg_meta = v2.finngen_candidates()
    ogw, og_meta = v2.opengwas_candidates(gcat, fgn, labels)
    cands = pd.concat([gcat, fgn, ogw], ignore_index=True)
    cands["qualifies"] = (cands["case_control"] & cands["european"] & cands["full_sumstats"]
                          & (cands["neff"].fillna(0) >= v2.NEFF_FLOOR))
    tiers = v2.indication_tiers(set(hyp["diseaseId"]), cands, dis)

    pair_sel = {}
    for g, d, u, dc, iv in (hyp[["gene", "diseaseId", "ukbppp", "decode", "interval"]]
                            .drop_duplicates(["gene", "diseaseId"]).itertuples(index=False)):
        srcs = [s for s, f in zip(SOURCE_ORDER, (u, dc, iv)) if f]
        if not srcs:
            continue
        t = tiers.loc[d]
        pools = {"t1": t["tier1"], "t2": t["tier2_descendant"], "t3": t["tier3_ancestor"]}
        pair_sel[(g, d)] = select_outcome(srcs, pools, cands)
    key = list(zip(hyp["gene"], hyp["diseaseId"]))
    sel = [pair_sel.get(k, {}) for k in key]
    hyp["sel_source"] = [s.get("source") for s in sel]
    hyp["sel_tier"] = [s.get("tier") for s in sel]
    hyp["sel_overlap"] = [s.get("overlap") for s in sel]
    hyp["sel_candidate_source"] = [cands.at[s["idx"], "source"] if s.get("idx") is not None else None for s in sel]
    hyp["tier1_unknown_elsewhere"] = [bool(s.get("tier1_unknown_available_other_source")) for s in sel]

    # Held-out keys: V5/V5.1 labels -> OT IDs via frozen_candidates_v34 names (offline).
    frozen = pd.read_csv(FROZEN_CANDIDATES, usecols=["disease", "ot_disease_id"]).drop_duplicates()
    name2id: dict[str, list] = {}
    for i, n in zip(dis["id"], dis["name"]):
        name2id.setdefault(str(n).lower(), []).append(i)
    label_map, label_map_problems = {}, {}
    for lab, name in zip(frozen["disease"], frozen["ot_disease_id"]):
        if name.lower() in MS_OVERRIDE:
            label_map[v2.norm_old(lab)] = MS_OVERRIDE[name.lower()]
            continue
        hits = name2id.get(name.lower(), [])
        if len(hits) == 1:
            label_map[v2.norm_old(lab)] = hits[0]
        else:
            label_map_problems[lab] = hits
    seen = pd.concat([pd.read_csv(v2.CLASSIFICATION_V5, usecols=["gene", "disease"]),
                      pd.read_csv(v2.CLASSIFICATION_V5_1, usecols=["gene", "disease"])])
    seen["diseaseId"] = seen["disease"].map(lambda x: label_map.get(v2.norm_old(x)))
    seen_keys = set(zip(seen["gene"], seen["diseaseId"]))
    pilot_ids = set(label_map.values())
    hyp["heldout"] = [k not in seen_keys for k in key]

    # Cross-check against the v2 count (held-out definition reproduced without coverage_counts.json).
    v2_json = json.loads(V2_JSON.read_text())
    pairs_all = hyp.drop_duplicates(["gene", "diseaseId"])
    v2_no_known = []
    for g, d, u, dc, iv in pairs_all[["gene", "diseaseId", "ukbppp", "decode", "interval"]].itertuples(index=False):
        srcs = [s for s, f in zip(SOURCE_ORDER, (u, dc, iv)) if f]
        t1 = tiers.at[d, "tier1"]
        v2_no_known.append(bool(t1) and v2.no_known_overlap(t1, cands, srcs))
    v2_final = pairs_all[pairs_all["pqtl"].values & pd.Series(v2_no_known, index=pairs_all.index).values
                         & pairs_all["heldout"].values]
    crosscheck = {
        "seen_keys_here": len(seen_keys),
        "seen_keys_v2": v2_json["seen_sets"]["keys"],
        "eligible_pairs_here": int(len(pairs_all)),
        "eligible_pairs_v2": v2_json["counts"]["pairs"]["eligible"]["pairs"],
        "v2_final_heldout_reproduced": {"pairs": int(len(v2_final)), "genes": int(v2_final["gene"].nunique()),
                                        "indications": int(v2_final["diseaseId"].nunique())},
        "v2_final_heldout_reported": v2_json["counts"]["pairs"]["final_heldout"],
        "label_map_problems": label_map_problems,
    }

    # ---- funnel (hypothesis level; pairs reported beside) -----------------------------------
    funnel = []
    cur = hyp

    def step(name: str, mask: pd.Series, reason: str) -> pd.DataFrame:
        nonlocal cur
        dropped = cur[~mask]
        cur = cur[mask]
        funnel.append({"step": name, "reason": reason, "dropped": block(dropped), "remaining": block(cur)})
        return cur

    funnel.append({"step": "phase2plus_hypotheses", "remaining": block(cur)})
    ind_reasons = dict(Counter(cur.loc[cur["ind_rule"] != "pass", "ind_rule"]))
    step("indication_rule", cur["ind_rule"] == "pass", json.dumps(ind_reasons))
    step("cis_pqtl_instrument", cur["pqtl"], "no cis-pQTL gene in UKB-PPP, deCODE ST02 or INTERVAL lists")
    step("selected_outcome_gwas_tier1", cur["sel_tier"] == 1,
         json.dumps(dict(Counter(cur.loc[cur["sel_tier"] != 1, "sel_tier"].map(lambda x: f"tier{int(x)}" if pd.notna(x) else "none")))))
    pre_overlap = cur.copy()
    step("selected_overlap_no", cur["sel_overlap"] == "no",
         json.dumps(dict(Counter(cur.loc[cur["sel_overlap"] != "no", "sel_overlap"]))))
    pool_s1 = cur.copy()

    # Related-indication rule, then held-out.
    pool_s1["related_keep"] = related_indication_rule(pool_s1, descendants, "drugs")
    pool_s1["related_keep_parent"] = related_indication_rule(pool_s1, descendants, "drug_parents")
    cur = pool_s1
    step("related_indication_rule", cur["related_keep"], "ontology ancestor whose Phase II+ drugs are all on retained descendants")
    after_related = cur.copy()
    step("heldout", cur["heldout"], "gene+indication key in classification_v5/v5_1")
    s1 = cur.copy()
    s1["shares_program"] = shares_program_within(s1, after_related)

    # S4: S1 plus overlap_unknown hypotheses.
    pool_s4 = pre_overlap[pre_overlap["sel_overlap"].isin(["no", "unknown"])].copy()
    pool_s4["related_keep"] = related_indication_rule(pool_s4, descendants, "drugs")
    unknown_added = pool_s4[(pool_s4["sel_overlap"] == "unknown") & pool_s4["related_keep"] & pool_s4["heldout"]]
    s4 = pd.concat([s1, unknown_added])
    s4_rerun = pool_s4[pool_s4["related_keep"] & pool_s4["heldout"]]

    # Alternatives recorded for the ambiguity list.
    s1_related_after_heldout = pool_s1[pool_s1["heldout"]].copy()
    s1_related_after_heldout = s1_related_after_heldout[
        related_indication_rule(s1_related_after_heldout, descendants, "drugs")]
    s1_parent_drug_ids = pool_s1[pool_s1["related_keep_parent"] & pool_s1["heldout"]]

    h1 = s1[s1["class"].isin(["aligned", "blocking"])]

    # Cross-tab of class inputs over S1 drug rows (table 4 of PREREG, no outcome).
    s1_keys = set(zip(s1["gene"], s1["diseaseId"], s1["direction"], s1["class"]))
    s1_rows = dt[[k in s1_keys for k in zip(dt["gene"], dt["diseaseId"], dt["direction"], dt["class"])]]
    xtab = (s1_rows.groupby(["drug_type", "actionType", "neutralizing_ba", "blood", "class"])
            .size().reset_index(name="drug_target_indication_rows"))

    counts = {
        "S1": {**block(s1), "pairs_per_gene": per_gene(s1),
               "hypotheses_per_gene": v2.distribution(s1.groupby("gene").size()),
               "by_class": class_block(s1),
               "shares_program": int(s1["shares_program"].sum()),
               "selected_instrument_source": dict(Counter(s1["sel_source"])),
               "selected_outcome_source": dict(Counter(s1["sel_candidate_source"])),
               "any_single_protein_moa": block(s1[s1["any_single_protein"]]),
               "pilot_indication_flag": int(s1["diseaseId"].isin(pilot_ids).sum()),
               "indication_id_prefixes": dict(Counter(d.split("_")[0] for d in s1["diseaseId"].unique()))},
        "S4": {**block(s4), "pairs_per_gene": per_gene(s4), "by_class": class_block(s4),
               "overlap_unknown_added": block(unknown_added),
               "alt_related_rule_rerun_on_no_plus_unknown_pool": block(s4_rerun)},
        "H1_eligible": {**block(h1), "pairs_per_gene": per_gene(h1),
                        "aligned": class_block(s1)["aligned"], "blocking": class_block(s1)["blocking"],
                        "other_excluded": class_block(s1)["other"],
                        "genes_in_both_classes": int(len(set(h1.loc[h1["class"] == "aligned", "gene"])
                                                         & set(h1.loc[h1["class"] == "blocking", "gene"])))},
        "alternatives": {
            "S1_related_rule_applied_after_heldout": block(s1_related_after_heldout),
            "S1_related_rule_with_parent_molecule_ids": block(s1_parent_drug_ids),
            "tier1_no_pool_before_related_and_heldout": block(pool_s1),
            "fallback_first_source_not_tier1_no_but_tier1_unknown_on_later_source":
                int(pre_overlap.loc[pre_overlap["sel_overlap"] == "yes", "tier1_unknown_elsewhere"].sum()),
        },
        "drug_gene_indication_with_multiple_class_direction": drug_gene_multi,
        "class_input_crosstab_S1": xtab.to_dict(orient="records"),
        "actionType_to_direction": {a: direction_of(a) for a in sorted(dt["actionType"].dropna().unique())},
    }

    hpa_head = {}
    try:
        h = requests.head(v2.HPA_URL, timeout=60, allow_redirects=True).headers
        hpa_head = {"last_modified": h.get("Last-Modified"), "content_length": h.get("Content-Length"),
                    "local_bytes": v2.HPA_FILE.stat().st_size}
    except requests.RequestException as e:
        hpa_head = {"error": str(e)}

    result = {
        "experiment": "08_mechanism_interaction feasibility v3: strict registered rules, outcome-blind",
        "date": run_ts[:10], "run_timestamp_utc": run_ts,
        "script": str(Path(__file__).resolve()), "script_sha256": v2.sha256(Path(__file__).resolve()),
        "v2_script_sha256": v2.sha256(V2_SCRIPT), "seed": None,
        "prereg": {"path": str(HERE.parents[1] / "PREREG.md"), "sha256": v2.sha256(HERE.parents[1] / "PREREG.md")},
        "blinding": __doc__.split("Blinding:")[1].split("Candidate mapping")[0].strip(),
        "definitions": {
            "hypothesis": "gene (HGNC symbol of OT MoA target) + OT diseaseId + intervention direction + mechanism class, over Phase II+ drug rows",
            "phase2plus": sorted(PHASE2_PLUS),
            "indication_rule": f"prefix in {sorted(ALLOWED_PREFIXES)}; in OT 26.09 disease index; ontology.isTherapeuticArea false; therapeuticAreas minus {sorted(NON_DISEASE_TAS)} non-empty",
            "blood_secreted": "HPA 'Secretome location' == 'Secreted to blood' for the MoA target Ensembl ID",
            "mechanism_class": "PREREG table, first matching row, per drug-target row; neutralizing BINDING AGENT (mechanismOfAction text contains inhibitor/antagonist/blocker/neutraliz/neutralis/sequestr) treated as INHIBITOR",
            "direction": "PREREG_direction_concordance.md §3.3 tiers on the (effective) actionType of each row",
            "outcome_selection": "per source in order UKB-PPP, deCODE, INTERVAL (those with the gene): first candidate by tier, overlap no<unknown<yes for that source's cohort, larger Neff, GWAS Catalog<FinnGen<OpenGWAS, lower accession; deCODE instrument excludes candidates with includes_iceland == yes; first source with tier1+no used, else first source's selection",
            "related_indication_rule": "within gene x direction x class, most specific first (fewest OT descendants); X dropped if every Phase II+ drugId of X is a drugId of a retained hypothesis on an OT descendant of X; applied after the overlap step and before the held-out filter",
            "heldout": "(gene, diseaseId) absent from classification_v5.csv and classification_v5_1.csv; labels mapped through frozen_candidates_v34.csv ot_disease_id names by exact name match in OT 26.09 (multiple sclerosis -> EFO_0803536)",
            "S4": "S1 plus hypotheses whose selected candidate is tier 1 with overlap unknown, retained by the related rule run on the tier-1 no+unknown pool, held out",
        },
        "releases": {"open_targets": v2.OT_RELEASE, "gwas_catalog": v2.GWASCAT_RELEASE,
                     "hpa": {"version": HPA_VERSION, "basis": HPA_VERSION_BASIS,
                             "zip_entry_date": list(zipfile.ZipFile(v2.HPA_FILE).getinfo("proteinatlas.tsv").date_time),
                             "live_download_head": hpa_head,
                             "secreted_to_blood_genes": len(hpa_blood)},
                     "finngen": fg_meta.get("finngen_projects")},
        "universe_meta": {"clinical_indication_rows": n_ci},
        "crosscheck_v2": crosscheck,
        "funnel": funnel,
        "counts": counts,
        "inputs": {
            "prereg": v2.input_record(HERE.parents[1] / "PREREG.md", "rules"),
            "v2_script": v2.input_record(V2_SCRIPT, "functions for candidates, tiers, instruments"),
            "v2_json": v2.input_record(V2_JSON, "seen_sets.keys, counts.pairs.eligible and final_heldout for cross-check only"),
            "frozen_candidates_v34": v2.input_record(FROZEN_CANDIDATES, "usecols disease, ot_disease_id"),
            "classification_v5": v2.input_record(v2.CLASSIFICATION_V5, "usecols gene, disease"),
            "classification_v5_1": v2.input_record(v2.CLASSIFICATION_V5_1, "usecols gene, disease"),
            "epigraphdb_catalog": v2.input_record(v2.EPIGRAPHDB_CATALOG, "usecols gene_symbol"),
            "ukbppp_cis_genes": v2.input_record(v2.UKBPPP_CIS, "UKB-PPP cis-pQTL gene list"),
            "olink_protein_map_3k": v2.input_record(v2.OLINK_MAP, "usecols HGNC.symbol, ensembl_id"),
            "decode_supplementary": v2.input_record(v2.DECODE_XLSX, "Ferkingstad 2021 ST01, ST02"),
            "sun2018_supplementary": v2.input_record(v2.SUN2018_XLSX, "sheet ST4 only", v2.SUN2018_URL),
            "hgnc": v2.input_record(v2.HGNC_FILE, "symbol, status, ensembl_gene_id, uniprot_ids", v2.HGNC_URL),
            "hpa": v2.input_record(v2.HPA_FILE, "Gene, Ensembl, Secretome location", v2.HPA_URL),
            "opengwas_gwasinfo": v2.input_record(v2.OPENGWAS_FILE, "full gwasinfo listing", v2.OPENGWAS_URL),
            **{f"ot_{k}": v2.input_record(OT_DIR / f"{k}.parquet", "", f"{v2.OT_BASE}/{p}") for k, p in v2.OT_FILES.items()},
            **{f"gwas_catalog_{fn}": v2.input_record(v2.GWASCAT_DIR / fn, "", f"{v2.GWASCAT_BASE}/{ep}")
               for fn, ep in v2.GWASCAT_ENDPOINTS.items()},
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    print(f"wrote {OUT_JSON}")
    print("crosscheck", json.dumps({k: v for k, v in crosscheck.items()}, default=str))
    for f in funnel:
        print(f"{f['step']:<36} remaining {f['remaining']}")
    for s in ("S1", "S4", "H1_eligible"):
        c = counts[s]
        print(f"{s:<12} hyp {c['hypotheses']:>6} pairs {c['pairs']:>6} genes {c['genes']:>5} indications {c['indications']:>5}")


def shares_program_within(s1: pd.DataFrame, retained: pd.DataFrame) -> pd.Series:
    """Flag S1 hypotheses sharing a drug with another retained hypothesis of the same gene."""
    flag = shares_program(retained)
    return flag.reindex(s1.index).fillna(False)


if __name__ == "__main__":
    main()
