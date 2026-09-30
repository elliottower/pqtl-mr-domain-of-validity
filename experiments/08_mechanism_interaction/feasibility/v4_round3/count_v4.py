# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
#     "pandas==2.2.3",
#     "pyarrow==21.0.0",
#     "requests==2.32.3",
#     "openpyxl==3.1.5",
# ]
# ///
"""Outcome-blind recount (v4) for experiment 08: the v3 strict count (v3_strict/count_strict.py)
with the round-3 adjudication rules (planning/perplexity_v8_design/01_ADJUDICATION.md, Round 3).

Changes from v3, and only these:
1. Instruments (S1): an instrument counts only if its assay maps to exactly one gene and one
   protein accession in its source. UKB-PPP: Olink assays in olink_protein_map_3k_v1.tsv with
   one HGNC symbol and one UniProt accession. deCODE: ST02 cis rows whose SeqId has, in ST01,
   exactly one Gene, one UniProt and at most one Ensembl ID. INTERVAL: Sun 2018 ST4 cis rows
   on a single-UniProt SOMAmer whose accession maps to exactly one HGNC symbol; the EpiGraphDB
   gene list carries no assay or accession and is not counted. The source-native inclusive
   lists (every gene an assay names, EpiGraphDB included) are a sensitivity; the v3 lists are
   run beside both as a bridge.
2. Mechanism rows: exact duplicate MoA rows are collapsed; class and direction are derived per
   row; per drug-target (canonical drug x gene), if the informative rows (direction not
   ambiguous) disagree in class or direction, the drug-target is other/ambiguous.
3. Program identity: the canonical parent molecule (drug_molecule.parentId where present)
   replaces the raw drugId for drug identity and the related-indication rule.
4. S12: within each gene, S1 hypotheses ordered by ontology depth (deepest first), then
   indication ID, then lowest program ID; kept if none of its programs is already kept.
5. Added counts: CNS hypotheses (therapeuticAreas contains MONDO_0005071 or MONDO_0002025) in
   S1 and in H1 by class; H1 restricted to HPA blood-secreted targets, by class.

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
- Instrument metadata: Olink map columns Assay, UniProt, UniProt2, HGNC.symbol, ensembl_id;
  deCODE ST01 SeqId, Gene, UniProt, Ensembl.Gene.ID, Included in analysis; ST02 gene, UniProt,
  SeqId, cis/trans; Sun 2018 ST4 UniProt, cis/trans; EpiGraphDB catalog gene_symbol only.
- The membership CSVs carry opaque integers only (no names, phase, outcome or evidence). IDs
  are ranks of sha256(ID_SALT + key) so their order carries no alphabetical information; the
  key is not written.

Candidate mapping, tiers, qualification and overlap flags are the functions of
count_all_indications.py, loaded from that file (PREREG: "mapped to ontology IDs exactly as
in count_all_indications.py").

Usage:
    cd ~/Documents/GitHub/pqtl-mr-domain-of-validity
    uv run experiments/08_mechanism_interaction/feasibility/v4_round3/count_v4.py
"""
import hashlib
import importlib.util
import io
import json
import re
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd
import pyarrow.parquet as pq
import requests

HERE = Path(__file__).resolve().parent
V2_DIR = HERE.parent / "v2_all_indications"
V2_SCRIPT = V2_DIR / "count_all_indications.py"
V2_JSON = V2_DIR / "coverage_all_indications.json"
V3_SCRIPT = HERE.parent / "v3_strict" / "count_strict.py"
ADJUDICATION = HERE.parents[3] / "planning" / "perplexity_v8_design" / "01_ADJUDICATION.md"
OUT_JSON = HERE / "coverage_v4.json"
OUT_H1 = HERE / "h1_membership.csv"
OUT_S1 = HERE / "s1_membership.csv"
ID_SALT = "exp08-v4-membership"
CNS_TAS = {"MONDO_0005071", "MONDO_0002025"}
UNIPROT_ACC = re.compile(r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$")

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


def resolve_drug_targets(rows: pd.DataFrame, informative_only: bool = True) -> pd.DataFrame:
    """One class and direction per (program, gene). rows carry program, gene, class, direction.
    Informative rows are those with a non-ambiguous direction; if they disagree in class or
    direction (or there are none), the drug-target is other/ambiguous."""
    out = []
    for (p, g), grp in rows.groupby(["program", "gene"]):
        pairs = set(zip(grp["class"], grp["direction"]))
        inf = {cd for cd in pairs if cd[1] != "ambiguous"} if informative_only else pairs
        if len(inf) == 1:
            c, d = next(iter(inf))
        else:
            c, d = "other", "ambiguous"
        out.append({"program": p, "gene": g, "class_dt": c, "direction_dt": d,
                    "n_distinct_class_direction": len(pairs), "n_informative": len(inf),
                    "conflict": len(inf) > 1})
    return pd.DataFrame(out)


def split_ids(x, sep: str = r"[ ,;|]+") -> list[str]:
    return [t for t in re.split(sep, str(x if x is not None else "")) if t and t not in ("None", "nan", "NA")]


def ukbppp_lists(uni2sym: dict) -> dict:
    """UKB-PPP cis list: strict (one symbol, one UniProt accession per Olink assay), v3
    (single-gene entries) and inclusive (every gene of every entry)."""
    raw = json.loads(v2.UKBPPP_CIS.read_text())
    olink = pd.read_csv(v2.OLINK_MAP, sep="\t", dtype=str,
                        usecols=["Assay", "UniProt", "UniProt2", "HGNC.symbol", "ensembl_id"])
    olink = olink[olink["Assay"].isin(raw)]
    strict_assays, dropped = set(), {}
    for assay, g in olink.groupby("Assay"):
        syms = set(g["HGNC.symbol"].dropna())
        accs = {a.split("-")[0] for u in g["UniProt"].dropna() for a in u.split("_")
                if UNIPROT_ACC.match(a.split("-")[0])}
        if not accs:  # UniProt holds a non-accession label (NTproBNP); UniProt2 carries it
            accs = {a.split("-")[0] for a in g["UniProt2"].dropna() if UNIPROT_ACC.match(a.split("-")[0])}
        if "_" not in assay and len(syms) == 1 and len(accs) == 1:
            strict_assays.add(assay)
        else:
            dropped[assay] = {"n_symbols": len(syms), "n_accessions": len(accs)}
    strict_rows = olink[olink["Assay"].isin(strict_assays)]
    v3_single = {x for x in raw if "_" not in x}
    return {
        "strict": {"sym": strict_assays | set(strict_rows["HGNC.symbol"]),
                   "ens": set(strict_rows["ensembl_id"].dropna())},
        "v3": {"sym": v3_single, "ens": set(olink.loc[olink["HGNC.symbol"].isin(v3_single), "ensembl_id"].dropna())},
        "inclusive": {"sym": {s for x in raw for s in x.split("_")} | set(olink["HGNC.symbol"].dropna()),
                      "ens": set(olink["ensembl_id"].dropna())},
        "meta": {"list_entries": len(raw), "entries_in_olink_map": int(olink["Assay"].nunique()),
                 "strict_assays": len(strict_assays), "dropped_assays": dropped,
                 "entries_absent_from_olink_map": sorted(set(raw) - set(olink["Assay"]))},
    }


def decode_lists() -> dict:
    """Ferkingstad 2021: strict keeps ST02 cis rows whose SeqId has one Gene, one UniProt and
    at most one Ensembl ID in ST01; v3 is the v2 loader; inclusive adds every ST01 gene."""
    wb = openpyxl.load_workbook(v2.DECODE_XLSX, read_only=True)
    st1 = wb["ST01"].iter_rows(min_row=3, values_only=True)
    h1 = [str(c) for c in next(st1)]
    i_seq, i_inc, i_ens = h1.index("SeqId"), h1.index("Included in\nanalysis"), h1.index("Ensembl.Gene.ID")
    i_gene, i_uni = h1.index("Gene"), h1.index("UniProt")
    seq = {}
    for r in st1:
        if r[i_seq] is None or str(r[i_inc]).strip() != "Yes":
            continue
        seq[str(r[i_seq])] = {"genes": set(split_ids(r[i_gene])), "unis": set(split_ids(r[i_uni])),
                              "ens": {e for e in split_ids(r[i_ens]) if e.startswith("ENSG")}}
    st2 = wb["ST02"].iter_rows(min_row=3, values_only=True)
    h2 = [str(c) for c in next(st2)]
    i_g, i_s, i_ct = h2.index("gene\n (prot.)"), h2.index("SeqId"), h2.index("cis/\ntrans")
    out = {k: {"sym": set(), "ens": set()} for k in ("strict", "v3", "inclusive")}
    cis_seqids, strict_seqids, multi_seqids = set(), set(), {}
    for r in st2:
        if r[i_ct] != "cis":
            continue
        s = str(r[i_s])
        info = seq.get(s, {"genes": set(), "unis": set(), "ens": set()})
        sym2 = {x for x in re.split(r"[ ,.;|]+", str(r[i_g])) if x and x != "NA"}
        cis_seqids.add(s)
        out["v3"]["sym"] |= sym2
        out["v3"]["ens"] |= info["ens"]
        out["inclusive"]["sym"] |= sym2 | info["genes"]
        out["inclusive"]["ens"] |= info["ens"]
        if len(info["genes"]) == 1 and len(info["unis"]) == 1 and len(info["ens"]) <= 1:
            strict_seqids.add(s)
            out["strict"]["sym"] |= sym2 | info["genes"]
            out["strict"]["ens"] |= info["ens"]
        else:
            multi_seqids[s] = {"n_genes": len(info["genes"]), "n_uniprot": len(info["unis"]),
                               "n_ensembl": len(info["ens"]), "in_st01_included": s in seq}
    out["meta"] = {"cis_seqids": len(cis_seqids), "strict_seqids": len(strict_seqids),
                   "dropped_seqids": multi_seqids}
    return out


def interval_lists(uni2sym: dict) -> dict:
    """Sun 2018 ST4: strict keeps single-UniProt cis SOMAmers whose accession maps to one HGNC
    symbol; v3 adds all symbols of single-UniProt rows plus the EpiGraphDB list; inclusive
    also takes every accession of multi-UniProt rows."""
    ws = openpyxl.load_workbook(v2.SUN2018_XLSX, read_only=True)[v2.SUN2018_SHEET]
    rows = ws.iter_rows(min_row=5, values_only=True)
    header = [str(c) if c is not None else "" for c in next(rows)]
    i_uni, i_ct = header.index("UniProt"), header.index("cis/ trans")
    strict, single, inclusive = set(), set(), set()
    n = Counter()
    for r in rows:
        if r[i_ct] is None or str(r[i_ct]).strip() != "cis":
            continue
        n["cis_rows"] += 1
        accs = [a for a in re.split(r"[ ,;|]+", str(r[i_uni] or "")) if a]
        syms_all = set().union(*[uni2sym.get(a, set()) for a in accs]) if accs else set()
        inclusive |= syms_all
        if len(accs) != 1:
            n["multi_uniprot_rows"] += 1
            continue
        s = uni2sym.get(accs[0], set())
        single |= s
        if len(s) == 1:
            strict |= s
        else:
            n["single_uniprot_rows_not_one_symbol"] += 1
    epi = pd.read_csv(v2.EPIGRAPHDB_CATALOG, usecols=["gene_symbol"])
    epi_genes = {g for s in epi["gene_symbol"].dropna().astype(str) for g in s.split(";")}
    return {"strict": {"sym": strict, "ens": set()},
            "v3": {"sym": single | epi_genes, "ens": set()},
            "inclusive": {"sym": inclusive | epi_genes, "ens": set()},
            "meta": {**dict(n), "st4_strict_genes": len(strict), "epigraphdb_genes": len(epi_genes),
                     "epigraphdb_genes_not_in_st4_strict": len(epi_genes - strict)}}


def ontology_depth(dis: pd.DataFrame) -> dict[str, int]:
    """Longest path from a root to each term along OT `parents`."""
    parents = {i: list(ps) if ps is not None else [] for i, ps in zip(dis["id"], dis["parents"])}
    depth: dict[str, int] = {}

    def d(i: str) -> int:
        if i not in depth:
            ps = [p for p in parents[i] if p in parents]
            depth[i] = 0 if not ps else 1 + max(d(p) for p in ps)
        return depth[i]

    for i in parents:
        d(i)
    return depth


def s12_one_per_program(s1: pd.DataFrame, depth: dict) -> pd.Series:
    """True = kept. Within gene: deepest indication first, then indication ID, then lowest
    program ID; a hypothesis is kept only if none of its programs is already kept."""
    keep = pd.Series(False, index=s1.index)
    for _, g in s1.groupby("gene"):
        order = sorted(g.index, key=lambda i: (-depth.get(g.at[i, "diseaseId"], 0), g.at[i, "diseaseId"],
                                               min(g.at[i, "programs"])))
        used: set = set()
        for i in order:
            if not (g.at[i, "programs"] & used):
                keep[i] = True
                used |= g.at[i, "programs"]
    return keep


def opaque_ids(keys) -> dict:
    ordered = sorted(set(keys), key=lambda k: hashlib.sha256(f"{ID_SALT}:{k}".encode()).hexdigest())
    return {k: n + 1 for n, k in enumerate(ordered)}


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
        program_count = Counter(p for s in g["programs"] for p in s)
        for i in g.index:
            flag[i] = any(program_count[p] > 1 for p in g.at[i, "programs"])
    return flag


def shares_program_within(s1: pd.DataFrame, retained: pd.DataFrame) -> pd.Series:
    """Flag S1 hypotheses sharing a program with another retained hypothesis of the same gene."""
    return shares_program(retained).reindex(s1.index).fillna(False)


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


def hkey(df: pd.DataFrame) -> list[tuple]:
    return list(zip(df["gene"], df["diseaseId"], df["direction"], df["class"]))


def build_hypotheses(dt: pd.DataFrame, resolved: pd.DataFrame) -> pd.DataFrame:
    """Hypotheses (gene + indication + direction + class) from drug-target-resolved rows."""
    d = dt.merge(resolved[["program", "gene", "class_dt", "direction_dt"]], on=["program", "gene"], how="left")
    d = d.rename(columns={"class_dt": "class", "direction_dt": "direction"})
    return (d.groupby(["gene", "diseaseId", "direction", "class"])
            .agg(ensembl=("ensembl", "first"), programs=("program", lambda s: frozenset(s)),
                 drugs=("drugId", lambda s: frozenset(s)),
                 any_single_protein=("targetType", lambda s: bool((s == "single protein").any())))
            .reset_index())


def run_variant(hyp: pd.DataFrame, lists: dict, tiers: pd.DataFrame, cands: pd.DataFrame,
                descendants: dict) -> dict:
    """Instrument flags from `lists`, joint outcome selection, funnel, S1, S4, H1."""
    h = hyp.copy()
    for s in SOURCE_ORDER:
        h[s] = h["gene"].isin(lists[s]["sym"]) | h["ensembl"].isin(lists[s]["ens"])
    h["pqtl"] = h[SOURCE_ORDER].any(axis=1)

    pair_sel = {}
    for g, d, u, dc, iv in (h[["gene", "diseaseId", "ukbppp", "decode", "interval"]]
                            .drop_duplicates(["gene", "diseaseId"]).itertuples(index=False)):
        srcs = [s for s, f in zip(SOURCE_ORDER, (u, dc, iv)) if f]
        if not srcs:
            continue
        t = tiers.loc[d]
        pools = {"t1": t["tier1"], "t2": t["tier2_descendant"], "t3": t["tier3_ancestor"]}
        pair_sel[(g, d)] = select_outcome(srcs, pools, cands)
    sel = [pair_sel.get(k, {}) for k in zip(h["gene"], h["diseaseId"])]
    h["sel_source"] = [s.get("source") for s in sel]
    h["sel_tier"] = [s.get("tier") for s in sel]
    h["sel_overlap"] = [s.get("overlap") for s in sel]
    h["sel_candidate_source"] = [cands.at[s["idx"], "source"] if s.get("idx") is not None else None for s in sel]
    h["tier1_unknown_elsewhere"] = [bool(s.get("tier1_unknown_available_other_source")) for s in sel]

    funnel = []
    cur = h

    def step(name: str, mask: pd.Series, reason: str) -> None:
        nonlocal cur
        funnel.append({"step": name, "reason": reason, "dropped": block(cur[~mask]), "remaining": block(cur[mask])})
        cur = cur[mask]

    funnel.append({"step": "phase2plus_hypotheses", "remaining": block(cur)})
    step("indication_rule", cur["ind_rule"] == "pass",
         json.dumps(dict(Counter(cur.loc[cur["ind_rule"] != "pass", "ind_rule"]))))
    step("cis_pqtl_instrument", cur["pqtl"], "no gene on the variant's UKB-PPP, deCODE or INTERVAL list")
    step("selected_outcome_gwas_tier1", cur["sel_tier"] == 1,
         json.dumps(dict(Counter(cur.loc[cur["sel_tier"] != 1, "sel_tier"].map(
             lambda x: f"tier{int(x)}" if pd.notna(x) else "none")))))
    pre_overlap = cur.copy()
    step("selected_overlap_no", cur["sel_overlap"] == "no",
         json.dumps(dict(Counter(cur.loc[cur["sel_overlap"] != "no", "sel_overlap"]))))
    pool_s1 = cur.copy()
    pool_s1["related_keep"] = related_indication_rule(pool_s1, descendants, "programs")
    cur = pool_s1
    step("related_indication_rule", cur["related_keep"],
         "ontology ancestor whose Phase II+ programs (canonical parent molecules) are all on retained descendants")
    after_related = cur.copy()
    step("heldout", cur["heldout"], "gene+indication key in classification_v5/v5_1")
    s1 = cur.copy()
    s1["shares_program"] = shares_program_within(s1, after_related)

    pool_s4 = pre_overlap[pre_overlap["sel_overlap"].isin(["no", "unknown"])].copy()
    pool_s4["related_keep"] = related_indication_rule(pool_s4, descendants, "programs")
    unknown_added = pool_s4[(pool_s4["sel_overlap"] == "unknown") & pool_s4["related_keep"] & pool_s4["heldout"]]
    s4 = pd.concat([s1, unknown_added])

    alt_raw = pool_s1[related_indication_rule(pool_s1, descendants, "drugs") & pool_s1["heldout"]]
    alt_after = pool_s1[pool_s1["heldout"]].copy()
    alt_after = alt_after[related_indication_rule(alt_after, descendants, "programs")]
    return {"h": h, "funnel": funnel, "s1": s1, "s4": s4, "unknown_added": unknown_added,
            "h1": s1[s1["class"].isin(["aligned", "blocking"])],
            "alternatives": {
                "S1_related_rule_with_raw_drugIds": block(alt_raw),
                "S1_related_rule_applied_after_heldout": block(alt_after),
                "S4_related_rule_rerun_on_no_plus_unknown_pool": block(
                    pool_s4[pool_s4["related_keep"] & pool_s4["heldout"]]),
                "fallback_first_source_not_tier1_no_but_tier1_unknown_on_later_source":
                    int(pre_overlap.loc[pre_overlap["sel_overlap"] == "yes", "tier1_unknown_elsewhere"].sum())}}


def summarize(r: dict, depth: dict) -> dict:
    s1, s4, h1 = r["s1"], r["s4"], r["h1"]
    al, bl = h1[h1["class"] == "aligned"], h1[h1["class"] == "blocking"]
    s12 = s1[s12_one_per_program(s1, depth)]
    s12_h1 = s12[s12["class"].isin(["aligned", "blocking"])]
    sec = h1[h1["secreted"]]
    sec_al, sec_bl = sec[sec["class"] == "aligned"], sec[sec["class"] == "blocking"]
    cns_al, cns_bl = al[al["cns"]], bl[bl["cns"]]

    def crossover(a: pd.DataFrame, b: pd.DataFrame) -> int:
        return len(set(a["gene"]) & set(b["gene"]))

    return {
        "S1": {**block(s1), "pairs_per_gene": per_gene(s1),
               "hypotheses_per_gene": v2.distribution(s1.groupby("gene").size()),
               "by_class": class_block(s1), "shares_program": int(s1["shares_program"].sum()),
               "selected_instrument_source": dict(Counter(s1["sel_source"])),
               "selected_outcome_source": dict(Counter(s1["sel_candidate_source"])),
               "any_single_protein_moa": block(s1[s1["any_single_protein"]]),
               "pilot_indication_flag": int(s1["pilot"].sum()),
               "secreted_target": block(s1[s1["secreted"]]),
               "indication_id_prefixes": dict(Counter(d.split("_")[0] for d in s1["diseaseId"].unique()))},
        "S4": {**block(s4), "by_class": class_block(s4), "overlap_unknown_added": block(r["unknown_added"])},
        "H1_eligible": {**block(h1), "pairs_per_gene": per_gene(h1),
                        "aligned": class_block(s1)["aligned"], "blocking": class_block(s1)["blocking"],
                        "other_excluded": class_block(s1)["other"],
                        "crossover_genes": crossover(al, bl)},
        "CNS": {"S1": {**block(s1[s1["cns"]]), "by_class": class_block(s1[s1["cns"]])},
                "S1_non_cns": block(s1[~s1["cns"]]),
                "H1": {**block(h1[h1["cns"]]), "aligned": block(cns_al), "blocking": block(cns_bl),
                       "crossover_genes": crossover(cns_al, cns_bl)},
                "H1_non_cns": {"aligned": block(al[~al["cns"]]), "blocking": block(bl[~bl["cns"]])}},
        "H1_secreted_only": {**block(sec), "aligned": block(sec_al), "blocking": block(sec_bl),
                             "crossover_genes": crossover(sec_al, sec_bl),
                             "aligned_not_secreted_check": int((~al["secreted"]).sum())},
        "S12_one_per_program": {**block(s12), "by_class": class_block(s12), "H1_classes": block(s12_h1),
                                "crossover_genes": crossover(s12_h1[s12_h1["class"] == "aligned"],
                                                             s12_h1[s12_h1["class"] == "blocking"])},
        "alternatives": r["alternatives"],
        "funnel": r["funnel"],
    }


def lost_by_source(frm: dict, to: dict, gene_lost_src: dict) -> dict:
    """Hypotheses in variant `frm` and not in `to`, grouped by the sources a gene lost."""
    out = {}
    for level in ("s1", "h1"):
        a, b = frm[level], to[level]
        ka, kb = set(hkey(a)), set(hkey(b))
        lost = a[[k not in kb for k in hkey(a)]]
        gained = b[[k not in ka for k in hkey(b)]]
        attr = Counter("+".join(gene_lost_src.get(g, ())) or "no_source_lost" for g in lost["gene"])
        genes_attr = Counter("+".join(gene_lost_src.get(g, ())) or "no_source_lost" for g in lost["gene"].unique())
        out[level.upper()] = {
            "lost": {**block(lost), "hypotheses_by_sources_lost": dict(attr), "genes_by_sources_lost": dict(genes_attr),
                     "genes_absent_entirely": len(set(a["gene"]) - set(b["gene"])),
                     "by_class": {c: int((lost["class"] == c).sum()) for c in ("aligned", "blocking", "other")}},
            "gained": {**block(gained), "hypotheses_by_class": dict(Counter(gained["class"]))},
        }
    return out


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
    n_moa_rows = len(moa)
    moa = moa.drop_duplicates()  # rule 2: exact duplicate MoA rows collapse
    n_moa_dedup = len(moa)
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
    # rule 3: canonical parent molecule (OT drug_molecule has no parent-of-parent chains; checked below)
    dt["program"] = [parent.get(d) if isinstance(parent.get(d), str) else d for d in dt["drugId"]]
    parent_chains = sum(1 for p in mol["parentId"].dropna() if isinstance(parent.get(p), str))

    # Localization (HPA secreted to blood, by the target's Ensembl ID).
    hpa = pd.read_csv(io.BytesIO(zipfile.ZipFile(v2.HPA_FILE).read("proteinatlas.tsv")), sep="\t",
                      usecols=["Gene", "Ensembl", "Secretome location"], dtype=str)
    hpa_blood = set(hpa.loc[hpa["Secretome location"].fillna("") == "Secreted to blood", "Ensembl"])
    dt["blood"] = dt["ensembl"].isin(hpa_blood)

    # Class and direction per row, then per drug-target (rule 2).
    dt["neutralizing_ba"] = [is_neutralizing(a, m) for a, m in zip(dt["actionType"], dt["mechanismOfAction"])]
    dt["action_eff"] = [effective_action(a, m) for a, m in zip(dt["actionType"], dt["mechanismOfAction"])]
    cr = [mechanism_class(t, a, b) for t, a, b in zip(dt["drug_type"], dt["action_eff"], dt["blood"])]
    dt["class"] = [c for c, _ in cr]
    dt["direction"] = dt["action_eff"].map(direction_of)
    dt_rows = dt[["program", "gene", "drug_type", "actionType", "mechanismOfAction", "class",
                  "direction"]].drop_duplicates()
    resolved = resolve_drug_targets(dt_rows)
    resolved_all = resolve_drug_targets(dt_rows, informative_only=False)
    conflicted = set(zip(resolved.loc[resolved["conflict"], "program"], resolved.loc[resolved["conflict"], "gene"]))
    dt = dt.rename(columns={"class": "class_row", "direction": "direction_row"})

    hyp = build_hypotheses(dt, resolved)
    hyp_all_rows = build_hypotheses(dt, resolved_all)

    # Indication rule, ontology, CNS.
    dis = pq.read_table(OT_DIR / "disease.parquet",
                        columns=["id", "name", "exactSynonyms", "parents", "children", "ancestors",
                                 "descendants", "obsoleteTerms", "therapeuticAreas", "ontology"]).to_pandas()
    ind_status = indication_rule(dis)
    descendants = {i: set(d) if d is not None else set() for i, d in zip(dis["id"], dis["descendants"])}
    depth = ontology_depth(dis)
    cns_ids = {i for i, tas in zip(dis["id"], dis["therapeuticAreas"]) if tas is not None and CNS_TAS & set(tas)}

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

    for h in (hyp, hyp_all_rows):
        h["ind_rule"] = h["diseaseId"].map(ind_status).fillna("not_in_ot_disease_index")
        h["heldout"] = [k not in seen_keys for k in zip(h["gene"], h["diseaseId"])]
        h["pilot"] = h["diseaseId"].isin(pilot_ids)
        h["cns"] = h["diseaseId"].isin(cns_ids)
        h["secreted"] = h["ensembl"].isin(hpa_blood)

    # Instrument lists (rule 1).
    ukb, dec, itv = ukbppp_lists(uni2sym), decode_lists(), interval_lists(uni2sym)
    lists = {v: {"ukbppp": ukb[v], "decode": dec[v], "interval": itv[v]} for v in ("strict", "v3", "inclusive")}

    # Outcome-GWAS candidates (count_all_indications.py functions).
    labels = v2.label_index(dis)
    gcat, gcat_raw = v2.gwas_catalog_candidates()
    fgn, fg_meta = v2.finngen_candidates()
    ogw, og_meta = v2.opengwas_candidates(gcat, fgn, labels)
    cands = pd.concat([gcat, fgn, ogw], ignore_index=True)
    cands["qualifies"] = (cands["case_control"] & cands["european"] & cands["full_sumstats"]
                          & (cands["neff"].fillna(0) >= v2.NEFF_FLOOR))
    tiers = v2.indication_tiers(set(hyp["diseaseId"]), cands, dis)

    runs = {v: run_variant(hyp, lists[v], tiers, cands, descendants) for v in ("strict", "v3", "inclusive")}
    run_all_rows = run_variant(hyp_all_rows, lists["strict"], tiers, cands, descendants)
    strict = runs["strict"]

    # Cross-check against the v2 count with the v3 lists (pair level; rules 2-4 do not touch pairs).
    v2_json = json.loads(V2_JSON.read_text())
    hv3 = runs["v3"]["h"]
    pairs_all = hv3.drop_duplicates(["gene", "diseaseId"])
    v2_no_known = []
    for g, d, u, dc, iv in pairs_all[["gene", "diseaseId", "ukbppp", "decode", "interval"]].itertuples(index=False):
        srcs = [s for s, f in zip(SOURCE_ORDER, (u, dc, iv)) if f]
        t1 = tiers.at[d, "tier1"]
        v2_no_known.append(bool(t1) and v2.no_known_overlap(t1, cands, srcs))
    v2_final = pairs_all[pairs_all["pqtl"].values & pd.Series(v2_no_known, index=pairs_all.index).values
                         & pairs_all["heldout"].values]
    crosscheck = {
        "seen_keys_here": len(seen_keys), "seen_keys_v2": v2_json["seen_sets"]["keys"],
        "eligible_pairs_here": int(len(pairs_all)), "eligible_pairs_v2": v2_json["counts"]["pairs"]["eligible"]["pairs"],
        "v2_final_heldout_reproduced": {"pairs": int(len(v2_final)), "genes": int(v2_final["gene"].nunique()),
                                        "indications": int(v2_final["diseaseId"].nunique())},
        "v2_final_heldout_reported": v2_json["counts"]["pairs"]["final_heldout"],
        "label_map_problems": label_map_problems,
    }

    # Genes and hypotheses lost by source.
    def lost_sources(frm: str) -> dict:
        hs, hf = runs["strict"]["h"], runs[frm]["h"]
        g_s = hs.drop_duplicates("gene").set_index("gene")[SOURCE_ORDER]
        g_f = hf.drop_duplicates("gene").set_index("gene")[SOURCE_ORDER]
        return {g: tuple(s for s in SOURCE_ORDER if g_f.at[g, s] and not g_s.at[g, s]) for g in g_f.index}

    universe_genes = set(hyp.loc[hyp["ind_rule"] == "pass", "gene"])
    strict_pqtl = runs["strict"]["h"].drop_duplicates("gene").set_index("gene")["pqtl"].to_dict()
    list_level = {}
    for frm in ("inclusive", "v3"):
        gl = lost_sources(frm)
        list_level[frm] = {s: {"list_symbols_removed": len(lists[frm][s]["sym"] - lists["strict"][s]["sym"]),
                               "indication_rule_genes_losing_source": sum(1 for g in universe_genes if s in gl.get(g, ())),
                               "of_which_left_with_no_strict_source": sum(
                                   1 for g in universe_genes if s in gl.get(g, ()) and not strict_pqtl[g])}
                           for s in SOURCE_ORDER}
    lost = {f"{frm}_to_strict": lost_by_source(runs[frm], strict, lost_sources(frm)) for frm in ("inclusive", "v3")}

    # Conflicting mechanism rows (rule 2).
    s1_keys = set(hkey(strict["s1"]))
    res_idx = resolved.set_index(["program", "gene"])
    dt_res = dt.join(res_idx[["class_dt", "direction_dt", "conflict"]], on=["program", "gene"])
    in_s1 = [k in s1_keys for k in zip(dt_res["gene"], dt_res["diseaseId"], dt_res["direction_dt"], dt_res["class_dt"])]
    s1_dt = set(zip(dt_res.loc[in_s1, "program"], dt_res.loc[in_s1, "gene"]))
    conflict_counts = {
        "moa_rows_after_explode": int(n_moa_rows), "exact_duplicate_moa_rows_collapsed": int(n_moa_rows - n_moa_dedup),
        "drug_targets_phase2plus_universe": int(len(resolved)),
        "drug_targets_with_more_than_one_distinct_class_direction_row": int((resolved["n_distinct_class_direction"] > 1).sum()),
        "drug_targets_conflicted_informative_rows": int(resolved["conflict"].sum()),
        "drug_targets_no_informative_row": int((resolved["n_informative"] == 0).sum()),
        "drug_targets_conflicted_if_all_rows_count": int(resolved_all["conflict"].sum()),
        "drug_targets_in_strict_S1": len(s1_dt),
        "conflicted_drug_targets_in_strict_S1": len(s1_dt & conflicted),
        "drug_target_indication_rows_conflicted_universe": int(dt_res["conflict"].sum()),
        "drug_target_indication_rows_conflicted_strict_S1": int(dt_res.loc[in_s1, "conflict"].sum()),
        "conflicted_by_row_classes": dict(Counter(
            "|".join(sorted({f"{c}/{d}" for c, d in zip(g["class"], g["direction"]) if d != "ambiguous"}))
            for (p, gn), g in dt_rows.groupby(["program", "gene"]) if (p, gn) in conflicted)),
    }
    program_counts = {"raw_drugIds": int(dt["drugId"].nunique()), "programs": int(dt["program"].nunique()),
                      "drugIds_mapped_to_a_parent": int((dt["drugId"] != dt["program"]).groupby(dt["drugId"]).any().sum()),
                      "parent_of_parent_chains_in_drug_molecule": parent_chains}

    # Cross-tab of class inputs over strict S1 drug rows (no outcome).
    xtab = (dt_res.loc[in_s1].groupby(["drug_type", "actionType", "neutralizing_ba", "blood", "class_row",
                                       "class_dt", "conflict"])
            .size().reset_index(name="drug_target_indication_rows"))

    # Membership files (strict S1, opaque IDs).
    s1 = strict["s1"]
    gid = opaque_ids(s1["gene"])
    iid = opaque_ids(s1["diseaseId"])
    pid = opaque_ids(p for ps in s1["programs"] for p in ps)
    hid = opaque_ids("|".join(k) for k in hkey(s1))
    mem = pd.DataFrame({
        "hypothesis_id": ["|".join(k) for k in hkey(s1)],
        "gene_id": s1["gene"].map(gid).values,
        "indication_id": s1["diseaseId"].map(iid).values,
        "class": s1["class"].values,
        "cns": s1["cns"].astype(int).values,
        "secreted": s1["secreted"].astype(int).values,
        "drug_program_ids": [";".join(str(x) for x in sorted(pid[p] for p in ps)) for ps in s1["programs"]],
    })
    mem["hypothesis_id"] = mem["hypothesis_id"].map(hid)
    mem = mem.sort_values("hypothesis_id").reset_index(drop=True)
    mem.to_csv(OUT_S1, index=False)
    mem[mem["class"].isin(["aligned", "blocking"])].to_csv(OUT_H1, index=False)

    hpa_head = {}
    try:
        hh = requests.head(v2.HPA_URL, timeout=60, allow_redirects=True).headers
        hpa_head = {"last_modified": hh.get("Last-Modified"), "content_length": hh.get("Content-Length"),
                    "local_bytes": v2.HPA_FILE.stat().st_size}
    except requests.RequestException as e:
        hpa_head = {"error": str(e)}

    counts = {v: summarize(runs[v], depth) for v in runs}
    result = {
        "experiment": "08_mechanism_interaction feasibility v4: round-3 rules, outcome-blind",
        "date": run_ts[:10], "run_timestamp_utc": run_ts,
        "script": str(Path(__file__).resolve()), "script_sha256": v2.sha256(Path(__file__).resolve()),
        "v2_script_sha256": v2.sha256(V2_SCRIPT), "v3_script_sha256": v2.sha256(V3_SCRIPT), "seed": None,
        "prereg": {"path": str(HERE.parents[1] / "PREREG.md"), "sha256": v2.sha256(HERE.parents[1] / "PREREG.md")},
        "blinding": __doc__.split("Blinding:")[1].split("Candidate mapping")[0].strip(),
        "definitions": {
            "hypothesis": "gene (HGNC symbol of OT MoA target) + OT diseaseId + intervention direction + mechanism class, over Phase II+ drug rows",
            "phase2plus": sorted(PHASE2_PLUS),
            "indication_rule": f"prefix in {sorted(ALLOWED_PREFIXES)}; in OT 26.09 disease index; ontology.isTherapeuticArea false; therapeuticAreas minus {sorted(NON_DISEASE_TAS)} non-empty",
            "blood_secreted": "HPA 'Secretome location' == 'Secreted to blood' for the MoA target Ensembl ID",
            "mechanism_class_row": "PREREG table, first matching row, per MoA row; neutralizing BINDING AGENT treated as INHIBITOR",
            "mechanism_class_drug_target": "exact duplicate exploded MoA rows (mechanismOfAction, actionType, drug, target, targetType) collapsed; per (canonical program, gene) the distinct (class, direction) of rows with non-ambiguous direction; one value -> that class and direction; none or several -> other/ambiguous",
            "program": "drug_molecule.parentId of the clinical_indication drugId where present, else the drugId",
            "instruments_strict": "UKB-PPP: list entry is an Olink assay with one HGNC symbol and one UniProt accession (isoform suffix stripped; UniProt2 used where UniProt holds a label); deCODE: ST02 cis rows whose SeqId has one Gene, one UniProt, <=1 Ensembl ID in ST01 (included in analysis); INTERVAL: ST4 cis rows with one UniProt mapping to one HGNC symbol; EpiGraphDB list not counted",
            "instruments_inclusive": "UKB-PPP: every gene of every list entry (multi-gene assays split); deCODE: every ST02 cis gene and every ST01 gene/Ensembl ID of its SeqId; INTERVAL: every HGNC symbol of every UniProt on ST4 cis rows, plus the EpiGraphDB gene list",
            "instruments_v3": "the lists of count_all_indications.py / v3_strict",
            "outcome_selection": "as v3: per source in order UKB-PPP, deCODE, INTERVAL; tier, overlap no<unknown<yes, larger Neff, GWAS Catalog<FinnGen<OpenGWAS, lower accession; deCODE instrument excludes includes_iceland == yes",
            "related_indication_rule": "within gene x direction x class, most specific first (fewest OT descendants); X dropped if every program of X is a program of a retained hypothesis on an OT descendant of X; after the overlap step, before held-out",
            "heldout": "(gene, diseaseId) absent from classification_v5.csv and classification_v5_1.csv (labels mapped as in v3)",
            "S4": "S1 plus tier-1 overlap-unknown hypotheses retained by the related rule on the no+unknown pool, held out",
            "S12": "within gene, S1 hypotheses by ontology depth descending (longest path to a root along OT parents), then indication ID, then lowest program ID; kept if none of its programs is already kept",
            "cns": f"OT therapeuticAreas of the indication contains one of {sorted(CNS_TAS)}",
            "H1_secreted_only": "H1 hypotheses whose target is HPA blood-secreted (aligned hypotheses are blood-secreted by the class rule)",
            "crossover_genes": "genes with at least one aligned and one blocking hypothesis in the set",
            "membership_ids": f"opaque integers: rank of sha256('{ID_SALT}:' + key); hypothesis key gene|diseaseId|direction|class; shared by both CSVs",
        },
        "releases": {"open_targets": v2.OT_RELEASE, "gwas_catalog": v2.GWASCAT_RELEASE,
                     "hpa": {"version": HPA_VERSION, "basis": HPA_VERSION_BASIS,
                             "zip_entry_date": list(zipfile.ZipFile(v2.HPA_FILE).getinfo("proteinatlas.tsv").date_time),
                             "live_download_head": hpa_head, "secreted_to_blood_genes": len(hpa_blood)},
                     "finngen": fg_meta.get("finngen_projects")},
        "universe_meta": {"clinical_indication_rows": n_ci, "phase2plus_hypotheses": int(len(hyp))},
        "crosscheck_v2": crosscheck,
        "instrument_lists": {
            "sizes": {v: {s: {"symbols": len(lists[v][s]["sym"]), "ensembl": len(lists[v][s]["ens"])}
                          for s in SOURCE_ORDER} for v in lists},
            "ukbppp": ukb["meta"], "decode": dec["meta"], "interval": itv["meta"],
            "genes_lost_list_level": list_level},
        "mechanism_rows": conflict_counts,
        "programs": program_counts,
        "counts": counts,
        "lost_by_source": lost,
        "alternative_conflict_rule_all_rows": {"S1": block(run_all_rows["s1"]),
                                               "H1": {**block(run_all_rows["h1"]),
                                                      "by_class": class_block(run_all_rows["s1"])}},
        "class_input_crosstab_strict_S1": xtab.to_dict(orient="records"),
        "actionType_to_direction": {a: direction_of(a) for a in sorted(dt["actionType"].dropna().unique())},
        "outputs": {"h1_membership": v2.input_record(OUT_H1, "H1 (strict) membership, opaque IDs"),
                    "s1_membership": v2.input_record(OUT_S1, "S1 (strict) membership, opaque IDs")},
        "inputs": {
            "prereg": v2.input_record(HERE.parents[1] / "PREREG.md", "rules"),
            "adjudication": v2.input_record(ADJUDICATION, "Round 3 rules"),
            "v2_script": v2.input_record(V2_SCRIPT, "functions for candidates, tiers, symbol maps"),
            "v3_script": v2.input_record(V3_SCRIPT, "copied and edited into this script"),
            "v2_json": v2.input_record(V2_JSON, "seen_sets.keys, counts.pairs.eligible and final_heldout for cross-check only"),
            "frozen_candidates_v34": v2.input_record(FROZEN_CANDIDATES, "usecols disease, ot_disease_id"),
            "classification_v5": v2.input_record(v2.CLASSIFICATION_V5, "usecols gene, disease"),
            "classification_v5_1": v2.input_record(v2.CLASSIFICATION_V5_1, "usecols gene, disease"),
            "epigraphdb_catalog": v2.input_record(v2.EPIGRAPHDB_CATALOG, "usecols gene_symbol"),
            "ukbppp_cis_genes": v2.input_record(v2.UKBPPP_CIS, "UKB-PPP cis-pQTL gene list"),
            "olink_protein_map_3k": v2.input_record(v2.OLINK_MAP, "usecols Assay, UniProt, UniProt2, HGNC.symbol, ensembl_id"),
            "decode_supplementary": v2.input_record(v2.DECODE_XLSX, "Ferkingstad 2021 ST01 (ids), ST02 (gene, SeqId, cis/trans)"),
            "sun2018_supplementary": v2.input_record(v2.SUN2018_XLSX, "sheet ST4 only (UniProt, cis/trans)", v2.SUN2018_URL),
            "hgnc": v2.input_record(v2.HGNC_FILE, "symbol, status, ensembl_gene_id, uniprot_ids", v2.HGNC_URL),
            "hpa": v2.input_record(v2.HPA_FILE, "Gene, Ensembl, Secretome location", v2.HPA_URL),
            "opengwas_gwasinfo": v2.input_record(v2.OPENGWAS_FILE, "full gwasinfo listing", v2.OPENGWAS_URL),
            **{f"ot_{k}": v2.input_record(OT_DIR / f"{k}.parquet", "", f"{v2.OT_BASE}/{p}") for k, p in v2.OT_FILES.items()},
            **{f"gwas_catalog_{fn}": v2.input_record(v2.GWASCAT_DIR / fn, "", f"{v2.GWASCAT_BASE}/{ep}")
               for fn, ep in v2.GWASCAT_ENDPOINTS.items()},
        },
    }
    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    print(f"wrote {OUT_JSON}, {OUT_H1}, {OUT_S1}")
    print("crosscheck", json.dumps(crosscheck, default=str))
    for v, c in counts.items():
        h1c = c["H1_eligible"]
        print(f"[{v}] S1 {c['S1']['hypotheses']} hyp / {c['S1']['genes']} genes; S4 {c['S4']['hypotheses']}; "
              f"H1 aligned {h1c['aligned']['hypotheses']} ({h1c['aligned']['genes']} genes), "
              f"blocking {h1c['blocking']['hypotheses']} ({h1c['blocking']['genes']} genes), "
              f"crossover {h1c['crossover_genes']}")


if __name__ == "__main__":
    main()
