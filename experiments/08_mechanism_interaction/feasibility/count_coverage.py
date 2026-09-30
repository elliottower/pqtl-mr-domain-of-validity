# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "pandas==2.2.3",
#     "requests==2.32.3",
#     "openpyxl==3.1.5",
# ]
# ///
"""Outcome-blind feasibility count for experiment 08 (mechanism x MR interaction).

Counts, for the Phase III (and Phase II-only) drug-target-indication universe over the
24 indications of the V3.4 frozen candidate set, how many gene-indication pairs can be
instrumented by each cis-QTL source, how many have a post-V5.1 outcome GWAS, how many
are new relative to the V5/V5.1 analysed sets, and how many genes lack a frozen
mechanism assignment.

No trial outcome, adjudication result, MR estimate or p-value is read. Files that mix
identifiers with such columns are loaded with `usecols` restricted to identifier columns.

Usage:
    cd ~/Documents/GitHub/pqtl-mr-domain-of-validity
    uv run experiments/08_mechanism_interaction/feasibility/count_coverage.py
"""
import hashlib
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
INPUTS = HERE / "inputs"
OT_CACHE = INPUTS / "ot_cache"
OUT_JSON = HERE / "coverage_counts.json"

TRANSPORT = Path.home() / "Documents" / "GitHub" / "transport-wrapper" / "DRUGS_EXPANDED"

FROZEN_CANDIDATES = REPO / "data" / "frozen_candidates_v34.csv"
EPIGRAPHDB_CATALOG = REPO / "zenodo_export" / "data" / "epigraphdb" / "v34_mr_catalog.csv"
CLASSIFICATION_V5 = REPO / "results" / "v5" / "classification_v5.csv"
CLASSIFICATION_V5_1 = REPO / "results" / "v5_1" / "classification_v5_1.csv"
DISEASE_GWAS_PY = REPO / "code" / "disease_gwas.py"
PRESPEC_V6 = REPO / "protocol" / "v6" / "PRESPEC_V6.md"
V7_BLIND = REPO / "protocol" / "v7" / "gene_tier_assignments_blind.md"
UKBPPP_CIS = TRANSPORT / "pqtl_sources" / "ukbppp_cis_pqtl_genes.json"
UKBPPP_ALL = TRANSPORT / "pqtl_sources" / "ukbppp_all_genes.json"
RESCUE_FUNNEL = TRANSPORT / "pqtl_sources" / "rescue_funnel.json"
DECODE_XLSX = INPUTS / "ferkingstad2021_MOESM4_ESM.xlsx"
DECODE_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-021-00978-w/"
              "MediaObjects/41588_2021_978_MOESM4_ESM.xlsx")
GWASINFO_CACHE = INPUTS / "opengwas_gwasinfo_subset.json"
EQTLGEN_FILE = INPUTS / "eqtlgen_cis_eQTLsFDR0.05_2019-12-11.txt.gz"
EQTLGEN_URL = ("https://download.gcc.rug.nl/downloads/eqtlgen/cis-eqtl/"
               "2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz")
EQTLGEN_SUMMARY = INPUTS / "eqtlgen_cis_gene_summary.csv"
GWS = 5e-8
ENV2 = REPO / ".env2"

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"
OPENGWAS = "https://api.opengwas.io/api"

# Stage strings -> numeric phase, identical to candidate_selection_v34.py PHASE_MAP.
PHASE_MAP = {
    "PRECLINICAL": -1, "IND": 0, "EARLY_PHASE_1": 0.5,
    "PHASE_1": 1, "PHASE_1_2": 1.5, "PHASE_2": 2,
    "PHASE_2_3": 2.5, "PHASE_3": 3, "PREAPPROVAL": 3.5,
    "APPROVAL": 4, "UNKNOWN": -1,
}

# The OT search returns no exact-name hit for "multiple sclerosis". The V3.3 disease
# table in candidate_selection_v34.py records MONDO_0005301, which OT 26.09 lists as an
# obsolete term of EFO_0803536 (children: relapsing-remitting and chronic progressive MS).
OT_ID_OVERRIDES = {"multiple sclerosis": "EFO_0803536"}

# Six further indications named in pqtl_sources/non_ukb_gwas_sources.json that are not
# among the 24 frozen indications. Counted as a sensitivity universe only.
EXTRA_INDICATIONS = {
    "Coronary artery disease": "MONDO_0005010",
    "Type 2 diabetes": "MONDO_0005148",
    "Celiac disease": "MONDO_0005130",
    "Atopic dermatitis": "MONDO_0004980",
    "IgA nephropathy": "MONDO_0005342",
    "Melanoma": "MONDO_0005105",
}

V6_SECTION_TO_V7 = {"A1": "D", "A2": "C", "A3": "C", "B1": "B", "B2": "B"}

OT_SEARCH_QUERY = """
query S($q: String!) {
  search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 10}) {
    hits { id name }
  }
}"""

OT_DRUG_QUERY = """
query DiseaseDrugs($diseaseId: String!) {
  meta { dataVersion { year month } }
  disease(efoId: $diseaseId) {
    id
    name
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id
          name
          drugType
          mechanismsOfAction {
            rows {
              actionType
              targets { id approvedSymbol }
            }
          }
        }
      }
    }
  }
}"""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def input_record(path: Path, note: str) -> dict:
    return {"path": str(path), "sha256": sha256(path), "note": note}


def ot_post(query: str, variables: dict) -> dict:
    for attempt in range(4):
        r = requests.post(OT_GRAPHQL, json={"query": query, "variables": variables}, timeout=120)
        if r.status_code == 200:
            return r.json()
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"Open Targets request failed with status {r.status_code}")


def resolve_ot_id(ot_name: str) -> str:
    if ot_name.lower() in OT_ID_OVERRIDES:
        return OT_ID_OVERRIDES[ot_name.lower()]
    hits = ot_post(OT_SEARCH_QUERY, {"q": ot_name})["data"]["search"]["hits"]
    exact = [h["id"] for h in hits if h["name"].lower() == ot_name.lower()]
    if not exact:
        raise ValueError(f"No exact Open Targets disease match for '{ot_name}'")
    return exact[0]


def fetch_ot_drugs(disease_id: str) -> dict:
    """Fetch drug-target rows for one disease; cached per disease (checkpoint per unit)."""
    cache = OT_CACHE / f"{disease_id}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    data = ot_post(OT_DRUG_QUERY, {"diseaseId": disease_id})
    cache.write_text(json.dumps(data))
    time.sleep(0.3)
    return data


def ot_rows(label: str, disease_id: str) -> tuple[list[dict], str]:
    data = fetch_ot_drugs(disease_id)
    version = data["data"]["meta"]["dataVersion"]
    d = data["data"]["disease"]
    out = []
    for row in (d.get("drugAndClinicalCandidates") or {}).get("rows", []):
        drug = row.get("drug") or {}
        stage = row.get("maxClinicalStage") or "UNKNOWN"
        for moa in ((drug.get("mechanismsOfAction") or {}).get("rows") or []):
            for t in moa.get("targets") or []:
                out.append({
                    "disease": label, "disease_id": disease_id,
                    "drug_id": drug.get("id", ""), "drug_type": drug.get("drugType", ""),
                    "stage": stage, "phase": PHASE_MAP.get(stage),
                    "action_type": moa.get("actionType", ""),
                    "ensembl": t.get("id", ""), "gene": t.get("approvedSymbol", ""),
                })
    return out, f"{version['year']}.{version['month']}"


def load_disease_gwas() -> dict:
    namespace: dict = {}
    exec(compile(DISEASE_GWAS_PY.read_text(), str(DISEASE_GWAS_PY), "exec"), namespace)
    return {"DISEASE_GWAS": namespace["DISEASE_GWAS"],
            "EFFECTIVE_N": namespace["EFFECTIVE_N"],
            "FLOOR": namespace["EFFECTIVE_N_FLOOR"]}


def load_v6_map() -> dict[str, str]:
    """Parse the frozen V6 gene -> subcategory tables (PRESPEC_V6.md section 5.4)."""
    genes: dict[str, str] = {}
    section = None
    in_54 = False
    for line in PRESPEC_V6.read_text().splitlines():
        if line.startswith("### 5.4"):
            in_54 = True
            continue
        if in_54 and line.startswith("### ") and not line.startswith("### 5.4"):
            break
        if not in_54:
            continue
        m = re.match(r"#### (A1|A2|A3|B1|B2)\.", line)
        if m:
            section = m.group(1)
            continue
        if line.startswith("#### "):
            section = None
            continue
        row = re.match(r"\| ([A-Z0-9]+) \|", line)
        if section and row and row.group(1) != "Gene":
            genes[row.group(1)] = section
    return genes


def load_v7_blind() -> dict[str, str]:
    genes = {}
    for line in V7_BLIND.read_text().splitlines():
        m = re.match(r"\| \*\*([A-Z0-9]+)\*\* \|.*\| \*\*([DCB])\*\* \|", line)
        if m:
            genes[m.group(1)] = m.group(2)
    return genes


def load_decode() -> dict:
    if not DECODE_XLSX.exists():
        r = requests.get(DECODE_URL, timeout=600)
        r.raise_for_status()
        DECODE_XLSX.write_bytes(r.content)
    wb = openpyxl.load_workbook(DECODE_XLSX, read_only=True)
    st1 = wb["ST01"].iter_rows(min_row=3, values_only=True)
    h1 = [str(c) for c in next(st1)]
    i_seq, i_gene, i_inc, i_ens = (h1.index("SeqId"), h1.index("Gene"),
                                   h1.index("Included in\nanalysis"), h1.index("Ensembl.Gene.ID"))
    seq_ens: dict[str, set] = {}
    measured_sym, measured_ens = set(), set()
    for r in st1:
        if r[i_seq] is None or str(r[i_inc]).strip() != "Yes":
            continue
        syms = {s for s in re.split(r"[ ,.;|]+", str(r[i_gene])) if s and s != "NA"}
        ens = {e for e in re.split(r"[ ,;|]+", str(r[i_ens])) if e.startswith("ENSG")}
        seq_ens[str(r[i_seq])] = ens
        measured_sym |= syms
        measured_ens |= ens
    st2 = wb["ST02"].iter_rows(min_row=3, values_only=True)
    h2 = [str(c) for c in next(st2)]
    i_g, i_s, i_ct = h2.index("gene\n (prot.)"), h2.index("SeqId"), h2.index("cis/\ntrans")
    cis_sym, cis_ens, n_cis_rows = set(), set(), 0
    for r in st2:
        if r[i_ct] != "cis":
            continue
        n_cis_rows += 1
        cis_sym |= {s for s in re.split(r"[ ,.;|]+", str(r[i_g])) if s and s != "NA"}
        cis_ens |= seq_ens.get(str(r[i_s]), set())
    return {"measured_sym": measured_sym, "measured_ens": measured_ens,
            "cis_sym": cis_sym, "cis_ens": cis_ens, "n_cis_rows": n_cis_rows,
            "n_aptamers_included": len(seq_ens)}


def read_token() -> str:
    for line in ENV2.read_text().splitlines():
        if line.startswith("OPEN_GWAS_TOKEN="):
            return line.split("=", 1)[1].strip().strip('"').strip("'")
    raise ValueError("OPEN_GWAS_TOKEN not found in .env2")


def load_gwasinfo(outcome_ids: list[str]) -> dict:
    """eqtl-a-* datasets and the outcome accessions, from one gwasinfo listing."""
    if GWASINFO_CACHE.exists():
        return json.loads(GWASINFO_CACHE.read_text())
    headers = {"Authorization": f"Bearer {read_token()}", "X-API-SOURCE": "exp08-feasibility/0.1"}
    r = requests.get(f"{OPENGWAS}/gwasinfo", headers=headers, timeout=600)
    result = {"status_code": r.status_code, "retrieved": datetime.now(timezone.utc).isoformat()}
    if r.status_code == 200:
        data = r.json()
        rows = data.values() if isinstance(data, dict) else data
        keep = ("id", "trait", "sample_size", "ncase", "ncontrol", "consortium", "author",
                "year", "nsnp", "build", "population")
        result["n_datasets_total"] = len(list(rows))
        rows = data.values() if isinstance(data, dict) else data
        result["eqtl_a"] = [{k: x.get(k) for k in keep} for x in rows
                            if str(x.get("id", "")).startswith("eqtl-a-")]
        rows = data.values() if isinstance(data, dict) else data
        result["outcomes"] = {x["id"]: {k: x.get(k) for k in keep} for x in rows
                              if x.get("id") in set(outcome_ids)}
        GWASINFO_CACHE.write_text(json.dumps(result))
    return result


def load_eqtlgen() -> pd.DataFrame:
    """Per-gene summary of the eQTLGen full-release cis-eQTL file (FDR < 0.05 rows)."""
    if EQTLGEN_SUMMARY.exists():
        return pd.read_csv(EQTLGEN_SUMMARY)
    if not EQTLGEN_FILE.exists():
        r = requests.get(EQTLGEN_URL, timeout=3600)
        r.raise_for_status()
        EQTLGEN_FILE.write_bytes(r.content)
    parts = []
    for chunk in pd.read_csv(EQTLGEN_FILE, sep="\t", chunksize=2_000_000,
                             usecols=["Pvalue", "Gene", "GeneSymbol", "NrSamples"]):
        parts.append(chunk.groupby(["Gene", "GeneSymbol"], as_index=False)
                     .agg(min_p=("Pvalue", "min"), max_n=("NrSamples", "max")))
    summary = (pd.concat(parts).groupby(["Gene", "GeneSymbol"], as_index=False)
               .agg(min_p=("min_p", "min"), max_n=("max_n", "max")))
    summary.to_csv(EQTLGEN_SUMMARY, index=False)
    return summary


def norm(s: str) -> str:
    return s.lower().replace("'", "").replace("-", " ").strip()


def main() -> None:
    OT_CACHE.mkdir(parents=True, exist_ok=True)
    run_ts = datetime.now(timezone.utc).isoformat()

    # ---- indications: the 24 in the frozen V3.4 candidate set --------------------------
    frozen = pd.read_csv(FROZEN_CANDIDATES,
                         usecols=["gene", "disease", "ot_disease_id", "mechanism_class"])
    indications = (frozen[["disease", "ot_disease_id"]].drop_duplicates()
                   .sort_values("disease").itertuples(index=False))
    disease_ids = {d: resolve_ot_id(name) for d, name in indications}

    rows, versions = [], set()
    for label, did in sorted(disease_ids.items()):
        r, v = ot_rows(label, did)
        rows.extend(r)
        versions.add(v)
    extra_rows = []
    for label, did in sorted(EXTRA_INDICATIONS.items()):
        r, v = ot_rows(label, did)
        extra_rows.extend(r)
        versions.add(v)

    ot = pd.DataFrame(rows)
    ot_extra = pd.DataFrame(extra_rows)
    unknown_stages = sorted(set(ot["stage"]) - set(PHASE_MAP))
    ot = ot[ot["gene"].astype(bool)]
    ot_extra = ot_extra[ot_extra["gene"].astype(bool)]

    def universes(df: pd.DataFrame) -> dict[str, pd.DataFrame]:
        mx = (df.dropna(subset=["phase"]).groupby(["gene", "disease"], as_index=False)
              .agg(max_phase=("phase", "max"), ensembl=("ensembl", "first")))
        return {"phase3": mx[mx["max_phase"] >= 3].reset_index(drop=True),
                "phase2_only": mx[(mx["max_phase"] >= 2) & (mx["max_phase"] < 3)]
                .reset_index(drop=True)}

    U = universes(ot)
    U_extra = universes(ot_extra)

    # ---- instrument sources ------------------------------------------------------------
    epi = pd.read_csv(EPIGRAPHDB_CATALOG, usecols=["gene_symbol"])
    epi_genes = {g for s in epi["gene_symbol"].dropna().astype(str) for g in s.split(";")}

    ukb_cis_raw = json.loads(UKBPPP_CIS.read_text())
    ukb_all_raw = json.loads(UKBPPP_ALL.read_text())
    ukb_single = {g for g in ukb_cis_raw if "_" not in g}
    ukb_multi_components = {p for g in ukb_cis_raw if "_" in g for p in g.split("_")}
    ukb_measured = {p for g in ukb_all_raw for p in g.split("_")}

    decode = load_decode()

    dg = load_disease_gwas()
    DISEASE_GWAS = dg["DISEASE_GWAS"]
    outcome_ids = sorted({v for v in DISEASE_GWAS.values() if v})
    gwasinfo = load_gwasinfo(outcome_ids)
    opengwas_eqtl_ens = {x["id"].replace("eqtl-a-", "") for x in gwasinfo.get("eqtl_a", [])}
    eqtl_n = [x.get("sample_size") for x in gwasinfo.get("eqtl_a", []) if x.get("sample_size")]
    eqg = load_eqtlgen()
    eqtl_fdr_ens = set(eqg["Gene"])
    eqtl_gws_ens = set(eqg.loc[eqg["min_p"] < GWS, "Gene"])

    # ---- seen pairs (identifier columns only) ------------------------------------------
    seen_v5 = pd.read_csv(CLASSIFICATION_V5, usecols=["gene", "disease"])
    seen_v51 = pd.read_csv(CLASSIFICATION_V5_1, usecols=["gene", "disease"])
    key = lambda g, d: (g, norm(d))
    seen_v5_keys = {key(g, d) for g, d in seen_v5.itertuples(index=False)}
    seen_v51_keys = {key(g, d) for g, d in seen_v51.itertuples(index=False)}
    frozen_keys = {key(g, d) for g, d in frozen[["gene", "disease"]].itertuples(index=False)}

    # ---- mechanism maps ------------------------------------------------------------------
    binary_map = (frozen.groupby("gene")["mechanism_class"]
                  .agg(lambda s: ";".join(sorted(set(s)))).to_dict())
    v6 = load_v6_map()
    v7_blind = load_v7_blind()
    v7_map = {g: V6_SECTION_TO_V7[s] for g, s in v6.items()} | v7_blind
    mixed_genes = {g for g, m in binary_map.items() if m == "mixed"}

    def annotate(df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["epigraphdb"] = df["gene"].isin(epi_genes)
        df["ukbppp"] = df["gene"].isin(ukb_single)
        df["ukbppp_multigene_assay_only"] = (~df["ukbppp"]) & df["gene"].isin(ukb_multi_components)
        df["decode"] = df["gene"].isin(decode["cis_sym"]) | df["ensembl"].isin(decode["cis_ens"])
        df["pqtl_union"] = df["epigraphdb"] | df["ukbppp"] | df["decode"]
        df["eqtlgen_fdr05"] = df["ensembl"].isin(eqtl_fdr_ens)
        df["eqtlgen"] = df["ensembl"].isin(eqtl_gws_ens)
        df["opengwas_eqtl_a"] = df["ensembl"].isin(opengwas_eqtl_ens)
        df["eqtl_only"] = df["eqtlgen"] & ~df["pqtl_union"]
        df["any_qtl"] = df["pqtl_union"] | df["eqtlgen"]
        df["gwas_id"] = df["disease"].map(lambda d: DISEASE_GWAS.get(d))
        df["has_outcome_gwas"] = df["gwas_id"].notna()
        eff = df["disease"].map(lambda d: dg["EFFECTIVE_N"].get(d))
        df["above_floor"] = df["has_outcome_gwas"] & (eff >= dg["FLOOR"])
        k = [key(g, d) for g, d in zip(df["gene"], df["disease"])]
        df["seen_v5"] = [x in seen_v5_keys for x in k]
        df["seen_v5_1"] = [x in seen_v51_keys for x in k]
        df["seen_any"] = df["seen_v5"] | df["seen_v5_1"]
        df["in_frozen_v34"] = [x in frozen_keys for x in k]
        df["has_binary_mech"] = df["gene"].isin(binary_map)
        df["has_v7_tier"] = df["gene"].isin(v7_map)
        return df

    sources = ["epigraphdb", "ukbppp", "decode", "pqtl_union", "eqtlgen_fdr05", "eqtlgen",
               "eqtl_only", "opengwas_eqtl_a", "any_qtl"]

    def summarise(df: pd.DataFrame) -> dict:
        out = {"n_genes": int(df["gene"].nunique()), "n_pairs": int(len(df)),
               "n_indications": int(df["disease"].nunique()),
               "ukbppp_multigene_assay_only": {
                   "genes": int(df.loc[df["ukbppp_multigene_assay_only"], "gene"].nunique()),
                   "pairs": int(df["ukbppp_multigene_assay_only"].sum())},
               "by_source": {}}
        for s in sources:
            sub = df[df[s]]
            gw = sub[sub["has_outcome_gwas"]]
            gwf = sub[sub["above_floor"]]
            new = gw[~gw["seen_any"]]
            genes_gw = set(gw["gene"])
            out["by_source"][s] = {
                "genes": int(sub["gene"].nunique()),
                "pairs": int(len(sub)),
                "pairs_with_outcome_gwas": int(len(gw)),
                "genes_with_outcome_gwas": int(len(genes_gw)),
                "pairs_with_outcome_gwas_above_neff_floor": int(len(gwf)),
                "pairs_with_outcome_gwas_seen_v5": int(gw["seen_v5"].sum()),
                "pairs_with_outcome_gwas_seen_v5_1": int(gw["seen_v5_1"].sum()),
                "pairs_with_outcome_gwas_seen_either": int(gw["seen_any"].sum()),
                "pairs_with_outcome_gwas_new": int(len(new)),
                "pairs_with_outcome_gwas_new_above_neff_floor": int(new["above_floor"].sum()),
                "new_pairs_whose_gene_lacks_binary_mechanism": int((~new["has_binary_mech"]).sum()),
                "new_pairs_whose_gene_lacks_v7_tier": int((~new["has_v7_tier"]
                                                           & ~new["gene"].isin(mixed_genes)).sum()),
                "genes_with_outcome_gwas_new_pairs": int(new["gene"].nunique()),
                "pairs_with_outcome_gwas_in_frozen_v34": int(gw["in_frozen_v34"].sum()),
                "genes_lacking_binary_mechanism": sorted(genes_gw - set(binary_map)),
                "n_genes_lacking_binary_mechanism": len(genes_gw - set(binary_map)),
                "genes_lacking_v7_tier": sorted(genes_gw - set(v7_map) - mixed_genes),
                "n_genes_lacking_v7_tier": len(genes_gw - set(v7_map) - mixed_genes),
                "pairs_whose_gene_lacks_v7_tier": int((~gw["has_v7_tier"]
                                                       & ~gw["gene"].isin(mixed_genes)).sum()),
                "pairs_by_indication": {k: int(v) for k, v in
                                        gw.groupby("disease").size().sort_index().items()},
            }
        return out

    A = {k: annotate(v) for k, v in U.items()}
    A_extra = {k: annotate(v) for k, v in U_extra.items()}

    combined = {}
    both = pd.concat([A["phase3"], A["phase2_only"]])
    for s in sources:
        g = set(both.loc[both[s] & both["has_outcome_gwas"], "gene"])
        combined[s] = {"genes_with_outcome_gwas": len(g),
                       "n_genes_lacking_binary_mechanism": len(g - set(binary_map)),
                       "n_genes_lacking_v7_tier": len(g - set(v7_map) - mixed_genes)}

    # Reproduction check against the July funnel (576 / 37 / 108 / 112).
    p3 = A["phase3"]
    funnel_check = {
        "july_funnel": json.loads(RESCUE_FUNNEL.read_text()),
        "rebuilt_phase3_genes": int(p3["gene"].nunique()),
        "rebuilt_epigraphdb_genes": int(p3.loc[p3["epigraphdb"], "gene"].nunique()),
        "rebuilt_ukbppp_genes": int(p3.loc[p3["ukbppp"], "gene"].nunique()),
        "rebuilt_epi_or_ukb_genes": int(p3.loc[p3["epigraphdb"] | p3["ukbppp"], "gene"].nunique()),
        "frozen_v34_pairs_present_in_rebuilt_phase3": int(p3["in_frozen_v34"].sum()),
        "frozen_v34_pairs_total": len(frozen_keys),
    }
    missing_frozen = sorted(f"{g}|{d}" for g, d in frozen[["gene", "disease"]].itertuples(index=False)
                            if key(g, d) not in set(zip(p3["gene"], p3["disease"].map(norm))))
    funnel_check["frozen_v34_pairs_absent_from_rebuilt_phase3"] = missing_frozen

    result = {
        "experiment": "08_mechanism_interaction feasibility: outcome-blind instrument coverage",
        "run_timestamp_utc": run_ts,
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "seed": None,
        "outcome_blindness": ("No outcome, adjudication, MR estimate or p-value column read. "
                              "Identifier columns only via usecols."),
        "definitions": {
            "universe": ("Open Targets drugAndClinicalCandidates for the 24 indications of "
                         "data/frozen_candidates_v34.csv (disease IDs resolved from its "
                         "ot_disease_id names), each drug's mechanismsOfAction targets."),
            "pair": "(target gene approvedSymbol, indication label as in frozen_candidates_v34)",
            "phase3": "max maxClinicalStage over drugs for the pair >= PHASE_3 (3, 3.5, 4)",
            "phase2_only": "max maxClinicalStage over drugs in {PHASE_2, PHASE_2_3}; no Phase III+ drug",
            "epigraphdb": "gene in gene_symbol of v34_mr_catalog.csv (INTERVAL cis, as used)",
            "ukbppp": "gene is a single-gene assay in ukbppp_cis_pqtl_genes.json (1,954 entries)",
            "decode": ("gene (symbol, or Ensembl ID via ST01 SeqId) has a cis row in "
                       "Ferkingstad 2021 ST02 (the paper's reported pQTL set; threshold as in the paper)"),
            "eqtlgen_fdr05": ("gene Ensembl ID has >=1 row in the eQTLGen full-release cis-eQTL "
                              "file (FDR < 0.05, 2019-12-11)"),
            "eqtlgen": "as eqtlgen_fdr05, with the gene's best cis-eQTL P < 5e-8",
            "eqtl_only": "eqtlgen true and pqtl_union false",
            "opengwas_eqtl_a": ("gene Ensembl ID has an eqtl-a-* dataset in the OpenGWAS gwasinfo "
                                "listing; zero if the listing request failed"),
            "any_qtl": "pqtl_union or eqtlgen",
            "has_outcome_gwas": "indication has a non-None accession in code/disease_gwas.py (post V5.1)",
            "above_neff_floor": "effective N >= EFFECTIVE_N_FLOOR in code/disease_gwas.py",
            "seen": "gene+indication key present in classification_v5.csv or classification_v5_1.csv",
            "binary_mechanism": "gene has mechanism_class in frozen_candidates_v34.csv",
            "v7_tier": ("gene in PRESPEC_V6 section 5.4 A1-B2 tables (mapped A1->D, A2/A3->C, "
                        "B1/B2->B) or in protocol/v7/gene_tier_assignments_blind.md; "
                        "mixed genes (APP, IGF1R) are not counted as lacking"),
        },
        "open_targets": {"data_versions": sorted(versions), "disease_ids": disease_ids,
                         "extra_disease_ids": EXTRA_INDICATIONS,
                         "unknown_stage_strings": unknown_stages,
                         "n_drug_target_rows": int(len(ot))},
        "source_sizes": {
            "epigraphdb_genes": len(epi_genes),
            "ukbppp_cis_entries": len(ukb_cis_raw),
            "ukbppp_cis_single_gene": len(ukb_single),
            "ukbppp_measured_genes": len(ukb_measured),
            "decode_aptamers_included": decode["n_aptamers_included"],
            "decode_measured_genes": len(decode["measured_sym"]),
            "decode_cis_pqtl_rows": decode["n_cis_rows"],
            "decode_cis_genes_symbol": len(decode["cis_sym"]),
            "eqtlgen_opengwas": {
                "gwasinfo_status": gwasinfo.get("status_code"),
                "gwasinfo_retrieved": gwasinfo.get("retrieved"),
                "n_datasets_total": gwasinfo.get("n_datasets_total"),
                "n_eqtl_a_datasets": len(gwasinfo.get("eqtl_a", [])),
                "sample_size_min": min(eqtl_n) if eqtl_n else None,
                "sample_size_max": max(eqtl_n) if eqtl_n else None,
                "example": (gwasinfo.get("eqtl_a") or [None])[0],
            },
            "eqtlgen_full_release": {
                "genes_fdr05": len(eqtl_fdr_ens),
                "genes_p_lt_5e-8": len(eqtl_gws_ens),
                "max_nrsamples": int(eqg["max_n"].max()),
                "median_gene_max_nrsamples": float(eqg["max_n"].median()),
            },
        },
        "outcome_gwas": {
            "mapping": DISEASE_GWAS,
            "effective_n": dg["EFFECTIVE_N"],
            "floor": dg["FLOOR"],
            "accessions_found_in_gwasinfo": sorted(gwasinfo.get("outcomes", {}).keys()),
            "accessions_missing_from_gwasinfo": sorted(set(outcome_ids)
                                                       - set(gwasinfo.get("outcomes", {}))),
        },
        "seen_sets": {"classification_v5_rows": len(seen_v5), "classification_v5_keys": len(seen_v5_keys),
                      "classification_v5_1_rows": len(seen_v51),
                      "classification_v5_1_keys": len(seen_v51_keys),
                      "union_keys": len(seen_v5_keys | seen_v51_keys)},
        "mechanism_maps": {"binary_genes": len(binary_map), "v6_genes": len(v6),
                           "v7_blind_genes": len(v7_blind), "v7_total_genes": len(v7_map),
                           "mixed_genes": sorted(mixed_genes)},
        "funnel_reproduction": funnel_check,
        "counts": {
            "phase3": summarise(A["phase3"]),
            "phase2_only": summarise(A["phase2_only"]),
            "phase3_or_phase2_genes_needing_assignment": combined,
        },
        "sensitivity_extra_six_indications": {
            "note": ("Indications in non_ukb_gwas_sources.json outside the frozen 24; none has "
                     "an accession in code/disease_gwas.py, so pairs_with_outcome_gwas is 0 "
                     "by construction."),
            "phase3": summarise(A_extra["phase3"]),
            "phase2_only": summarise(A_extra["phase2_only"]),
        },
        "inputs": {
            "frozen_candidates_v34": input_record(FROZEN_CANDIDATES, "usecols gene, disease, ot_disease_id, mechanism_class"),
            "epigraphdb_catalog": input_record(EPIGRAPHDB_CATALOG, "usecols gene_symbol"),
            "classification_v5": input_record(CLASSIFICATION_V5, "usecols gene, disease"),
            "classification_v5_1": input_record(CLASSIFICATION_V5_1, "usecols gene, disease"),
            "disease_gwas": input_record(DISEASE_GWAS_PY, "accessions and effective N"),
            "prespec_v6": input_record(PRESPEC_V6, "section 5.4 gene tables"),
            "v7_blind": input_record(V7_BLIND, "13 blind tier assignments"),
            "ukbppp_cis_genes": input_record(UKBPPP_CIS, "1,954 UKB-PPP cis-pQTL genes"),
            "ukbppp_all_genes": input_record(UKBPPP_ALL, "2,922 UKB-PPP measured assays"),
            "rescue_funnel": input_record(RESCUE_FUNNEL, "July funnel counts for reproduction"),
            "decode_supplementary": input_record(DECODE_XLSX, f"Ferkingstad 2021 Supplementary Tables, {DECODE_URL}"),
            "opengwas_gwasinfo_subset": (input_record(GWASINFO_CACHE, "eqtl-a-* rows and outcome accessions")
                                         if GWASINFO_CACHE.exists() else
                                         {"path": str(GWASINFO_CACHE), "sha256": None,
                                          "note": f"not written; gwasinfo returned HTTP {gwasinfo.get('status_code')}"}),
            "eqtlgen_cis_fdr05": input_record(EQTLGEN_FILE, f"eQTLGen full release, {EQTLGEN_URL}"),
            "eqtlgen_gene_summary": input_record(EQTLGEN_SUMMARY, "per-gene min P and max NrSamples derived from the file above"),
            "open_targets_cache": {p.name: sha256(p) for p in sorted(OT_CACHE.glob("*.json"))},
        },
    }

    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    for phase in ("phase3", "phase2_only"):
        c = result["counts"][phase]
        print(f"\n{phase}: {c['n_genes']} genes, {c['n_pairs']} pairs, {c['n_indications']} indications")
        print(f"  {'source':<12}{'genes':>7}{'pairs':>7}{'+gwas':>7}{'>=floor':>8}{'seen':>6}{'new':>6}"
              f"{'noBin':>7}{'noV7':>6}")
        for s, v in c["by_source"].items():
            print(f"  {s:<12}{v['genes']:>7}{v['pairs']:>7}{v['pairs_with_outcome_gwas']:>7}"
                  f"{v['pairs_with_outcome_gwas_above_neff_floor']:>8}"
                  f"{v['pairs_with_outcome_gwas_seen_either']:>6}{v['pairs_with_outcome_gwas_new']:>6}"
                  f"{v['n_genes_lacking_binary_mechanism']:>7}{v['n_genes_lacking_v7_tier']:>6}")
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
