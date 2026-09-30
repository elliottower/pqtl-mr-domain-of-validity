# /// script
# requires-python = ">=3.10,<3.14"
# dependencies = [
#     "pandas==2.2.3",
#     "pyarrow==21.0.0",
#     "requests==2.32.3",
#     "openpyxl==3.1.5",
# ]
# ///
"""Outcome-blind feasibility count for experiment 08, extended to every indication.

Unit: target gene x indication (Open Targets disease ID). A pair is eligible if at least one
drug acting on the gene reached Phase II or later for the indication (Open Targets 26.09
`clinical_indication`, `maxClinicalStage`). Phase is used for eligibility only and is
reported only as the overall II vs III+ marginal of the final held-out set.

Blinding:
- `clinical_indication` is read with columns drugId, diseaseId, maxClinicalStage only;
  `clinicalReportIds` and `clinical_report` are never read.
- `drug_molecule` is read with columns id, drugType, parentId only.
- No MR, colocalization or disease-association statistic is computed or read. From Sun et
  al. 2018 only sheet ST4 is read; ST14-ST16 are never opened. No Eldjarn 2023 sheet is read.
- `classification_v5*.csv` are read with columns gene, disease only.
- Per-pair phase is never written to disk or printed.

Usage (the OpenGWAS token is needed only if inputs/opengwas/gwasinfo_all.json is absent):
    cd ~/Documents/GitHub/pqtl-mr-domain-of-validity
    uv run experiments/08_mechanism_interaction/feasibility/v2_all_indications/count_all_indications.py
"""
import hashlib
import io
import json
import os
import re
import time
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import openpyxl
import pandas as pd
import pyarrow.parquet as pq
import requests

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
INPUTS = HERE / "inputs"
OUT_JSON = HERE / "coverage_all_indications.json"
OUT_PAIRS = HERE / "final_heldout_pairs.csv"
OUT_CANDIDATES = HERE / "outcome_gwas_candidates.csv"

FEAS_V1 = HERE.parent
TRANSPORT = Path.home() / "Documents" / "GitHub" / "transport-wrapper" / "DRUGS_EXPANDED"

# ---- inputs already on disk (as in count_coverage.py) -----------------------------------
EARLIER_JSON = FEAS_V1 / "coverage_counts.json"
CLASSIFICATION_V5 = REPO / "results" / "v5" / "classification_v5.csv"
CLASSIFICATION_V5_1 = REPO / "results" / "v5_1" / "classification_v5_1.csv"
EPIGRAPHDB_CATALOG = REPO / "zenodo_export" / "data" / "epigraphdb" / "v34_mr_catalog.csv"
UKBPPP_CIS = TRANSPORT / "pqtl_sources" / "ukbppp_cis_pqtl_genes.json"
OLINK_MAP = REPO / "planning" / "perplexity_v8_design" / "F_data_readmes" / "olink_protein_map_3k_v1.tsv"
DECODE_XLSX = FEAS_V1 / "inputs" / "ferkingstad2021_MOESM4_ESM.xlsx"

# ---- downloaded inputs -------------------------------------------------------------------
OT_RELEASE = "26.09"
OT_BASE = f"https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/{OT_RELEASE}/output"
OT_DIR = INPUTS / f"ot_{OT_RELEASE}"
OT_FILES = {
    "clinical_indication": "clinical_indication/00000000.parquet",
    "drug_mechanism_of_action": "drug_mechanism_of_action/00000000.parquet",
    "drug_molecule": "drug_molecule/00000000.parquet",
    "disease": "disease/00000000.parquet",
    "study": "study/part-00000-90d9a8f0-e00a-492e-b961-40b56a21dbd0-c000.zstd.parquet",
}
GWASCAT_DIR = INPUTS / "gwas_catalog"
# REST download endpoints (the FTP mirror serves the same files at ~25 KB/s). The release
# date is in the served filename and recorded in GWASCAT_RELEASE.
GWASCAT_BASE = "https://www.ebi.ac.uk/gwas/api/search/downloads"
GWASCAT_RELEASE = "r2026-09-13"
GWASCAT_ENDPOINTS = {f"gwas-catalog-v1.0.3.1-{p.replace('_', '-')}-{GWASCAT_RELEASE}.tsv": f"{p}/v1.0.3.1"
                     for p in ("studies", "ancestries", "unpublished_studies", "unpublished_ancestries")}
GWASCAT_FILES = list(GWASCAT_ENDPOINTS)
HGNC_FILE = INPUTS / "hgnc" / "hgnc_complete_set.txt"
HGNC_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HPA_FILE = INPUTS / "hpa" / "proteinatlas.tsv.zip"
HPA_URL = "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"
UNIPROT_DIR = INPUTS / "uniprot"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
UNIPROT_QUERIES = {"KW-0964": "(organism_id:9606) AND (reviewed:true) AND (keyword:KW-0964)",
                   "SL-0243": "(organism_id:9606) AND (reviewed:true) AND (cc_scl_term:SL-0243)"}
OPENGWAS_FILE = INPUTS / "opengwas" / "gwasinfo_all.json"
OPENGWAS_META = INPUTS / "opengwas" / "gwasinfo_all.retrieval.json"
OPENGWAS_URL = "https://api.opengwas.io/api/gwasinfo"
SUN2018_XLSX = INPUTS / "interval" / "sun2018_MOESM4_ESM.xlsx"
SUN2018_URL = ("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/"
               "MediaObjects/41586_2018_175_MOESM4_ESM.xlsx")
SUN2018_SHEET = "ST4 - pQTL summary"

# Identical to count_coverage.py / candidate_selection_v34.py.
PHASE_MAP = {
    "PRECLINICAL": -1, "IND": 0, "EARLY_PHASE_1": 0.5,
    "PHASE_1": 1, "PHASE_1_2": 1.5, "PHASE_2": 2,
    "PHASE_2_3": 2.5, "PHASE_3": 3, "PREAPPROVAL": 3.5,
    "APPROVAL": 4, "UNKNOWN": -1,
}
NEFF_FLOOR = 2000
KNOWN_DRUG_TYPES = {"Antibody", "Antibody drug conjugate", "Cell", "Enzyme", "Gene",
                    "Oligonucleotide", "Oligosaccharide", "Protein", "Small molecule",
                    "Vaccine component"}
UKB_PATTERN = re.compile(r"\bUKB\b|UKBB|UK ?Biobank|UKBiobank", re.IGNORECASE)
ICELAND_PATTERN = re.compile(r"deCODE|Iceland", re.IGNORECASE)
INTERVAL_PATTERN = re.compile(r"\bINTERVAL\b")


# ---- utilities ---------------------------------------------------------------------------

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def input_record(path: Path, note: str, url: str | None = None) -> dict:
    rec = {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size, "note": note}
    if url:
        rec["url"] = url
    return rec


def fetch(url: str, dest: Path, attempts: int = 30) -> None:
    """Download with HTTP Range resume until the local size equals Content-Length."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    total = int(requests.head(url, timeout=120, allow_redirects=True).headers.get("Content-Length", -1))
    for _ in range(attempts):
        have = dest.stat().st_size if dest.exists() else 0
        if total >= 0 and have == total:
            return
        if total >= 0 and have > total:
            raise RuntimeError(f"{dest} is larger than the remote file; delete it and rerun")
        headers = {"Range": f"bytes={have}-"} if have else {}
        try:
            with requests.get(url, headers=headers, stream=True, timeout=300) as r:
                if have and r.status_code != 206:
                    raise RuntimeError(f"server ignored Range for {url} (HTTP {r.status_code})")
                r.raise_for_status()
                with open(dest, "ab") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
        except requests.RequestException:
            time.sleep(10)
        if total < 0:
            return
    raise RuntimeError(f"download of {url} incomplete after {attempts} attempts")


def fetch_whole(url: str, dest: Path, expect_name: str, attempts: int = 5) -> None:
    """Download a generated file in one piece; keep it only if the size matches and the
    served filename carries the pinned release."""
    if dest.exists():
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    for _ in range(attempts):
        try:
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                served = r.headers.get("Content-Disposition", "")
                if expect_name not in served:
                    raise RuntimeError(f"{url} serves '{served}', expected {expect_name}; "
                                       "the catalog release changed")
                total = int(r.headers.get("Content-Length", -1))
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(1 << 20):
                        f.write(chunk)
            if total < 0 or tmp.stat().st_size == total:
                tmp.rename(dest)
                return
        except requests.RequestException:
            time.sleep(10)
    raise RuntimeError(f"download of {url} incomplete after {attempts} attempts")


def curie(uri: str) -> str:
    return uri.strip().rstrip("/").split("/")[-1].replace(":", "_")


def norm_label(s: str) -> str:
    s = str(s).lower().replace("'", "").replace("’", "")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return s.strip()


def norm_old(s: str) -> str:
    """Key normalization of count_coverage.py, used for the seen-set labels."""
    return s.lower().replace("'", "").replace("-", " ").strip()


def effective_n(ncase: float, ncontrol: float) -> float | None:
    if not ncase or not ncontrol or ncase <= 0 or ncontrol <= 0:
        return None
    return 4.0 / (1.0 / ncase + 1.0 / ncontrol)


def to_num(x) -> float | None:
    try:
        v = float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return None
    return v if v == v else None


# ---- step 1: universe --------------------------------------------------------------------

def download_all() -> dict:
    status = {}
    for name, rel in OT_FILES.items():
        fetch(f"{OT_BASE}/{rel}", OT_DIR / f"{name}.parquet")
    for fn, ep in GWASCAT_ENDPOINTS.items():
        fetch_whole(f"{GWASCAT_BASE}/{ep}", GWASCAT_DIR / fn, expect_name=fn)
    fetch(HGNC_URL, HGNC_FILE)
    fetch(HPA_URL, HPA_FILE)
    fetch(SUN2018_URL, SUN2018_XLSX)
    status["uniprot"] = {k: fetch_uniprot(k, q) for k, q in UNIPROT_QUERIES.items()}
    status["opengwas"] = fetch_opengwas()
    return status


def fetch_uniprot(tag: str, query: str) -> dict:
    """Paged UniProt search, checkpointed per page; complete when rows == X-Total-Results."""
    out = UNIPROT_DIR / f"secreted_{tag}.tsv"
    meta = UNIPROT_DIR / f"secreted_{tag}.meta.json"
    if out.exists() and meta.exists():
        return json.loads(meta.read_text())
    part = UNIPROT_DIR / f"secreted_{tag}.partial.jsonl"
    UNIPROT_DIR.mkdir(parents=True, exist_ok=True)
    pages = [json.loads(line) for line in part.read_text().splitlines()] if part.exists() else []
    url = pages[-1]["next"] if pages else None
    params = None if url else {"query": query, "fields": "accession,gene_primary",
                               "format": "tsv", "size": 500}
    url = url or UNIPROT_SEARCH
    headers: dict = {}
    while url:
        r = requests.get(url, params=params, timeout=300)
        r.raise_for_status()
        headers = dict(r.headers)
        lines = r.text.splitlines()
        nxt = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("Link", ""))
        page = {"rows": lines[1:], "header": lines[0], "next": nxt.group(1) if nxt else None,
                "total": int(r.headers["X-Total-Results"]),
                "release": r.headers.get("X-UniProt-Release"),
                "release_date": r.headers.get("X-UniProt-Release-Date")}
        with open(part, "a") as f:
            f.write(json.dumps(page) + "\n")
        pages.append(page)
        url, params = page["next"], None
    rows = [row for p in pages for row in p["rows"]]
    if len(rows) != pages[-1]["total"]:
        raise RuntimeError(f"UniProt {tag}: {len(rows)} rows, expected {pages[-1]['total']}")
    out.write_text("\n".join([pages[0]["header"]] + rows) + "\n")
    info = {"query": query, "rows": len(rows), "release": pages[-1]["release"],
            "release_date": pages[-1]["release_date"],
            "retrieved": datetime.now(timezone.utc).isoformat()}
    meta.write_text(json.dumps(info, indent=2))
    part.unlink()
    return info


def fetch_opengwas() -> dict:
    if OPENGWAS_FILE.exists():
        return json.loads(OPENGWAS_META.read_text()) if OPENGWAS_META.exists() else {"status": "cached"}
    token = os.environ.get("OPENGWAS", "")
    if not token:
        return {"status": "not retrieved: OPENGWAS token absent from environment"}
    r = requests.get(OPENGWAS_URL, timeout=900,
                     headers={"Authorization": f"Bearer {token}",
                              "X-API-SOURCE": "exp08-feasibility-v2/0.1"})
    info = {"url": OPENGWAS_URL, "method": "GET", "status": r.status_code,
            "retrieved": datetime.now(timezone.utc).isoformat()}
    if r.status_code == 200:
        OPENGWAS_FILE.parent.mkdir(parents=True, exist_ok=True)
        OPENGWAS_FILE.write_bytes(r.content)
        OPENGWAS_META.write_text(json.dumps(info))
    return info


def load_symbol_maps() -> tuple[dict, dict, dict]:
    hg = pd.read_csv(HGNC_FILE, sep="\t", dtype=str,
                     usecols=["symbol", "status", "ensembl_gene_id", "uniprot_ids"])
    hg = hg[hg["status"] == "Approved"]
    ens2sym = dict(zip(hg["ensembl_gene_id"].dropna(), hg.loc[hg["ensembl_gene_id"].notna(), "symbol"]))
    sym2ens = {s: e for e, s in ens2sym.items()}
    uni2sym: dict[str, set] = {}
    for s, u in zip(hg["symbol"], hg["uniprot_ids"].fillna("")):
        for acc in u.split("|"):
            if acc:
                uni2sym.setdefault(acc, set()).add(s)
    return ens2sym, sym2ens, uni2sym


def build_universe(ens2sym: dict) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ci = pq.read_table(OT_DIR / "clinical_indication.parquet",
                       columns=["drugId", "diseaseId", "maxClinicalStage"]).to_pandas()
    unknown_stages = sorted(set(ci["maxClinicalStage"].dropna()) - set(PHASE_MAP))
    ci["phase"] = ci["maxClinicalStage"].map(PHASE_MAP)

    moa = pq.read_table(OT_DIR / "drug_mechanism_of_action.parquet",
                        columns=["actionType", "chemblIds", "targets", "targetType"]).to_pandas()
    moa = moa.explode("chemblIds").explode("targets").dropna(subset=["chemblIds", "targets"])
    moa = moa.rename(columns={"chemblIds": "moa_drug", "targets": "ensembl"})

    mol = pq.read_table(OT_DIR / "drug_molecule.parquet",
                        columns=["id", "drugType", "parentId"]).to_pandas()
    parent = dict(zip(mol["id"], mol["parentId"]))
    dtype = dict(zip(mol["id"], mol["drugType"]))

    # A drug's mechanisms: its own rows; if none, its parent's rows.
    own = set(moa["moa_drug"])
    ci["moa_drug"] = [d if d in own else (parent.get(d) if parent.get(d) in own else None)
                      for d in ci["drugId"]]
    n_ci_no_moa = int(ci["moa_drug"].isna().sum())
    dt = ci.dropna(subset=["moa_drug"]).merge(moa, on="moa_drug", how="inner")
    dt["gene"] = dt["ensembl"].map(ens2sym)
    n_unmapped_ens = int(dt.loc[dt["gene"].isna(), "ensembl"].nunique())
    dt = dt.dropna(subset=["gene"])
    dt["drug_type"] = dt["drugId"].map(dtype).fillna(dt["moa_drug"].map(dtype))
    dt["single_protein"] = dt["targetType"] == "single protein"

    elig_rows = dt[dt["phase"] >= 2]
    pairs = (elig_rows.groupby(["gene", "diseaseId"], as_index=False)
             .agg(ensembl=("ensembl", "first"), max_phase=("phase", "max"),
                  any_single_protein=("single_protein", "max")))
    meta = {"clinical_indication_rows": int(len(ci)),
            "clinical_indication_rows_without_mechanism": n_ci_no_moa,
            "drug_target_indication_rows": int(len(dt)),
            "unmapped_target_ensembl_ids": n_unmapped_ens,
            "unknown_stage_strings": unknown_stages}
    return pairs, elig_rows, meta


# ---- step 2: instruments -----------------------------------------------------------------

def load_decode() -> dict:
    """Ferkingstad 2021 ST01/ST02, as in count_coverage.py."""
    wb = openpyxl.load_workbook(DECODE_XLSX, read_only=True)
    st1 = wb["ST01"].iter_rows(min_row=3, values_only=True)
    h1 = [str(c) for c in next(st1)]
    i_seq, i_inc, i_ens = h1.index("SeqId"), h1.index("Included in\nanalysis"), h1.index("Ensembl.Gene.ID")
    seq_ens: dict[str, set] = {}
    for r in st1:
        if r[i_seq] is None or str(r[i_inc]).strip() != "Yes":
            continue
        seq_ens[str(r[i_seq])] = {e for e in re.split(r"[ ,;|]+", str(r[i_ens])) if e.startswith("ENSG")}
    st2 = wb["ST02"].iter_rows(min_row=3, values_only=True)
    h2 = [str(c) for c in next(st2)]
    i_g, i_s, i_ct = h2.index("gene\n (prot.)"), h2.index("SeqId"), h2.index("cis/\ntrans")
    cis_sym, cis_ens = set(), set()
    for r in st2:
        if r[i_ct] != "cis":
            continue
        cis_sym |= {s for s in re.split(r"[ ,.;|]+", str(r[i_g])) if s and s != "NA"}
        cis_ens |= seq_ens.get(str(r[i_s]), set())
    return {"sym": cis_sym, "ens": cis_ens}


def load_interval_st4(uni2sym: dict) -> dict:
    """Sun 2018 ST4 cis rows on single-UniProt SOMAmers (study-wide threshold of the paper)."""
    ws = openpyxl.load_workbook(SUN2018_XLSX, read_only=True)[SUN2018_SHEET]
    rows = ws.iter_rows(min_row=5, values_only=True)
    header = [str(c) if c is not None else "" for c in next(rows)]
    i_uni, i_ct = header.index("UniProt"), header.index("cis/ trans")
    syms, n_cis, n_multi, n_unmapped = set(), 0, 0, 0
    for r in rows:
        if r[i_ct] is None or str(r[i_ct]).strip() != "cis":
            continue
        n_cis += 1
        accs = [a for a in re.split(r"[ ,;|]+", str(r[i_uni] or "")) if a]
        if len(accs) != 1:
            n_multi += 1
            continue
        s = uni2sym.get(accs[0])
        if not s:
            n_unmapped += 1
            continue
        syms |= s
    return {"sym": syms, "n_cis_rows": n_cis, "n_multi_uniprot_rows": n_multi,
            "n_unmapped_rows": n_unmapped}


def load_instruments(uni2sym: dict) -> dict:
    ukb_cis_raw = json.loads(UKBPPP_CIS.read_text())
    ukb_single = {g for g in ukb_cis_raw if "_" not in g}
    olink = pd.read_csv(OLINK_MAP, sep="\t", usecols=["HGNC.symbol", "ensembl_id"], dtype=str)
    ukb_ens = set(olink.loc[olink["HGNC.symbol"].isin(ukb_single), "ensembl_id"].dropna())
    epi = pd.read_csv(EPIGRAPHDB_CATALOG, usecols=["gene_symbol"])
    epi_genes = {g for s in epi["gene_symbol"].dropna().astype(str) for g in s.split(";")}
    return {"ukbppp": {"sym": ukb_single, "ens": ukb_ens},
            "decode": load_decode(),
            "interval_st4": load_interval_st4(uni2sym),
            "interval_epigraphdb": {"sym": epi_genes, "ens": set()}}


# ---- step 3: outcome GWAS candidates -------------------------------------------------------

def load_disease() -> pd.DataFrame:
    return pq.read_table(OT_DIR / "disease.parquet",
                         columns=["id", "name", "exactSynonyms", "parents", "children",
                                  "ancestors", "descendants", "obsoleteTerms"]).to_pandas()


def label_index(dis: pd.DataFrame) -> dict[str, set]:
    idx: dict[str, set] = {}
    for did, name, syn in zip(dis["id"], dis["name"], dis["exactSynonyms"]):
        for lab in [name, *(list(syn) if syn is not None else [])]:
            idx.setdefault(norm_label(lab), set()).add(did)
    return idx


def parse_cases_controls(text: str) -> tuple[float, float]:
    ncase = ncontrol = 0.0
    for num, word in re.findall(r"(\d[\d,]*)\s[^\d]*?\b(cases?|controls?)\b", str(text)):
        v = float(num.replace(",", ""))
        if word.startswith("case"):
            ncase += v
        else:
            ncontrol += v
    return ncase, ncontrol


def overlap_flag(pattern: re.Pattern, text: str, cohort_listed: bool) -> str:
    if pattern.search(text):
        return "yes"
    return "no_per_cohort_list" if cohort_listed else "unknown"


def gwas_catalog_candidates() -> tuple[pd.DataFrame, dict]:
    frames = []
    for studies_fn, anc_fn, kind in [(GWASCAT_FILES[0], GWASCAT_FILES[1], "published"),
                                     (GWASCAT_FILES[2], GWASCAT_FILES[3], "unpublished")]:
        st = pd.read_csv(GWASCAT_DIR / studies_fn, sep="\t", dtype=str, quoting=3,
                         on_bad_lines="skip")
        an = pd.read_csv(GWASCAT_DIR / anc_fn, sep="\t", dtype=str, quoting=3, on_bad_lines="skip")
        st.columns = [c.strip().upper() for c in st.columns]
        an.columns = [c.strip().upper() for c in an.columns]
        frames.append((st, an, kind))
    rows, raw_counts = [], {}
    for st, an, kind in frames:
        raw_counts[kind] = int(len(st))
        acc_col = "STUDY ACCESSION"
        uri_col = "MAPPED_TRAIT_URI" if "MAPPED_TRAIT_URI" in st.columns else next(
            c for c in st.columns if "TRAIT" in c and "URI" in c and "BACKGROUND" not in c)
        full_col = next((c for c in st.columns if c.startswith("FULL SUMMARY")), None)
        loc_col = next((c for c in st.columns if "SUMMARY STATS LOCATION" in c), None)
        size_col = next((c for c in st.columns if c.startswith("INITIAL SAMPLE")), None)
        trait_col = next(c for c in st.columns if c in ("DISEASE/TRAIT", "TRAIT"))
        cohort_col = "COHORT" if "COHORT" in st.columns else None
        an_init = an[an["STAGE"].str.lower().eq("initial")] if "STAGE" in an.columns else an
        anc = an_init.groupby("STUDY ACCESSION").agg(
            ancestries=("BROAD ANCESTRAL CATEGORY", lambda s: "|".join(sorted(set(s.dropna())))),
            countries=("COUNTRY OF RECRUITMENT", lambda s: "|".join(sorted(set(s.dropna())))),
            ncase_anc=("NUMBER OF CASES", lambda s: sum(to_num(x) or 0 for x in s)),
            ncontrol_anc=("NUMBER OF CONTROLS", lambda s: sum(to_num(x) or 0 for x in s)),
        )
        for r in st.itertuples(index=False):
            rec = dict(zip(st.columns, r))
            acc = rec.get(acc_col)
            ids = sorted({curie(u) for u in str(rec.get(uri_col) or "").split(",") if u.strip()
                          and u.strip().lower() != "nan"})
            a = anc.loc[acc] if acc in anc.index else None
            ancestries = a["ancestries"] if a is not None else ""
            ncase = a["ncase_anc"] if a is not None else 0
            ncontrol = a["ncontrol_anc"] if a is not None else 0
            if not (ncase and ncontrol):
                ncase, ncontrol = parse_cases_controls(rec.get(size_col, ""))
            full = str(rec.get(full_col, "")).strip().lower() == "yes" if full_col else True
            if kind == "unpublished":
                full = full or bool(str(rec.get(loc_col, "") or "").strip() not in ("", "nan"))
            cohort = str(rec.get(cohort_col, "") or "") if cohort_col else ""
            cohort = "" if cohort.lower() == "nan" else cohort
            text = " ".join([cohort, str(rec.get(size_col, "")), a["countries"] if a is not None else ""])
            rows.append({
                "source": "gwas_catalog", "accession": acc, "trait": rec.get(trait_col),
                "mapped_ids": ids, "mapping_method": "gwas_catalog_mapped_trait_uri",
                "n_mapped_ids": len(ids), "ncase": ncase, "ncontrol": ncontrol,
                "neff": effective_n(ncase, ncontrol),
                "european": ancestries == "European",
                "ancestries": ancestries, "full_sumstats": full,
                "case_control": bool(ncase and ncontrol),
                "includes_ukb": overlap_flag(UKB_PATTERN, text, bool(cohort)),
                "includes_iceland": overlap_flag(ICELAND_PATTERN, text, bool(cohort)),
                "includes_interval": overlap_flag(INTERVAL_PATTERN, text, bool(cohort)),
                "cohorts": cohort,
            })
    return pd.DataFrame(rows), raw_counts


def finngen_candidates() -> tuple[pd.DataFrame, dict]:
    cols = [c for c in ["studyId", "projectId", "studyType", "traitFromSource",
                        "traitFromSourceMappedIds", "diseaseIds", "nCases", "nControls",
                        "cohorts", "hasSumstats", "ldPopulationStructure"]
            if c in pq.read_schema(OT_DIR / "study.parquet").names]
    st = pq.read_table(OT_DIR / "study.parquet", columns=cols).to_pandas()
    fg = st[st["studyId"].str.upper().str.startswith("FINNGEN")].copy()
    releases = sorted(set(fg["projectId"].dropna())) if "projectId" in fg.columns else []
    rows = []
    for r in fg.itertuples(index=False):
        rec = dict(zip(fg.columns, r))
        ids = set()
        for c in ("diseaseIds", "traitFromSourceMappedIds"):
            v = rec.get(c)
            if v is not None and not (isinstance(v, float)):
                ids |= {curie(x) for x in v}
        ncase, ncontrol = to_num(rec.get("nCases")) or 0, to_num(rec.get("nControls")) or 0
        rows.append({
            "source": "finngen", "accession": rec["studyId"], "trait": rec.get("traitFromSource"),
            "mapped_ids": sorted(ids), "mapping_method": "open_targets_study_index_finngen",
            "n_mapped_ids": len(ids), "ncase": ncase, "ncontrol": ncontrol,
            "neff": effective_n(ncase, ncontrol), "european": True, "ancestries": "Finnish",
            "full_sumstats": bool(rec.get("hasSumstats", True)),
            "case_control": bool(ncase and ncontrol),
            "includes_ukb": "no_per_cohort_list", "includes_iceland": "no_per_cohort_list",
            "includes_interval": "no_per_cohort_list", "cohorts": "FinnGen",
        })
    return pd.DataFrame(rows), {"finngen_studies_in_ot_study_index": int(len(fg)),
                                "finngen_projects": releases,
                                "ot_study_columns_read": cols}


def opengwas_candidates(gcat: pd.DataFrame, fgn: pd.DataFrame, labels: dict) -> tuple[pd.DataFrame, dict]:
    if not OPENGWAS_FILE.exists():
        return pd.DataFrame(), {"status": "gwasinfo listing not available"}
    info = json.loads(OPENGWAS_FILE.read_text())
    gcat_ids = dict(zip(gcat["accession"], gcat["mapped_ids"]))
    fg_code = {}
    for acc, ids in zip(fgn["accession"], fgn["mapped_ids"]):
        code = re.sub(r"^FINNGEN_R\d+_", "", acc, flags=re.IGNORECASE)
        fg_code.setdefault(code.upper(), set()).update(ids)
    skip_prefix = ("eqtl-a-", "prot-", "met-", "ubm-")
    rows, methods = [], Counter()
    for gid, x in info.items():
        if gid.startswith(skip_prefix):
            continue
        ids, method = set(), "unmapped"
        if gid.startswith("ebi-a-") and gid[6:] in gcat_ids:
            ids, method = set(gcat_ids[gid[6:]]), "gwas_catalog_accession"
        elif gid.startswith("finn-b-") and gid[7:].upper() in fg_code:
            ids, method = fg_code[gid[7:].upper()], "open_targets_finngen_endpoint_code"
        ont = str(x.get("ontology") or "")
        if not ids and ont not in ("", "NA", "None", "nan"):
            ids, method = {curie(o) for o in re.split(r"[;,| ]+", ont) if o}, "opengwas_ontology_field"
        if not ids and norm_label(x.get("trait", "")) in labels:
            ids, method = set(labels[norm_label(x["trait"])]), "exact_label_or_synonym"
        methods[method] += 1
        ncase, ncontrol = to_num(x.get("ncase")) or 0, to_num(x.get("ncontrol")) or 0
        text = " ".join(str(x.get(k, "")) for k in ("consortium", "author", "note"))
        ukb = "yes" if (gid.startswith("ukb-") or UKB_PATTERN.search(text)) else (
            "no_per_cohort_list" if gid.startswith("finn-b-") else "unknown")
        rows.append({
            "source": "opengwas", "accession": gid, "trait": x.get("trait"),
            "mapped_ids": sorted(ids), "mapping_method": method, "n_mapped_ids": len(ids),
            "ncase": ncase, "ncontrol": ncontrol, "neff": effective_n(ncase, ncontrol),
            "european": str(x.get("population")) == "European",
            "ancestries": str(x.get("population")), "full_sumstats": True,
            "case_control": bool(ncase and ncontrol),
            "includes_ukb": ukb,
            "includes_iceland": "yes" if ICELAND_PATTERN.search(text) else (
                "no_per_cohort_list" if gid.startswith("finn-b-") else "unknown"),
            "includes_interval": "yes" if INTERVAL_PATTERN.search(text) else (
                "no_per_cohort_list" if gid.startswith("finn-b-") else "unknown"),
            "cohorts": str(x.get("consortium", "")),
        })
    return pd.DataFrame(rows), {"datasets_total": len(info),
                                "datasets_considered": len(rows),
                                "mapping_methods": dict(methods)}


def indication_tiers(indications: set, cands: pd.DataFrame, dis: pd.DataFrame) -> pd.DataFrame:
    """Per indication, qualifying candidates by tier (1 exact, 2 descendant, 3 ancestor)."""
    q = cands[cands["qualifies"]]
    by_id: dict[str, list[int]] = {}
    for i, ids in zip(q.index, q["mapped_ids"]):
        for d in ids:
            by_id.setdefault(d, []).append(i)
    rel = dis.set_index("id")
    out = []
    for d in sorted(indications):
        if d in rel.index:
            r = rel.loc[d]
            desc = set(r["descendants"]) if r["descendants"] is not None else set()
            anc = set(r["ancestors"]) if r["ancestors"] is not None else set()
            par = set(r["parents"]) if r["parents"] is not None else set()
            chi = set(r["children"]) if r["children"] is not None else set()
            obs = set(r["obsoleteTerms"]) if r["obsoleteTerms"] is not None else set()
        else:
            desc = anc = par = chi = obs = set()
        exact_ids = {d} | {curie(o) for o in obs}
        t1 = sorted({i for x in exact_ids for i in by_id.get(x, [])})
        t2 = sorted({i for x in desc for i in by_id.get(x, [])} - set(t1))
        t2c = sorted({i for x in chi for i in by_id.get(x, [])} - set(t1))
        t3 = sorted({i for x in anc for i in by_id.get(x, [])} - set(t1))
        t3p = sorted({i for x in par for i in by_id.get(x, [])} - set(t1))
        out.append({"diseaseId": d, "in_ot_disease_index": d in rel.index,
                    "tier1": t1, "tier2_descendant": t2, "tier2_child": t2c,
                    "tier3_ancestor": t3, "tier3_parent": t3p})
    return pd.DataFrame(out).set_index("diseaseId")


def no_known_overlap(cand_idx: list[int], cands: pd.DataFrame, sources: list[str]) -> bool:
    """True if some (instrument source, candidate) combination has no known sample overlap."""
    flag = {"ukbppp": "includes_ukb", "decode": "includes_iceland",
            "interval": "includes_interval"}
    for s in sources:
        col = flag[s]
        for i in cand_idx:
            if cands.at[i, col] != "yes":
                return True
    return False


def distribution(counts: pd.Series) -> dict:
    if counts.empty:
        return {}
    bins = pd.cut(counts, [0, 1, 2, 5, 10, 20, 10**6], labels=["1", "2", "3-5", "6-10", "11-20", ">20"])
    return {"n_genes": int(len(counts)), "min": int(counts.min()), "median": float(counts.median()),
            "mean": round(float(counts.mean()), 3), "p90": float(counts.quantile(0.9)),
            "max": int(counts.max()),
            "histogram": {str(k): int(v) for k, v in bins.value_counts(sort=False).items()}}


def main() -> None:
    run_ts = datetime.now(timezone.utc).isoformat()
    fetch_status = download_all()

    ens2sym, sym2ens, uni2sym = load_symbol_maps()
    pairs, elig_rows, uni_meta = build_universe(ens2sym)
    dis = load_disease()
    dis_name = dict(zip(dis["id"], dis["name"]))
    ont = pq.read_table(OT_DIR / "disease.parquet", columns=["id", "ontology"]).to_pandas()
    ta_flag = {i: bool(o["isTherapeuticArea"]) for i, o in zip(ont["id"], ont["ontology"])}
    leaf_flag = {i: bool(o["leaf"]) for i, o in zip(ont["id"], ont["ontology"])}

    # ---- step 2: instrument flags ---------------------------------------------------------
    inst = load_instruments(uni2sym)
    for s in ("ukbppp", "decode", "interval_st4", "interval_epigraphdb"):
        pairs[s] = pairs["gene"].isin(inst[s]["sym"]) | pairs["ensembl"].isin(inst[s].get("ens", set()))
    pairs["interval"] = pairs["interval_st4"] | pairs["interval_epigraphdb"]
    pairs["pqtl"] = pairs["ukbppp"] | pairs["decode"] | pairs["interval"]

    def sources_of(r) -> list[str]:
        return [s for s in ("ukbppp", "decode", "interval") if r[s]]

    # ---- step 3: outcome GWAS candidates ---------------------------------------------------
    labels = label_index(dis)
    gcat, gcat_raw = gwas_catalog_candidates()
    fgn, fg_meta = finngen_candidates()
    ogw, og_meta = opengwas_candidates(gcat, fgn, labels)
    cands = pd.concat([gcat, fgn, ogw], ignore_index=True)
    cands["qualifies"] = (cands["case_control"] & cands["european"] & cands["full_sumstats"]
                          & (cands["neff"].fillna(0) >= NEFF_FLOOR))

    indications = set(pairs["diseaseId"])
    tiers = indication_tiers(indications, cands, dis)
    pairs["tier1_idx"] = pairs["diseaseId"].map(tiers["tier1"])
    pairs["has_tier1"] = pairs["tier1_idx"].map(len) > 0
    has_t2 = pairs["diseaseId"].map(tiers["tier2_descendant"]).map(len) > 0
    has_t3 = pairs["diseaseId"].map(tiers["tier3_ancestor"]).map(len) > 0
    has_t2c = pairs["diseaseId"].map(tiers["tier2_child"]).map(len) > 0
    has_t3p = pairs["diseaseId"].map(tiers["tier3_parent"]).map(len) > 0
    pairs["t1_no_overlap"] = [bool(r["has_tier1"]) and no_known_overlap(r["tier1_idx"], cands, sources_of(r))
                              for _, r in pairs.iterrows()]

    # ---- held-out keys ----------------------------------------------------------------------
    earlier = json.loads(EARLIER_JSON.read_text())
    label2id = {norm_old(k): v for k, v in earlier["open_targets"]["disease_ids"].items()}
    seen = pd.concat([pd.read_csv(CLASSIFICATION_V5, usecols=["gene", "disease"]),
                      pd.read_csv(CLASSIFICATION_V5_1, usecols=["gene", "disease"])])
    seen["diseaseId"] = seen["disease"].map(lambda d: label2id.get(norm_old(d)))
    unmapped_seen = sorted(set(seen.loc[seen["diseaseId"].isna(), "disease"]))
    seen_keys = set(zip(seen["gene"], seen["diseaseId"]))
    pairs["heldout"] = [(g, d) not in seen_keys for g, d in zip(pairs["gene"], pairs["diseaseId"])]
    # Conservative variant: also drop pairs whose indication is an ancestor or descendant of a
    # seen indication for the same gene.
    rel = dis.set_index("id")
    related: dict[str, set] = {}
    for d in set(seen["diseaseId"].dropna()):
        if d in rel.index:
            r = rel.loc[d]
            related[d] = ({d} | set(r["descendants"] if r["descendants"] is not None else [])
                          | set(r["ancestors"] if r["ancestors"] is not None else []))
        else:
            related[d] = {d}
    seen_gene_related = {}
    for g, d in seen_keys:
        seen_gene_related.setdefault(g, set()).update(related.get(d, {d}))
    pairs["heldout_strict"] = [d not in seen_gene_related.get(g, set())
                               for g, d in zip(pairs["gene"], pairs["diseaseId"])]

    # ---- reproduction of the earlier 24-indication count (sum over phases only) ------------
    ids24 = set(earlier["open_targets"]["disease_ids"].values())
    p24 = pairs[pairs["diseaseId"].isin(ids24)]
    earlier_total = earlier["counts"]["phase3"]["n_pairs"] + earlier["counts"]["phase2_only"]["n_pairs"]
    repro = {"earlier_24_indication_eligible_pairs": earlier_total,
             "bulk_24_indication_eligible_pairs": int(len(p24)),
             "earlier_24_pqtl_union_pairs": (earlier["counts"]["phase3"]["by_source"]["pqtl_union"]["pairs"]
                                             + earlier["counts"]["phase2_only"]["by_source"]["pqtl_union"]["pairs"]),
             "bulk_24_pqtl_pairs_same_sources_as_earlier": int(
                 (p24["ukbppp"] | p24["decode"] | p24["interval_epigraphdb"]).sum())}

    # ---- funnel ------------------------------------------------------------------------------
    P = pairs
    inst_ = P[P["pqtl"]]
    t1 = inst_[inst_["has_tier1"]]
    t1no = t1[t1["t1_no_overlap"]]
    final = t1no[t1no["heldout"]]
    final_strict = t1no[t1no["heldout_strict"]]
    ind_all = pd.Series(sorted(indications))
    ind_t1 = set(tiers.index[tiers["tier1"].map(len) > 0])
    ind_t2 = set(tiers.index[tiers["tier2_descendant"].map(len) > 0])
    ind_t3 = set(tiers.index[tiers["tier3_ancestor"].map(len) > 0])
    ind_t2c = set(tiers.index[tiers["tier2_child"].map(len) > 0])
    ind_t3p = set(tiers.index[tiers["tier3_parent"].map(len) > 0])

    def block(df: pd.DataFrame) -> dict:
        return {"pairs": int(len(df)), "genes": int(df["gene"].nunique()),
                "indications": int(df["diseaseId"].nunique())}

    # Final set, looser candidate definitions (no phase).
    def final_with(mask_candidate: pd.Series) -> dict:
        d = P[P["pqtl"] & mask_candidate & P["heldout"]]
        return block(d)

    # ---- localization (genes only) ------------------------------------------------------------
    hpa = pd.read_csv(io.BytesIO(zipfile.ZipFile(HPA_FILE).read("proteinatlas.tsv")), sep="\t",
                      usecols=["Gene", "Ensembl", "Secretome location"], dtype=str)
    hpa_blood = set(hpa.loc[hpa["Secretome location"].fillna("").str.contains("Secreted to blood"), "Ensembl"])
    hpa_any_secreted = set(hpa.loc[hpa["Secretome location"].notna(), "Ensembl"])
    uni_sec: dict[str, set] = {}
    for tag in UNIPROT_QUERIES:
        t = pd.read_csv(UNIPROT_DIR / f"secreted_{tag}.tsv", sep="\t", dtype=str)
        uni_sec[tag] = set(t["Gene Names (primary)"].dropna())
    uni_secreted = uni_sec["KW-0964"] | uni_sec["SL-0243"]

    def localization(genes: set) -> dict:
        ens = {g: sym2ens.get(g) for g in genes}
        blood = {g for g in genes if ens[g] in hpa_blood}
        uni = genes & uni_secreted
        return {"n_genes": len(genes), "hpa_secreted_to_blood": len(blood),
                "hpa_any_secretome_location": len({g for g in genes if ens[g] in hpa_any_secreted}),
                "uniprot_secreted": len(uni), "uniprot_kw0964": len(genes & uni_sec["KW-0964"]),
                "uniprot_sl0243": len(genes & uni_sec["SL-0243"]),
                "both": len(blood & uni), "hpa_blood_only": len(blood - uni),
                "uniprot_only": len(uni - blood), "neither": len(genes - blood - uni)}

    # ---- mechanism-rule input availability (final set; totals only) ---------------------------
    fk = set(zip(final["gene"], final["diseaseId"]))
    fr = elig_rows[[k in fk for k in zip(elig_rows["gene"], elig_rows["diseaseId"])]]
    fr = fr.assign(has_type=fr["drug_type"].isin(KNOWN_DRUG_TYPES), has_action=fr["actionType"].notna())
    per_pair = fr.groupby(["gene", "diseaseId"]).agg(
        all_inputs=("has_type", lambda s: bool(s.all())), any_type=("has_type", "max"),
        all_action=("has_action", "min"))
    per_pair["all_inputs"] = per_pair["all_inputs"] & per_pair["all_action"]
    genes_all = per_pair.groupby(level="gene")["all_inputs"].min()
    hg_uni = pd.read_csv(HGNC_FILE, sep="\t", dtype=str, usecols=["symbol", "uniprot_ids"])
    has_uni = set(hg_uni.loc[hg_uni["uniprot_ids"].notna(), "symbol"])
    mech_inputs = {
        "definition": ("pair classifiable if every Phase II+ drug-target row for the pair has a "
                       "known ChEMBL molecule type (Open Targets drugType not Unknown) and an "
                       "action type; gene classifiable if all its final pairs are; localization "
                       "input = gene has a UniProt accession in HGNC"),
        "pairs_classifiable": int(per_pair["all_inputs"].sum()),
        "pairs_with_any_known_molecule_type": int(per_pair["any_type"].sum()),
        "genes_classifiable": int(genes_all.sum()),
        "genes_with_uniprot_accession": len(set(final["gene"]) & has_uni),
        "genes_classifiable_and_uniprot": int(sum(1 for g, ok in genes_all.items() if ok and g in has_uni)),
        "mechanism_classes_computed": False,
    }

    phase_marginal = {"max_phase_II": int(((final["max_phase"] >= 2) & (final["max_phase"] < 3)).sum()),
                      "max_phase_III_or_later": int((final["max_phase"] >= 3).sum())}

    cand_q = cands[cands["qualifies"]]
    by_source = {}
    for s, g in cands.groupby("source"):
        by_source[s] = {"rows": int(len(g)), "qualifying": int(g["qualifies"].sum()),
                        "qualifying_mapped": int((g["qualifies"] & (g["n_mapped_ids"] > 0)).sum()),
                        "qualifying_includes_ukb": dict(Counter(g.loc[g["qualifies"], "includes_ukb"])),
                        "qualifying_includes_iceland": dict(Counter(g.loc[g["qualifies"], "includes_iceland"]))}

    counts = {
        "indications": {
            "total_with_eligible_pair": len(ind_all),
            "with_tier1_candidate": len(ind_t1),
            "tier2_or_3_only": len((ind_t2 | ind_t3) - ind_t1),
            "tier2_descendant_only": len(ind_t2 - ind_t1),
            "tier3_ancestor_only_no_tier2": len(ind_t3 - ind_t1 - ind_t2),
            "tier2_child_or_tier3_parent_only": len((ind_t2c | ind_t3p) - ind_t1),
            "no_candidate_any_tier": len(set(ind_all) - ind_t1 - ind_t2 - ind_t3),
            "not_in_ot_disease_index": int((~tiers["in_ot_disease_index"]).sum()),
            "id_prefixes": dict(Counter(d.split("_")[0] for d in ind_all)),
        },
        "pairs": {
            "eligible": block(P),
            "eligible_any_single_protein_moa": block(P[P["any_single_protein"]]),
            "pqtl_instrumentable": block(inst_),
            "pqtl_by_source": {s: block(P[P[s]]) for s in
                               ("ukbppp", "decode", "interval_st4", "interval_epigraphdb", "interval")},
            "pqtl_interval_only": block(P[P["interval"] & ~P["ukbppp"] & ~P["decode"]]),
            "pqtl_with_tier1": block(t1),
            "pqtl_with_tier1_no_known_overlap": block(t1no),
            "final_heldout": block(final),
            "final_heldout_strict_ontology": block(final_strict),
            "pqtl_heldout_any_candidate_status": block(inst_[inst_["heldout"]]),
            "pqtl_with_tier1_heldout_ignoring_overlap": block(t1[t1["heldout"]]),
            "final_heldout_if_tier1_or_2": final_with(P["has_tier1"] | has_t2),
            "final_heldout_if_tier1_or_child_or_parent": final_with(P["has_tier1"] | has_t2c | has_t3p),
            "final_heldout_if_any_tier": final_with(P["has_tier1"] | has_t2 | has_t3),
            "seen_keys_present_in_eligible_universe": int((~P["heldout"]).sum()),
        },
        "final_heldout_pairs_per_gene": distribution(final.groupby("gene").size()),
        "final_heldout_indication_granularity": {
            "pairs_by_id_prefix": dict(Counter(d.split("_")[0] for d in final["diseaseId"])),
            "pairs_on_therapeutic_area_terms": int(final["diseaseId"].map(ta_flag).fillna(False).sum()),
            "pairs_on_leaf_terms": int(final["diseaseId"].map(leaf_flag).fillna(False).sum()),
            "indications_on_therapeutic_area_terms": int(pd.Series(final["diseaseId"].unique()).map(ta_flag).fillna(False).sum()),
            "indications_on_leaf_terms": int(pd.Series(final["diseaseId"].unique()).map(leaf_flag).fillna(False).sum()),
        },
        "final_heldout_max_phase_marginal": phase_marginal,
        "localization_genes": {"final_heldout": localization(set(final["gene"])),
                               "pqtl_instrumentable_eligible": localization(set(inst_["gene"]))},
        "mechanism_rule_inputs_final_heldout": mech_inputs,
    }

    tier_rows = {d: {"name": dis_name.get(d), **{k: [cands.at[i, "accession"] for i in v]
                                                   for k, v in r.items() if k != "in_ot_disease_index"}}
                 for d, r in tiers.iterrows() if len(r["tier1"]) or len(r["tier2_descendant"])}

    result = {
        "experiment": "08_mechanism_interaction feasibility v2: every indication, outcome-blind",
        "run_timestamp_utc": run_ts,
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "seed": None,
        "outcome_blindness": (
            "clinical_indication read with drugId, diseaseId, maxClinicalStage only; "
            "clinical_report never read; drug_molecule read with id, drugType, parentId only; "
            "classification_v5*.csv read with gene, disease only; Sun 2018 sheet ST4 only; no "
            "Eldjarn 2023 sheet read; no MR, colocalization or disease-association statistic "
            "computed. Per-pair phase is not written; phase appears only as the final held-out "
            "II vs III+ marginal."),
        "definitions": {
            "pair": "(HGNC approved symbol of an Open Targets MoA target, Open Targets diseaseId)",
            "eligible": "max over drugs of PHASE_MAP[maxClinicalStage] >= 2 (PHASE_2, PHASE_2_3, PHASE_3, PREAPPROVAL, APPROVAL)",
            "drug_mechanisms": "drug's own drug_mechanism_of_action rows; if none, its parent molecule's rows; all targets of each row",
            "ukbppp": "gene symbol single-gene entry in ukbppp_cis_pqtl_genes.json, or Ensembl ID of that symbol in olink_protein_map_3k_v1.tsv",
            "decode": "Ferkingstad 2021 ST02 cis rows (symbol, or Ensembl via ST01 SeqId), as in count_coverage.py",
            "interval_st4": "Sun 2018 ST4 cis rows on single-UniProt SOMAmers, UniProt -> HGNC symbol",
            "interval_epigraphdb": "gene_symbol of v34_mr_catalog.csv, as in count_coverage.py (V3.4-restricted list)",
            "pqtl_instrumentable": "ukbppp or decode or interval_st4 or interval_epigraphdb",
            "candidate_qualifies": f"case-control (ncase, ncontrol > 0), European only, full summary statistics, effective N 4/(1/ncase+1/ncontrol) >= {NEFF_FLOOR}",
            "tier1": "candidate mapped ID equals the indication ID (or an obsolete term OT lists for it)",
            "tier2": "candidate mapped ID is a descendant of the indication (tier2_child: direct child)",
            "tier3": "candidate mapped ID is an ancestor of the indication (tier3_parent: direct parent)",
            "overlap_flags": "yes if metadata names UK Biobank / Iceland or deCODE / INTERVAL; no_per_cohort_list if the source lists cohorts and none matches (FinnGen by construction); otherwise unknown",
            "no_known_overlap": "some pQTL source available for the gene and some tier-1 candidate whose flag for that source (UKB-PPP: includes_ukb; deCODE: includes_iceland; INTERVAL: includes_interval) is not yes",
            "heldout": "(gene, diseaseId) absent from classification_v5.csv and classification_v5_1.csv, disease labels mapped to IDs via coverage_counts.json disease_ids",
            "heldout_strict_ontology": "also drops pairs whose indication is an ancestor or descendant of a seen indication for the same gene",
            "gwas_catalog_mapping": "MAPPED_TRAIT_URI; ancestry from the ancestries file (initial stage all European); ncase/ncontrol from the ancestries file, else parsed from INITIAL SAMPLE SIZE",
            "opengwas_mapping": "ebi-a-* via the GWAS Catalog accession; finn-b-* via the Open Targets FinnGen endpoint code; else the ontology field; else exact match of trait label to an OT disease name or exact synonym",
            "finngen": "FinnGen studies in the Open Targets 26.09 study index with their OT disease mappings; the FinnGen manifest carries no ontology terms",
        },
        "releases": {
            "open_targets": OT_RELEASE,
            "gwas_catalog": "v1.0.3.1 files from releases/latest (see inputs sha256)",
            "finngen": fg_meta.get("finngen_projects"),
            "uniprot": {k: v.get("release") for k, v in fetch_status["uniprot"].items()},
            "hpa": zipfile.ZipFile(HPA_FILE).getinfo("proteinatlas.tsv").date_time,
            "opengwas": fetch_status["opengwas"],
        },
        "universe_meta": uni_meta,
        "instrument_source_sizes": {
            "ukbppp_single_gene": len(inst["ukbppp"]["sym"]),
            "decode_cis_symbols": len(inst["decode"]["sym"]),
            "interval_st4": {k: (len(v) if isinstance(v, set) else v) for k, v in inst["interval_st4"].items()},
            "interval_epigraphdb": len(inst["interval_epigraphdb"]["sym"]),
        },
        "outcome_gwas_sources": {"gwas_catalog_raw_rows": gcat_raw, "finngen": fg_meta,
                                 "opengwas": og_meta, "by_source": by_source,
                                 "qualifying_total": int(len(cand_q))},
        "reproduction_24_indications": repro,
        "seen_sets": {"keys": len(seen_keys), "unmapped_labels": unmapped_seen},
        "counts": counts,
        "indications_with_tier1_or_tier2_candidates": tier_rows,
        "outputs": {"final_heldout_pairs_csv": str(OUT_PAIRS),
                    "outcome_gwas_candidates_csv": str(OUT_CANDIDATES)},
    }

    cand_out = cands[cands["qualifies"] & (cands["n_mapped_ids"] > 0)].copy()
    cand_out["mapped_ids"] = cand_out["mapped_ids"].map(";".join)
    cand_out.to_csv(OUT_CANDIDATES, index=False)
    fo = final[["gene", "ensembl", "diseaseId", "ukbppp", "decode", "interval"]].copy()
    fo["disease_name"] = fo["diseaseId"].map(dis_name)
    fo["n_tier1_candidates"] = final["tier1_idx"].map(len)
    fo.sort_values(["diseaseId", "gene"]).to_csv(OUT_PAIRS, index=False)

    result["inputs"] = {
        "classification_v5": input_record(CLASSIFICATION_V5, "usecols gene, disease"),
        "classification_v5_1": input_record(CLASSIFICATION_V5_1, "usecols gene, disease"),
        "earlier_coverage_counts": input_record(EARLIER_JSON, "disease_ids label map and 24-indication totals"),
        "epigraphdb_catalog": input_record(EPIGRAPHDB_CATALOG, "usecols gene_symbol"),
        "ukbppp_cis_genes": input_record(UKBPPP_CIS, "UKB-PPP cis-pQTL gene list"),
        "olink_protein_map_3k": input_record(OLINK_MAP, "usecols HGNC.symbol, ensembl_id"),
        "decode_supplementary": input_record(DECODE_XLSX, "Ferkingstad 2021 ST01, ST02"),
        "sun2018_supplementary": input_record(SUN2018_XLSX, "sheet ST4 only", SUN2018_URL),
        "hgnc": input_record(HGNC_FILE, "symbol, status, ensembl_gene_id, uniprot_ids", HGNC_URL),
        "hpa": input_record(HPA_FILE, "Gene, Ensembl, Secretome location", HPA_URL),
        "opengwas_gwasinfo": input_record(OPENGWAS_FILE, "full gwasinfo listing", OPENGWAS_URL)
        if OPENGWAS_FILE.exists() else {"path": str(OPENGWAS_FILE), "sha256": None},
        **{f"ot_{k}": input_record(OT_DIR / f"{k}.parquet", "", f"{OT_BASE}/{v}") for k, v in OT_FILES.items()},
        **{f"gwas_catalog_{fn}": input_record(GWASCAT_DIR / fn, "", f"{GWASCAT_BASE}/{ep}")
           for fn, ep in GWASCAT_ENDPOINTS.items()},
        **{f"uniprot_{k}": input_record(UNIPROT_DIR / f"secreted_{k}.tsv", q, UNIPROT_SEARCH)
           for k, q in UNIPROT_QUERIES.items()},
    }
    OUT_JSON.write_text(json.dumps(result, indent=2, default=str))

    c = counts
    print(f"indications with an eligible pair       {c['indications']['total_with_eligible_pair']}")
    print(f"  with tier-1 candidate                 {c['indications']['with_tier1_candidate']}")
    print(f"  tier-2/3 only                         {c['indications']['tier2_or_3_only']}")
    for k in ("eligible", "pqtl_instrumentable", "pqtl_with_tier1", "pqtl_with_tier1_no_known_overlap",
              "final_heldout", "final_heldout_strict_ontology"):
        b = c["pairs"][k]
        print(f"{k:<40}{b['pairs']:>7} pairs {b['genes']:>6} genes {b['indications']:>6} indications")
    print(f"final held-out max phase: II {phase_marginal['max_phase_II']}, "
          f"III+ {phase_marginal['max_phase_III_or_later']}")
    print(f"\nWrote {OUT_JSON}")


if __name__ == "__main__":
    main()
