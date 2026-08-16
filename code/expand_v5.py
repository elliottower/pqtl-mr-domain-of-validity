"""
V5 Expansion: Recover missing pQTL MR results for abundance-modulating pairs.

Four recovery steps:
  1. Add outcome GWAS for JIA and neuroblastoma
  2. Resolve non-rs variant IDs via ENSEMBL variant recoder
  3. LD proxy lookup for remaining SNP-GWAS mismatches
  4. deCODE pQTL source for EpiGraphDB gaps (Ferkingstad et al. 2021)

Reads frozen_candidates_v34.csv + classification_v34.csv, outputs expanded
classification to results/v5/.

Usage:
    cd ~/Documents/GitHub/pqtl-mr-domain-of-validity
    export $(grep GWAS_KEY ~/Documents/GitHub/causal-inference-neuro-epidemiology/.env)
    uv run python code/expand_v5.py
"""
import gzip
import json
import os
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from scipy import stats
from scipy.stats import binomtest
from tqdm import tqdm

ROOT = Path(__file__).parent.parent
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "results" / "v5"
OUT_DIR.mkdir(parents=True, exist_ok=True)

API = "https://api.opengwas.io/api"
ENSEMBL_API = "https://rest.ensembl.org"
ENSEMBL_GRCH37 = "https://grch37.rest.ensembl.org"

DISEASE_GWAS = {
    "Alzheimers disease": "ieu-b-2",
    "Amyotrophic lateral sclerosis": "ebi-a-GCST90027163",
    "Anorexia nervosa": "ieu-b-61",
    "Autism spectrum disorder": "ieu-b-87",
    "Bipolar disorder": "ieu-b-41",
    "Chronic kidney disease": "ieu-b-4874",
    "Crohns disease": "ieu-a-12",
    "Glioma": "ieu-b-4987",
    "Hypercholesterolemia": "ieu-a-300",
    "Inflammatory bowel disease": "ieu-a-31",
    "Juvenile idiopathic arthritis": "ebi-a-GCST90018873",
    "Lung cancer": "ieu-a-984",
    "Major depressive disorder": "ieu-b-102",
    "Melanoma": "ieu-a-62",
    "Multiple sclerosis": "ieu-b-18",
    "Myocardial infarction": "ieu-a-798",
    "Neuroblastoma": "ieu-a-816",
    "Ovarian cancer": "ieu-a-1120",
    "Pancreatic cancer": "ieu-b-4866",
    "Parkinsons disease": "ieu-b-7",
    "Rheumatoid arthritis": "ieu-a-833",
    "Schizophrenia": "ieu-b-5102",
    "Systemic lupus erythematosus": "ieu-a-1073",
    "Thyroid cancer": "ieu-a-1082",
    "Ulcerative colitis": "ieu-a-970",
}

GWAS_CACHE_FILE = OUT_DIR / "outcome_gwas_cache_v5.json"
VARIANT_MAP_FILE = OUT_DIR / "variant_id_map_v5.json"
LD_PROXY_FILE = OUT_DIR / "ld_proxy_map_v5.json"
DECODE_CACHE_FILE = OUT_DIR / "decode_cis_pqtl_v5.json"

DECODE_DIR = DATA_DIR / "decode"

# GRCh38 gene coordinates from ENSEMBL (matching deCODE coordinate system)
DECODE_GENES = {
    "PCSK9": {
        "chrom": "chr1", "start": 55039445, "end": 55064852,
        "files": ["5231_79_PCSK9_PCSK9.txt.gz"],
        "diseases": ["Hypercholesterolemia", "Myocardial infarction"],
    },
    "GHR": {
        "chrom": "chr5", "start": 42423439, "end": 42721878,
        "files": ["2948_58_GHR_Growth_hormone_receptor.txt.gz"],
        "diseases": ["Chronic kidney disease", "Juvenile idiopathic arthritis"],
    },
    "CSF2RB": {
        "chrom": "chr22", "start": 36913590, "end": 36940454,
        "files": [
            "10512_13_CSF2RB_IL3RB.txt.gz",
            "11137_43_CSF2RB_IL3RB.txt.gz",
        ],
        "diseases": ["Crohns disease"],
    },
    "CD86": {
        "chrom": "chr3", "start": 122055356, "end": 122121557,
        "diseases": ["Chronic kidney disease", "Crohns disease",
                     "Juvenile idiopathic arthritis", "Rheumatoid arthritis",
                     "Ulcerative colitis"],
    },
    "PDCD1": {
        "chrom": "chr2", "start": 241849884, "end": 241858894,
        "diseases": ["Lung cancer", "Ovarian cancer"],
    },
    "MMP13": {
        "chrom": "chr11", "start": 102942994, "end": 102955740,
        "diseases": ["Alzheimers disease"],
    },
    "MIF": {
        "chrom": "chr22", "start": 23894372, "end": 23895227,
        "diseases": ["Rheumatoid arthritis"],
    },
    "TNF": {
        "chrom": "chr6", "start": 31575558, "end": 31578336,
        "diseases": ["Crohns disease", "Rheumatoid arthritis",
                     "Ulcerative colitis"],
    },
    "ITGB7": {
        "chrom": "chr12", "start": 53191318, "end": 53207310,
        "diseases": ["Crohns disease", "Multiple sclerosis",
                     "Ulcerative colitis"],
    },
    "PLA2G10": {
        "chrom": "chr16", "start": 14672382, "end": 14694493,
        "diseases": ["Myocardial infarction", "Rheumatoid arthritis"],
    },
    "SERPINC1": {
        "chrom": "chr1", "start": 173903414, "end": 173917543,
        "diseases": ["Chronic kidney disease", "Myocardial infarction"],
    },
    "TNFSF13B": {
        "chrom": "chr13", "start": 108251240, "end": 108308484,
        "diseases": ["Rheumatoid arthritis",
                     "Systemic lupus erythematosus"],
    },
    "PTH1R": {
        "chrom": "chr3", "start": 46877706, "end": 46903799,
        "diseases": ["Anorexia nervosa"],
    },
    "CSF3R": {
        "chrom": "chr1", "start": 36465669, "end": 36483322,
        "diseases": ["Neuroblastoma"],
    },
}
CIS_WINDOW = 1_000_000


def find_decode_files(gene):
    """Find deCODE summary stats files for a gene in the decode directory."""
    if not DECODE_DIR.exists():
        return []
    matches = []
    gene_upper = gene.upper()
    for f in sorted(DECODE_DIR.glob("*.txt.gz")):
        parts = f.stem.replace(".txt", "").split("_")
        gene_parts = [p.upper() for p in parts[2:]]
        if gene_upper in gene_parts:
            matches.append(f.name)
    return matches


def extract_top_cis_pqtl(gene, gene_info):
    """Extract the top cis-pQTL from deCODE summary stats (streaming, ~900 MB/file)."""
    chrom = gene_info["chrom"]
    cis_start = gene_info["start"] - CIS_WINDOW
    cis_end = gene_info["end"] + CIS_WINDOW
    best = None

    files = gene_info.get("files") or find_decode_files(gene)
    if not files:
        print(f"    No deCODE files found for {gene}")
        return None

    for fname in files:
        fpath = DECODE_DIR / fname
        if not fpath.exists():
            print(f"    WARNING: {fpath} not found, skipping")
            continue

        print(f"    Scanning {fname} for cis-pQTLs in {chrom}:{cis_start}-{cis_end}...")
        n_cis = 0
        with gzip.open(fpath, "rt") as f:
            header = f.readline().strip().split("\t")
            col = {name: i for i, name in enumerate(header)}

            for line in tqdm(f, desc=f"    {fname}", unit=" variants", mininterval=5):
                fields = line.strip().split("\t")
                if fields[col["Chrom"]] != chrom:
                    continue
                pos = int(fields[col["Pos"]])
                if pos < cis_start or pos > cis_end:
                    continue
                n_cis += 1

                pval = float(fields[col["Pval"]])
                if pval >= 5e-8:
                    continue

                rsid_field = fields[col["rsids"]]
                if rsid_field == "NA" or not rsid_field.startswith("rs"):
                    continue

                beta = float(fields[col["Beta"]])
                se = float(fields[col["SE"]])

                if best is None or pval < best["pval"]:
                    best = {
                        "rsid": rsid_field,
                        "chrom": chrom,
                        "pos": pos,
                        "effect_allele": fields[col["effectAllele"]],
                        "other_allele": fields[col["otherAllele"]],
                        "beta": beta,
                        "se": se,
                        "pval": pval,
                        "n": int(fields[col["N"]]),
                        "maf": float(fields[col["ImpMAF"]]),
                        "source_file": fname,
                    }

        print(f"    {fname}: {n_cis} cis variants, best p={best['pval']:.2e} ({best['rsid']})" if best else f"    {fname}: {n_cis} cis variants, no genome-wide significant cis-pQTL")

    return best


def _gwas_headers() -> dict:
    token = os.environ.get("GWAS_KEY", "")
    h = {"X-API-SOURCE": "v5-expansion/0.1"}
    if token:
        h["Authorization"] = f"Bearer {token}"
    return h


def compute_wald_ratio(beta_exposure, se_exposure, beta_outcome, se_outcome):
    wald_beta = beta_outcome / beta_exposure
    wald_se = abs(se_outcome / beta_exposure)
    z = wald_beta / wald_se if wald_se > 0 else 0.0
    p = float(2 * stats.norm.sf(abs(z)))
    return {
        "beta": float(wald_beta),
        "se": float(wald_se),
        "p": p,
        "ci_lower": float(wald_beta - 1.96 * wald_se),
        "ci_upper": float(wald_beta + 1.96 * wald_se),
        "mr_supports_causal": p < 0.05,
    }


def compute_metrics(predictions, labels):
    tp = sum(p == "SUCCESS" and l == "SUCCESS" for p, l in zip(predictions, labels))
    tn = sum(p == "FAILURE" and l == "FAILURE" for p, l in zip(predictions, labels))
    fp = sum(p == "SUCCESS" and l == "FAILURE" for p, l in zip(predictions, labels))
    fn = sum(p == "FAILURE" and l == "SUCCESS" for p, l in zip(predictions, labels))
    sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    ba = (sens + spec) / 2 if not (np.isnan(sens) or np.isnan(spec)) else float("nan")
    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = float((tp * tn - fp * fn) / denom) if denom > 0 else float("nan")
    return {
        "n": len(predictions), "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sensitivity": sens, "specificity": spec, "balanced_accuracy": ba, "mcc": mcc,
    }


def bootstrap_ba(predictions, labels, n_boot=10000):
    rng = np.random.default_rng()
    bas = []
    n = len(predictions)
    preds_arr = np.array(predictions)
    labs_arr = np.array(labels)
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        p_boot = preds_arr[idx]
        l_boot = labs_arr[idx]
        tp = ((p_boot == "SUCCESS") & (l_boot == "SUCCESS")).sum()
        fn = ((p_boot == "FAILURE") & (l_boot == "SUCCESS")).sum()
        tn = ((p_boot == "FAILURE") & (l_boot == "FAILURE")).sum()
        fp = ((p_boot == "SUCCESS") & (l_boot == "FAILURE")).sum()
        sens = tp / (tp + fn) if (tp + fn) > 0 else np.nan
        spec = tn / (tn + fp) if (tn + fp) > 0 else np.nan
        if not (np.isnan(sens) or np.isnan(spec)):
            bas.append((sens + spec) / 2)
    bas = np.array(bas)
    return {
        "ba_mean": float(np.mean(bas)),
        "ba_ci_lo": float(np.percentile(bas, 5)),
        "ba_ci_hi": float(np.percentile(bas, 95)),
        "n_boot": n_boot,
    }


def load_cache(path):
    if path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(cache, path):
    with open(path, "w") as f:
        json.dump(cache, f, indent=2)


def is_non_rs_variant(rsid):
    if pd.isna(rsid) or not rsid:
        return False
    return not rsid.startswith("rs")


def parse_chrpos_variant(variant_id):
    """Parse chr:pos_ref_alt format to (chr, pos, ref, alt)."""
    parts = variant_id.split("_")
    if len(parts) >= 2:
        chrpos = parts[0]
        if ":" in chrpos:
            chrom, pos = chrpos.split(":", 1)
            ref = parts[1] if len(parts) > 1 else None
            alt = parts[2] if len(parts) > 2 else None
            return chrom, int(pos), ref, alt
    return None, None, None, None


def resolve_variant_ensembl(variant_id, max_retries=3):
    """Use ENSEMBL GRCh37 variant recoder to find rs-number for chr:pos variant."""
    chrom, pos, ref, alt = parse_chrpos_variant(variant_id)
    if chrom is None:
        return None

    for attempt in range(max_retries):
        try:
            url = f"{ENSEMBL_GRCH37}/vep/human/region/{chrom}:{pos}-{pos}/{alt}"
            r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list) and len(data) > 0:
                    rsid = data[0].get("id")
                    if rsid and rsid.startswith("rs"):
                        return rsid
                    colocated = data[0].get("colocated_variants", [])
                    for cv in colocated:
                        cv_id = cv.get("id", "")
                        if cv_id.startswith("rs"):
                            return cv_id
            elif r.status_code == 429:
                time.sleep(2 * (attempt + 1))
                continue
            return None
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            time.sleep(2 * (attempt + 1))
    return None


def find_ld_proxy(rsid, gwas_id, gwas_cache, max_retries=3):
    """Find LD proxy SNP that exists in the outcome GWAS."""
    for attempt in range(max_retries):
        try:
            url = (
                f"{ENSEMBL_GRCH37}/ld/human/{rsid}/1000GENOMES:phase_3:EUR"
                f"?r2=0.6&window_size=500"
            )
            r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list):
                    proxies = sorted(data, key=lambda x: float(x.get("r2", 0)), reverse=True)
                    proxy_rsids = []
                    for p in proxies[:20]:
                        other = p.get("variation2") if p.get("variation1") == rsid else p.get("variation1")
                        if other and other.startswith("rs"):
                            proxy_rsids.append((other, float(p.get("r2", 0))))

                    if proxy_rsids:
                        rsids_to_check = [p[0] for p in proxy_rsids]
                        associations = query_outcome_batch(rsids_to_check, gwas_id)
                        for proxy_rs, r2 in proxy_rsids:
                            if proxy_rs in associations:
                                assoc = associations[proxy_rs]
                                if assoc.get("beta") is not None and assoc.get("se") is not None:
                                    return {
                                        "proxy_rsid": proxy_rs,
                                        "r2": r2,
                                        "beta": float(assoc["beta"]),
                                        "se": float(assoc["se"]),
                                    }
            elif r.status_code == 429:
                time.sleep(5 * (attempt + 1))
                continue
            return None
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            time.sleep(3 * (attempt + 1))
    return None


def query_outcome_batch(rsids, gwas_id, max_retries=3):
    results = {}
    for attempt in range(max_retries):
        try:
            r = requests.post(
                f"{API}/associations",
                headers=_gwas_headers(),
                data={"variant": rsids, "id": gwas_id},
                timeout=180,
            )
            if r.status_code == 200:
                data = r.json()
                if isinstance(data, list):
                    for row in data:
                        snp = row.get("rsid") or row.get("name")
                        if snp:
                            results[snp] = row
                return results
            if r.status_code == 429:
                time.sleep(10 * (attempt + 1))
                continue
            return results
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError):
            time.sleep(5 * (attempt + 1))
    return results


def query_single_snp(rsid, gwas_id):
    results = query_outcome_batch([rsid], gwas_id)
    return results.get(rsid)


def run_expansion():
    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] V5 Expansion: Recovering missing MR results")
    print(f"{'='*70}")

    if not os.environ.get("GWAS_KEY"):
        print("WARNING: GWAS_KEY not set. Set it with:")
        print("  export $(grep GWAS_KEY ~/Documents/GitHub/causal-inference-neuro-epidemiology/.env)")
        print()

    candidates = pd.read_csv(DATA_DIR / "frozen_candidates_v34.csv")
    existing_clf = pd.read_csv(DATA_DIR / "classification_v34.csv")
    classified_ids = set(existing_clf["pair_id"])

    missing = candidates[~candidates["pair_id"].isin(classified_ids)].copy()
    abund_missing = missing[missing["mechanism_class"] == "abundance_modulating"]
    act_missing = missing[missing["mechanism_class"] == "activity_blocking"]

    print(f"Total candidates: {len(candidates)}")
    print(f"Already classified: {len(existing_clf)}")
    print(f"Missing total: {len(missing)}")
    print(f"  abundance_modulating: {len(abund_missing)}")
    print(f"  activity_blocking: {len(act_missing)}")
    print()

    gwas_cache = load_cache(GWAS_CACHE_FILE)
    variant_map = load_cache(VARIANT_MAP_FILE)
    ld_proxy_map = load_cache(LD_PROXY_FILE)

    recovered = {}
    recovery_log = []

    # ================================================================
    # STEP 1: New GWAS IDs (JIA, neuroblastoma)
    # ================================================================
    print(f"--- Step 1: New outcome GWAS IDs ---")
    step1_diseases = {"Juvenile idiopathic arthritis", "Neuroblastoma"}
    step1_pairs = missing[missing["disease"].isin(step1_diseases)]
    step1_with_rsid = step1_pairs[step1_pairs["sentinel_rsid"].notna()]

    print(f"  Pairs with new GWAS: {len(step1_pairs)} ({len(step1_with_rsid)} with sentinel RSIDs)")

    for _, row in tqdm(step1_with_rsid.iterrows(), total=len(step1_with_rsid), desc="Step 1"):
        pid = row["pair_id"]
        gwas_id = DISEASE_GWAS[row["disease"]]
        rsid = row["sentinel_rsid"]

        if is_non_rs_variant(rsid):
            continue

        cache_key = f"{rsid}_{gwas_id}"
        if cache_key not in gwas_cache:
            result = query_single_snp(rsid, gwas_id)
            if result and result.get("beta") is not None:
                gwas_cache[cache_key] = {
                    "beta": float(result["beta"]),
                    "se": float(result["se"]),
                }
            else:
                gwas_cache[cache_key] = None
            save_cache(gwas_cache, GWAS_CACHE_FILE)
            time.sleep(0.5)

        cached = gwas_cache.get(cache_key)
        if cached and pd.notna(row["sentinel_beta"]) and pd.notna(row["sentinel_se"]):
            mr = compute_wald_ratio(
                float(row["sentinel_beta"]), float(row["sentinel_se"]),
                cached["beta"], cached["se"],
            )
            mr["rsid"] = rsid
            mr["source"] = "v5_step1_new_gwas"
            recovered[pid] = mr
            recovery_log.append({
                "pair_id": pid, "step": 1, "method": "new_gwas",
                "gwas_id": gwas_id, "rsid": rsid, "mr_p": mr["p"],
            })

    print(f"  Recovered: {len([r for r in recovery_log if r['step']==1])}")

    # ================================================================
    # STEP 2: Resolve non-rs variant IDs
    # ================================================================
    print(f"\n--- Step 2: Resolve non-rs variant IDs ---")
    non_rs_pairs = missing[
        missing["sentinel_rsid"].apply(lambda x: is_non_rs_variant(x))
        & ~missing["pair_id"].isin(recovered)
    ]
    print(f"  Pairs with non-rs variants: {len(non_rs_pairs)}")

    unique_non_rs = non_rs_pairs["sentinel_rsid"].unique()
    print(f"  Unique non-rs variants to resolve: {len(unique_non_rs)}")

    for variant_id in tqdm(unique_non_rs, desc="Step 2: Resolving"):
        if variant_id in variant_map:
            continue
        resolved = resolve_variant_ensembl(variant_id)
        variant_map[variant_id] = resolved
        save_cache(variant_map, VARIANT_MAP_FILE)
        time.sleep(0.3)
        if resolved:
            print(f"    {variant_id} -> {resolved}")
        else:
            print(f"    {variant_id} -> NO RS-NUMBER FOUND")

    for _, row in non_rs_pairs.iterrows():
        pid = row["pair_id"]
        if pid in recovered:
            continue

        gwas_id = DISEASE_GWAS.get(row["disease"])
        if gwas_id is None:
            continue

        original_rsid = row["sentinel_rsid"]
        resolved_rsid = variant_map.get(original_rsid)
        if not resolved_rsid:
            continue

        cache_key = f"{resolved_rsid}_{gwas_id}"
        if cache_key not in gwas_cache:
            result = query_single_snp(resolved_rsid, gwas_id)
            if result and result.get("beta") is not None:
                gwas_cache[cache_key] = {
                    "beta": float(result["beta"]),
                    "se": float(result["se"]),
                }
            else:
                gwas_cache[cache_key] = None
            save_cache(gwas_cache, GWAS_CACHE_FILE)
            time.sleep(0.5)

        cached = gwas_cache.get(cache_key)
        if cached and pd.notna(row["sentinel_beta"]) and pd.notna(row["sentinel_se"]):
            mr = compute_wald_ratio(
                float(row["sentinel_beta"]), float(row["sentinel_se"]),
                cached["beta"], cached["se"],
            )
            mr["rsid"] = resolved_rsid
            mr["source"] = "v5_step2_variant_resolved"
            mr["original_variant"] = original_rsid
            recovered[pid] = mr
            recovery_log.append({
                "pair_id": pid, "step": 2, "method": "variant_resolved",
                "original_variant": original_rsid, "resolved_rsid": resolved_rsid,
                "gwas_id": gwas_id, "mr_p": mr["p"],
            })

    print(f"  Recovered: {len([r for r in recovery_log if r['step']==2])}")

    # ================================================================
    # STEP 3: LD proxy lookup
    # ================================================================
    print(f"\n--- Step 3: LD proxy lookup ---")
    still_missing = missing[
        ~missing["pair_id"].isin(recovered)
        & ~missing["pair_id"].isin(step1_pairs[step1_pairs["sentinel_rsid"].isna()]["pair_id"])
        & missing["sentinel_rsid"].notna()
    ]

    rs_missing = still_missing[still_missing["sentinel_rsid"].apply(lambda x: str(x).startswith("rs"))]
    non_rs_resolved = still_missing[
        still_missing["sentinel_rsid"].apply(lambda x: is_non_rs_variant(x))
    ]

    proxy_candidates = []
    for _, row in rs_missing.iterrows():
        gwas_id = DISEASE_GWAS.get(row["disease"])
        if gwas_id:
            proxy_candidates.append(row)

    for _, row in non_rs_resolved.iterrows():
        resolved = variant_map.get(row["sentinel_rsid"])
        if resolved:
            gwas_id = DISEASE_GWAS.get(row["disease"])
            if gwas_id:
                proxy_candidates.append(row)

    print(f"  Pairs eligible for LD proxy: {len(proxy_candidates)}")

    for row in tqdm(proxy_candidates, desc="Step 3: LD proxy"):
        pid = row["pair_id"] if isinstance(row, pd.Series) else row.get("pair_id")
        if pid in recovered:
            continue

        rsid = row["sentinel_rsid"] if isinstance(row, pd.Series) else row.get("sentinel_rsid")
        if is_non_rs_variant(rsid):
            rsid = variant_map.get(rsid)
        if not rsid:
            continue

        gwas_id = DISEASE_GWAS.get(row["disease"] if isinstance(row, pd.Series) else row.get("disease"))
        if not gwas_id:
            continue

        proxy_key = f"{rsid}_{gwas_id}"
        if proxy_key in ld_proxy_map:
            proxy_info = ld_proxy_map[proxy_key]
        else:
            proxy_info = find_ld_proxy(rsid, gwas_id, gwas_cache)
            ld_proxy_map[proxy_key] = proxy_info
            save_cache(ld_proxy_map, LD_PROXY_FILE)
            time.sleep(1.0)

        if proxy_info and pd.notna(row["sentinel_beta"]) and pd.notna(row["sentinel_se"]):
            mr = compute_wald_ratio(
                float(row["sentinel_beta"]), float(row["sentinel_se"]),
                proxy_info["beta"], proxy_info["se"],
            )
            mr["rsid"] = proxy_info["proxy_rsid"]
            mr["source"] = "v5_step3_ld_proxy"
            mr["proxy_r2"] = proxy_info["r2"]
            mr["original_rsid"] = rsid
            recovered[pid] = mr
            recovery_log.append({
                "pair_id": pid, "step": 3, "method": "ld_proxy",
                "original_rsid": rsid, "proxy_rsid": proxy_info["proxy_rsid"],
                "r2": proxy_info["r2"], "gwas_id": gwas_id, "mr_p": mr["p"],
            })

    print(f"  Recovered: {len([r for r in recovery_log if r['step']==3])}")

    # ================================================================
    # STEP 4: deCODE cis-pQTL instruments
    # ================================================================
    print(f"\n--- Step 4: deCODE cis-pQTL instruments (Ferkingstad et al. 2021) ---")
    decode_cache = load_cache(DECODE_CACHE_FILE)

    for gene, ginfo in DECODE_GENES.items():
        if gene in decode_cache:
            pqtl = decode_cache[gene]
            if pqtl:
                print(f"  {gene}: cached sentinel {pqtl['rsid']} (p={pqtl['pval']:.2e})")
            else:
                print(f"  {gene}: cached — no genome-wide significant cis-pQTL")
        else:
            print(f"  Extracting cis-pQTL for {gene}...")
            pqtl = extract_top_cis_pqtl(gene, ginfo)
            decode_cache[gene] = pqtl
            save_cache(decode_cache, DECODE_CACHE_FILE)
            if pqtl:
                print(f"  {gene}: sentinel {pqtl['rsid']} (p={pqtl['pval']:.2e}, beta={pqtl['beta']:.4f})")
            else:
                print(f"  {gene}: NO genome-wide significant cis-pQTL found")

        if not pqtl:
            continue

        for disease in ginfo["diseases"]:
            pid = f"{gene}|{disease}"
            if pid in recovered or pid in classified_ids:
                continue

            gwas_id = DISEASE_GWAS.get(disease)
            if not gwas_id:
                print(f"    {pid}: no GWAS ID for {disease}, skipping")
                continue

            cache_key = f"{pqtl['rsid']}_{gwas_id}"
            if cache_key not in gwas_cache:
                result = query_single_snp(pqtl["rsid"], gwas_id)
                if result and result.get("beta") is not None:
                    gwas_cache[cache_key] = {
                        "beta": float(result["beta"]),
                        "se": float(result["se"]),
                    }
                else:
                    gwas_cache[cache_key] = None
                save_cache(gwas_cache, GWAS_CACHE_FILE)
                time.sleep(0.5)

            cached = gwas_cache.get(cache_key)
            if cached:
                mr = compute_wald_ratio(
                    pqtl["beta"], pqtl["se"],
                    cached["beta"], cached["se"],
                )
                mr["rsid"] = pqtl["rsid"]
                mr["source"] = "v5_step4_decode"
                mr["decode_file"] = pqtl["source_file"]
                mr["decode_pval"] = pqtl["pval"]
                recovered[pid] = mr
                recovery_log.append({
                    "pair_id": pid, "step": 4, "method": "decode_cis_pqtl",
                    "rsid": pqtl["rsid"], "gwas_id": gwas_id,
                    "decode_pval": pqtl["pval"], "mr_p": mr["p"],
                })
                print(f"    {pid}: MR p={mr['p']:.4e} (decode sentinel {pqtl['rsid']})")
            else:
                print(f"    {pid}: sentinel {pqtl['rsid']} not in {gwas_id}")
                # Try LD proxy for the deCODE sentinel
                proxy_key = f"{pqtl['rsid']}_{gwas_id}"
                if proxy_key not in ld_proxy_map:
                    proxy_info = find_ld_proxy(pqtl["rsid"], gwas_id, gwas_cache)
                    ld_proxy_map[proxy_key] = proxy_info
                    save_cache(ld_proxy_map, LD_PROXY_FILE)
                    time.sleep(1.0)
                else:
                    proxy_info = ld_proxy_map[proxy_key]

                if proxy_info:
                    mr = compute_wald_ratio(
                        pqtl["beta"], pqtl["se"],
                        proxy_info["beta"], proxy_info["se"],
                    )
                    mr["rsid"] = proxy_info["proxy_rsid"]
                    mr["source"] = "v5_step4_decode_ld_proxy"
                    mr["decode_file"] = pqtl["source_file"]
                    mr["decode_pval"] = pqtl["pval"]
                    mr["proxy_r2"] = proxy_info["r2"]
                    mr["original_rsid"] = pqtl["rsid"]
                    recovered[pid] = mr
                    recovery_log.append({
                        "pair_id": pid, "step": 4, "method": "decode_ld_proxy",
                        "original_rsid": pqtl["rsid"],
                        "proxy_rsid": proxy_info["proxy_rsid"],
                        "r2": proxy_info["r2"],
                        "gwas_id": gwas_id, "mr_p": mr["p"],
                    })
                    print(f"    {pid}: MR p={mr['p']:.4e} (LD proxy {proxy_info['proxy_rsid']}, r2={proxy_info['r2']:.2f})")
                else:
                    print(f"    {pid}: no LD proxy found either")

    print(f"  Recovered: {len([r for r in recovery_log if r['step']==4])}")

    # ================================================================
    # RESULTS
    # ================================================================
    print(f"\n{'='*70}")
    print(f"  EXPANSION RESULTS")
    print(f"{'='*70}")

    total_recovered = len(recovered)
    print(f"\nTotal recovered: {total_recovered}")
    for step in [1, 2, 3, 4]:
        n = len([r for r in recovery_log if r["step"] == step])
        print(f"  Step {step}: {n}")

    by_mech = {}
    for pid, mr in recovered.items():
        row = candidates[candidates["pair_id"] == pid].iloc[0]
        mech = row["mechanism_class"]
        by_mech.setdefault(mech, []).append(pid)
    for mech, pids in by_mech.items():
        print(f"  {mech}: {len(pids)} recovered")

    # Load adjudicated outcomes
    adj_path = DATA_DIR / "adjudicated_v34.csv"
    adj = pd.read_csv(adj_path)
    outcome_map = dict(zip(adj["pair_id"], adj["outcome"]))

    # Build expanded classification
    new_rows = []
    for pid, mr in recovered.items():
        row = candidates[candidates["pair_id"] == pid].iloc[0]
        outcome = outcome_map.get(pid, "")
        prediction = "SUCCESS" if mr["mr_supports_causal"] else "FAILURE"
        new_rows.append({
            "pair_id": pid,
            "gene": row["gene"],
            "disease": row["disease"],
            "outcome": outcome,
            "mechanism_class": row["mechanism_class"],
            "disease_area": row["disease_area"],
            "instrument_source": row["instrument_source"],
            "sample_overlap": row["sample_overlap"],
            "mr_beta": mr["beta"],
            "mr_se": mr["se"],
            "mr_p": mr["p"],
            "mr_causal": mr["mr_supports_causal"],
            "mr_source": mr["source"],
            "mr_rsid": mr["rsid"],
            "prediction": prediction,
            "correct": prediction == outcome,
        })

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        expanded = pd.concat([existing_clf, new_df], ignore_index=True)
    else:
        expanded = existing_clf.copy()

    expanded.to_csv(OUT_DIR / "classification_v5.csv", index=False)
    print(f"\nWrote expanded classification: {OUT_DIR / 'classification_v5.csv'}")
    print(f"  V3.4 pairs: {len(existing_clf)}")
    print(f"  New pairs:  {len(new_rows)}")
    print(f"  Total:      {len(expanded)}")

    # Filter to evaluable (have outcome)
    evaluable = expanded[expanded["outcome"].isin(["SUCCESS", "FAILURE"])].copy()

    if len(evaluable) == 0:
        print("  No evaluable pairs with outcomes — cannot compute metrics")
        save_cache({"recovery_log": recovery_log, "total_recovered": total_recovered}, OUT_DIR / "expansion_log_v5.json")
        return

    # ================================================================
    # RE-RUN ANALYSES on expanded set
    # ================================================================
    preds = evaluable["prediction"].tolist()
    labs = evaluable["outcome"].tolist()

    print(f"\n{'='*70}")
    print(f"  EXPANDED ANALYSES (V5)")
    print(f"{'='*70}")

    metrics = compute_metrics(preds, labs)
    boot = bootstrap_ba(preds, labs)

    print(f"\n  [1] PRIMARY: Pooled BA (n={metrics['n']})")
    print(f"      BA: {metrics['balanced_accuracy']:.3f} [{boot['ba_ci_lo']:.3f}, {boot['ba_ci_hi']:.3f}]")
    print(f"      TP={metrics['tp']} TN={metrics['tn']} FP={metrics['fp']} FN={metrics['fn']}")

    print(f"\n  [2] CO-PRIMARY: Mechanism stratification")
    mech_results = {}
    for mech in ["abundance_modulating", "activity_blocking"]:
        sub = evaluable[evaluable["mechanism_class"] == mech]
        if len(sub) < 2:
            continue
        s_preds = sub["prediction"].tolist()
        s_labs = sub["outcome"].tolist()
        s_metrics = compute_metrics(s_preds, s_labs)
        s_boot = bootstrap_ba(s_preds, s_labs)
        mech_results[mech] = {"metrics": s_metrics, "bootstrap": s_boot}
        ci_status = "ABOVE 0.50" if s_boot["ba_ci_lo"] > 0.50 else "includes 0.50"
        print(f"    {mech} (n={s_metrics['n']}): BA={s_metrics['balanced_accuracy']:.3f} "
              f"[{s_boot['ba_ci_lo']:.3f}, {s_boot['ba_ci_hi']:.3f}] -- {ci_status}")

    # Comparison with V3.4
    print(f"\n  [3] V3.4 vs V5 comparison")
    v34_abund = existing_clf[existing_clf["mechanism_class"] == "abundance_modulating"]
    v34_abund = v34_abund[v34_abund["outcome"].isin(["SUCCESS", "FAILURE"])]
    if len(v34_abund) > 0:
        v34_m = compute_metrics(v34_abund["prediction"].tolist(), v34_abund["outcome"].tolist())
        v34_b = bootstrap_ba(v34_abund["prediction"].tolist(), v34_abund["outcome"].tolist())
        print(f"    V3.4 abundance: n={v34_m['n']}, BA={v34_m['balanced_accuracy']:.3f} "
              f"[{v34_b['ba_ci_lo']:.3f}, {v34_b['ba_ci_hi']:.3f}]")
        if "abundance_modulating" in mech_results:
            v5_m = mech_results["abundance_modulating"]["metrics"]
            v5_b = mech_results["abundance_modulating"]["bootstrap"]
            delta = v5_m["balanced_accuracy"] - v34_m["balanced_accuracy"]
            print(f"    V5 abundance:   n={v5_m['n']}, BA={v5_m['balanced_accuracy']:.3f} "
                  f"[{v5_b['ba_ci_lo']:.3f}, {v5_b['ba_ci_hi']:.3f}]")
            print(f"    Delta BA: {delta:+.3f}, Delta n: +{v5_m['n'] - v34_m['n']}")

    # Save everything
    results = {
        "timestamp": ts,
        "version": "V5",
        "total_recovered": total_recovered,
        "recovery_by_step": {
            str(step): len([r for r in recovery_log if r["step"] == step])
            for step in [1, 2, 3, 4]
        },
        "n_v34": len(existing_clf),
        "n_v5": len(expanded),
        "n_evaluable": len(evaluable),
        "primary_metrics": metrics,
        "primary_bootstrap": boot,
        "mechanism_strat": {
            mech: {"metrics": mr["metrics"], "bootstrap": mr["bootstrap"]}
            for mech, mr in mech_results.items()
        },
        "recovery_log": recovery_log,
    }

    with open(OUT_DIR / "evaluation_v5.json", "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nWrote evaluation: {OUT_DIR / 'evaluation_v5.json'}")

    ts2 = datetime.now(timezone.utc).isoformat()
    print(f"\n[{ts2}] V5 expansion complete")


if __name__ == "__main__":
    run_expansion()
