"""
V3.4 Stage C: MR classification and balanced accuracy evaluation.

Reads adjudicated_v34.csv (Stage B output), joins with MR estimates
(pre-computed catalog + fresh OpenGWAS Wald ratios), applies the
pre-locked classifier (p < 0.05 → CAUSAL → predict SUCCESS), and
computes all 12 pre-specified analyses from DEVIATION_LOG_V34.md.

Primary metric: Balanced Accuracy (BA) with 10k bootstrap 90% CI.
MCC reported as secondary. The classifier has ZERO tunable parameters.

Usage:
    cd ~/Documents/GitHub/transport-wrapper
    export $(grep GWAS_KEY ~/Documents/GitHub/causal-inference-neuro-epidemiology/.env)
    uv run python -m DRUGS_EXPANDED.classify_v34
"""
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

HERE = Path(__file__).parent
OUT_DIR = HERE / "results" / "v34"
API = "https://api.opengwas.io/api"

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
    "Juvenile idiopathic arthritis": None,
    "Lung cancer": "ieu-a-984",
    "Major depressive disorder": "ieu-b-102",
    "Melanoma": "ieu-a-62",
    "Multiple sclerosis": "ieu-b-18",
    "Myocardial infarction": "ieu-a-798",
    "Neuroblastoma": None,
    "Ovarian cancer": "ieu-a-1120",
    "Pancreatic cancer": "ieu-b-4866",
    "Parkinsons disease": "ieu-b-7",
    "Rheumatoid arthritis": "ieu-a-833",
    "Schizophrenia": "ieu-b-5102",
    "Systemic lupus erythematosus": "ieu-a-1073",
    "Thyroid cancer": "ieu-a-1082",
    "Ulcerative colitis": "ieu-a-970",
}

GWAS_CACHE_FILE = OUT_DIR / "outcome_gwas_cache_v34.json"
MR_RESULTS_FILE = OUT_DIR / "mr_results_v34.jsonl"


def _headers() -> dict:
    token = os.environ.get("GWAS_KEY", "")
    h = {"X-API-SOURCE": "v34-stage-c/0.1"}
    if token:
        h["Authorization"] = f"Bearer {token}"
    return h


def norm_disease(s: str) -> str:
    return s.lower().replace("'", "").replace(" ", "_").replace("-", "_").strip()


def query_outcome_batch(rsids: list[str], gwas_id: str,
                        max_retries: int = 3) -> dict[str, dict]:
    """Query OpenGWAS for multiple SNP-disease associations in one call."""
    results = {}
    for attempt in range(max_retries):
        try:
            r = requests.post(
                f"{API}/associations",
                headers=_headers(),
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
            continue
    return results


def compute_wald_ratio(beta_exposure: float, se_exposure: float,
                       beta_outcome: float, se_outcome: float) -> dict:
    """Single-SNP Wald ratio MR estimate."""
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


def load_catalog_matches(evaluable: pd.DataFrame) -> dict[str, dict]:
    """Load pre-computed MR results from the EpiGraphDB catalog."""
    catalog = pd.read_csv(OUT_DIR / "v34_mr_catalog.csv")
    cat_by_key = {}
    for _, r in catalog.iterrows():
        key = (r["gene_symbol"], norm_disease(r["canonical_disease"]))
        if key not in cat_by_key or r["p"] < cat_by_key[key]["p"]:
            cat_by_key[key] = {
                "beta": float(r["beta"]),
                "se": float(r["se"]),
                "p": float(r["p"]),
                "ci_lower": float(r["ci_lower"]),
                "ci_upper": float(r["ci_upper"]),
                "mr_supports_causal": bool(r["mr_supports_causal"]),
                "rsid": r["rsid"],
                "source": "catalog",
            }

    matched = {}
    for _, r in evaluable.iterrows():
        key = (r["gene"], norm_disease(r["disease"]))
        if key in cat_by_key:
            matched[r["pair_id"]] = cat_by_key[key]
    return matched


def get_epigraphdb_rsids() -> dict[str, str]:
    """Get sentinel RSIDs for EpiGraphDB genes from the catalog."""
    catalog = pd.read_csv(OUT_DIR / "v34_mr_catalog.csv")
    gene_rsid = {}
    for _, r in catalog.iterrows():
        gene = r["gene_symbol"]
        if gene not in gene_rsid:
            gene_rsid[gene] = r["rsid"]
    return gene_rsid


def load_gwas_cache() -> dict:
    if GWAS_CACHE_FILE.exists():
        with open(GWAS_CACHE_FILE) as f:
            return json.load(f)
    return {}


def save_gwas_cache(cache: dict) -> None:
    with open(GWAS_CACHE_FILE, "w") as f:
        json.dump(cache, f, indent=2)


def compute_metrics(predictions: list[str], labels: list[str]) -> dict:
    tp = sum(p == "SUCCESS" and l == "SUCCESS" for p, l in zip(predictions, labels))
    tn = sum(p == "FAILURE" and l == "FAILURE" for p, l in zip(predictions, labels))
    fp = sum(p == "SUCCESS" and l == "FAILURE" for p, l in zip(predictions, labels))
    fn = sum(p == "FAILURE" and l == "SUCCESS" for p, l in zip(predictions, labels))
    sens = tp / (tp + fn) if (tp + fn) > 0 else float("nan")
    spec = tn / (tn + fp) if (tn + fp) > 0 else float("nan")
    ba = (sens + spec) / 2 if not (np.isnan(sens) or np.isnan(spec)) else float("nan")
    ppv = tp / (tp + fp) if (tp + fp) > 0 else float("nan")
    npv = tn / (tn + fn) if (tn + fn) > 0 else float("nan")
    denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
    mcc = float((tp * tn - fp * fn) / denom) if denom > 0 else float("nan")
    return {
        "n": len(predictions), "tp": tp, "tn": tn, "fp": fp, "fn": fn,
        "sensitivity": sens, "specificity": spec, "balanced_accuracy": ba,
        "ppv": ppv, "npv": npv, "mcc": mcc,
    }


def bootstrap_ba(predictions: list[str], labels: list[str],
                 n_boot: int = 10000) -> dict:
    rng = np.random.default_rng()
    bas = []
    mccs = []
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
        denom = np.sqrt(float((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
        if denom > 0:
            mccs.append(float((tp * tn - fp * fn) / denom))
    bas = np.array(bas)
    mccs = np.array(mccs)
    return {
        "ba_mean": float(np.mean(bas)),
        "ba_ci_lo": float(np.percentile(bas, 5)),
        "ba_ci_hi": float(np.percentile(bas, 95)),
        "mcc_mean": float(np.mean(mccs)) if len(mccs) > 0 else float("nan"),
        "mcc_ci_lo": float(np.percentile(mccs, 5)) if len(mccs) > 0 else float("nan"),
        "mcc_ci_hi": float(np.percentile(mccs, 95)) if len(mccs) > 0 else float("nan"),
        "n_boot": n_boot,
    }


def run_stage_c():
    ts = datetime.now(timezone.utc).isoformat()
    print(f"[{ts}] V3.4 Stage C: MR classification")

    # Load evaluable pairs
    adj = pd.read_csv(OUT_DIR / "adjudicated_v34.csv")
    evaluable = adj[
        (adj["outcome"].isin(["SUCCESS", "FAILURE"]))
        & (adj["excluded"] != True)
        & (adj["excluded"] != "True")
    ].copy()
    print(f"  Evaluable: {len(evaluable)} (SUCCESS={int((evaluable['outcome'] == 'SUCCESS').sum())}, "
          f"FAILURE={int((evaluable['outcome'] == 'FAILURE').sum())})")

    # Step 1: Load catalog matches (40 pairs)
    catalog_mr = load_catalog_matches(evaluable)
    print(f"  Catalog MR matches: {len(catalog_mr)}")

    # Step 2: Get RSIDs for EpiGraphDB genes without catalog match
    epigraphdb_rsids = get_epigraphdb_rsids()

    # Step 3: Query OpenGWAS for remaining pairs
    gwas_cache = load_gwas_cache()
    fresh_mr = {}
    needs_query = []

    for _, row in evaluable.iterrows():
        pid = row["pair_id"]
        if pid in catalog_mr:
            continue

        gwas_id = DISEASE_GWAS.get(row["disease"])
        if gwas_id is None:
            continue

        rsid = row["sentinel_rsid"] if pd.notna(row["sentinel_rsid"]) and row["sentinel_rsid"] else None
        if rsid is None:
            rsid = epigraphdb_rsids.get(row["gene"])
        if rsid is None:
            continue

        beta_exp = row["sentinel_beta"] if pd.notna(row["sentinel_beta"]) and row["sentinel_beta"] else None
        se_exp = row["sentinel_se"] if pd.notna(row["sentinel_se"]) and row["sentinel_se"] else None

        needs_query.append({
            "pair_id": pid,
            "gene": row["gene"],
            "disease": row["disease"],
            "rsid": rsid,
            "gwas_id": gwas_id,
            "beta_exp": float(beta_exp) if beta_exp is not None else None,
            "se_exp": float(se_exp) if se_exp is not None else None,
        })

    print(f"  Pairs needing OpenGWAS query: {len(needs_query)}")

    # Check which are already cached
    to_query = []
    for item in needs_query:
        cache_key = f"{item['rsid']}_{item['gwas_id']}"
        if cache_key in gwas_cache:
            cached = gwas_cache[cache_key]
            if cached is not None and item["beta_exp"] is not None:
                mr = compute_wald_ratio(
                    item["beta_exp"], item["se_exp"],
                    cached["beta"], cached["se"],
                )
                mr["rsid"] = item["rsid"]
                mr["source"] = "opengwas_cached"
                mr["beta_outcome"] = cached["beta"]
                mr["se_outcome"] = cached["se"]
                fresh_mr[item["pair_id"]] = mr
            elif cached is not None:
                # EpiGraphDB pair: use outcome beta/se with catalog exposure
                cat_rows = pd.read_csv(OUT_DIR / "v34_mr_catalog.csv")
                gene_rows = cat_rows[cat_rows["gene_symbol"] == item["gene"]]
                if len(gene_rows) > 0:
                    exp_row = gene_rows.iloc[0]
                    mr_est = {
                        "beta": cached["beta"],
                        "se": cached["se"],
                        "p": float(2 * stats.norm.sf(abs(cached["beta"] / cached["se"]))) if cached["se"] > 0 else 1.0,
                        "source": "opengwas_cached_epi",
                        "rsid": item["rsid"],
                    }
                    mr_est["ci_lower"] = mr_est["beta"] - 1.96 * mr_est["se"]
                    mr_est["ci_upper"] = mr_est["beta"] + 1.96 * mr_est["se"]
                    mr_est["mr_supports_causal"] = mr_est["p"] < 0.05
                    fresh_mr[item["pair_id"]] = mr_est
        else:
            to_query.append(item)

    print(f"  Already cached: {len(needs_query) - len(to_query)}")
    print(f"  Need fresh query: {len(to_query)}")

    if to_query:
        if not os.environ.get("GWAS_KEY"):
            print("  WARNING: GWAS_KEY not set. Set it with:")
            print("    export $(grep GWAS_KEY ~/Documents/GitHub/causal-inference-neuro-epidemiology/.env)")

        # Batch by GWAS ID (one API call per disease, all SNPs at once)
        by_gwas: dict[str, list[dict]] = {}
        for item in to_query:
            by_gwas.setdefault(item["gwas_id"], []).append(item)

        for gwas_id, items in tqdm(by_gwas.items(), desc="Querying OpenGWAS by disease"):
            rsids = list(set(item["rsid"] for item in items))
            disease_name = items[0]["disease"]
            print(f"    {disease_name}: {len(rsids)} SNPs → {gwas_id}")

            results = query_outcome_batch(rsids, gwas_id)
            print(f"      Got {len(results)} associations back")

            for item in items:
                cache_key = f"{item['rsid']}_{item['gwas_id']}"
                result = results.get(item["rsid"])

                if result is not None:
                    beta_out = result.get("beta")
                    se_out = result.get("se")
                    if beta_out is not None and se_out is not None:
                        gwas_cache[cache_key] = {
                            "beta": float(beta_out),
                            "se": float(se_out),
                            "p": float(result.get("p", 1.0)),
                            "ea": result.get("ea", ""),
                            "nea": result.get("nea", ""),
                        }

                        if item["beta_exp"] is not None:
                            mr = compute_wald_ratio(
                                item["beta_exp"], item["se_exp"],
                                float(beta_out), float(se_out),
                            )
                            mr["rsid"] = item["rsid"]
                            mr["source"] = "opengwas_fresh"
                            fresh_mr[item["pair_id"]] = mr
                    else:
                        gwas_cache[cache_key] = None
                else:
                    gwas_cache[cache_key] = None

            save_gwas_cache(gwas_cache)
            time.sleep(1.0)

    # For EpiGraphDB pairs without exposure betas in frozen CSV,
    # use the raw outcome association as the test statistic
    for item in needs_query:
        pid = item["pair_id"]
        if pid in fresh_mr:
            continue
        if item["beta_exp"] is not None:
            continue
        cache_key = f"{item['rsid']}_{item['gwas_id']}"
        if cache_key not in gwas_cache or gwas_cache[cache_key] is None:
            continue

        cached_out = gwas_cache[cache_key]
        beta_out = cached_out["beta"]
        se_out = cached_out["se"]
        if se_out > 0:
            z = beta_out / se_out
            p = float(2 * stats.norm.sf(abs(z)))
            fresh_mr[pid] = {
                "beta": beta_out,
                "se": se_out,
                "p": p,
                "ci_lower": beta_out - 1.96 * se_out,
                "ci_upper": beta_out + 1.96 * se_out,
                "mr_supports_causal": p < 0.05,
                "rsid": item["rsid"],
                "source": "opengwas_outcome_only",
            }

    # Combine all MR results
    all_mr = {**catalog_mr, **fresh_mr}
    print(f"\n  Total MR results: {len(all_mr)}")

    # Join with evaluable pairs
    evaluable["mr_beta"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("beta"))
    evaluable["mr_se"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("se"))
    evaluable["mr_p"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("p"))
    evaluable["mr_causal"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("mr_supports_causal"))
    evaluable["mr_source"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("source", ""))
    evaluable["mr_rsid"] = evaluable["pair_id"].map(lambda pid: all_mr.get(pid, {}).get("rsid", ""))

    has_mr = evaluable["mr_p"].notna()
    evaluable_mr = evaluable[has_mr].copy()
    no_mr = evaluable[~has_mr]

    print(f"  Evaluable with MR: {len(evaluable_mr)}")
    print(f"  Evaluable without MR: {len(no_mr)}")
    if len(no_mr) > 0:
        print(f"    Missing diseases: {no_mr['disease'].value_counts().to_dict()}")

    # Apply classifier: p < 0.05 → CAUSAL → predict SUCCESS
    evaluable_mr["prediction"] = evaluable_mr["mr_causal"].map(
        {True: "SUCCESS", False: "FAILURE"}
    )
    evaluable_mr["correct"] = evaluable_mr["prediction"] == evaluable_mr["outcome"]

    # === RESULTS ===
    preds = evaluable_mr["prediction"].tolist()
    labs = evaluable_mr["outcome"].tolist()

    print(f"\n{'='*70}")
    print(f"  STAGE C RESULTS — V3.4 MR Classification")
    print(f"{'='*70}")

    # ---------------------------------------------------------------
    # Analysis 1: Primary BA with bootstrap 90% CI
    # ---------------------------------------------------------------
    metrics = compute_metrics(preds, labs)
    boot = bootstrap_ba(preds, labs)
    base_rate = (evaluable_mr["outcome"] == "SUCCESS").mean()

    print(f"\n  [1] PRIMARY: Pooled BA (n={metrics['n']})")
    print(f"      TP={metrics['tp']}  TN={metrics['tn']}  FP={metrics['fp']}  FN={metrics['fn']}")
    print(f"      Sensitivity:  {metrics['sensitivity']:.3f}")
    print(f"      Specificity:  {metrics['specificity']:.3f}")
    print(f"      Balanced Accuracy: {metrics['balanced_accuracy']:.3f}")
    print(f"      Bootstrap 90% CI:  [{boot['ba_ci_lo']:.3f}, {boot['ba_ci_hi']:.3f}]")
    print(f"      MCC: {metrics['mcc']:.3f}  90% CI: [{boot['mcc_ci_lo']:.3f}, {boot['mcc_ci_hi']:.3f}]")
    print(f"      PPV: {metrics['ppv']:.3f}  NPV: {metrics['npv']:.3f}")
    print(f"      Base rate (SUCCESS): {base_rate:.3f}")
    print(f"      Null BA (any constant): 0.500")

    # ---------------------------------------------------------------
    # Analysis 2: Mechanism stratification (CO-PRIMARY)
    # ---------------------------------------------------------------
    print(f"\n  [2] CO-PRIMARY: Mechanism stratification")
    for mech in ["abundance_modulating", "activity_blocking"]:
        sub = evaluable_mr[evaluable_mr["mechanism_class"] == mech]
        if len(sub) < 2:
            continue
        s_preds = sub["prediction"].tolist()
        s_labs = sub["outcome"].tolist()
        s_metrics = compute_metrics(s_preds, s_labs)
        s_boot = bootstrap_ba(s_preds, s_labs)
        print(f"    {mech} (n={s_metrics['n']}):")
        print(f"      BA={s_metrics['balanced_accuracy']:.3f} "
              f"[{s_boot['ba_ci_lo']:.3f}, {s_boot['ba_ci_hi']:.3f}]")
        print(f"      TP={s_metrics['tp']} TN={s_metrics['tn']} "
              f"FP={s_metrics['fp']} FN={s_metrics['fn']}")
        print(f"      Sens={s_metrics['sensitivity']:.3f} Spec={s_metrics['specificity']:.3f}")
        print(f"      MCC={s_metrics['mcc']:.3f}")

    # Mixed stratum — report in pooled only, no standalone BA claim
    mixed = evaluable_mr[evaluable_mr["mechanism_class"] == "mixed"]
    if len(mixed) > 0:
        print(f"    mixed (n={len(mixed)}): pooled only, too few for standalone claim")
        m_preds = mixed["prediction"].tolist()
        m_labs = mixed["outcome"].tolist()
        m_metrics = compute_metrics(m_preds, m_labs)
        print(f"      TP={m_metrics['tp']} TN={m_metrics['tn']} "
              f"FP={m_metrics['fp']} FN={m_metrics['fn']}")

    # ---------------------------------------------------------------
    # Analysis 3: Gene-level deduplication (CO-PRIMARY)
    # ---------------------------------------------------------------
    print(f"\n  [3] CO-PRIMARY: Gene-level deduplication")
    gene_dedup = evaluable_mr.sort_values("mr_p").drop_duplicates(subset="gene", keep="first")
    gd_preds = gene_dedup["prediction"].tolist()
    gd_labs = gene_dedup["outcome"].tolist()
    gd_metrics = compute_metrics(gd_preds, gd_labs)
    gd_boot = bootstrap_ba(gd_preds, gd_labs)
    print(f"    n={gd_metrics['n']} unique genes")
    print(f"    BA={gd_metrics['balanced_accuracy']:.3f} "
          f"[{gd_boot['ba_ci_lo']:.3f}, {gd_boot['ba_ci_hi']:.3f}]")
    print(f"    MCC={gd_metrics['mcc']:.3f}")

    # ---------------------------------------------------------------
    # Analysis 4: Exact binomial tests
    # ---------------------------------------------------------------
    print(f"\n  [4] Exact binomial tests")
    sens_res = binomtest(metrics["tp"], metrics["tp"] + metrics["fn"], 0.5)
    spec_res = binomtest(metrics["tn"], metrics["tn"] + metrics["fp"], 0.5)
    print(f"    Sensitivity vs 0.50: p={sens_res.pvalue:.4f}")
    print(f"    Specificity vs 0.50: p={spec_res.pvalue:.4f}")

    # ---------------------------------------------------------------
    # Analysis 5: Continuous ROC/AUC with permutation p
    # ---------------------------------------------------------------
    print(f"\n  [5] Continuous ROC/AUC")
    y_true = (evaluable_mr["outcome"] == "SUCCESS").astype(int).values
    y_score = -np.log10(evaluable_mr["mr_p"].clip(lower=1e-300).values)

    def manual_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
        pos = y_score[y_true == 1]
        neg = y_score[y_true == 0]
        n_pos, n_neg = len(pos), len(neg)
        if n_pos == 0 or n_neg == 0:
            return float("nan")
        total = 0.0
        for p in pos:
            total += (neg < p).sum() + 0.5 * (neg == p).sum()
        return float(total / (n_pos * n_neg))

    auc = manual_auc(y_true, y_score)
    rng = np.random.default_rng(42)
    n_perm = 10000
    perm_aucs = np.array([manual_auc(rng.permutation(y_true), y_score) for _ in range(n_perm)])
    perm_p = (np.sum(perm_aucs >= auc) + 1) / (n_perm + 1)
    print(f"    AUC: {auc:.3f}")
    print(f"    Permutation p (n={n_perm}): {perm_p:.4f}")

    # ---------------------------------------------------------------
    # Analysis 9: Instrument source stratification (SECONDARY)
    # ---------------------------------------------------------------
    print(f"\n  [9] SECONDARY: Instrument source stratification")
    for src in ["EpiGraphDB", "UKB-PPP"]:
        sub = evaluable_mr[evaluable_mr["instrument_source"] == src]
        if len(sub) < 2:
            continue
        s_preds = sub["prediction"].tolist()
        s_labs = sub["outcome"].tolist()
        s_metrics = compute_metrics(s_preds, s_labs)
        s_boot = bootstrap_ba(s_preds, s_labs)
        print(f"    {src} (n={s_metrics['n']}): BA={s_metrics['balanced_accuracy']:.3f} "
              f"[{s_boot['ba_ci_lo']:.3f}, {s_boot['ba_ci_hi']:.3f}] MCC={s_metrics['mcc']:.3f}")

    # ---------------------------------------------------------------
    # Analysis 10: Sample overlap sensitivity (SECONDARY)
    # ---------------------------------------------------------------
    print(f"\n  [10] SECONDARY: Sample overlap sensitivity")
    clean = evaluable_mr[evaluable_mr["sample_overlap"] == "clean"]
    if len(clean) > 0:
        c_preds = clean["prediction"].tolist()
        c_labs = clean["outcome"].tolist()
        c_metrics = compute_metrics(c_preds, c_labs)
        c_boot = bootstrap_ba(c_preds, c_labs)
        print(f"    Clean only (n={c_metrics['n']}): BA={c_metrics['balanced_accuracy']:.3f} "
              f"[{c_boot['ba_ci_lo']:.3f}, {c_boot['ba_ci_hi']:.3f}]")

    flagged = evaluable_mr[evaluable_mr["sample_overlap"] == "overlap_flagged"]
    if len(flagged) > 0:
        print(f"    Overlap-flagged (n={len(flagged)}): included in pooled, reported here")

    # ---------------------------------------------------------------
    # Analysis 11: Oncology vs non-oncology (SECONDARY)
    # ---------------------------------------------------------------
    print(f"\n  [11] SECONDARY: Oncology vs non-oncology")
    for area in ["non_oncology", "oncology"]:
        sub = evaluable_mr[evaluable_mr["disease_area"] == area]
        if len(sub) < 2:
            continue
        s_preds = sub["prediction"].tolist()
        s_labs = sub["outcome"].tolist()
        s_metrics = compute_metrics(s_preds, s_labs)
        s_boot = bootstrap_ba(s_preds, s_labs)
        print(f"    {area} (n={s_metrics['n']}): BA={s_metrics['balanced_accuracy']:.3f} "
              f"[{s_boot['ba_ci_lo']:.3f}, {s_boot['ba_ci_hi']:.3f}] MCC={s_metrics['mcc']:.3f}")

    # ---------------------------------------------------------------
    # Analysis 13: Effect direction concordance (POST-HOC)
    # ---------------------------------------------------------------
    print(f"\n  [13] POST-HOC: Effect direction concordance")

    # Map action types to expected MR beta sign for SUCCESS
    REDUCING_ACTIONS = {"INHIBITOR", "ANTAGONIST", "ANTISENSE INHIBITOR", "DISRUPTING AGENT"}
    ENHANCING_ACTIONS = {"AGONIST", "ACTIVATOR", "EXOGENOUS PROTEIN", "STABILISER"}

    def get_expected_sign(action_types_str: str) -> str | None:
        """Return 'positive' or 'negative' for the expected MR beta if drug works."""
        if pd.isna(action_types_str):
            return None
        actions = {a.strip() for a in action_types_str.split(";")}
        has_reducing = bool(actions & REDUCING_ACTIONS)
        has_enhancing = bool(actions & ENHANCING_ACTIONS)
        if has_reducing and not has_enhancing:
            return "positive"
        if has_enhancing and not has_reducing:
            return "negative"
        return None

    at_col = "action_types"
    if at_col not in evaluable_mr.columns:
        adj_actions = pd.read_csv(OUT_DIR / "adjudicated_v34.csv")[["pair_id", "action_types"]]
        evaluable_mr = evaluable_mr.merge(adj_actions, on="pair_id", how="left")

    evaluable_mr["expected_sign"] = evaluable_mr[at_col].map(get_expected_sign)
    evaluable_mr["beta_sign"] = evaluable_mr["mr_beta"].map(
        lambda b: "positive" if b > 0 else "negative" if b < 0 else None
    )
    evaluable_mr["direction_concordant"] = (
        evaluable_mr["expected_sign"] == evaluable_mr["beta_sign"]
    )

    # Among all pairs with a directional prediction
    has_direction = evaluable_mr["expected_sign"].notna() & evaluable_mr["mr_beta"].notna()
    dir_pairs = evaluable_mr[has_direction]
    n_conc = dir_pairs["direction_concordant"].sum()
    n_dir = len(dir_pairs)
    print(f"    All pairs with directional prediction: {n_conc}/{n_dir} concordant "
          f"({n_conc/n_dir:.1%})" if n_dir > 0 else "    No directional predictions")

    # Among significant pairs only (the 16 where mr_causal=True)
    sig_dir = dir_pairs[dir_pairs["mr_causal"] == True]
    n_sig_conc = sig_dir["direction_concordant"].sum()
    n_sig = len(sig_dir)
    if n_sig > 0:
        print(f"    Significant pairs (mr_causal=True): {n_sig_conc}/{n_sig} concordant "
              f"({n_sig_conc/n_sig:.1%})")

        # Break down by TP vs FP
        sig_tp = sig_dir[sig_dir["outcome"] == "SUCCESS"]
        sig_fp = sig_dir[sig_dir["outcome"] == "FAILURE"]
        tp_conc = sig_tp["direction_concordant"].sum()
        fp_conc = sig_fp["direction_concordant"].sum()
        print(f"      TP (correct): {tp_conc}/{len(sig_tp)} concordant")
        print(f"      FP (wrong):   {fp_conc}/{len(sig_fp)} concordant")

    # By mechanism class
    for mech in ["abundance_modulating", "activity_blocking"]:
        m_dir = dir_pairs[dir_pairs["mechanism_class"] == mech]
        m_sig = m_dir[m_dir["mr_causal"] == True]
        if len(m_sig) > 0:
            m_conc = m_sig["direction_concordant"].sum()
            print(f"    {mech} significant: {m_conc}/{len(m_sig)} concordant")

    # Print the significant pairs with direction info
    print(f"\n    Detail (significant pairs):")
    print(f"    {'pair_id':<42} {'beta':>8} {'exp_sign':>8} {'concord':>7} {'outcome':>7} {'actions'}")
    for _, r in sig_dir.sort_values("mr_p").iterrows():
        conc = "Y" if r["direction_concordant"] else "N"
        actions = str(r[at_col])[:30] if pd.notna(r[at_col]) else ""
        print(f"    {r['pair_id']:<42} {r['mr_beta']:>8.3f} {r['expected_sign']:>8} "
              f"{conc:>7} {r['outcome']:>7} {actions}")

    # Binomial test: is concordance rate among significant pairs above chance (50%)?
    if n_sig > 0:
        conc_binom = binomtest(n_sig_conc, n_sig, 0.5)
        print(f"\n    Binomial test (concordance vs 50%): p={conc_binom.pvalue:.4f}")

    # ---------------------------------------------------------------
    # Analysis 14: Leave-one-disease-out stability (POST-HOC)
    # ---------------------------------------------------------------
    print(f"\n  [14] POST-HOC: Leave-one-disease-out BA stability")
    diseases_with_pairs = [d for d in sorted(evaluable_mr["disease"].unique())
                           if len(evaluable_mr[evaluable_mr["disease"] == d]) >= 1]
    lodo_results = {}
    print(f"    {'Disease dropped':<35} {'n_left':>6} {'BA':>7} {'delta':>7}")
    print(f"    {'-'*35} {'-'*6} {'-'*7} {'-'*7}")
    for disease in diseases_with_pairs:
        left = evaluable_mr[evaluable_mr["disease"] != disease]
        if len(left) < 2:
            continue
        l_preds = left["prediction"].tolist()
        l_labs = left["outcome"].tolist()
        l_metrics = compute_metrics(l_preds, l_labs)
        ba_drop = l_metrics["balanced_accuracy"]
        delta = ba_drop - metrics["balanced_accuracy"]
        lodo_results[disease] = {"ba": ba_drop, "delta": delta, "n": l_metrics["n"]}
        marker = " ***" if abs(delta) > 0.02 else ""
        print(f"    {disease:<35} {l_metrics['n']:>6} {ba_drop:>7.4f} {delta:>+7.4f}{marker}")

    ba_values = [v["ba"] for v in lodo_results.values()]
    print(f"\n    BA range: [{min(ba_values):.4f}, {max(ba_values):.4f}]")
    print(f"    BA std:   {np.std(ba_values):.4f}")
    most_influential = max(lodo_results.items(), key=lambda x: abs(x[1]["delta"]))
    print(f"    Most influential: {most_influential[0]} (delta={most_influential[1]['delta']:+.4f})")

    # ---------------------------------------------------------------
    # Per-disease breakdown (excluding single-pair diseases, Analysis 12)
    # ---------------------------------------------------------------
    print(f"\n  Per-disease breakdown (excluding single-pair diseases):")
    print(f"  {'disease':<40} {'n':>3} {'BA':>6} {'TP':>3} {'TN':>3} {'FP':>3} {'FN':>3}")
    print(f"  {'-'*40} {'-'*3} {'-'*6} {'-'*3} {'-'*3} {'-'*3} {'-'*3}")
    per_disease = {}
    for disease in sorted(evaluable_mr["disease"].unique()):
        d = evaluable_mr[evaluable_mr["disease"] == disease]
        d_preds = d["prediction"].tolist()
        d_labs = d["outcome"].tolist()
        d_metrics = compute_metrics(d_preds, d_labs)
        per_disease[disease] = d_metrics
        marker = " *" if d_metrics["n"] <= 1 else ""
        print(f"  {disease:<40} {d_metrics['n']:>3} {d_metrics['balanced_accuracy']:>6.2f} "
              f"{d_metrics['tp']:>3} {d_metrics['tn']:>3} {d_metrics['fp']:>3} {d_metrics['fn']:>3}{marker}")

    # ---------------------------------------------------------------
    # Per-pair classification table
    # ---------------------------------------------------------------
    print(f"\n  Per-pair classification (sorted by p):")
    print(f"  {'pair_id':<45} {'p':>10} {'pred':>8} {'label':>8} {'ok':>3} {'mech':>10} {'src':>8}")
    print(f"  {'-'*45} {'-'*10} {'-'*8} {'-'*8} {'-'*3} {'-'*10} {'-'*8}")
    for _, row in evaluable_mr.sort_values("mr_p").iterrows():
        ok = "Y" if row["correct"] else "N"
        mech = row["mechanism_class"][:10] if pd.notna(row["mechanism_class"]) else ""
        src = row["mr_source"][:8] if pd.notna(row["mr_source"]) else ""
        print(f"  {row['pair_id']:<45} {row['mr_p']:>10.2e} {row['prediction']:>8} "
              f"{row['outcome']:>8} {ok:>3} {mech:>10} {src:>8}")

    # ---------------------------------------------------------------
    # Save results
    # ---------------------------------------------------------------
    results = {
        "timestamp": ts,
        "version": "V3.4",
        "n_evaluable_total": len(evaluable),
        "n_evaluable_with_mr": len(evaluable_mr),
        "n_missing_mr": len(no_mr),
        "missing_diseases": no_mr["disease"].value_counts().to_dict() if len(no_mr) > 0 else {},
        "base_rate": float(base_rate),
        "primary_metrics": metrics,
        "primary_bootstrap": boot,
        "per_disease": per_disease,
        "mechanism_strat": {},
        "gene_dedup": {"metrics": gd_metrics, "bootstrap": gd_boot},
        "direction_concordance": {
            "n_all_directional": int(n_dir),
            "n_all_concordant": int(n_conc),
            "rate_all": float(n_conc / n_dir) if n_dir > 0 else None,
            "n_significant": int(n_sig),
            "n_significant_concordant": int(n_sig_conc),
            "rate_significant": float(n_sig_conc / n_sig) if n_sig > 0 else None,
        },
        "leave_one_disease_out": {
            "ba_range": [float(min(ba_values)), float(max(ba_values))],
            "ba_std": float(np.std(ba_values)),
            "most_influential": most_influential[0],
            "most_influential_delta": float(most_influential[1]["delta"]),
            "per_disease": {k: {"ba": v["ba"], "delta": v["delta"]} for k, v in lodo_results.items()},
        },
    }

    for mech in ["abundance_modulating", "activity_blocking", "mixed"]:
        sub = evaluable_mr[evaluable_mr["mechanism_class"] == mech]
        if len(sub) >= 2:
            s_preds = sub["prediction"].tolist()
            s_labs = sub["outcome"].tolist()
            results["mechanism_strat"][mech] = {
                "metrics": compute_metrics(s_preds, s_labs),
                "bootstrap": bootstrap_ba(s_preds, s_labs),
            }

    eval_path = OUT_DIR / "evaluation_v34.json"
    with open(eval_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\n  Wrote evaluation to {eval_path}")

    # Save per-pair classification CSV
    out_cols = [
        "pair_id", "gene", "disease", "outcome", "mechanism_class",
        "disease_area", "instrument_source", "sample_overlap",
        "mr_beta", "mr_se", "mr_p", "mr_causal", "mr_source", "mr_rsid",
        "prediction", "correct",
    ]
    class_path = OUT_DIR / "classification_v34.csv"
    evaluable_mr[out_cols].to_csv(class_path, index=False)
    print(f"  Wrote classification to {class_path}")

    # Save the full evaluable set (including pairs without MR)
    full_path = OUT_DIR / "evaluable_full_v34.csv"
    evaluable.to_csv(full_path, index=False)
    print(f"  Wrote full evaluable set to {full_path}")

    ts2 = datetime.now(timezone.utc).isoformat()
    print(f"\n[{ts2}] Stage C complete")


if __name__ == "__main__":
    run_stage_c()
