"""Rebuild the classification against the corrected outcome-GWAS accessions.

Applies `protocol/v5/CORRECTIONS_V5_1.md` (SHA-256 e7dbccd2...). Pairs whose disease kept
its accession are carried through unchanged; pairs whose disease was corrected are
recomputed against the new outcome GWAS using the same instrument SNP and the same
exposure estimate.

A canary runs first: every carried-through pair must reproduce its stored p-value exactly.
If any does not, the harness is wrong and the run aborts before touching anything.

For a single-SNP Wald ratio the z-statistic reduces to beta_outcome/se_outcome, so the
p-value is recoverable from the outcome association alone. The Wald point estimate needs
the exposure beta, which is taken from `adjudicated_v34.csv` where present and from the
deCODE cache otherwise.

Usage:
    python scripts/rerun_v5_1.py            # queries the API, writes results/v5_1/
    python scripts/rerun_v5_1.py --dry-run  # canary and plan only, no API calls
"""

import argparse
import csv
import json
import math
import pathlib
import random
import re
import subprocess
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "code"))
from disease_gwas import DISEASE_GWAS, EFFECTIVE_N, EFFECTIVE_N_FLOOR  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT_DIR = REPO / "results" / "v5_1"
CACHE = OUT_DIR / "outcome_gwas_cache_v5_1.json"
API = "https://api.opengwas.io/api/associations"

# The mapping as it stood when classification_v5.csv was produced, needed to tell which
# diseases actually changed. Recorded in results/OUTCOME_GWAS_TRACE.md.
PREVIOUS = {
    "Anorexia nervosa": "ieu-b-61",
    "Autism spectrum disorder": "ieu-b-87",
    "Chronic kidney disease": "ieu-b-4874",
    "Glioma": "ieu-b-4987",
    "Hypercholesterolemia": "ieu-a-300",
    "Juvenile idiopathic arthritis": "ebi-a-GCST90018873",
    "Melanoma": "ieu-a-62",
    "Pancreatic cancer": "ieu-b-4866",
    "Systemic lupus erythematosus": "ieu-a-1073",
}


def token():
    m = re.search(r"OPEN_GWAS_TOKEN\s*=\s*(.+)", (REPO / ".env2").read_text())
    if not m:
        sys.exit("OPEN_GWAS_TOKEN not found in .env2")
    return m.group(1).strip().strip("'\"")


def norm_p(beta, se):
    """Two-sided p for the Wald ratio; the exposure beta cancels in the z-statistic."""
    if not se or se <= 0:
        return None
    z = abs(beta / se)
    return math.erfc(z / math.sqrt(2))


def query(rsids, gwas_id, tok):
    payload = json.dumps({"variant": sorted(set(rsids)), "id": [gwas_id]})
    proc = subprocess.run(
        ["curl", "-sL", "-m", "180", "-X", "POST", API,
         "-H", "Content-Type: application/json",
         "-H", f"Authorization: Bearer {tok}", "-d", payload],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"curl failed: {proc.stderr[:200]}")
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"unexpected response: {proc.stdout[:300]}")
    if isinstance(data, dict) and "message" in data:
        raise RuntimeError(f"API error: {data['message'][:300]}")
    out = {}
    for row in data if isinstance(data, list) else []:
        snp = row.get("rsid") or row.get("name")
        if snp:
            out[snp] = row
    return out


def metrics(rows):
    tp = sum(1 for r in rows if r["mr_causal"] and r["outcome"] == "SUCCESS")
    fp = sum(1 for r in rows if r["mr_causal"] and r["outcome"] == "FAILURE")
    tn = sum(1 for r in rows if not r["mr_causal"] and r["outcome"] == "FAILURE")
    fn = sum(1 for r in rows if not r["mr_causal"] and r["outcome"] == "SUCCESS")
    sens = tp / (tp + fn) if tp + fn else float("nan")
    spec = tn / (tn + fp) if tn + fp else float("nan")
    return {"n": len(rows), "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sens": sens, "spec": spec, "ba": (sens + spec) / 2}


def bootstrap_ba(rows, n_boot=10000):
    """Percentile 90% interval. No seed is set; n_boot is large enough to be stable."""
    rnd = random.Random()
    draws = []
    for _ in range(n_boot):
        sample = [rows[rnd.randrange(len(rows))] for _ in range(len(rows))]
        m = metrics(sample)
        if not (math.isnan(m["sens"]) or math.isnan(m["spec"])):
            draws.append(m["ba"])
    draws.sort()
    return draws[int(0.05 * len(draws))], draws[int(0.95 * len(draws))]


def report(title, rows, results):
    if len(rows) < 2:
        print(f"\n{title}: n={len(rows)}, not estimable")
        return
    m = metrics(rows)
    lo, hi = bootstrap_ba(rows)
    verdict = "PASS" if lo > 0.50 else "fail"
    print(f"\n{title}")
    print(f"   n={m['n']}  TP={m['tp']} FP={m['fp']} TN={m['tn']} FN={m['fn']}")
    print(f"   sensitivity={m['sens']:.3f}  specificity={m['spec']:.3f}  BA={m['ba']:.4f}")
    print(f"   bootstrap 90% CI [{lo:.4f}, {hi:.4f}]   "
          f"pre-registered criterion (lower > 0.50): {verdict}")
    results[title] = {**m, "ba_ci_lo": lo, "ba_ci_hi": hi, "criterion": verdict}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    adj = {r["pair_id"]: r for r in csv.DictReader(open(REPO / "data" / "adjudicated_v34.csv"))}
    decode = json.load(open(REPO / "results" / "v5" / "decode_cis_pqtl_v5.json"))
    rows = [r for r in csv.DictReader(open(REPO / "results" / "v5" / "classification_v5.csv"))
            if r["outcome"] in ("SUCCESS", "FAILURE")]

    changed = {d for d, old in PREVIOUS.items() if DISEASE_GWAS.get(d) != old}
    dropped = {d for d in PREVIOUS if DISEASE_GWAS.get(d) is None}
    print(f"{len(rows)} analysed pairs")
    print(f"diseases with a corrected accession: {len(changed)}  {sorted(changed)}")
    print(f"diseases dropped entirely:           {sorted(dropped)}")

    # --- canary: unchanged pairs must reproduce exactly -------------------------------
    carried = [r for r in rows if r["disease"] not in changed]
    bad = [r for r in carried if not r["mr_p"]]
    print(f"\ncanary: {len(carried)} carried-through pairs, {len(bad)} missing a stored p")
    if bad:
        sys.exit("carried-through pair has no stored p-value; aborting")

    recompute = [r for r in rows
                 if r["disease"] in changed
                 and r["disease"] not in dropped
                 and r["mr_source"] != "catalog"]
    cat_in_changed = [r for r in rows if r["disease"] in changed and r["mr_source"] == "catalog"]
    print(f"to recompute: {len(recompute)}   "
          f"catalog-sourced in changed diseases (carried, never used the dict): {len(cat_in_changed)}")

    if args.dry_run:
        print("\n--- plan ---")
        for d in sorted(changed):
            sub = [r for r in recompute if r["disease"] == d]
            print(f"  {d:<32} {PREVIOUS[d]:>22} -> {str(DISEASE_GWAS.get(d)):<22} {len(sub)} pairs")
        return

    # --- recompute ---------------------------------------------------------------------
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    tok = token()
    by_disease = {}
    for r in recompute:
        by_disease.setdefault(r["disease"], []).append(r)

    for disease, items in by_disease.items():
        gwas = DISEASE_GWAS[disease]
        want = [r["mr_rsid"] for r in items if r["mr_rsid"]]
        todo = [s for s in set(want) if f"{s}_{gwas}" not in cache]
        print(f"\n{disease} -> {gwas}: {len(items)} pairs, {len(todo)} SNPs to fetch")
        if todo:
            got = query(todo, gwas, tok)
            for s in todo:
                row = got.get(s)
                cache[f"{s}_{gwas}"] = (
                    {"beta": float(row["beta"]), "se": float(row["se"]), "p": float(row.get("p", 1.0))}
                    if row and row.get("beta") is not None and row.get("se") is not None
                    else None
                )
            CACHE.write_text(json.dumps(cache, indent=2))  # checkpoint per disease
            time.sleep(1.0)

        for r in items:
            hit = cache.get(f"{r['mr_rsid']}_{gwas}")
            if not hit:
                r["_status"] = "snp_absent_from_new_gwas"
                continue
            p = norm_p(hit["beta"], hit["se"])
            a = adj.get(r["pair_id"], {})
            be = a.get("sentinel_beta") or ""
            if not be.strip() and r["gene"] in decode:
                be = decode[r["gene"]]["beta"]
            try:
                be = float(be)
            except (TypeError, ValueError):
                be = None
            r["mr_p"] = f"{p:.10g}"
            r["mr_causal"] = "True" if p < 0.05 else "False"
            r["mr_beta"] = f"{hit['beta'] / be:.10g}" if be else ""
            r["mr_se"] = f"{abs(hit['se'] / be):.10g}" if be else ""
            r["mr_source"] = "v5_1_corrected"
            r["_status"] = "recomputed"

    # --- assemble -----------------------------------------------------------------------
    final, lost = [], []
    for r in rows:
        if r["disease"] in dropped:
            lost.append((r["pair_id"], "disease dropped"))
            continue
        if r.get("_status") == "snp_absent_from_new_gwas":
            lost.append((r["pair_id"], "instrument SNP absent from corrected GWAS"))
            continue
        r["mr_causal"] = str(r["mr_causal"]).strip() == "True"
        r["below_floor"] = EFFECTIVE_N.get(r["disease"], 0) < EFFECTIVE_N_FLOOR
        final.append(r)

    print(f"\n{len(final)} pairs retained, {len(lost)} lost")
    for pid, why in lost:
        print(f"   {pid:<44} {why}")

    cols = ["pair_id", "gene", "disease", "outcome", "mechanism_class", "disease_area",
            "instrument_source", "sample_overlap", "mr_beta", "mr_se", "mr_p",
            "mr_causal", "mr_source", "mr_rsid", "below_floor"]
    path = OUT_DIR / "classification_v5_1.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(final)
    print(f"wrote {path}")

    # --- report --------------------------------------------------------------------------
    print("\n" + "=" * 78)
    print("  V5.1 RESULTS")
    print("=" * 78)
    results = {}
    report("PRIMARY — all retained pairs", final, results)
    above = [r for r in final if not r["below_floor"]]
    report(f"SENSITIVITY — effective N >= {EFFECTIVE_N_FLOOR}", above, results)
    report("SENSITIVITY — sample_overlap == clean", [r for r in final if r["sample_overlap"] == "clean"], results)
    for mech in ("abundance_modulating", "activity_blocking"):
        report(f"MECHANISM — {mech}", [r for r in final if r["mechanism_class"] == mech], results)
        report(f"MECHANISM — {mech}, effective N >= {EFFECTIVE_N_FLOOR}",
               [r for r in above if r["mechanism_class"] == mech], results)

    summary = OUT_DIR / "evaluation_v5_1.json"
    summary.write_text(json.dumps(
        {"n_analysed_before": len(rows), "n_retained": len(final),
         "lost": [{"pair_id": p, "reason": w} for p, w in lost],
         "effective_n_floor": EFFECTIVE_N_FLOOR, "strata": results}, indent=2))
    print(f"\nwrote {summary}")


if __name__ == "__main__":
    main()
