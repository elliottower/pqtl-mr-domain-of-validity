"""Query Genebass LoF burden for all V3.4 gene-disease pairs (v2, fixed).

Fixes from Perplexity review:
- One canonical phenocode per disease, chosen by N not p-value
- Exact phenocode matching, no substring
- icd_first_occurrence handled correctly
- Continuous signal collected (beta, z-score), not just thresholded
- NO classification — instrument collection only
"""

import csv
import gzip
import json
import math
import time
import urllib.request

CLASSIFICATION = "data/classification_v34.csv"
CANONICAL_MAP = "data/canonical_disease_mapping.json"
ENSEMBL_CACHE = "data/ensembl_cache.json"
GENEBASS_CACHE = "data/genebass_raw_cache.json"
OUTPUT = "results/genebass_lof_results.json"

ENSEMBL_REST = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
GENEBASS_API = "https://main.genebass.org/api/phewas/{ensembl_id}"


def load_json(path):
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def save_json(data, path):
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def get_ensembl_id(gene_symbol: str, cache: dict) -> str | None:
    if gene_symbol in cache:
        return cache[gene_symbol]
    url = ENSEMBL_REST.format(gene=gene_symbol)
    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read())
            ens_id = data.get("id")
            cache[gene_symbol] = ens_id
            return ens_id
    except Exception as e:
        print(f"  Ensembl lookup failed for {gene_symbol}: {e}")
        cache[gene_symbol] = None
        return None


def query_genebass(ensembl_id: str, cache: dict) -> list[dict]:
    if ensembl_id in cache:
        return cache[ensembl_id]
    url = GENEBASS_API.format(ensembl_id=ensembl_id)
    try:
        req = urllib.request.Request(url)
        req.add_header("Accept-Encoding", "gzip")
        with urllib.request.urlopen(req, timeout=60) as resp:
            raw = resp.read()
            if raw[:2] == b'\x1f\x8b':
                raw = gzip.decompress(raw)
            data = json.loads(raw)
            results = data.get("phewas", []) if isinstance(data, dict) else data
            plof = [r for r in results if r.get("annotation") == "pLoF"]
            cache[ensembl_id] = plof
            return plof
    except Exception as e:
        print(f"  Genebass query failed for {ensembl_id}: {e}")
        cache[ensembl_id] = []
        return []


def match_canonical(plof_results: list[dict], disease: str, canonical_map: dict) -> dict | None:
    mapping = canonical_map.get(disease)
    if not mapping or mapping.get("phenocode") is None:
        return None

    target_phenocode = mapping["phenocode"]
    target_trait_type = mapping["trait_type"]
    target_coding = mapping.get("coding", "")

    for r in plof_results:
        phenocode_match = r.get("phenocode", "") == target_phenocode
        trait_match = r.get("trait_type", "") == target_trait_type
        coding_match = r.get("coding", "") == target_coding

        if phenocode_match and trait_match and coding_match:
            return r

    return None


def compute_z_score(beta: float, p: float) -> float | None:
    if p is None or p <= 0 or p >= 1 or beta is None:
        return None
    from scipy.stats import norm
    z_unsigned = norm.ppf(1 - p / 2)
    return z_unsigned if beta >= 0 else -z_unsigned


def main():
    with open(CLASSIFICATION) as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    with open(CANONICAL_MAP) as f:
        canonical_map = json.load(f)

    ens_cache = load_json(ENSEMBL_CACHE)
    gb_cache = load_json(GENEBASS_CACHE)

    genes = sorted(set(r["gene"] for r in rows))
    diseases = sorted(set(r["disease"] for r in rows))
    print(f"Genes: {len(genes)}, Diseases: {len(diseases)}")
    print(f"Gene-disease pairs in V3.4: {len(rows)}")

    mappable = {d for d in diseases if canonical_map.get(d, {}).get("phenocode") is not None}
    unmappable = sorted(set(diseases) - mappable)
    print(f"Mappable diseases: {len(mappable)}")
    print(f"Unmappable (structural missingness): {unmappable}")

    genes_queried = set()
    for row in rows:
        gene = row["gene"]
        if gene not in genes_queried:
            print(f"[{len(genes_queried)+1}/{len(genes)}] {gene}...", end=" ", flush=True)
            ens_id = get_ensembl_id(gene, ens_cache)
            if ens_id and ens_id not in gb_cache:
                time.sleep(0.3)
                query_genebass(ens_id, gb_cache)
            n_plof = len(gb_cache.get(ens_id, [])) if ens_id else 0
            print(f"pLoF={n_plof}")
            genes_queried.add(gene)

            if len(genes_queried) % 10 == 0:
                save_json(ens_cache, ENSEMBL_CACHE)
                save_json(gb_cache, GENEBASS_CACHE)

    save_json(ens_cache, ENSEMBL_CACHE)
    save_json(gb_cache, GENEBASS_CACHE)

    pair_results = {}
    for row in rows:
        gene = row["gene"]
        disease = row["disease"]
        pair_id = row["pair_id"]
        ens_id = ens_cache.get(gene)

        plof_all = gb_cache.get(ens_id, []) if ens_id else []
        match = match_canonical(plof_all, disease, canonical_map)

        lof_beta = float(match["BETA_Burden"]) if match and match.get("BETA_Burden") is not None else None
        lof_p = float(match["Pvalue_Burden"]) if match and match.get("Pvalue_Burden") is not None else None
        lof_z = compute_z_score(lof_beta, lof_p) if lof_beta is not None and lof_p is not None else None

        pair_results[pair_id] = {
            "gene": gene,
            "disease": disease,
            "ensembl_id": ens_id,
            "mechanism_class": row["mechanism_class"],
            "v34_outcome": row["outcome"],
            "v34_mr_beta": float(row["mr_beta"]) if row.get("mr_beta") else None,
            "v34_mr_p": float(row["mr_p"]) if row.get("mr_p") else None,
            "lof_phenocode": match.get("phenocode") if match else None,
            "lof_trait_type": match.get("trait_type") if match else None,
            "lof_total_variants": match.get("total_variants") if match else None,
            "lof_beta_burden": lof_beta,
            "lof_p_burden": lof_p,
            "lof_z_burden": lof_z,
            "lof_missing_reason": (
                canonical_map.get(disease, {}).get("rationale", "unknown")
                if match is None and canonical_map.get(disease, {}).get("phenocode") is None
                else ("no_genebass_result" if match is None else None)
            ),
        }

    has_lof = sum(1 for r in pair_results.values() if r["lof_p_burden"] is not None)
    has_sig = sum(1 for r in pair_results.values() if r["lof_p_burden"] is not None and r["lof_p_burden"] < 0.05)
    structural_miss = sum(1 for r in pair_results.values() if r["lof_missing_reason"] and "structural" in r["lof_missing_reason"].lower())

    print(f"\n=== SUMMARY (v2, canonical phenocodes) ===")
    print(f"Total V3.4 pairs: {len(pair_results)}")
    print(f"Pairs with LoF burden result: {has_lof}")
    print(f"Pairs with sig LoF (p<0.05): {has_sig}")
    print(f"Structural missingness (pediatric/absent): {structural_miss}")
    print(f"Other missing (no Genebass result): {len(pair_results) - has_lof - structural_miss}")

    by_mechanism = {}
    for r in pair_results.values():
        mech = r["mechanism_class"]
        by_mechanism.setdefault(mech, {"total": 0, "has_lof": 0, "sig_lof": 0, "missing": 0})
        by_mechanism[mech]["total"] += 1
        if r["lof_p_burden"] is not None:
            by_mechanism[mech]["has_lof"] += 1
            if r["lof_p_burden"] < 0.05:
                by_mechanism[mech]["sig_lof"] += 1
        else:
            by_mechanism[mech]["missing"] += 1

    print(f"\nBy mechanism:")
    for mech, c in sorted(by_mechanism.items()):
        print(f"  {mech:25s}: {c['has_lof']}/{c['total']} with LoF ({c['missing']} missing), {c['sig_lof']} sig")

    if has_lof > 0:
        betas = [r["lof_beta_burden"] for r in pair_results.values() if r["lof_beta_burden"] is not None]
        zs = [r["lof_z_burden"] for r in pair_results.values() if r["lof_z_burden"] is not None]
        print(f"\nContinuous signal stats (for pre-registration, NOT classification):")
        print(f"  Beta range: [{min(betas):.4f}, {max(betas):.4f}]")
        print(f"  Z-score range: [{min(zs):.2f}, {max(zs):.2f}]")
        print(f"  Pairs with |z| > 1.96: {sum(1 for z in zs if abs(z) > 1.96)}")

    by_disease_missing = {}
    for r in pair_results.values():
        if r["lof_p_burden"] is None:
            by_disease_missing.setdefault(r["disease"], []).append(r["lof_missing_reason"])
    if by_disease_missing:
        print(f"\nMissing by disease:")
        for d, reasons in sorted(by_disease_missing.items()):
            print(f"  {d}: {len(reasons)} pairs — {reasons[0][:80]}")

    save_json(pair_results, OUTPUT)
    print(f"\nSaved to {OUTPUT}")


if __name__ == "__main__":
    main()
