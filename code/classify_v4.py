"""V4 Classification: LoF Burden + Evidence Routing Analyses.

Runs all 6 pre-registered analyses from PRESPEC_V4.md.
Uses frozen instrument data from genebass_lof_results_v2.json.
Outcomes carried from V3.4 adjudication.
"""

import json
import math
import random

from scipy.stats import binomtest, norm, mannwhitneyu

LOF_RESULTS = "results/genebass_lof_results.json"
V34_EVAL = "results/results_v4.json"
OUTPUT = "results/results_v4.json"

random.seed(None)


def load_pairs():
    with open(LOF_RESULTS) as f:
        all_pairs = json.load(f)
    has_lof = {k: v for k, v in all_pairs.items() if v["lof_p_burden"] is not None}
    return all_pairs, has_lof


def signed_z(beta, p):
    if beta is None or p is None or p <= 0 or p >= 1:
        return None
    z = norm.ppf(1 - p / 2)
    return z if beta >= 0 else -z


def pqtl_z(pair):
    beta = pair.get("v34_mr_beta")
    p = pair.get("v34_mr_p")
    if beta is None or p is None:
        return None
    return signed_z(beta, p)


def analysis_1_lof_auc(pairs):
    """Primary: LoF Burden AUC with permutation test."""
    print("=== Analysis 1: LoF Burden AUC ===")

    outcomes = []
    scores = []
    for r in pairs.values():
        z = r.get("lof_z_burden")
        if z is None:
            continue
        outcome = 1 if r["v34_outcome"] == "SUCCESS" else 0
        outcomes.append(outcome)
        scores.append(z)

    n = len(outcomes)
    n_success = sum(outcomes)
    n_failure = n - n_success
    print(f"  N={n}, SUCCESS={n_success}, FAILURE={n_failure}")

    stat, p_mw = mannwhitneyu(
        [s for s, o in zip(scores, outcomes) if o == 1],
        [s for s, o in zip(scores, outcomes) if o == 0],
        alternative="two-sided"
    )
    auc = stat / (n_success * n_failure)
    print(f"  AUC = {auc:.3f} (Mann-Whitney p = {p_mw:.4f})")

    n_perm = 10000
    perm_aucs = []
    for _ in range(n_perm):
        perm_outcomes = outcomes[:]
        random.shuffle(perm_outcomes)
        s1 = [s for s, o in zip(scores, perm_outcomes) if o == 1]
        s0 = [s for s, o in zip(scores, perm_outcomes) if o == 0]
        if len(s1) == 0 or len(s0) == 0:
            continue
        u, _ = mannwhitneyu(s1, s0, alternative="two-sided")
        perm_aucs.append(u / (len(s1) * len(s0)))

    perm_p = sum(1 for a in perm_aucs if abs(a - 0.5) >= abs(auc - 0.5)) / len(perm_aucs)
    print(f"  Permutation p = {perm_p:.4f} (n_perm={n_perm})")

    return {
        "n": n, "n_success": n_success, "n_failure": n_failure,
        "auc": auc, "mannwhitney_p": p_mw,
        "permutation_p": perm_p, "n_permutations": n_perm
    }


def analysis_2_concordance(pairs):
    """Primary: Directional concordance between pQTL and LoF."""
    print("\n=== Analysis 2: Directional Concordance (pQTL vs LoF) ===")

    concordant = 0
    discordant = 0
    by_mechanism = {}
    by_outcome = {}

    lof_neg = 0
    lof_pos = 0

    for r in pairs.values():
        lof_z = r.get("lof_z_burden")
        mr_beta = r.get("v34_mr_beta")
        if lof_z is None or mr_beta is None or lof_z == 0 or mr_beta == 0:
            continue

        lof_sign = 1 if lof_z > 0 else -1
        mr_sign = 1 if mr_beta > 0 else -1
        agree = lof_sign == mr_sign

        if lof_z > 0:
            lof_pos += 1
        else:
            lof_neg += 1

        if agree:
            concordant += 1
        else:
            discordant += 1

        mech = r["mechanism_class"]
        by_mechanism.setdefault(mech, {"concordant": 0, "discordant": 0})
        by_mechanism[mech]["concordant" if agree else "discordant"] += 1

        outcome = r["v34_outcome"]
        by_outcome.setdefault(outcome, {"concordant": 0, "discordant": 0})
        by_outcome[outcome]["concordant" if agree else "discordant"] += 1

    total = concordant + discordant
    rate = concordant / total if total > 0 else 0
    lof_base_rate = lof_pos / (lof_pos + lof_neg) if (lof_pos + lof_neg) > 0 else 0.5

    binom = binomtest(concordant, total, 0.5, alternative="two-sided")

    print(f"  Total pairs with both instruments: {total}")
    print(f"  Concordant: {concordant} ({rate:.1%})")
    print(f"  Discordant: {discordant}")
    print(f"  Binomial p (vs 50%): {binom.pvalue:.4f}")
    print(f"  LoF directional base rate: {lof_pos} pos / {lof_neg} neg ({lof_base_rate:.1%} positive)")

    binom_lof_base = binomtest(concordant, total, lof_base_rate, alternative="two-sided")
    print(f"  Binomial p (vs LoF base rate {lof_base_rate:.1%}): {binom_lof_base.pvalue:.4f}")

    print(f"\n  By mechanism:")
    for mech, c in sorted(by_mechanism.items()):
        t = c["concordant"] + c["discordant"]
        r = c["concordant"] / t if t > 0 else 0
        print(f"    {mech:25s}: {c['concordant']}/{t} ({r:.1%})")

    print(f"\n  By outcome:")
    for outcome, c in sorted(by_outcome.items()):
        t = c["concordant"] + c["discordant"]
        r = c["concordant"] / t if t > 0 else 0
        print(f"    {outcome:10s}: {c['concordant']}/{t} ({r:.1%})")

    return {
        "total": total, "concordant": concordant, "discordant": discordant,
        "concordance_rate": rate,
        "binomial_p_vs_50pct": binom.pvalue,
        "lof_directional_base_rate": lof_base_rate,
        "binomial_p_vs_lof_base": binom_lof_base.pvalue,
        "by_mechanism": by_mechanism,
        "by_outcome": by_outcome,
    }


def analysis_3_depth2(pairs):
    """Primary: Depth-2 licensing."""
    print("\n=== Analysis 3: Depth-2 Licensing ===")

    depth2 = []
    all_with_both = []

    for pair_id, r in pairs.items():
        lof_z = r.get("lof_z_burden")
        lof_beta = r.get("lof_beta_burden")
        lof_p = r.get("lof_p_burden")
        mr_beta = r.get("v34_mr_beta")
        mr_p = r.get("v34_mr_p")

        if lof_z is None or mr_beta is None:
            continue

        outcome = 1 if r["v34_outcome"] == "SUCCESS" else 0
        lof_sign = 1 if lof_beta > 0 else -1
        mr_sign = 1 if mr_beta > 0 else -1
        agree = lof_sign == mr_sign
        either_sig = (lof_p < 0.05) or (mr_p < 0.05)

        all_with_both.append({"outcome": outcome, "pair_id": pair_id, "agree": agree, "either_sig": either_sig})

        if agree and either_sig:
            depth2.append({"outcome": outcome, "pair_id": pair_id, "gene": r["gene"], "disease": r["disease"], "mechanism": r["mechanism_class"]})

    n_all = len(all_with_both)
    n_d2 = len(depth2)
    print(f"  Pairs with both instruments: {n_all}")
    print(f"  Depth-2 licensed: {n_d2}")

    if n_d2 > 0:
        d2_success = sum(1 for d in depth2 if d["outcome"] == 1)
        d2_failure = n_d2 - d2_success
        d2_ba = (d2_success / max(sum(1 for d in depth2 if d["outcome"] == 1), 1) +
                 (d2_failure - sum(1 for d in depth2 if d["outcome"] == 0)) / max(1, 1)) if n_d2 > 1 else 0

        d2_outcomes = [d["outcome"] for d in depth2]
        d2_pred_success_rate = sum(d2_outcomes) / len(d2_outcomes)

        all_success_rate = sum(d["outcome"] for d in all_with_both) / n_all

        print(f"  Depth-2 SUCCESS rate: {d2_success}/{n_d2} ({d2_pred_success_rate:.1%})")
        print(f"  All-pairs SUCCESS rate: {all_success_rate:.1%}")
        print(f"  Enrichment: {d2_pred_success_rate / all_success_rate:.2f}x" if all_success_rate > 0 else "  Enrichment: N/A")

        print(f"\n  Depth-2 licensed pairs:")
        for d in depth2:
            label = "SUCCESS" if d["outcome"] == 1 else "FAILURE"
            print(f"    {d['gene']:12s} | {d['disease']:30s} | {d['mechanism']:25s} | {label}")

        d2_scores = []
        d2_labels = []
        for d in depth2:
            r = pairs[d["pair_id"]]
            combined_z = (r["lof_z_burden"] + pqtl_z(r)) / 2
            d2_scores.append(combined_z)
            d2_labels.append(d["outcome"])

        if sum(d2_labels) > 0 and sum(d2_labels) < len(d2_labels):
            s1 = [s for s, o in zip(d2_scores, d2_labels) if o == 1]
            s0 = [s for s, o in zip(d2_scores, d2_labels) if o == 0]
            u, p_mw = mannwhitneyu(s1, s0, alternative="two-sided")
            d2_auc = u / (len(s1) * len(s0))
            print(f"\n  Depth-2 AUC (combined z): {d2_auc:.3f} (MW p={p_mw:.4f})")
        else:
            d2_auc = None
            print(f"\n  Depth-2 AUC: cannot compute (all same outcome)")
    else:
        d2_success = 0
        d2_pred_success_rate = 0
        d2_auc = None

    return {
        "n_all": n_all, "n_depth2": n_d2,
        "depth2_success": d2_success,
        "depth2_success_rate": d2_pred_success_rate if n_d2 > 0 else None,
        "all_success_rate": sum(d["outcome"] for d in all_with_both) / n_all if n_all > 0 else None,
        "depth2_auc": d2_auc,
        "depth2_pairs": [{"gene": d["gene"], "disease": d["disease"], "mechanism": d["mechanism"], "outcome": "SUCCESS" if d["outcome"] == 1 else "FAILURE"} for d in depth2],
    }


def analysis_4_thresholded(pairs):
    """Secondary: Thresholded LoF classifier (p<0.05 → SUCCESS)."""
    print("\n=== Analysis 4: Thresholded LoF Classifier ===")

    tp = fp = tn = fn = 0
    for r in pairs.values():
        if r.get("lof_p_burden") is None:
            continue
        pred = r["lof_p_burden"] < 0.05
        actual = r["v34_outcome"] == "SUCCESS"
        if pred and actual: tp += 1
        elif pred and not actual: fp += 1
        elif not pred and actual: fn += 1
        else: tn += 1

    n = tp + fp + tn + fn
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0
    ba = (sens + spec) / 2
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0

    print(f"  N={n}, TP={tp}, FP={fp}, TN={tn}, FN={fn}")
    print(f"  Sensitivity: {sens:.3f}")
    print(f"  Specificity: {spec:.3f}")
    print(f"  BA: {ba:.3f}")
    print(f"  PPV: {ppv:.3f}")

    return {"n": n, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
            "sensitivity": sens, "specificity": spec, "ba": ba, "ppv": ppv}


def analysis_5_mechanism_auc(pairs):
    """Secondary: Mechanism-stratified LoF AUC."""
    print("\n=== Analysis 5: Mechanism-Stratified LoF AUC ===")

    results = {}
    for mech in ["activity_blocking", "abundance_modulating"]:
        subset = {k: v for k, v in pairs.items() if v["mechanism_class"] == mech and v.get("lof_z_burden") is not None}
        outcomes = [1 if r["v34_outcome"] == "SUCCESS" else 0 for r in subset.values()]
        scores = [r["lof_z_burden"] for r in subset.values()]

        n = len(outcomes)
        n_s = sum(outcomes)
        n_f = n - n_s

        if n_s > 0 and n_f > 0:
            s1 = [s for s, o in zip(scores, outcomes) if o == 1]
            s0 = [s for s, o in zip(scores, outcomes) if o == 0]
            u, p_mw = mannwhitneyu(s1, s0, alternative="two-sided")
            auc = u / (n_s * n_f)
            print(f"  {mech}: N={n} (S={n_s}, F={n_f}), AUC={auc:.3f} (MW p={p_mw:.4f})")
            results[mech] = {"n": n, "n_success": n_s, "n_failure": n_f, "auc": auc, "mannwhitney_p": p_mw}
        else:
            print(f"  {mech}: N={n} (S={n_s}, F={n_f}), AUC=N/A (single class)")
            results[mech] = {"n": n, "n_success": n_s, "n_failure": n_f, "auc": None}

    return results


def analysis_6_combined(pairs):
    """Exploratory: Combined instrument score (avg of pQTL z + LoF z)."""
    print("\n=== Analysis 6: Combined Instrument Score ===")

    outcomes = []
    combined_scores = []
    pqtl_scores = []
    lof_scores = []

    for r in pairs.values():
        lof_z = r.get("lof_z_burden")
        pz = pqtl_z(r)
        if lof_z is None or pz is None:
            continue
        outcome = 1 if r["v34_outcome"] == "SUCCESS" else 0
        outcomes.append(outcome)
        combined_scores.append((lof_z + pz) / 2)
        pqtl_scores.append(pz)
        lof_scores.append(lof_z)

    n = len(outcomes)
    n_s = sum(outcomes)
    n_f = n - n_s
    print(f"  N={n} (S={n_s}, F={n_f})")

    def compute_auc(scores, outcomes):
        s1 = [s for s, o in zip(scores, outcomes) if o == 1]
        s0 = [s for s, o in zip(scores, outcomes) if o == 0]
        if len(s1) == 0 or len(s0) == 0:
            return None, None
        u, p = mannwhitneyu(s1, s0, alternative="two-sided")
        return u / (len(s1) * len(s0)), p

    auc_combined, p_combined = compute_auc(combined_scores, outcomes)
    auc_pqtl, p_pqtl = compute_auc(pqtl_scores, outcomes)
    auc_lof, p_lof = compute_auc(lof_scores, outcomes)

    print(f"  pQTL-only AUC:  {auc_pqtl:.3f} (p={p_pqtl:.4f})" if auc_pqtl else "  pQTL-only AUC: N/A")
    print(f"  LoF-only AUC:   {auc_lof:.3f} (p={p_lof:.4f})" if auc_lof else "  LoF-only AUC: N/A")
    print(f"  Combined AUC:   {auc_combined:.3f} (p={p_combined:.4f})" if auc_combined else "  Combined AUC: N/A")

    if auc_combined and auc_pqtl:
        delta = auc_combined - auc_pqtl
        print(f"  Delta (combined - pQTL): {delta:+.3f}")

    return {
        "n": n, "n_success": n_s, "n_failure": n_f,
        "auc_pqtl_only": auc_pqtl, "p_pqtl": p_pqtl,
        "auc_lof_only": auc_lof, "p_lof": p_lof,
        "auc_combined": auc_combined, "p_combined": p_combined,
    }


def main():
    all_pairs, has_lof = load_pairs()
    print(f"Total pairs: {len(all_pairs)}")
    print(f"Pairs with LoF data: {len(has_lof)}")
    print(f"Structural missingness: {len(all_pairs) - len(has_lof)}")
    print()

    results = {}
    results["analysis_1_lof_auc"] = analysis_1_lof_auc(has_lof)
    results["analysis_2_concordance"] = analysis_2_concordance(has_lof)
    results["analysis_3_depth2"] = analysis_3_depth2(has_lof)
    results["analysis_4_thresholded"] = analysis_4_thresholded(has_lof)
    results["analysis_5_mechanism_auc"] = analysis_5_mechanism_auc(has_lof)
    results["analysis_6_combined"] = analysis_6_combined(has_lof)

    results["meta"] = {
        "total_pairs": len(all_pairs),
        "pairs_with_lof": len(has_lof),
        "structural_missingness": len(all_pairs) - len(has_lof),
        "prespec": "PRESPEC_V4.md",
        "instrument_data": "genebass_lof_results_v2.json",
    }

    with open(OUTPUT, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\n\nSaved all results to {OUTPUT}")


if __name__ == "__main__":
    main()
