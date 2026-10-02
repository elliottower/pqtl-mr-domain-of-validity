"""The outcome-blind counts the plan registers, and the same counts taken from stage A's output.

PREREG §Sample size states the counts of the round-3 feasibility run
(`feasibility/v4_round3/count_v4.py`, `coverage_v4.json`), which "applies the rules in this
document". Stage A implements the same rules on the same pinned inputs, so its output must give the
same numbers. `run_stage_a.py` writes both sides to `registered_count_check.json` and refuses to
seal the output (no MANIFEST.tsv) when any differs. Stage A reads no outcome, so a difference is
investigated and the stage rerun under a new run token.

Each quantity is computed as count_v4.py computes it (`summarize`): hypotheses are rows, pairs are
distinct (gene symbol, indication ID), genes are distinct gene symbols.
"""
import pandas as pd

# PREREG.md §Sample size, second paragraph, in the order the paragraph states them.
REGISTERED_COUNTS: dict[str, int] = {
    "s1_hypotheses": 5064,                           # "S1 has 5,064 held-out hypotheses"
    "s1_gene_indication_pairs": 4778,                # "(4,778 gene-indication pairs,"
    "s1_genes": 384,                                 # "384 genes,"
    "s1_indications": 493,                           # "493 indications)"
    "h1_hypotheses": 3253,                           # "The H1 model has 3,253:"
    "h1_aligned_hypotheses": 582,                    # "582 abundance-aligned"
    "h1_aligned_genes": 69,                          # "in 69 genes"
    "h1_blocking_hypotheses": 2671,                  # "2,671 function-blocking"
    "h1_blocking_genes": 170,                        # "in 170 genes"
    "h1_genes_in_both_classes": 11,                  # "with 11 genes in both classes"
    "s1_other_hypotheses": 1811,                     # "1,811 are other"
    "h4_neuro_psych_hypotheses": 698,                # "698 hypotheses"
    "h4_neuro_psych_genes": 168,                     # "in 168 genes on neurologic or psychiatric indications"
    "h4_other_indication_hypotheses": 4366,          # "4,366 on other indications"
    "s4_hypotheses": 5355,                           # "S4 has 5,355"
    "h1_blood_secreted_aligned_hypotheses": 582,     # "Restricted to blood-secreted targets on both sides (S19),
    "h1_blood_secreted_aligned_genes": 69,           # H1 has 582 aligned hypotheses in 69 genes
    "h1_blood_secreted_blocking_hypotheses": 408,    # and 408 blocking
    "h1_blood_secreted_blocking_genes": 30,          # in 30 genes"
}
SOURCE = "PREREG.md §Sample size (counts of feasibility/v4_round3/count_v4.py, coverage_v4.json)"


class RegisteredCountMismatch(RuntimeError):
    """Stage A's output does not give the outcome-blind counts the plan registers."""


def obtained_counts(hyp: pd.DataFrame) -> dict[str, int]:
    """The registered quantities from `hypotheses.csv` rows (HypothesisRow columns).

    S1 is `in_s1 & heldout`; H1 its aligned and blocking rows; the neurologic-or-psychiatric split
    is `neuro_psych`; S4 is `in_s4 & heldout`. The last four quantities are count_v4.py's
    `H1_secreted_only`: H1 rows whose target is HPA blood-secreted (`blood_secreted_hpa`). That is
    the count the plan's sentence states; the analysis set S19 is formed in stage D from `s19_arm`
    and is narrower."""
    s1 = hyp[hyp["in_s1"] & hyp["heldout"]]
    s4 = hyp[hyp["in_s4"] & hyp["heldout"]]
    aligned, blocking = s1[s1["mechanism_class"] == "aligned"], s1[s1["mechanism_class"] == "blocking"]
    neuro = s1[s1["neuro_psych"]]
    blood_aligned, blood_blocking = aligned[aligned["blood_secreted_hpa"]], blocking[blocking["blood_secreted_hpa"]]
    counts = {
        "s1_hypotheses": len(s1),
        "s1_gene_indication_pairs": len(s1[["gene_symbol", "indication_id"]].drop_duplicates()),
        "s1_genes": s1["gene_symbol"].nunique(),
        "s1_indications": s1["indication_id"].nunique(),
        "h1_hypotheses": len(aligned) + len(blocking),
        "h1_aligned_hypotheses": len(aligned),
        "h1_aligned_genes": aligned["gene_symbol"].nunique(),
        "h1_blocking_hypotheses": len(blocking),
        "h1_blocking_genes": blocking["gene_symbol"].nunique(),
        "h1_genes_in_both_classes": len(set(aligned["gene_symbol"]) & set(blocking["gene_symbol"])),
        "s1_other_hypotheses": int((s1["mechanism_class"] == "other").sum()),
        "h4_neuro_psych_hypotheses": len(neuro),
        "h4_neuro_psych_genes": neuro["gene_symbol"].nunique(),
        "h4_other_indication_hypotheses": len(s1) - len(neuro),
        "s4_hypotheses": len(s4),
        "h1_blood_secreted_aligned_hypotheses": len(blood_aligned),
        "h1_blood_secreted_aligned_genes": blood_aligned["gene_symbol"].nunique(),
        "h1_blood_secreted_blocking_hypotheses": len(blood_blocking),
        "h1_blood_secreted_blocking_genes": blood_blocking["gene_symbol"].nunique(),
    }
    return {k: int(v) for k, v in counts.items()}


def count_check(hyp: pd.DataFrame, registered: dict[str, int]) -> dict:
    """The record written to registered_count_check.json: per quantity the registered and the
    obtained value and whether they agree, and `all_match`."""
    obtained = obtained_counts(hyp)
    if set(obtained) != set(registered):
        raise RegisteredCountMismatch(f"registered quantities {sorted(set(registered) ^ set(obtained))} have no counterpart")
    quantities = {k: {"registered": int(registered[k]), "obtained": obtained[k], "match": obtained[k] == int(registered[k])}
                  for k in registered}
    return {"source": SOURCE, "quantities": quantities, "all_match": all(q["match"] for q in quantities.values())}
