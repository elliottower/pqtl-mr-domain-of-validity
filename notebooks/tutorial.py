"""Tutorial: Using pqtl-validity to assess MR informativeness by drug mechanism.

This script demonstrates the pqtl-validity library on three cases:
1. ACE inhibitors — activity-blocking drugs where MR nulls are expected
2. TNF neutralization — soluble ligand mechanism where MR should work
3. Batch scoring: evaluate your own drug-target pairs against the taxonomy

Run as: python notebooks/tutorial.py
Or install and run: pip install -e . && python notebooks/tutorial.py
"""

import pqtl_validity as pv

# %% Case 1: ACE inhibitors — MR null is the correct answer to the wrong question
#
# ACE inhibitors block enzymatic activity. The cis-pQTL measures ACE protein
# abundance, which the drug does not change. MR should return null, and it does.

result = pv.classify_pair("ACE", mr_p=0.72)
print("=== Case 1: ACE (enzyme inhibitor) ===")
print(result)
print(f"  Mechanism: {result.mechanism.value}")
print(f"  MR informative expected? {result.mr_informative_expected}")
print(f"  MR significant? {result.mr_significant}")
print(f"  Prediction: {result.prediction}")
print()
print("  Interpretation: The MR null for ACE is not evidence against the target.")
print("  The instrument measures abundance; the drug acts on enzymatic function.")
print()


# %% Case 2: TNF — MR works because the drug depletes the protein
#
# Anti-TNF biologics (infliximab, adalimumab) neutralize soluble TNF.
# cis-pQTL MR should detect the causal effect because the mechanism IS
# abundance reduction.

result_tnf = pv.classify_pair("TNF", mr_p=0.003)
print("=== Case 2: TNF (soluble ligand neutralization) ===")
print(result_tnf)
print(f"  Mechanism: {result_tnf.mechanism.value}")
print(f"  MR informative expected? {result_tnf.mr_informative_expected}")
print(f"  MR significant? {result_tnf.mr_significant}")
print(f"  Prediction: {result_tnf.prediction}")
print()
print("  Interpretation: MR correctly detects the causal effect for TNF because")
print("  the drug's mechanism (neutralizing soluble ligand) is what the pQTL measures.")
print()


# %% Case 3: EGFR — kinase inhibitor, same story as ACE

result_egfr = pv.classify_pair("EGFR", mr_p=0.45)
print("=== Case 3: EGFR (receptor/kinase antagonist) ===")
print(result_egfr)
print(f"  Mechanism: {result_egfr.mechanism.value}")
print(f"  MR informative expected? {result_egfr.mr_informative_expected}")
print()


# %% Case 4: Unknown gene — not in the pre-registered taxonomy

result_unknown = pv.classify_pair("BRCA1")
print("=== Case 4: BRCA1 (not in taxonomy) ===")
print(result_unknown)
print(f"  Mechanism: {result_unknown.mechanism}")
print()


# %% Case 5: Taxonomy overview

s = pv.summary()
print("=== Taxonomy Summary ===")
print(f"  Total classified genes: {s.total}")
print(f"  A1 (soluble ligand):     {s.A1}")
print(f"  A2 (surface target):     {s.A2}")
print(f"  A3 (indirect abundance): {s.A3}")
print(f"  B1 (enzyme inhibitor):   {s.B1}")
print(f"  B2 (receptor/kinase):    {s.B2}")
print(f"  Mixed:                   {s.mixed}")
print(f"  Abundance-mediated (MR informative): {s.abundance_mediated}")
print(f"  Activity-blocking (MR uninformative): {s.activity_blocking}")
print()


# %% Case 6: Batch scoring with known outcomes
#
# Score a batch of pairs where we know drug approval outcome, using the same
# screen the paper evaluates (predict SUCCESS when MR is significant).
# This is an illustrative subset, not a reproduction: every outcome here is
# SUCCESS, so specificity is undefined and reported as zero. For the published
# figures see scripts/deviation_log_gaps.py, which runs the full tables.

genes = ["TNF", "VEGFA", "IL6R", "ACE", "EGFR", "DPP4", "PARP1", "KDR"]
outcomes = ["SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS", "SUCCESS"]
mr_ps = [0.003, 0.01, 0.008, 0.72, 0.45, 0.88, 0.61, 0.53]

scoring = pv.score_pairs(genes, outcomes, mr_ps=mr_ps, seed=42)
print("=== Case 6: Batch Scoring ===")
print(scoring)
print(f"  Balanced accuracy: {scoring.balanced_accuracy:.3f}")
print(f"  90% CI: [{scoring.ba_ci_low:.3f}, {scoring.ba_ci_high:.3f}]")
print(f"  Sensitivity: {scoring.sensitivity:.3f}")
print(f"  Specificity: {scoring.specificity:.3f}")
print()


# %% Case 7: Bring your own data
#
# To classify your own drug-target pairs:
#
#   import pqtl_validity as pv
#
#   # Single pair
#   result = pv.classify_pair("YOUR_GENE", mr_p=your_pvalue)
#   print(result.prediction)  # "INFORMATIVE" or "UNINFORMATIVE"
#
#   # Check if a gene is in the taxonomy
#   mech = pv.classify_gene("YOUR_GENE")
#   if mech is None:
#       print("Gene not in pre-registered taxonomy")
#
#   # Get the full gene list
#   all_assignments = pv.all_genes()
#   print(f"{len(all_assignments)} genes classified")
