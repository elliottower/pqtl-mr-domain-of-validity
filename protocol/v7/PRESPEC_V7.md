# Pre-Specification V7: Three-Level Mechanism Stratification

**Author:** Elliot Tower
**Date:** 2026-08-16
**Status:** FROZEN 2026-08-16 (hash in PRESPEC_V7_sha256.txt)
**Depends on:** V3.4, V5 (instrument expansion), V6 (five-way taxonomy, superseded here)

---

## 0. What was known when this was written

This document is written with knowledge of prior results. It is not blind, and the
provenance is in the repository, so there is no point pretending otherwise.

**Known at time of writing:**

- Stratum-level balanced accuracies from the V3.4 and V5 datasets, for the binary
  abundance/activity split and for the V6 five-way subcategories.
- That a test of the binary stratum difference, run 2026-08-16, did not reach
  conventional significance.
- That eight outcome-GWAS accessions in `code/classify_v34.py` are wrong or
  unresolvable, affecting roughly 20% of analysed rows.

**Not known, and the reason this document can still do its job:**

The corrected dataset does not exist. Every accession error is unfixed at the time of
writing. No stratum estimate on corrected data has been computed by anyone. What is
frozen here — the stratification, the test, the decision rules — is therefore fixed
before the numbers it will be applied to exist.

**The principle this rests on.** Pre-registration is valid when the *outcome on the
frozen dataset is unknown*, not when the author's prior reasoning has been erased.
Registrations routinely follow pilot data, failed specifications, and abandoned
analyses — that is the normal case, not an exception, and it is why deviation logs
exist. What would invalidate this document is knowing how the corrected data comes out.
Nobody does, because it has not been produced.

**What this document does not claim.** It does not claim blind predictions. Predictions
in §4 are stated directionally only, because an interval prediction from an author who
has seen adjacent results is worthless whatever label is attached to it. The load-bearing
content is §2 and §3: what the strata are, and what test is run.

---

## 1. Why the stratification is changing

V6 specified a five-way taxonomy: A1 soluble ligand neutralisation, A2 surface target
depletion, A3 indirect abundance modulation, B1 enzyme inhibitor, B2 receptor/kinase
antagonist. V6 also defined a three-way collapse in its §5.3 but treated it as a
convenience for small strata rather than as the primary contrast.

V7 makes the three-way collapse primary and retires the five-way. The reason is that the
five-way encodes distinctions the mechanism argument asserts but does not derive.

**B1 and B2 make the identical prediction.** An enzyme active-site inhibitor and a
receptor antagonist both leave circulating protein abundance unchanged while disrupting
what the protein does. A cis-pQTL is blind to both by exactly the same route. Nothing in
the theory says catalytic inhibition and receptor blockade should differ in how invisible
they are to an abundance instrument. Two bins are doing one bin's work.

**A2 and A3 are not separated by a derived argument.** A2 says the assay reads a shed
soluble form while the drug engages the membrane-bound form. A3 says the drug changes
abundance through an intermediate. Both describe a measured quantity that is correlated
with, but not identical to, the drug's target. Which should attenuate more is not
predicted by the theory — one could argue it either way, and nobody has.

**The distinction the theory does derive is between identity and correlation.** Either
the pQTL measures the same molecular species the drug acts on, or it measures something
correlated with it, or it measures a quantity the drug leaves untouched. Those three make
genuinely different predictions. Finer cuts do not.

**Why not the binary either.** The abundance/activity split lumps identity and
correlation together. Direct engagement is the only case where the genetic perturbation
and the pharmacological perturbation are the *same quantity*; correlation introduces an
attenuating step. Collapsing those discards a prediction the theory does make.

Three levels is the resolution the mechanism argument supports. Not two, not five.

---

## 2. The three strata

Assignment is on the drug's relationship to the circulating molecular species the
cis-pQTL instruments.

**D — Direct.** The drug binds, neutralises, replaces, or is the same circulating
molecule the pQTL measures. Monoclonal antibodies against soluble cytokines, ligand
traps, recombinant protein replacement. Genetic variation in abundance and
pharmacological change in abundance are the same quantity.
*Maps from V6: A1.*

**C — Correlated.** The drug acts on a species correlated with, but not identical to,
what the assay measures. The assay reads a shed soluble form while the drug engages the
membrane-bound form; an antisense oligonucleotide reduces transcript rather than binding
protein; an inhibitor of a secreted enzyme where level and activity covary.
*Maps from V6: A2, A3.*

**B — Function-blocking.** The drug leaves circulating abundance unchanged and disrupts
function. Enzyme active-site inhibitors, receptor antagonists, kinase inhibitors.
*Maps from V6: B1, B2.*

**Ordering.** D > C > B in expected MR informativeness. The ordering is fixed here and is
not revisited.

### 2.1 Assignment of genes not covered by the V6 frozen map

The V6 gene-to-subcategory assignment was frozen on the V3.4 gene set. The V5 expansion
introduced 13 genes with no assignment: CD86, EPHA2, EPHA4, EPHB4, ERBB4, GHR, IL17F,
IL2RA, ITGB5, MIF, MUC16, PCSK9, TNFSF13B.

These are assigned to D, C, or B by a rater who has not seen, for any of these genes, the
trial outcome, the MR estimate, the MR p-value, or any stratum result. The rater sees the
gene symbol, the drug or drug class, and the mechanism of action only. Assignments are
recorded and committed before the corrected dataset is produced.

Genes whose mechanism cannot be determined from the drug label are assigned UNCLASSIFIED
and excluded from the stratified analysis, never guessed into a bin.

---

## 3. The test

**Primary — Firth penalised logistic interaction with an ordered mechanism tier.**

```
success ~ mr_evidence + tier + mr_evidence × tier
```

where `mr_evidence` is −log₁₀(p_MR) as a continuous predictor, `tier` is coded linearly
as D = 1, C = 2, B = 3, and `success` is the binary adjudicated Phase III outcome.
**The interaction coefficient is the hypothesis**: it tests whether the association
between genetic evidence strength and trial success weakens as the drug's action moves
away from the species the instrument measures. A negative interaction in the direction
D → B supports the mechanism account.

**One-sided.** The direction is derived from the mechanism argument in §1 and fixed
before any data, so the test is one-sided at α = 0.05. A positive interaction, however
large, does not support the account and is reported as a failure of P1.

**Firth penalisation.** Cell counts within tier are expected to be small, and standard
maximum likelihood is unstable or non-convergent under separation. Firth's penalised
likelihood is specified in advance for this reason — it is a numerical-stability choice,
not a power choice, and applies regardless of what the data show.

### 3.1 Why continuous evidence here and a binary threshold elsewhere

These are two instruments answering two questions, both fixed in advance.

**The binary threshold (p < 0.05) is retained for the decision-rule analysis** — the
positive and negative predictive values a practitioner uses when deciding what a
significant or null MR result means for a target. That question is inherently
dichotomous: an analyst has a result in hand and must act. That analysis was
pre-specified in V3 and is unchanged by this document.

**Continuous evidence is used for the mechanism interaction** because the question is
different: does the *strength* of the genetic-evidence-to-success association scale with
mechanistic alignment? Thresholding discards the information that question needs. Only
pairs crossing p < 0.05 carry any signal under a binary coding, which is 19 of 157 —
roughly an 88% loss of informative observations for this hypothesis specifically.

The choice is also supported externally and from a different direction: Ravarani et al.
(2026), at 11,482 target-indication pairs in Phase II, found that binary MR significance
does not enrich for success while continuous MR-derived features do. Different dataset,
different phase, different group, same structural conclusion.

Neither instrument is selected on the basis of results. Both are specified here, before
the corrected dataset exists, with the reason for each stated.

**Why this rather than a trend test on proportions.** A Cochran–Armitage or
Jonckheere–Terpstra trend test asks whether the *success rate* differs across ordered
groups. That is a fact about which drugs occupy which tier, not about whether the
classifier discriminates better in some tiers than others. The claim under test concerns
classifier performance, which is a property of the 2×2 within each tier. Only the
interaction expresses that. Cochran–Armitage on the main effect is reported as a
descriptive secondary so that a reader can see whether tiers differ in base rate, which
is a distinct and also interesting question.

**Clustering.** Stratum membership is a property of the gene, not the pair, and genes
recur across diseases. Inference is therefore gene-clustered: the interaction p-value is
obtained by permuting tier labels **at the gene level** over 20,000 iterations, not from
the asymptotic standard error and not by row-level shuffling. Row-level permutation is
anti-conservative here because it breaks the gene-tier dependency that actually exists.
Cluster-robust standard errors by gene are reported alongside as a check.

**Reported alongside, for every stratum:** n, the full 2×2 confusion table, sensitivity,
specificity, PPV, NPV, and balanced accuracy with a 90% bootstrap interval over 10,000
resamples, resampled at gene level.

**Also reported:** the pairwise D vs B contrast, since it is the largest expected gap and
the one comparable to the previously reported binary split.

**Not reported as primary:** any single stratum in isolation. No stratum becomes the
headline on its own, whatever it shows.

---

## 4. Prediction

Directional only.

**P1.** Balanced accuracy is ordered D > C > B.

**P2.** The D versus B contrast is the largest of the three pairwise gaps.

**P3.** B does not differ detectably from chance.

No interval, point estimate, or significance level is predicted. See §0 for why.

**Falsifier stated in advance:** if B exceeds D, the mechanism account is wrong, and this
document says so before the test is run.

---

## 5. Interpretation, fixed before scoring

**Trend significant, ordering as predicted.** Licensed: that MR informativeness varies
with how directly the drug engages the species the instrument measures, on this set at
this sample size. Not licensed: any claim about a specific stratum's absolute performance,
or extrapolation beyond cis-pQTL instruments in Phase III.

**Trend not significant, ordering as predicted.** Licensed: the ordering is consistent
with the mechanism account and the study is underpowered to establish it. Report the
sample size that would be required. This is the most likely outcome and it is not a
failure of the study; it is a measurement of what the available instrument universe can
support.

**Trend not significant, ordering not as predicted.** The mechanism account is
unsupported on this data. Report it plainly. Do not re-bin, do not add strata, do not
switch metric.

### 5.1 Power, computed before freezing

Analytic Wald power for the interaction coefficient, assuming the published marginal
structure (157 pairs, base success 0.306), a main effect of 0.45 log-OR per unit of
−log₁₀(p), and one-sided α = 0.05. Script: `scripts/power_v7_interaction.py`. No outcome
data, tier assignment, or result file enters this calculation.

| Design | Interaction 0.10 | 0.20 | 0.30 | 0.45 | MDE at 80% |
|---|---|---|---|---|---|
| 3 tiers, ~52/tier | 26% | 62% | 85% | 96% | **0.27** |
| 2 tiers, ~78/tier | 18% | 42% | 70% | 93% | 0.35 |

**Three tiers is chosen, and the power calculation agrees with the mechanistic argument
in §1.** Under linear tier coding, three levels span two units of the predictor while two
span one, so the interaction column carries more variance and its coefficient is more
precisely estimated. The familiar "fewer groups, more power" heuristic applies to
unordered comparisons of separate group means, not to a linear-coded ordinal term.

**What the study can and cannot detect.** The mechanism account predicts complete
attenuation: the evidence-success association is present in D and absent in B. Against an
assumed main effect of 0.45, complete attenuation over two tier steps corresponds to an
interaction near 0.225, where power is roughly 65–70%. **The design is therefore powered
at about two-in-three odds to detect complete attenuation and is not powered to detect
partial attenuation.** This is stated here, before the test is run, so that a null result
is read as an underpowered test of a strong prediction rather than as evidence against
the mechanism account.

The ceiling is structural. Only ~211 Phase III drug-target pairs in Open Targets have any
cis-pQTL instrument, because most drug targets are not measurable in plasma. No design
choice recovers power that the instrument universe does not contain.

**Caveat on the calculation.** The Wald approximation above assumes an exponential
evidence distribution and applies a fixed adjustment for collinearity among the four
design columns. It is a design aid, not a precise operating characteristic, and the
realised interval will differ.

A second caveat applies to the two-versus-three-tier row specifically. The comparison
assumes the same log-OR per unit tier step under both parameterisations, which is a
simplifying assumption rather than a guarantee that the two designs test equivalent
hypotheses. Collapsing to D versus C+B merges two theoretical steps into one, and the
true per-step effect under that collapsed coding need not equal the per-step effect under
three levels. The three-tier choice does not rest on this row alone — it follows from §1,
and it retains the C tier's information rather than discarding it.

### 5.2 Minimum stratum size

Any tier with fewer than **20 pairs** is reported with its counts and excluded from the
trend estimate, with the exclusion and the reason stated. Twenty is roughly 40% of the
balanced allocation of 52; below that the tier contributes too little to the interaction
column to move the estimate, while adding an extreme value to the linear coding. If a
tier is excluded, the trend is estimated on the remaining tiers and the reduced design's
MDE is recomputed and reported alongside.

---

## 6. What is frozen

The three strata and their ordering (§2). The V6-to-V7 mapping (§2). The blind assignment
procedure for the 13 unassigned genes (§2.1). The test, its clustering level, and its
iteration count (§3). The reporting set (§3). The interpretation rules (§5).

**No element above may change after the corrected dataset exists.** Any change is a new
protocol version with its own disclosure, and does not replace this one.

---

## 7. What this cannot establish

That the three-level scheme is the correct resolution — only that it is the resolution the
mechanism argument supports, stated before the corrected data existed.

That mechanism causes the difference. Stratum membership is confounded with drug modality,
target biology, and disease area, and this design cannot separate them.

Anything about instrument classes other than cis-pQTLs, or trial phases other than III.
