# Direction-Concordance Classifier: Disclosed Exploratory Analysis

**Author:** Elliot Tower
**Date:** 2026-08-16
**Status:** EXPLORATORY — **NOT PRE-REGISTERED**
**Depends on:** V3.4 (Zenodo Version 2), V5 (pQTL Instrument Expansion), V6 (Expanded Mechanism Taxonomy)

---

## 0. Disclosure — read before anything else

**This analysis was run before this document existed. It cannot be called pre-registered, and it will not be called pre-registered anywhere in the paper, the supplement, the Zenodo record, or any cover letter.**

The specific facts:

1. On **2026-08-16**, earlier the same day this document was written, an assistant ran an unblinded scratch computation of the direction-concordance classifier described in §3.
2. The numeric results of that run — balanced accuracies, reclassification counts, or both — were **displayed in a working session**.
3. That run happened **before any protocol for it was written**. No hypothesis, classifier definition, stratum list, or decision rule existed on paper at the time it executed.
4. The author of this document has **not seen those results**. The study author states that he did not read them.

Point 4 does not repair points 1 through 3. The results were computed and rendered. They existed in a context that the study author could have read, and the analysis design in §3 was described to the protocol author by someone who had already been in that context. Blinding that depends on a person choosing not to look is not the same guarantee as blinding that depends on the numbers not existing. Every prediction in §5 is therefore stated under weaker protection than the V1–V6 predictions, and must be reported as such.

**The filename `PREREG_direction_concordance.md` is a misnomer** retained only so that existing path references do not break. This document is a *disclosed exploratory protocol*. Do not cite it as a pre-registration.

**Hashing does not fix this.** This document may be SHA-256 hashed for tamper-evidence, so that the predictions in §5 can be shown to predate the *reported* scoring run. A hash proves only when the text was frozen. It does not and cannot establish that the analysis was designed before results existed, because it was not. Any hash file for this document must carry the sentence: *"Hash establishes document integrity only. The analysis was run unblinded on 2026-08-16 before this protocol was written."*

**Required reporting language.** Wherever this analysis appears, it is introduced as: *"A post-hoc direction-concordance analysis, run once in unblinded form before its protocol was written (see experiments/PREREG_direction_concordance.md §0), and reported here with that history disclosed."* This sentence is not optional and is not to be shortened.

**Contrast with the frozen work.** V1–V3.4 were hashed before outcome adjudication. V4, V5, and V6 were hashed before their respective results were computed. This analysis has none of that protection. It sits in `experiments/`, outside `protocol/`, for exactly that reason. The frozen classifier (MR p < 0.05 → SUCCESS) remains the primary result of the study and is not revised by anything below.

---

## 1. Rationale

### 1.1 What the frozen classifier ignores

The primary classifier is zero-parameter: `mr_p < 0.05` → predict SUCCESS, otherwise predict FAILURE. It reads the p-value and discards the effect estimate. A pair whose MR analysis returns a significant estimate pointing *opposite* to the direction the drug's mechanism requires is currently scored as a prediction of SUCCESS.

That case is substantively strange. Consider an inhibitor of a target whose cis-pQTL MR estimate says, with p < 0.05, that higher circulating levels of that protein are *protective*. The genetic evidence, read at face value, argues that lowering the protein should make the disease worse. The frozen classifier counts this as genetic support for developing the inhibitor. The proposed analysis asks what happens when that case is instead scored as a prediction of FAILURE.

### 1.2 Why the frozen classifier omits direction

The omission was deliberate and is defensible. The classifier was specified in V3.1 to have zero tunable parameters, so that Stage C could not tune anything against outcome labels (see `protocol/DEVIATION_LOG.md`, "MR classifier threshold"). Adding a direction requirement adds a second input — the drug's mechanism direction — that must be derived from annotation rather than from the MR result. Annotation introduces a judgement call, and judgement calls are what the outcome-blind freeze was designed to exclude. A simpler classifier is also a stronger claim: if a rule with no free parameters and no curation beyond mechanism-class stratification carries signal, the signal is hard to attribute to analyst choices.

The cost of that simplicity is the case in §1.1. This analysis measures the size of that cost. It does not propose replacing the frozen classifier.

### 1.3 Why direction is the contested variable in this literature

Ravarani et al. (2026, medRxiv doi 10.64898/2026.02.19.26346536) evaluated MR against clinical trial outcomes at scale: 11,482 target–indication pairs with Phase II outcomes, built on the Minikel dataset. They found that MR statistical significance alone does not enrich for Phase II success, and rescued prediction only with an MR-informed XGBoost classifier.

They name the absence of mechanism-of-action annotation as their central limitation. Their resource "does not encode the directionality of pharmacologic modulation, i.e., whether a given drug acts as an inhibitor or activator of its target," so they "were unable to ... align genetic estimates with the mechanism of action of individual compounds."

The parent study already supplies half of what they say is missing: a frozen mechanism-class annotation (abundance-modulating vs activity-blocking), which is what produced the domain-of-validity result. It does not supply the other half — the *sign* alignment between the genetic estimate and the pharmacologic direction. The `action_types` field in `data/adjudicated_v34.csv` (INHIBITOR, ANTAGONIST, AGONIST, ACTIVATOR, MODULATOR, and others) makes that alignment computable on this set. Whether the alignment helps is an open empirical question, and it is the one this analysis addresses.

### 1.4 A caution carried forward from V3.4

`action_types` is not `mechanism_class`, and the parent protocol says so explicitly: "The `action_type` field from Open Targets (INHIBITOR, ANTAGONIST, etc.) does NOT determine the classification. 69 pairs have `action_type=INHIBITOR` but are classified `abundance_modulating`" (`protocol/DEVIATION_LOG.md`, Change 4). The two fields answer different questions. `mechanism_class` asks whether the drug's effect runs through circulating protein level. `action_types` asks whether the drug turns the target up or down. This analysis uses the second, and any result must be read with the annotation quality of that second field in mind.

---

## 2. Analysed set

The analysed set is exactly the set of pairs scored by the parent primary analysis: **n = 138** drug–target–indication pairs with an MR result and an adjudicated outcome in {SUCCESS, FAILURE}.

| Stratum | n | Source column |
|---|---|---|
| Pooled | 138 | all analysed pairs |
| Activity-blocking | 76 | `mechanism_class == activity_blocking` |
| Abundance-modulating | 59 | `mechanism_class == abundance_modulating` |
| Mixed | 3 | `mechanism_class == mixed` |

Mixed pairs are retained in pooled and excluded from the two stratified analyses, matching the V6 convention. No pair is added, removed, re-adjudicated, or re-instrumented. Outcomes and MR estimates are carried unchanged from the frozen files. The only pairs whose membership changes are those dropped under Sensitivity Variant B (§3.6), which re-scores both classifiers on the reduced set so the comparison stays like-for-like.

---

## 3. Classifier definition

The rule must be specified tightly enough that two implementers produce identical predictions for all 138 pairs. Each step below is deterministic.

### 3.1 Sign convention for `mr_beta`

`mr_beta` is the Wald ratio `beta_outcome / beta_exposure` (`code/classify_v34.py:111`): the effect of a one-unit increase in genetically predicted circulating protein on the outcome GWAS trait.

The concordance rule assumes each outcome GWAS is coded so that **a positive beta means more disease**. Under that coding:

- `mr_beta > 0` → the target protein is **harmful** (raises disease risk)
- `mr_beta < 0` → the target protein is **protective** (lowers disease risk)

**Required audit before scoring.** Not every outcome GWAS is guaranteed to be coded in the risk-increasing direction. A trait such as eGFR is coded so that higher values mean *better* kidney function, which inverts the interpretation of the sign. Before the classifier is applied, produce a table listing every disease in the analysed set, its outcome GWAS source, its trait definition, and a boolean `flip_sign`. Apply the flip where the source trait is coded in the protective direction. Save the table to `results/direction_concordance/outcome_trait_direction_audit.csv` and report it in the supplement.

Any disease whose trait coding cannot be determined from source metadata is assigned `flip_sign = UNKNOWN`, and every pair for that disease is routed to AMBIGUOUS (§3.4), never guessed.

### 3.2 Parsing `action_types`

`action_types` in `data/adjudicated_v34.csv` may hold multiple tokens per pair, because `n_drugs` can exceed 1. Parse as follows:

1. If the field is empty, null, or whitespace-only, the pair has **no tokens**.
2. Otherwise split the field on its delimiter, strip surrounding whitespace from each element, uppercase, and collapse internal runs of whitespace to a single space.
3. Deduplicate. Token order is discarded.

### 3.3 Token → direction class

Each token is assigned exactly one of BLOCKING, ACTIVATING, or AMBIGUOUS, by evaluating the three tiers **in order** and stopping at the first match.

**Tier 1 — exact-match overrides** (checked first, because these tokens contain substrings that would otherwise misclassify them):

| Token (exact, after normalisation) | Class |
|---|---|
| `INVERSE AGONIST` | BLOCKING |
| `NEGATIVE ALLOSTERIC MODULATOR` | BLOCKING |
| `POSITIVE ALLOSTERIC MODULATOR` | ACTIVATING |
| `PARTIAL AGONIST` | ACTIVATING |

**Tier 2 — BLOCKING substrings.** The token is BLOCKING if it contains any of: `INHIBITOR`, `ANTAGONIST`, `BLOCKER`, `DEGRADER`, `DISRUPT`, `SUPPRESSOR`.

**Tier 3 — ACTIVATING substrings.** The token is ACTIVATING if it contains any of: `AGONIST`, `ACTIVATOR`, `OPENER`, `STABILISER`, `STABILIZER`, `RELEASING AGENT`.

**Tier 4 — everything else is AMBIGUOUS.** This deliberately catches bare `MODULATOR`, `BINDING AGENT`, `OTHER`, `UNKNOWN`, and any token not enumerated above, including tokens that appear in the data but were not anticipated when this table was written. Unlisted tokens are never assigned a direction by inference. The full list of tokens observed in the data, with the class each received, is reported in the supplement so that Tier 4 catches can be audited.

### 3.4 Token set → pair direction requirement

For each pair:

1. Drop AMBIGUOUS tokens from the pair's token set.
2. If at least one directional token remains **and all remaining directional tokens agree**, the pair's requirement is that direction.
   - All BLOCKING → requirement is **POSITIVE** (drug reduces target function, so the target must be harmful, so a concordant MR estimate has `mr_beta > 0`).
   - All ACTIVATING → requirement is **NEGATIVE** (drug increases target function, so the target must be protective, so a concordant MR estimate has `mr_beta < 0`).
3. Otherwise — no directional tokens remain, or the remaining directional tokens disagree, or the pair has no tokens at all, or the pair's disease has `flip_sign = UNKNOWN` — the requirement is **AMBIGUOUS**.

### 3.5 The classifier

Let `concordant(pair)` be:

- requirement POSITIVE and `mr_beta > 0` → TRUE
- requirement NEGATIVE and `mr_beta < 0` → TRUE
- requirement AMBIGUOUS → determined by the variant in §3.6
- otherwise → FALSE

`mr_beta` exactly equal to 0, missing, or NaN is **not concordant**. Report the count of such pairs; it is expected to be zero.

**Direction-concordance classifier:** predict SUCCESS if and only if `mr_p < 0.05` AND `concordant(pair)`. Otherwise predict FAILURE.

### 3.6 Handling AMBIGUOUS pairs — one primary and two sensitivity variants

The treatment of AMBIGUOUS pairs can move the answer on its own, so all three variants are fixed here and all three are reported.

- **Primary — PASS-THROUGH.** AMBIGUOUS pairs keep the frozen classifier's prediction; the direction requirement is not applied to them. This isolates the effect of direction to the pairs where direction is actually defined, and is the variant reported in the main text.
- **Sensitivity A — STRICT.** An AMBIGUOUS pair with `mr_p < 0.05` is predicted FAILURE. Reading: absent directional evidence, there is no licence to predict success.
- **Sensitivity B — DROP.** AMBIGUOUS pairs are removed from the analysed set entirely, and *both* the frozen and the direction classifier are re-scored on the reduced set.

### 3.7 Arithmetic consequences, stated in advance

The rule can only move predictions from SUCCESS to FAILURE. No pair moves from FAILURE to SUCCESS. Therefore **sensitivity is non-increasing and specificity is non-decreasing**, by construction, in every stratum and every variant. Any implementation that produces an increase in sensitivity has a bug.

Let `P` be the number of SUCCESS pairs and `N` the number of FAILURE pairs in a stratum. Let `a` be the number of reclassified pairs whose outcome is SUCCESS (true positives destroyed) and `b` the number whose outcome is FAILURE (false positives removed). Then

```
ΔBA = (b/N − a/P) / 2
```

Two consequences follow, and both should be reported alongside the result:

1. **The null expectation is zero change.** If direction concordance carries no information, discordant pairs are a random subset of the significant pairs, so `a : b` matches the success-to-failure ratio among significant pairs, `b/N ≈ a/P`, and `ΔBA ≈ 0`. A finding of "no change" is therefore the *expected* result under the null and is not evidence of a null.
2. **The rule is asymmetrically fragile.** The parent set is roughly one-third SUCCESS, so `N/P ≈ 2`. Each removed false positive buys about `+1/(2N)`; each destroyed true positive costs about `−1/(2P)`, roughly twice as much. Improvement requires the discordant pairs to be enriched for FAILURE by better than 2:1.

---

## 4. Analysis plan

Conventions match the parent study throughout.

**Strata.** Pooled (n = 138), activity-blocking (n = 76), abundance-modulating (n = 59). Mixed (n = 3) is pooled-only.

**Metrics, reported for both classifiers in every stratum.** Balanced accuracy, sensitivity, specificity, PPV, and the full 2×2 confusion table with raw counts.

**Confidence intervals.** 10,000-iteration bootstrap, 90% percentile interval, resampling pairs within stratum — identical to V3.4/V5/V6.

**Difference in balanced accuracy.** `ΔBA = BA(direction) − BA(frozen)`, computed on **paired** bootstrap resamples: within each of the 10,000 iterations, score both classifiers on the same resampled pairs and record the difference. Report the point estimate and the 90% percentile interval of the difference distribution. The two classifiers are scored on identical rows, so the paired interval is the only correct interval for the comparison; the marginal intervals for the two classifiers will overlap heavily and must not be used to judge the difference.

**Reclassification accounting.** For each stratum and variant, report: the number of pairs with `mr_p < 0.05`; how many of those are concordant, discordant, and AMBIGUOUS; and the SUCCESS/FAILURE split of the reclassified pairs (`a` and `b` from §3.7). These counts are the primary output. The BA comparison is a summary of them.

**Fragility check.** Leave-one-out over reclassified pairs: drop each reclassified pair in turn, recompute `ΔBA`, and report the range. Also report the single largest contribution any one pair makes to `ΔBA`.

**Annotation audit.** Report the token-to-class table actually realised in the data (§3.3), the count of pairs in each requirement category, and the outcome-trait direction audit (§3.1).

**No confirmatory testing.** A permutation test on `ΔBA` may be reported for descriptive purposes. Its p-value has no nominal coverage here, because the analysis was already run once unblinded, and it must be labelled descriptive wherever it appears. No falsification criterion is attached to this analysis, and no result from it can satisfy or fail the study's pre-declared falsification criterion.

**Outputs.** All results to `results/direction_concordance/`, including per-pair predictions under both classifiers so that reclassifications can be inspected individually.

---

## 5. Predictions

Stated before any result is consulted, by an author who has not seen the 2026-08-16 unblinded run. Reference values are the frozen V3.4/V5 figures quoted in `protocol/v5/PRESPEC_V5.md` and `protocol/v6/PRESPEC_V6.md`: activity-blocking BA = 0.511, abundance-modulating BA = 0.590 with 90% CI [0.496, 0.691].

**Scoring note added 2026-08-16, after the predictions were written and before scoring.** P1's first clause is checkable against figures already published in the parent manuscript rather than against the new run: the number of pairs currently predicted SUCCESS is TP + FP = 16 pooled, not the 20–45 predicted. P1 is therefore already missed on that clause, and it is scored as missed. The predictions are left exactly as written — they are not revised to match. Readers should note that P2 and P4 derive their magnitude intervals from the inflated volume assumed in P1, so those intervals are calibrated to a scale roughly three times too large and should be read as directional predictions rather than as interval forecasts.

### P1 — Reclassification volume (pooled)

**Prediction:** between 20 and 45 of the 138 pairs currently carry a SUCCESS prediction, and between 6 and 18 of those are discordant and get reclassified under the PASS-THROUGH variant.

**Reasoning:** The frozen classifier fires only on `mr_p < 0.05`. Under a null of no genetic signal, about 5% of pairs would clear that threshold; the observed rate is higher, since the strata carry some real signal and instruments were selected as genome-wide-significant cis-pQTLs, but it stays well below half. Among those that fire, sign is close to a coin flip in the activity-blocking half of the set and should lean concordant in the abundance-modulating half, giving a discordant fraction somewhat under 50%.

### P2 — Pooled balanced accuracy: **small increase**

**Prediction:** `ΔBA` pooled is positive and small, most likely in **[0.00, +0.04]**, with the paired 90% CI including zero.

**Reasoning:** The pooled set is a mixture. Roughly 55% of it is activity-blocking, where §5.3 predicts direction is uninformative and reclassification removes a base-rate mix of true and false positives for no net gain. The remaining abundance-modulating pairs should contribute a real but small gain (§5.4). Mixing a small positive with a zero yields a smaller positive. The 2:1 asymmetry in §3.7 caps the plausible upside: with `P = 40` pooled, losing even two true positives costs about 0.025 and eats most of what four removed false positives buy.

### P3 — Activity-blocking: **no material change; stays at chance**

**Prediction:** `ΔBA` in the activity-blocking stratum is within **±0.03** of zero, the resulting BA stays within **[0.48, 0.55]**, and its 90% CI continues to include 0.50. Direction concordance does not lift this stratum above chance.

**Reasoning:** This is the sharpest prediction in the document and follows directly from the parent study's causal account. The cis-pQTL instrument moves circulating protein *level*. An active-site inhibitor abolishes catalytic *function* while leaving level intact. The two act on orthogonal axes, so the MR estimate for these pairs carries no information about what the drug will do. An estimate that carries no information in its magnitude carries none in its sign either. Filtering noise on the sign of noise leaves noise. There is a second, mechanical reason: activity-blocking targets are almost entirely INHIBITOR and ANTAGONIST, so the requirement collapses to `mr_beta > 0` for nearly every pair, and among significant-but-uninformative estimates that is a coin flip that discards roughly half of them at base-rate composition.

### P4 — Abundance-modulating: **increase, bounded small**

**Prediction:** `ΔBA` in the abundance-modulating stratum is positive, most likely in **[+0.01, +0.06]**, putting the point estimate in roughly **[0.60, 0.65]**. I predict it does not exceed 0.66, and I predict its 90% CI lower bound still does not clear 0.50 cleanly.

**Reasoning:** Here the drug and the instrument share a causal axis — both act on the amount of circulating protein — so the MR estimate is measuring something close to what the drug does, and its sign should therefore be meaningful. A significant estimate pointing the wrong way is genuine counter-evidence for that drug. The bound on the size of the gain is arithmetic rather than biological: with `P ≈ 19` and `N ≈ 40` in this stratum, each removed false positive buys about +0.013 and each destroyed true positive costs about −0.026, so reaching +0.06 requires roughly five discordant failures reclassified with no true positive lost. With this few significant pairs, that is near the ceiling of what is mechanically available.

### P5 — Ordering across strata

**Prediction:** `ΔBA(abundance) > ΔBA(pooled) > ΔBA(activity)`, with `ΔBA(activity) ≈ 0`.

**Reasoning:** Direction should help exactly where mechanism-class stratification already showed the instrument to be informative, and nowhere else. The pooled value is a size-weighted blend of the two and must land between them. This ordering, rather than any single BA value, is the quantity worth looking at, because it is the prediction the domain-of-validity account actually makes.

### P6 — Sensitivity variants

**Prediction:** STRICT (Sensitivity A) produces a *smaller* `ΔBA` than PASS-THROUGH in the abundance-modulating stratum, and may turn it negative. DROP (Sensitivity B) lands between the two.

**Reasoning:** STRICT reclassifies AMBIGUOUS-annotated significant pairs to FAILURE on no directional evidence. That is a random cut with respect to outcome, so it removes true and false positives at base rate, and given the 2:1 asymmetry a base-rate cut is mildly BA-negative rather than neutral.

**If P3 and P5 both fail — that is, if direction concordance lifts activity-blocking above chance while doing nothing for abundance-modulating — the domain-of-validity account is in trouble and this document says so in advance.**

---

## 6. Interpretation rules, fixed before scoring

Fixed now so that the reading of the result is not chosen after seeing it. Every branch below is subject to the fragility rules in §7 and to the disclosure in §0: no branch licenses a confirmatory claim.

### 6.1 Direction improves performance

*Pooled and abundance-modulating both rise, activity-blocking flat (the predicted pattern).*

**Licensed:** that within the mechanism class where cis-pQTL MR was already informative, the sign of the estimate carries additional information beyond the p-value; that mechanism-of-action directionality — the annotation Ravarani et al. identify as missing — is recoverable and appears useful on this set; that a future pre-registered study should specify a direction-aware classifier from the outset.

**Not licensed:** replacing the frozen classifier or restating the paper's primary result. The rule was written with the frozen results already in hand, and was scored once unblinded before this protocol existed. Reporting the higher number as the study's headline would be selecting a classifier on the same data that scored it. The frozen number stays primary; this one appears as a disclosed secondary with §0 attached.

### 6.2 Direction improves only activity-blocking

**Licensed:** nothing yet. This contradicts the parent causal account, which predicts the instrument is blind on that stratum in magnitude and therefore in sign.

**Required:** treat as hypothesis-generating and report it as an anomaly rather than a finding. Before any claim, check three things: whether the gain survives leave-one-out (§7); whether it is driven by a single disease or a single gene; and whether it replicates on an independent set, such as the Phase III pairs excluded for sample overlap, the pre-declared Phase II tier, or the Ravarani Phase II set. A result of this shape reported without replication would be the exact failure mode this repository's protocol series exists to prevent.

### 6.3 Direction degrades performance

**Licensed:** that the direction rule as specified in §3, applied to `action_types` as annotated, does not improve prediction on this set.

**Not licensed:** the inference that direction is biologically uninformative. Degradation has at least three explanations that this analysis cannot separate: the sign mapping in §3.3 is wrong for some drug class (checkpoint antibodies and shed-ectodomain targets are the obvious candidates, since genetic liability to disease *incidence* need not align with drug effect on *established* disease); the `action_types` annotation is too noisy at the pair level, which V3.4 already warned about (§1.4); or the discordant-significant pairs really are enriched for successes, which would undercut the causal reading of the MR estimate itself.

**Required:** report the direction of the change and the reclassification table, then state the three explanations as unresolved. Do not run further variants hunting for a version that improves.

### 6.4 Direction changes nothing

Two distinct situations produce this, and they must be distinguished before anything is written.

- **Too few pairs reclassify** (below the §7 thresholds): the analysis is uninformative. Report the counts and stop. No conclusion about direction is available, in either direction.
- **Enough pairs reclassify and `ΔBA ≈ 0`**: the discordant pairs split SUCCESS/FAILURE at approximately the base rate. This licenses the statement that direction concordance adds no information beyond significance *on this set at this sample size*, with the reclassification counts reported so the reader can judge the power. Note §3.7: zero change is the null expectation, so this outcome is weak evidence and must not be written up as a positive finding that "the frozen classifier already captures direction."

### 6.5 Activity-blocking stays at chance either way

This is the most likely single outcome and the least informative. Fix the reading now, because it is the branch most open to over-claiming.

**Licensed:** an internal consistency check — the activity-blocking stratum behaves as the domain-of-validity account says it should under a modification designed to help it.

**Not licensed, and this is the important part:** treating it as independent confirmation of the parent finding. It re-uses the same 76 pairs, the same adjudicated outcomes, the same MR estimates, and the same mechanism annotation. It contributes no new data. Worse, the outcome is close to guaranteed by construction: applying a filter to a chance-level classifier yields a chance-level classifier unless the filter itself is informative, and if the filter were informative on this stratum the parent account would already be wrong. The finding cannot fail in a way that would surprise the account, so its success does not support the account. It is reported as a null, in one sentence, in the supplement.

---

## 7. Stopping and fragility rules

### 7.0 Amendment, 2026-08-16 — thresholds corrected before scoring

The thresholds originally written here (8 pooled, 5 per stratum) were set by an author blinded to the analysed data, and are **not achievable**. The frozen confusion matrices already published in the parent manuscript fix the mechanical ceiling:

| Stratum | Pairs with `mr_p < 0.05` (= TP + FP, published) | Max possible reclassifications |
|---|---|---|
| Pooled | 8 + 8 = **16** | 16 |
| Activity-blocking | 2 + 4 = **6** | 6 |
| Abundance-modulating | 5 + 4 = **9** | 9 |

The original pooled threshold of 8 demanded that half of every significant pair be discordant, which is precisely the null expectation for a coin flip on sign — so the analysis would have been voided by chance about half the time. The activity-blocking threshold of 5 demanded 5 of a possible 6. Both would have declared the analysis uninterpretable by construction rather than by evidence.

These counts are taken from the parent manuscript's published tables, not from the 2026-08-16 unblinded run. This amendment is made before scoring and is logged in `protocol/DEVIATION_LOG.md`.

### 7.1 What replaces them

Given a ceiling of sixteen reclassifiable pairs, `ΔBA` cannot carry an interpretable magnitude under any threshold. It is quantised at roughly 0.005–0.011 per pair pooled and 0.013–0.026 per pair in the abundance stratum, so every observable value is one to three pairs. Choosing a cut-off would put a decimal point on that arithmetic rather than fix it.

**The reclassification table is therefore the result, and `ΔBA` is a descriptive summary of it.** Three rules replace the thresholds:

1. **`ΔBA` is never reported without its pair count in the same sentence.** The required form is: *"ΔBA = X, driven by a reclassified SUCCESS pairs and b reclassified FAILURE pairs out of n significant."* A bare `ΔBA` is not a permitted output of this analysis.

2. **Only directional statements are licensed, and only when stable.** The analysis may say direction concordance *helped*, *hurt*, or *did nothing* in a stratum, and only when the leave-one-out range for `ΔBA` in that stratum does not span zero. No magnitude claim, no comparison of effect sizes across strata beyond the ordinal prediction in P5, and no statement that a stratum "reached" or "approached" any particular BA value.

3. **Zero reclassifications is a reportable result.** If no significant pair is discordant in a stratum, that is the finding for that stratum — the sign already agreed everywhere the classifier fired — and it is reported as such rather than as a failed analysis.

**Single-pair dominance and sign-flip rules below are unchanged and now carry the fragility load.** They are data-driven rather than arbitrary, which is why they survive the amendment while the count gates do not.

**Single-pair dominance rule.** If any one reclassified pair accounts for more than half of the observed `ΔBA` in a stratum, label that stratum's result **single-pair-driven**, name the pair, and interpret it no further, whatever the threshold count says.

**Sign-flip rule.** If the leave-one-out range for `ΔBA` in a stratum spans zero, report the result as directionally unstable and do not describe it as an improvement or a degradation.

**No stopping for an unfavourable answer, and no iteration.** The analysis is run once, exactly as specified in §3, and **the result is reported regardless of direction** — improvement, degradation, or nothing. This applies to the primary PASS-THROUGH variant and to both sensitivity variants, all three of which are reported whatever they show. If the result is unfavourable to the direction hypothesis, it appears in the paper or supplement with the same prominence it would have had if favourable. No additional variant, threshold, token mapping, or stratum may be added after the numbers are seen. Any change to §3 after scoring is a new analysis, is logged in `protocol/DEVIATION_LOG.md` with its own disclosure, and does not replace this one.

---

## 8. What this analysis cannot establish

1. **That the frozen classifier should be revised.** The comparison is on the same 138 pairs that scored the frozen classifier, using a rule designed after those scores were known and run before this protocol existed. It measures a difference; it does not select a classifier.

2. **That direction-aware MR predicts Phase III success.** One instrument per pair, one significance threshold, no colocalisation, n = 138, a single mechanism annotation source. Any generalisation requires an independent set.

3. **A calibrated p-value or a confidence interval with nominal coverage.** The analysis was run once unblinded before its protocol was written. Every interval and p-value reported is descriptive.

4. **The contribution of direction separately from the quality of the annotation.** `action_types` supplies the direction requirement. A negative result is jointly a test of the direction hypothesis and of the annotation, and this design cannot decompose the two. Improving the annotation and re-running would be a different, and again post-hoc, analysis.

5. **A resolution of the Ravarani et al. limitation.** They analyse 11,482 Phase II pairs; this is 138 Phase III pairs. Phase II and Phase III outcomes differ in base rate, attrition mechanism, and adjudication. Agreement with them would not be replication, and disagreement would not be refutation.

6. **Anything causal about the drugs.** A reclassified pair is a changed prediction, not evidence about whether the compound would have worked, or about why the trial read out as it did.

7. **Independent support for the domain-of-validity finding.** No branch of §6 adds data. The analysis re-partitions pairs already scored. It can refine the account or embarrass it; it cannot corroborate it.

---

## 9. Integrity statement

This document specifies one exploratory analysis with a fixed classifier definition (§3), a fixed analysis plan (§4), directional predictions (§5), interpretation rules for every outcome branch (§6), and fragility thresholds (§7). The predictions were written by an author who has not seen the results of the 2026-08-16 unblinded run.

That protection is weaker than the protection on V1–V6, and the difference is not cosmetic. The V1–V6 analyses were hashed before their results existed. This one was run first. The disclosure in §0 travels with every report of this analysis, in full, unshortened, and is not to be softened in review response, cover letter, or supplement.

The frozen classifier (MR p < 0.05 → predict SUCCESS) and the V3.4 mechanism-stratified result remain the study's primary findings and are unaffected by anything in this document.
