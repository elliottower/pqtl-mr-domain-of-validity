# V7 Amendment A — mechanism tier assignment

**Amends:** `PRESPEC_V7.md` (SHA-256 `2b9173f8f0d3c2ca9a0a1d45be6cd99179fae1d381c0b9e2981cffb33fc4f84d`, frozen 2026-08-16T22:54:16Z)
**Written:** 2026-08-16
**Status:** amendment, not a new version. V7's hypothesis, tiers, ordering, test, and
power analysis are unchanged.

## 0. Standing

V7 §6 permits this. It fixes the protocol against the moment the corrected dataset
exists: *"No element above may change after the corrected dataset exists."* The corrected
dataset does not exist. The eight erroneous outcome-GWAS accessions are unfixed and no
V7 quantity has been computed.

What prompted the amendment was an audit of the *instruction* given to the blinded rater,
not of the results it returned. The instruction told the rater to score "the principal or
most advanced agent," which is two criteria that can disagree, and did not tell the rater
to use the drug the dataset pairs with the gene. The rater scored pegvisomant for GHR
where the dataset names somatropin, imalumab for MIF where the dataset names iguratimod,
and a soluble-protein construct for EPHB4 where the dataset has an unnamed inhibitor for
thyroid cancer. The blinded assignments produced under the defective instruction
(`gene_tier_assignments_blind.md`) are void and are not used.

## 1. Which drug is scored

**The tier is assigned from the drug whose Phase III outcome is the pair's label.** Scoring
a compound other than the one that generated the outcome is incoherent, since the outcome
belongs to a specific drug-indication pair.

Where `approved_drugs` names a drug for the pair, that drug is scored. Where it is empty,
the rater identifies the Phase III compound for that gene in that indication and scores it,
recording the compound by name. Where no compound can be identified, the pair is
UNCLASSIFIED and excluded from the stratified analysis, per V7 §2.1.

This supersedes the "principal or most advanced agent" language.

## 2. Agonists and activators are tier B

The taxonomy's axis is whether the drug changes the quantity the pQTL measures, which is
abundance. It is not the direction in which the drug pushes function. **An abundance
instrument is blind to the direction of functional change and sensitive only to whether
abundance moves.** An agonist and an antagonist are therefore symmetric with respect to
the instrument: both leave circulating abundance unchanged and alter function, which is
V7 §2's definition of B.

This covers 11 AGONIST and 3 ACTIVATOR pairs, 14 of 157. It is stated here rather than
left to rater judgment because V7 §2's examples for B are all inhibitory and a rater could
reasonably read agonism as falling outside the tier.

A consequence worth recording: GHR is B under either candidate compound. Pegvisomant
blocks the receptor and somatropin activates it; neither changes GHR abundance. The
discrepancy between the V3.4 deviation log's stated rationale for GHR, which cites
pegvisomant, and the dataset's `approved_drugs`, which names somatropin, does not affect
the tier.

## 3. Assignment stays gene-level, with one named exception

V7 §3 justifies gene-level permutation on tier being a gene property. Auditing all 37
genes that appear in more than one indication, 5 carry differing `action_types`. The
criterion applied: **a gene requires pair-level splitting only where its pairs name
pharmacologically heterogeneous drugs.** A difference in annotation vocabulary over an
identical drug class does not.

| gene | pairs and drugs | verdict |
|---|---|---|
| CD274 | lung: atezolizumab, durvalumab · ovarian: + avelumab | not split. Nested drug list, all anti-PD-L1 monoclonals. `OTHER` denotes avelumab's Fc effector function, not a different target species. |
| IL6R | JIA: tocilizumab · RA: + sarilumab, levilimab, satralizumab | not split. Nested drug list, all anti-IL6R monoclonals. `ANTAGONIST` is a synonym applied to the same class. |
| PCSK9 | hypercholesterolemia: inclisiran, evolocumab, alirocumab · MI: evolocumab | not split. The heterogeneous row is dropped under §5. |
| IFNAR1 | MS: interferon beta · SLE: anifrolumab | not split. Frozen at D in the V6 map; outside V7 §2.1's scope, which covers unassigned genes only. |
| **IL2RA** | MS: daclizumab · neuroblastoma, SLE: unnamed IL-2 agonist | **split.** Daclizumab blocks CD25; IL-2 activates it. Genuinely heterogeneous and in scope. |

**IL2RA is the sole pair-level exception.** Its MS row is assigned separately from its
neuroblastoma and SLE rows.

**Clustering is unchanged.** Permutation remains gene-level over 20,000 iterations, with
the IL2RA rows treated as two clusters. Cluster-robust standard errors by gene are
reported alongside, as in V7 §3. Moving the whole design to pair-level permutation would
change nothing for 90 of 91 genes and would replace a justified clustering story with a
weaker one.

## 4. Two raters, with agreement reported

V7 §2.1 specifies a single blinded rater. It is amended to **two independent raters,
blinded as V7 §2.1 requires**, assigning in parallel without contact.

The reason is that the dataset cannot do this job. 123 of 157 pairs (78%) carry only
`INHIBITOR`, which does not distinguish an antibody neutralising a circulating cytokine
(D) from a small molecule occupying an active site (B). Modality determines the tier and
modality is not a recorded field. The curated map is therefore analytic content rather
than bookkeeping, and its reliability is a property of the result.

- **Cohen's κ is computed on the three-level tier** across all genes assigned under §5 and
  reported in the amendment record and in the paper's Methods, **whatever its value**.
- **κ ≥ 0.8** is the threshold, matching the adjudication stop condition already in force
  in the V3.4 deviation log. Below 0.8, the assignment is treated as unreliable and the
  stratified analysis is reported with that stated.
- **Disagreements are resolved by a third blinded rater**, not by discussion between the
  two, which would destroy the independence the κ estimates.
- A **15% random sample of the V6 frozen map is re-rated blind** by the same procedure as
  a check on the inherited assignments, and its κ is reported separately. The V6 map is
  not revised on the basis of this check; a low κ is reported as a limitation.

## 5. The gene enumeration is corrected to 15

V7 §2.1 states that the V5 expansion "introduced 13 genes with no assignment" and
enumerates them. The analysed set contains 91 genes: 78 are in the V6 frozen map and 11 of
the 13 named genes appear. **AGER and CLU are in neither.** Fifteen genes lack an
assignment, not thirteen.

This is an error of counting in the frozen document, not a reinterpretation of its scope.
AGER and CLU are assigned by the same blinded procedure as the other thirteen. The
enumeration in V7 §2.1 is read as *all genes in the analysed set absent from the V6 frozen
map*, which is what it was intended to denote.

**PCSK9|Hypercholesterolemia is dropped** on the independent ground that it is circular:
the outcome GWAS is an LDL cholesterol GWAS and the Phase III endpoint is LDL lowering.
PCSK9|Myocardial infarction is retained. This drop is recorded here rather than in the
accession fix because it is a design decision, not a data correction.

## 6. What is unchanged

The three tiers and their ordering (V7 §2). The V6-to-V7 mapping (V7 §2). The Firth
penalised one-sided interaction on continuous evidence, its clustering level, and its
20,000 iterations (V7 §3). The two-instruments justification (V7 §3.1). The predictions
(V7 §4). The power analysis and the 20-pair minimum stratum size (V7 §5). The
interpretation rules and the limits in V7 §7.

## 7. Order of operations

Nothing below may be reordered.

1. This amendment is frozen and committed.
2. Two raters assign the 15 genes and the IL2RA split, blind. κ computed and recorded.
3. Assignments committed.
4. The eight erroneous outcome-GWAS accessions are fixed; the three withdrawn accessions
   are recovered or their pairs dropped.
5. The corrected dataset is built.
6. The analysis is run **once**. The binary decision-rule result and the V7 interaction
   are reported together.
