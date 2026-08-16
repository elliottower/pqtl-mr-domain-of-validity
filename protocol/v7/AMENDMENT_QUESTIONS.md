# V7 amendment — five questions before proceeding

Paused. V7 is frozen (SHA-256 `2b9173f8…`, 2026-08-16). Blinded gene assignment surfaced
problems with it. The corrected dataset still does not exist, so all of this is a method
fix rather than a reaction to results — but it needs deciding before the accession fix.

## What happened

A blinded agent assigned the 13 genes V5 added into the three tiers (D direct / C
correlated / B function-blocking). It returned D=6, C=3, B=4, no UNCLASSIFIED.

Checking those against the dataset's own drug annotation showed the instruction was
faulty. I told the agent to score "the principal or most advanced agent," which is two
criteria that can disagree, and never told it to use the drug the dataset actually pairs
with the gene. Consequences:

- **MIF** was scored on imalumab and ibudilast. The dataset's drug is **iguratimod**, an
  oral small-molecule DMARD. Neither compound the agent weighed.
- **EPHB4** was scored on sEphB4-HSA, a recombinant soluble protein. The dataset pairs
  EPHB4 with thyroid cancer and an unnamed INHIBITOR — almost certainly a multi-kinase
  inhibitor. sEphB4-HSA is not a thyroid cancer drug.
- **GHR** was scored on pegvisomant, a receptor *antagonist*. The dataset's drug is
  **somatropin**, an *agonist* — recombinant growth hormone.

## Four findings

**1. Tier is a property of the drug-target-indication triple, not the gene — but this
affects one pair in the analysed set.** V7 §2 asserts gene-level membership and §3
justifies gene-level permutation on that basis. Of 37 genes appearing in more than one
disease, 5 carry differing `action_types`. Applying a checkable criterion — a gene needs
splitting only if its pairs name *pharmacologically heterogeneous* drugs — two of the five
resolve as artifacts of annotation vocabulary:

| gene | pairs | drugs | verdict |
|---|---|---|---|
| CD274 | lung, ovarian | atezolizumab, durvalumab / + avelumab | nested superset, all anti-PD-L1 mAbs; `OTHER` is avelumab's Fc effector function | 
| IL6R | JIA, RA | tocilizumab / + sarilumab, levilimab, satralizumab | nested superset, all anti-IL6R mAbs; `ANTAGONIST` is a synonym for the same class |

Three are genuinely heterogeneous, and two of those resolve on independent grounds:

| gene | heterogeneity | status |
|---|---|---|
| PCSK9 | inclisiran (siRNA, C) vs evolocumab/alirocumab (mAb, D) | the mixed row is `Hypercholesterolemia`, already queued for dropping as circular; PCSK9 is then D via evolocumab/MI |
| IFNAR1 | interferon beta (agonist) vs anifrolumab (antagonist) | frozen at D in the V6 map; outside V7 §2.1's scope |
| IL2RA | daclizumab (inhibitor, MS) vs unnamed IL-2 agonist (neuroblastoma, SLE) | **open** — one of the 13 in scope |

PCSK9|Hypercholesterolemia is also the only pair in the dataset whose own drug list spans
modalities. The other four multi-type rows (APP, CD274, FOLR1, IL6R) are single-target.

**1b. V7 §2.1's gene enumeration is incomplete.** It states the V5 expansion "introduced
13 genes with no assignment" and lists them. The analysed set has 91 genes; 78 are in the
V6 frozen map and 11 of the 13 named genes appear. **AGER and CLU are in neither** — 15
genes lack an assignment, not 13.

**2. Agonists appear to belong in B.** An agonist activates the target without changing
its circulating abundance, so an abundance instrument is blind to it for the same reason
it is blind to an antagonist. Only whether abundance moves matters, not the direction of
functional change. 11 of 157 pairs (7%) are agonists; 3 more are activators.

**3. The dataset's annotation cannot assign tier for most pairs.** 123 of 157 (78%) carry
only `INHIBITOR`, which does not distinguish an antibody neutralising a circulating
cytokine (D) from a small molecule blocking an active site (B). Modality decides the tier
and modality is not recorded. Assignment therefore requires per-drug lookup, which is what
the V6 hand-curated map did.

**4. Only some pairs name a drug.** 8 of the 13 new genes have an empty `approved_drugs`
field, so the compound has to be identified externally before its modality can be judged.

## The five questions

**Q1. Are agonists correctly assigned to B?** The argument is that abundance is unchanged
so the instrument is blind, identically to an antagonist. Is there a reason an
abundance-based instrument should be differentially informative about agonism?

**Q2. Pair-level assignment, or gene-level with one documented exception?** After the
audit above, exactly one gene in the analysed set is both genuinely heterogeneous and
unresolved by other means: IL2RA. Is gene-level assignment with IL2RA split at the pair
level acceptable, or does the principle require pair-level assignment throughout even
though it changes nothing else?

**Q3. Does V7 §3's gene-level permutation survive?** It is justified by tier being a gene
property. With one gene split, is gene-level permutation still defensible — treating the
IL2RA rows as two clusters — or must inference move to pair-level permutation with
gene-clustered robust standard errors?

**Q6. AGER and CLU have no tier assignment and are named nowhere.** They are absent from
both the V6 frozen map and V7 §2.1's list of 13. Assign them by the same blinded procedure
as the 13, or set them UNCLASSIFIED and exclude? Assigning them is the more informative
choice and is still blind to outcome, but it means the frozen document's enumeration was
wrong and is being extended rather than followed.

**Q4. How should modality lookup be handled given 78% of pairs carry only `INHIBITOR`?**
A hand-curated modality map is unavoidable. Should it be built blind to outcomes by a
single rater, by two raters with agreement reported, or can `action_types` plus drug name
be combined into a mechanical rule for some subset? How should the residual judgment be
disclosed?

**Q5. Is this an amendment to V7 or a new V8?** It changes the assignment procedure and
possibly the clustering justification, both frozen. It is prompted by a blinded process
before the corrected data exists. Amend with a timestamped entry, or supersede?

## Not in question

The three tiers and their ordering. The mechanistic argument for three rather than five or
two. The Firth-penalised one-sided interaction on continuous evidence. The power analysis
and the 20-pair floor. The two-instruments justification in §3.1. None of that is affected
by these findings.

## Sequence, unchanged and still pending

Resolve the above → re-derive assignments → commit → fix the eight bad GWAS accessions →
drop PCSK9|Hypercholesterolemia → re-run once → report binary decision-rule numbers and
the V7 interaction together.
