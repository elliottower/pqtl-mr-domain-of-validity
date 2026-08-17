# What is this paper now — decision before drafting

Written 2026-08-17 after the mechanism dissociation failed its own test. Nothing is
rewritten until this is settled. The point of the document is to choose a claim, list what
each choice requires, and mark what is already evidenced against what still has to be run.

## Where the evidence actually stands

Verified today against `results/v5_1/classification_v5_1.csv`.

| statement | status |
|---|---|
| Pooled BA 0.5783, 90% CI [0.5223, 0.6364] excludes 0.50 | holds; gene-clustered CI [0.5138, 0.6446] still excludes, by 0.014 |
| abundance_modulating BA 0.6165, CI excludes 0.50 | holds; gene-clustered [0.5122, 0.7200], by 0.012 |
| activity_blocking did not reach the criterion | holds |
| **The two strata differ** | **fails.** Gene-clustered permutation p = 0.2144; bootstrap 90% CI on the gap [−0.0198, +0.2262] |
| Gap outside immune-mediated disease | +0.0148 |
| Gap within UKB-PPP instruments | −0.0349, sign reversed |
| Gap dropping PCSK9, IL6R, IL12B, IL2RA, IFNAR1 | −0.0531, abundance retains zero true positives |
| Gap reclassifying IL6R and IL2RA as receptor blockade | −0.0111 |
| Wald p equals the outcome-GWAS p at the sentinel | 134 of 139 pairs exactly; 5 catalog pairs use a different SE |
| Abundance true positives | 7 pairs, 5 genes, 2 double-counted across correlated outcomes; 4 of 7 from the provenance-unverifiable catalog |
| Applicability funnel | 576 Phase III target genes, 112 with any plasma cis-pQTL instrument |
| **Gap among pairs with identifiable outcome-GWAS provenance** | **+0.0400** (abundance 0.5808 n=38, activity 0.5408 n=62). The 40 catalog pairs split 27 abundance / 13 activity, so the signal concentrates in the subset nobody can verify |
| **Two other pre-specified stratifications give the same-size null** | disease area +0.0941, p = 0.281; instrument source +0.0757, p = 0.385; mechanism +0.1028, p = 0.220 |
| Oncology stratum | n = 49, BA exactly 0.5000, zero true positives, predominantly activity-blocking |

The mechanism claim is not recoverable by reframing. Three independent confounds each
remove it, and the restriction the critique itself proposes — non-immune, UKB-PPP only —
gives 0.5231 against 0.5083.

**Two findings settle it beyond the significance question.** First, the gap is largest in
exactly the pairs whose outcome GWAS cannot be named. Second, the paper selected the
largest of three pre-specified stratifications that are indistinguishable from one another
and individually null; pre-registering all three does not license reporting one.

**The design cannot answer the question at any reachable n.** Detecting a true gap of
0.103 at 80% power needs roughly 345–435 analysed pairs. The funnel puts the entire
universe of Phase III drug-target pairs with a plasma cis-pQTL instrument at 211. The
mechanism hypothesis is not merely untested here; it is untestable with existing data.
Stating that, with both numbers, is itself a contribution and belongs in the paper.

## Four candidate papers

**A. Applicability.** cis-pQTL MR cannot be applied to most drug targets, because most are
not assayed in plasma. Purely descriptive, rests on the funnel, no contested inference.
Unattackable and probably too thin on its own — it quantifies something the field already
assumes. Viable as a component, not as a paper.

**B. The negative result.** A pre-registered evaluation finds cis-pQTL MR does not
discriminate Phase III success, and the apparent mechanism dissociation is confounded with
disease area, with instrument source, and with five well-known loci. Everything needed is
computed. The pre-registration stops being decoration here and becomes load-bearing: a
registered prediction that failed is the one thing a reader cannot attribute to fishing.

**C. The reduced-form observation.** For a single-SNP Wald ratio with a first-order
standard error, the MR p-value is algebraically the outcome-GWAS p-value at the sentinel;
the protein measurement never enters the test. The field's showcase cis-pQTL MR successes
therefore reproduce under a plain cis-window GWAS lookup. Sharper and more novel than B,
and it explains B rather than merely reporting it.

**D. B with C as its mechanism.** The null, plus the structural reason the field believes
otherwise.

**Recommendation: D, led by a positive result rather than by the null.**

The one claim supported by analyses actually performed: the classifier fires 17 times in
143 and is right 10 of them, against a 30.8% base rate. PPV 0.588, a 1.9× lift. NPV 0.730
against a 69.2% failure base rate, a 1.05× lift. **A positive cis-pQTL MR result carries
information; a null carries almost none.** That asymmetry needs no between-stratum
contrast, it is what a translational reader can act on, and it lands on the field's
canonical 2× genetic-support anchor.

Around that: the applicability ceiling (A), the null with its confounds (B), and the
reduced-form explanation (C). B alone says "we failed to find this." D says "here is the
one thing the method does, here is why the field believes it does more, and here is why
the larger question cannot be settled with data that exist."

## What D requires before a word is written

**1. A novelty check on C, and it is a real risk.** The first-order-SE cancellation is
known in the MR literature; the question is whether anyone has published the empirical
consequence — that canonical cis-pQTL MR results reproduce without the protein. If that
is known, C collapses into a methods footnote and the paper falls back to B.
*Not started. Blocking. Do this first.*

**2. The incremental-value analysis, which is the paper's new primary.** Compare the
cis-pQTL MR classifier against a baseline using the same SNP and the outcome-GWAS p-value
alone, no proteomics. Report ΔAUC with a gene-clustered interval. The prediction is that
it is indistinguishable from zero. If it is, C is demonstrated rather than argued.
*Not run. Must be pre-registered first, since it is a new primary analysis.*

**3. Gene-clustered inference throughout, and a multiplicity rule.** The clustered primary
is already computed and still clears 0.50, by 0.014. Every stratum needs the same
treatment, and a pre-specified correction across the strata tested.
*Partly done.*

**4. The confound decomposition promoted to the main text.** Disease area, instrument
source, and the five-locus leave-out. Currently these live nowhere; they are the argument
now.
*Computed, needs presenting.*

**5. Colocalization.** Demanded by the methodological critique, and the abundance true
positives sit in gene-dense immune regions where LD contamination is most likely. Zheng
2020, the instrument source, ran this step and it was dropped here.
*Not run. Expensive. Decide whether the paper needs it or states its absence.*

**6. A decision on the mechanism hypothesis.** It does not vanish; it becomes a hypothesis
the data were underpowered to test, stated once with the interaction estimate and its
interval, and not in the title. Reporting it honestly is a strength in a null paper.

## What must come out of the current draft

The dissociation as a finding. "Domain of validity" as a frame. The word "uninformative"
anywhere. Title options (a) and (c) — (a) asserts mechanism-dependence, (c) asserts a null
from a failure to reject. The abstract's three concessions, which now misstate corrected
numbers that no longer need hedging. The SHA-256 material in the abstract and conclusions,
which belongs in Methods. The residual BMC and GigaScience section headings. The pLoF null
experiment, out of the abstract and into an appendix.

## What the framing gains

The pre-registration inverts. In a positive paper it read as compensation, which is what
the triage critique said. In a null paper it is the whole argument: a hypothesis registered
before the data existed, an audit that corrected eight wrong outcome-GWAS accessions
without the answer moving, and a registered prediction that then failed. That sequence is
hard to fake and harder to dismiss.

## Order

1. Novelty check on C. If it fails, fall back to B and revisit this document.
2. Pre-register the incremental-value analysis, freeze, commit.
3. Run it once.
4. Choose the title from what survives.
5. Then, and only then, draft.

---

## Verified in the review round, added after the four checks above

**The Data Availability statement is false as written.** It says all data are in the public
repository. `v34_mr_catalog.csv`, `outcome_gwas_cache_v34.json` and `evaluation_v34.json`
are not; they live in `transport-wrapper/DRUGS_EXPANDED/results/v34/`. Twenty-eight percent
of the analysed set cannot be regenerated from anything a reader can obtain. Deposit them.

**Both hand-assigned variables were assigned by the sole unblinded author.** Mechanism
class is the exposure of the primary analysis; the outcome label is the endpoint. This is
the κ question from earlier in the day, and it now has a second variable attached to it.
Reclassifying IL6R and IL2RA — both receptor-antagonist antibodies, which the Methods
already concede straddle the categories — moves the gap to −0.011.

**Straw nulls to delete.** Testing specificity 90/98 and sensitivity 8/40 against 0.50
establishes that a classifier firing on 12% of pairs fires rarely. "The asymmetry is
statistically real" does not follow from those tests.

**An analysis is asserted that does not exist.** The robustness section reports a
seven-pair overlap sensitivity as completed, with "analysis not shown," concluding
ΔBA < 0.005. The corrected dataset has one flagged pair and the analysis is not evaluable.

**Multiplicity.** Roughly 59 balanced-accuracy estimates and 25 p-values are reported with
no correction and no declared families. Pre-registration protects against selective
reporting, not against multiplicity: 25 registered tests under a global null still yield a
p < 0.05 with probability about 0.72.

**The sequence of looks.** The falsification criterion has been evaluated at n = 30, 50,
138, 157 and 143, with expansion triggered by failure to clear it, and no alpha-spending
function pre-specified. The manuscript presents this as a virtue. It needs stating as what
it is.

## Revised order

1. Novelty check on C. Blocking.
2. Deposit the three missing catalog files. Non-negotiable and independent of everything else.
3. Pre-register, as one document: the incremental-value analysis against the reduced-form
   baseline; the complete-provenance co-primary; gene-clustered inference throughout; a
   multiplicity rule; and the power/ceiling calculation. Freeze.
4. Run once.
5. Choose a title from what survives. Noun phrase, no verdict.
6. Then draft.

---

## Novelty and scope check — outcome

**C is not novel as algebra and must never be written as one.** The cancellation is a known
property of the single-instrument Wald ratio with a first-order standard error, and it is
why the cis-MR methods literature works with multi-SNP, correlated-instrument and cML
estimators rather than sentinel-SNP ratios. State it in Methods as known, with a citation.

The contribution is the empirical consequence: a corpus of canonical cis-pQTL MR
drug-target findings reproduces under a plain outcome-GWAS sentinel lookup, so the protein
contributes locus selection and not inference. That is an audit of practice, not a claim
about the estimator.

**A is partly pre-empted and needs a citation.** Zheng et al. report the complementary
coverage number — 682 of 1,002 instrumentable proteins overlap the druggable genome — and
name the inability to instrument the proteome as a limitation. The distinction here is the
denominator: not what fraction of instrumentable proteins are druggable, but what fraction
of actual Phase III targets are instrumentable, scored against trial outcomes. Keep the
funnel, cite Zheng as prior art, state the distinction in one sentence.

**B lands on the Nelson 2015 / King 2019 anchor, and that cuts both ways.** A 1.9× lift is
instantly legible. It also invites the question of whether cis-pQTL MR adds anything over
the cheaper genetic-support signal, particularly since King found the >2× lift where causal
genes are clear through coding variants, which is close to what a cis-pQTL sentinel picks.

### The B–C tension, and how it is resolved

If C holds, B's lift is attributable to outcome-GWAS association at a druggable locus
rather than to anything proteomic. Written as separate sections, the first referee to
notice reads a contradiction. Written as one argument, there is none:

> cis-pQTL MR's rule-in value is real and sits at the known genetic-support benchmark. It
> reduces to outcome-GWAS association at a pQTL-selected locus. The set of targets where it
> applies at all is small.

### The comparison: run it as an audit, not as an increment

Asking "does MR beat genetic support" at 17 positive calls is underpowered, and a null
reads as a failed test. Asking "what does the pQTL contribute" is the same computation, and
a null is the result.

The increment is not zero by construction, which is why the comparison is worth running.
The MR classifier tests one pre-specified SNP at p < 0.05; generic genetic support asks
whether any variant in the locus reaches genome-wide significance. The pQTL converts a
multiple-testing problem across the cis window into a single pre-specified test, which is a
real mechanism for increment.

**Baseline to pre-specify:** for each pair, the most strongly associated variant in the
same cis window of the same outcome GWAS, at genome-wide significance and separately at
p < 0.05. Compare against the sentinel-SNP classifier. Report discordant pairs with
McNemar rather than ΔAUC, which will not resolve at this n. Approximately 143 cis-window
queries.
