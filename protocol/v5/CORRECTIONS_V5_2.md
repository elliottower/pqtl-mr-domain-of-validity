# V5.2 — LD-proxy recovery applied symmetrically

**Written:** 2026-08-17, before the recovery pass is run. No V5.2 quantity has been
computed.
**Amends:** `protocol/v5/CORRECTIONS_V5_1.md` (SHA-256 `e7dbccd2…`)
**Scope:** one procedural change. No hypothesis, classifier, threshold, outcome label, or
accession changes.

## 1. The asymmetry being removed

V5 assembled its 157 pairs using an LD-proxy recovery step: where an instrument SNP was
absent from the outcome GWAS, a proxy in linkage disequilibrium was substituted. Nine of
the 157 pairs entered the analysed set that way (`mr_source = v5_step3_ld_proxy`).

The V5.1 correction pass did not apply that step. Eighteen pairs were therefore recorded
as lost on a criterion — presence of the exact sentinel SNP — that the original 157 were
never held to. The corrected accessions are older and more sparsely genotyped (PanScan1
2009, Hinks 2013, CKDGen 2015) than the accessions they replace, so the asymmetry falls
disproportionately on the corrected diseases.

Whether a pair survives is therefore currently determined in part by which run it went
through. That makes the 157-versus-138 comparison uninterpretable, independent of any
effect on the estimate.

## 2. The procedure

The V5 procedure is applied unchanged, to every pair in the analysed set, against whatever
outcome GWAS that pair's disease now maps to. Parameters are taken verbatim from
`find_ld_proxy` in `code/expand_v5.py` and are not re-tuned:

| parameter | value |
|---|---|
| service | ENSEMBL REST, GRCh37 assembly |
| reference panel | `1000GENOMES:phase_3:EUR` |
| r² threshold | 0.6 |
| window | 500 kb |
| candidates considered | top 20 by descending r² |
| selection | the highest-r² proxy with a non-null beta and se in the outcome GWAS |

The Wald ratio uses the proxy's outcome association with the original instrument's
exposure estimate, as in V5. A pair with no qualifying proxy is dropped, and the reason is
recorded per pair.

**Symmetry requirement.** The pass runs over all pairs, not only the 18 lost in V5.1. A
pair whose sentinel is present keeps its sentinel; the proxy search is attempted only on
absence. Pairs already carrying a V5 proxy retain it, since their outcome GWAS is
unchanged.

## 3. What is checked afterwards, and reported either way

Three checks, all specified here rather than after seeing the output.

**3.1 Recovery.** The number of the 18 lost pairs recovered, by disease. Reported whether
recovery is complete, partial, or nil.

**3.2 Residual dropout.** For pairs still unrecovered, the minor allele frequency and
genotyping provenance of the sentinel are reported against the distribution for recovered
pairs. If unrecovered instruments are systematically rarer, that is stated as a limitation
with the numbers behind it. If they are not, no limitation is claimed. The concern is not
written into the manuscript unless it survives this check.

**3.3 Overlap sensitivity.** Whether the five `overlap_flagged` CKD pairs remain coincident
with the dropped pairs. If they do, the sample-overlap sensitivity analysis is reported as
not evaluable, in those words, rather than as a clean result. A sensitivity analysis whose
comparison stratum is empty is not reported as if it functioned.

## 4. Unchanged

The classifier (`mr_p < 0.05` → SUCCESS). The falsification criterion (BA 90% bootstrap
lower bound above 0.50). The effective-N floor of 2,000 and its power justification. Every
accession in `code/disease_gwas.py`. The drop of PCSK9 | Hypercholesterolemia. The
mechanism classification.

## 5. Order

1. This document is frozen and committed.
2. The recovery pass runs over all pairs.
3. The three checks in §3 are computed.
4. The analysis runs **once**. Primary, effective-N sensitivity, overlap sensitivity if
   evaluable, and the mechanism strata.
5. Whatever comes out is reported, including a null recovery.
