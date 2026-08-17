# V5.2 — the three pre-specified checks

Computed after the recovery pass, per `protocol/v5/CORRECTIONS_V5_2.md` §3
(SHA-256 `40bf7441…`). Each was specified before the pass ran and is reported here whatever
it returned.

## 3.1 Recovery

**Five of the eighteen pairs recovered.** n moves from 138 to 143.

| pair | proxy | r² |
|---|---|---|
| DPP4 \| Chronic kidney disease | rs6733162 | 0.87 |
| GHR \| Chronic kidney disease | rs6451627 | 0.99 |
| KIT \| Pancreatic cancer | rs218237 | 1.00 |
| PDGFRB \| Pancreatic cancer | rs740750 | 0.99 |
| RET \| Pancreatic cancer | rs2744085 | 0.77 |

Thirteen remain unrecovered: Pancreatic cancer 6, Chronic kidney disease 4, Juvenile
idiopathic arthritis 3. The proxy search ran on all eighteen; thirteen returned no variant
at r² ≥ 0.6 carrying an association in the corrected outcome GWAS.

## 3.2 Residual dropout — the rarity concern is confirmed

Minor allele frequency of the sentinel instrument, from `sentinel_eaf` in
`adjudicated_v34.csv`:

| group | n with MAF recorded | median MAF | mean | range |
|---|---|---|---|---|
| unrecovered | 9 of 13 | **0.100** | 0.107 | 0.001 – 0.453 |
| recovered via proxy | 3 of 5 | 0.379 | 0.315 | 0.124 – 0.443 |
| all retained pairs | 89 of 143 | 0.215 | 0.247 | 0.001 – 0.492 |

Unrecovered instruments sit at roughly half the minor allele frequency of retained ones.
The pattern is consistent with array-era genotyping in the corrected outcome GWAS: PanScan1
(2009), Hinks (2013) and CKDGen (2015) are sparser than the accessions they replace, and a
rare variant is both less likely to be typed and less likely to have a common proxy at
r² ≥ 0.6.

**This survives the check and is therefore stated as a limitation**, with these numbers.
The analysed set is depleted of low-frequency instruments in three diseases, so the
corrected result speaks to common-variant instruments more than the uncorrected one did.

The comparison rests on 9 and 3 observations respectively, since MAF is not recorded for
every pair. The direction is consistent across all three groups but the group medians are
not precisely estimated.

## 3.3 Sample-overlap sensitivity — not evaluable

| stratum | pairs retained |
|---|---|
| clean | 142 |
| overlap_flagged | 1 (DPP4 \| Chronic kidney disease) |

Five of the six `overlap_flagged` pairs lost their instrument in the corrected CKD GWAS and
only DPP4 recovered a proxy. A comparison against a single-pair stratum estimates nothing.

**The sample-overlap sensitivity analysis is reported as not evaluable in the corrected
dataset**, in those words. The row labelled `SENSITIVITY — sample_overlap == clean` in the
run output is 142 of 143 pairs and is not a sensitivity analysis; it is the primary
analysis less one pair, and is not reported as evidence of anything.

## Result after the recovery pass

| stratum | n | BA | 90% CI | criterion |
|---|---|---|---|---|
| Primary | 143 | 0.5783 | [0.5223, 0.6364] | PASS |
| Effective N ≥ 2,000 | 134 | 0.5810 | [0.5243, 0.6415] | PASS |
| abundance_modulating | 65 | 0.6165 | [0.5309, 0.7040] | PASS |
| abundance_modulating, above floor | 64 | 0.6156 | [0.5304, 0.7042] | PASS |
| activity_blocking | 75 | 0.5136 | [0.4554, 0.5826] | fail |
| activity_blocking, above floor | 67 | 0.5147 | [0.4479, 0.5876] | fail |

Mechanism gap 0.1029, against 0.1028 before the recovery pass and 0.0955 before the
accession corrections.
