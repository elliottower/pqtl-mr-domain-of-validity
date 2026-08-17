# Outcome-GWAS provenance trace — full audit

Traced 2026-08-16, before any correction. Covers all 157 analysed pairs (`outcome ∈
{SUCCESS, FAILURE}`) in `results/v5/classification_v5.csv`. No data was modified.

The first audit checked the 23 accessions in `code/classify_v34.py` and found 8 problems.
That audit covered 62% of pairs. This trace covers all of them and finds four further
issues the first pass could not see.

## How an MR estimate reaches a pair

Three paths, none documented together until now.

| path | pairs | % | outcome GWAS comes from |
|---|---|---|---|
| `opengwas_cached` | 83 | 53% | `DISEASE_GWAS` in `classify_v34.py` |
| `catalog` | 40 | 25% | EpiGraphDB MR-EvE, **accession not recorded** |
| `opengwas_cached_epi` | 15 | 10% | `DISEASE_GWAS` in `classify_v34.py` |
| `v5_step3_ld_proxy` | 9 | 6% | `DISEASE_GWAS` in **`expand_v5.py`** (a second dict) |
| `v5_step4_decode` | 6 | 4% | same |
| `v5_step1_new_gwas` | 3 | 2% | same |
| `v5_step4_decode_ld_proxy` | 1 | 1% | same |

## A. Wrong trait (4 diseases, 12 corrupted pairs)

Confirmed against live OpenGWAS metadata.

| disease | accession | actual trait | corrupted pairs |
|---|---|---|---|
| Chronic kidney disease | `ieu-b-4874` | Bladder cancer | 6 of 8 |
| Systemic lupus erythematosus | `ieu-a-1073` | Copper (serum) | 4 of 4 |
| Autism spectrum disorder | `ieu-b-87` | Oral cavity and pharyngeal cancer | 1 of 1 |
| Hypercholesterolemia | `ieu-a-300` | LDL cholesterol | 1 of 1 |
| Melanoma | `ieu-a-62` | Waist circumference | 0 — no pairs use it |

Hypercholesterolemia is separately circular: the outcome GWAS is LDL cholesterol and the
Phase III endpoint for PCSK9 inhibitors is LDL lowering.

## B. Withdrawn accessions (3 diseases, 16 corrupted pairs)

| disease | accession | pairs | note |
|---|---|---|---|
| Pancreatic cancer | `ieu-b-4866` | 10 of 12 | 11 cache hits — data was retrieved before withdrawal |
| Anorexia nervosa | `ieu-b-61` | 3 of 3 | 4 cache hits |
| Glioma | `ieu-b-4987` | 3 of 3 | 3 cache hits |

These are not missing data. The estimates exist in `results/v5/outcome_gwas_cache_v5.json`
and were computed in July against a study that has since been removed from the index. What
is lost is the ability to verify what the study was.

## C. A second `DISEASE_GWAS` dict, never audited

`code/expand_v5.py` line 41 holds its own copy. It is identical to `classify_v34.py`'s
except for two entries the earlier dict left as `None`:

| disease | accession | actual trait | cases / controls | pre-specified source |
|---|---|---|---|---|
| Juvenile idiopathic arthritis | `ebi-a-GCST90018873` | Juvenile rheumatoid arthritis (Sakaue 2021) | **216 / 409,001** | `DEVIATION_LOG.md` names **Hinks 2013** |
| Neuroblastoma | `ieu-a-816` | Neuroblastoma (Capasso 2013) | 1,627 / 3,254 | log names COG; Capasso is COG-affiliated |

Two problems with the JIA entry. It is not the source the deviation log pre-specified, and
Sakaue 2021 is a Biobank Japan + UK Biobank + FinnGen meta-analysis, so its 409,001
controls include UK Biobank. Five of the six JIA pairs use UKB-PPP instruments and none is
flagged `sample_overlap`, which the log's own policy requires.

## D. 40 pairs (25%) have no recorded outcome GWAS

`v34_mr_catalog.csv` was absent from this repository. It is in the original project at
`~/Documents/GitHub/transport-wrapper/DRUGS_EXPANDED/results/v34/`, along with
`outcome_gwas_cache_v34.json` and `evaluation_v34.json`. All 40 catalog-sourced pairs
reproduce their stored p-values exactly from that file, so the numbers are sound.

The catalog records `source = EpiGraphDB_INTERVAL`, `method = Wald`, `n_snps = 1`, and
`rsid`. It does **not** record which outcome GWAS EpiGraphDB used. For 25% of the analysed
set, the outcome GWAS is therefore unidentifiable from any artifact in either repository.

This cuts both ways. Four pairs in the diseases listed in A and B came through the catalog
rather than the dict, so they never touched a wrong accession — which is why the corrupted
count is 28 and not 32.

## E. Power

Effective sample size, `4 / (1/cases + 1/controls)`, after applying the replacements
proposed in F.

| disease | pairs | eff. N |
|---|---|---|
| Juvenile idiopathic arthritis | 6 | **864** |
| Thyroid cancer | 6 | **1,036** |
| Glioma (proposed FinnGen) | 3 | **1,852** |
| Pancreatic cancer (proposed PanScan1) | 12 | **3,835** |
| Neuroblastoma | 3 | 4,339 |
| Anorexia nervosa | 3 | 10,605 |
| all remaining 17 diseases | 121 | 13,220 – 449,856 |

Thirty pairs (19%) sit in diseases below effective N = 5,000; fifteen (10%) below 2,000.
A single-SNP Wald ratio in those diseases will return null nearly regardless of the truth,
so their nulls are being counted as measurements. JIA at 864 is the weakest and was never
audited because it lives in the second dict.

## F. Proposed replacements, each matching a source named in `DEVIATION_LOG.md`

| disease | from | to | consortium | matches log |
|---|---|---|---|---|
| Autism spectrum disorder | `ieu-b-87` | `ieu-a-1185` | iPSYCH-PGC | yes |
| Chronic kidney disease | `ieu-b-4874` | `ieu-a-1102` | CKDGen (Pattaro 2015) | yes |
| Systemic lupus erythematosus | `ieu-a-1073` | `ebi-a-GCST003156` | Bentham 2015 | yes |
| Anorexia nervosa | `ieu-b-61` | `ieu-a-1186` | PGC-ED (Duncan 2017) | yes |
| Pancreatic cancer | `ieu-b-4866` | `ieu-a-822` | PanScan1 (Amundadottir 2009) | yes |
| Glioma | `ieu-b-4987` | `finn-b-C3_BRAIN` | FinnGen | **no** — GICC is absent from OpenGWAS |
| Melanoma | `ieu-a-62` | corrected | — | 0 pairs affected |
| Hypercholesterolemia | `ieu-a-300` | pair dropped | — | circular |
| Juvenile idiopathic arthritis | `ebi-a-GCST90018873` | **unresolved** | log names Hinks 2013 | **no** |

Glioma has no GWAS of that phenotype in OpenGWAS; the only "glioma" hit is `prot-a-1217`,
a protein assay. `finn-b-C3_BRAIN` is "malignant neoplasm of brain", broader than glioma,
and rests on the log's pre-declared FinnGen fallback.

## Damage summary

```
157  analysed pairs
 28  corrupted — queried a wrong or withdrawn accession   (22 FAILURE, 6 SUCCESS)
  3  of those 28 currently classified causal (p < 0.05)
  4  in affected diseases but catalog-sourced, so unaffected
  9  JIA + neuroblastoma, outcome GWAS never audited until now
 40  catalog-sourced, outcome GWAS unidentifiable
```

## Open questions

1. **Glioma** — FinnGen fallback at 464 cases, or drop the 3 pairs? Dropping is
   inconsistent with retaining thyroid cancer (1,036) and JIA (864), which are smaller.
2. **JIA** — Hinks 2013 was pre-specified but the code used Sakaue 2021 at 216 cases with
   UKB controls. Find Hinks in OpenGWAS, flag the 5 UKB-PPP pairs as overlap, or drop?
3. **Low power generally** — is a pre-specified effective-N sensitivity analysis the right
   answer to the 30 pairs below 5,000, rather than per-disease exclusions?
4. **The 40 catalog pairs** — is "EpiGraphDB MR-EvE, outcome GWAS not recorded upstream"
   an acceptable disclosure, or do they need re-computation against known accessions?
5. **Two dicts** — they should be one shared constant. Does that matter beyond hygiene,
   given they agree except on the two entries above?
