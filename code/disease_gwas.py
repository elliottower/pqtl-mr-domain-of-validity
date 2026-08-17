"""Outcome-GWAS accession per disease, and the effective-N floor.

Single source of truth. `classify_v34.py` and `expand_v5.py` previously held separate
copies of this mapping which diverged on two entries, so Juvenile idiopathic arthritis and
Neuroblastoma were supplied by a dict the accession audit never parsed.

Corrections applied per `protocol/v5/CORRECTIONS_V5_1.md`
(SHA-256 e7dbccd2367398bf70f80bcd4b378777f6821de65ed170d7635f5335e77910b5), each restoring
a source named in `protocol/DEVIATION_LOG.md`. Provenance for every entry is in
`results/OUTCOME_GWAS_TRACE.md`.
"""

DISEASE_GWAS = {
    "Alzheimers disease": "ieu-b-2",
    "Amyotrophic lateral sclerosis": "ebi-a-GCST90027163",
    # was ieu-b-61, withdrawn from the index -> PGC-ED, Duncan 2017
    "Anorexia nervosa": "ieu-a-1186",
    # was ieu-b-87, which is oral cavity and pharyngeal cancer -> iPSYCH-PGC
    "Autism spectrum disorder": "ieu-a-1185",
    "Bipolar disorder": "ieu-b-41",
    # was ieu-b-4874, which is bladder cancer -> CKDGen, Pattaro 2015
    "Chronic kidney disease": "ieu-a-1102",
    "Crohns disease": "ieu-a-12",
    # was ieu-b-4987, withdrawn. GICC is absent from OpenGWAS; FinnGen fallback per
    # DEVIATION_LOG.md. Phenotype is broader than glioma and falls below the N floor.
    "Glioma": "finn-b-C3_BRAIN",
    # dropped as circular: LDL cholesterol GWAS against an LDL-lowering endpoint
    "Hypercholesterolemia": None,
    "Inflammatory bowel disease": "ieu-a-31",
    # was ebi-a-GCST90018873 (Sakaue 2021, 216 cases, UK Biobank controls) -> Hinks 2013,
    # the pre-specified source, which pre-dates UK Biobank
    "Juvenile idiopathic arthritis": "ebi-a-GCST005528",
    "Lung cancer": "ieu-a-984",
    "Major depressive disorder": "ieu-b-102",
    # was ieu-a-62, which is waist circumference. GenoMEL is absent from OpenGWAS and the
    # best remaining option is finn-b-C3_MELANOMA_SKIN at 98 cases, far below the floor.
    # No pairs use this entry, so the wrong accession is removed rather than replaced.
    "Melanoma": None,
    "Multiple sclerosis": "ieu-b-18",
    "Myocardial infarction": "ieu-a-798",
    "Neuroblastoma": "ieu-a-816",
    "Ovarian cancer": "ieu-a-1120",
    # was ieu-b-4866, withdrawn from the index -> PanScan1, Amundadottir 2009
    "Pancreatic cancer": "ieu-a-822",
    "Parkinsons disease": "ieu-b-7",
    "Rheumatoid arthritis": "ieu-a-833",
    "Schizophrenia": "ieu-b-5102",
    # was ieu-a-1073, which is serum copper -> Bentham 2015
    "Systemic lupus erythematosus": "ebi-a-GCST003156",
    "Thyroid cancer": "ieu-a-1082",
    "Ulcerative colitis": "ieu-a-970",
}

# Effective sample size, 4 / (1/ncase + 1/ncontrol), from OpenGWAS metadata.
# Recorded here so the floor can be applied without a live API call.
EFFECTIVE_N = {
    "Alzheimers disease": 57693,
    "Amyotrophic lateral sclerosis": 87381,
    "Anorexia nervosa": 10605,
    "Autism spectrum disorder": 44368,
    "Bipolar disorder": 49367,
    "Chronic kidney disease": 44303,
    "Crohns disease": 46889,
    "Glioma": 1852,
    "Inflammatory bowel disease": 32372,
    "Juvenile idiopathic arthritis": 9266,
    "Lung cancer": 37301,
    "Major depressive disorder": 449856,
    "Multiple sclerosis": 112015,
    "Myocardial infarction": 130309,
    "Neuroblastoma": 4339,
    "Ovarian cancer": 62866,
    "Pancreatic cancer": 3835,
    "Parkinsons disease": 125300,
    "Rheumatoid arthritis": 58622,
    "Schizophrenia": 123451,
    "Systemic lupus erythematosus": 13220,
    "Thyroid cancer": 1036,
    "Ulcerative colitis": 39191,
}

# Pre-specified in CORRECTIONS_V5_1.md §3: the effective N at which a single-SNP Wald
# ratio has 80% power at alpha=0.05 two-sided to detect a per-allele OR of 1.10 at
# MAF 0.30. Diseases below it are reported in a sensitivity stratum, not excluded.
EFFECTIVE_N_FLOOR = 2000


def below_floor(disease: str) -> bool:
    """True when the disease's outcome GWAS cannot detect a typical effect."""
    n = EFFECTIVE_N.get(disease)
    return n is not None and n < EFFECTIVE_N_FLOOR
