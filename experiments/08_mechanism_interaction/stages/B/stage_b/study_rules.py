"""Study-level uncertainty rules: GWAS Catalog studies whose odds-ratio files are read in the p-value
form of coloc.abf (pvalue_coloc) whether or not their confidence limits pass the CI-versus-p check
(outcome_files.choose_uncertainty_mode; 2026-10-08 amendment).

Verma et al. 2024, "Diversity and scale: Genetic architecture of 2068 traits in the VA Million
Veteran Program", Science (PubMed 39024449). Binary traits were tested by a SAIGE-family pipeline
whose published p-value is a saddlepoint-calibrated score test; the odds ratio and its 95% limits
need not encode that test, and the published p is the association test of record. In the
pre-analysis validation of 2026-10-07, 40 of the study's 207 outcome files failed the rounding-aware
CI-versus-p check (85.3% to 94.8% of rows agreeing) and 167 passed; the disagreement is one-sided
(the published p smaller) and grows with the association's strength. Reading the whole study in one
mode keeps the uncertainty representation from depending on a file-level agreement threshold.

The accessions are the stage A outcome accessions whose GWAS Catalog study has pubmed_id 39024449
(https://www.ebi.ac.uk/gwas/rest/api/v2/studies/<accession>, retrieved 2026-10-08T00:28:08Z; the
map of all 374 GWAS Catalog outcome accessions has sha256
2f51e1feb751302844f82994b5cf8b0bedd75f6fade398af2ea702b2fa85c679, in
inputs/stage_b/diagnostics/20261007/gwas_catalog_study_pubmed.json). Every one is an odds-ratio
GWAS-SSF file with ci_lower, ci_upper, p_value and effect_allele_frequency.
"""
from pydantic import BaseModel, ConfigDict


class StudyPvalueForm(BaseModel):
    """One study read in pvalue_coloc: its PubMed ID, the reason, and its outcome accessions."""
    model_config = ConfigDict(frozen=True)

    pubmed_id: int
    reason: str
    accessions: frozenset[str]


VERMA_2024_MVP = StudyPvalueForm(
    pubmed_id=39024449,
    reason=("Verma et al. 2024 MVP (PubMed 39024449): the published p is a SAIGE saddlepoint-calibrated score test "
            "that the odds ratio and its limits need not encode; the whole study is read in the p-value form"),
    accessions=frozenset({
        "GCST90475170", "GCST90475241", "GCST90475272", "GCST90475553", "GCST90475558", "GCST90475567",
        "GCST90475568", "GCST90475575", "GCST90475577", "GCST90475582", "GCST90475583", "GCST90475591",
        "GCST90475593", "GCST90475595", "GCST90475598", "GCST90475603", "GCST90475606", "GCST90475608",
        "GCST90475609", "GCST90475611", "GCST90475612", "GCST90475614", "GCST90475615", "GCST90475617",
        "GCST90475618", "GCST90475647", "GCST90475658", "GCST90475661", "GCST90475667", "GCST90475676",
        "GCST90475686", "GCST90475689", "GCST90475696", "GCST90475715", "GCST90475718", "GCST90475724",
        "GCST90475727", "GCST90475731", "GCST90475740", "GCST90475748", "GCST90475759", "GCST90475765",
        "GCST90475779", "GCST90475781", "GCST90475796", "GCST90475802", "GCST90475807", "GCST90475812",
        "GCST90475817", "GCST90475821", "GCST90475824", "GCST90475825", "GCST90475826", "GCST90475829",
        "GCST90475830", "GCST90475832", "GCST90475844", "GCST90475859", "GCST90475869", "GCST90475872",
        "GCST90475876", "GCST90475877", "GCST90475884", "GCST90475888", "GCST90475895", "GCST90475897",
        "GCST90475909", "GCST90475911", "GCST90475914", "GCST90475922", "GCST90475924", "GCST90475929",
        "GCST90475932", "GCST90475936", "GCST90475949", "GCST90475963", "GCST90475967", "GCST90475981",
        "GCST90475983", "GCST90475984", "GCST90475988", "GCST90475989", "GCST90475992", "GCST90475993",
        "GCST90475994", "GCST90475996", "GCST90475997", "GCST90475999", "GCST90476020", "GCST90476021",
        "GCST90476027", "GCST90476028", "GCST90476029", "GCST90476035", "GCST90476039", "GCST90476045",
        "GCST90476069", "GCST90476070", "GCST90476078", "GCST90476084", "GCST90476086", "GCST90476088",
        "GCST90476089", "GCST90476118", "GCST90476119", "GCST90476123", "GCST90476126", "GCST90476136",
        "GCST90476139", "GCST90476153", "GCST90476155", "GCST90476160", "GCST90476161", "GCST90476165",
        "GCST90476171", "GCST90476173", "GCST90476178", "GCST90476182", "GCST90476183", "GCST90476190",
        "GCST90476200", "GCST90476202", "GCST90476204", "GCST90476207", "GCST90476238", "GCST90476242",
        "GCST90476245", "GCST90476247", "GCST90476258", "GCST90476259", "GCST90476272", "GCST90476274",
        "GCST90476289", "GCST90476745", "GCST90476951", "GCST90476965", "GCST90476971", "GCST90477023",
        "GCST90477037", "GCST90477044", "GCST90477048", "GCST90477056", "GCST90477078", "GCST90477122",
        "GCST90477164", "GCST90477166", "GCST90477167", "GCST90477168", "GCST90477170", "GCST90477173",
        "GCST90477174", "GCST90477187", "GCST90477203", "GCST90477211", "GCST90477217", "GCST90477218",
        "GCST90477247", "GCST90477252", "GCST90477281", "GCST90477292", "GCST90477334", "GCST90477367",
        "GCST90477374", "GCST90477545", "GCST90477551", "GCST90477582", "GCST90477616", "GCST90477644",
        "GCST90477735", "GCST90477776", "GCST90477778", "GCST90477893", "GCST90477896", "GCST90477980",
        "GCST90477996", "GCST90478054", "GCST90478057", "GCST90478110", "GCST90478121", "GCST90478130",
        "GCST90478153", "GCST90478161", "GCST90478256", "GCST90478260", "GCST90478266", "GCST90478353",
        "GCST90478363", "GCST90478419", "GCST90478541", "GCST90478551", "GCST90478614", "GCST90478640",
        "GCST90478649", "GCST90478795", "GCST90478973", "GCST90479124", "GCST90479178", "GCST90479186",
        "GCST90479202", "GCST90479210", "GCST90479216", "GCST90479312", "GCST90479317", "GCST90479399",
        "GCST90479424", "GCST90479433", "GCST90479471",
    }))
STUDY_PVALUE_FORM: tuple[StudyPvalueForm, ...] = (VERMA_2024_MVP,)


def study_pvalue_form(accession: str) -> StudyPvalueForm | None:
    """The study rule that reads `accession` in pvalue_coloc, or None."""
    return next((s for s in STUDY_PVALUE_FORM if accession in s.accessions), None)
