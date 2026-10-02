"""Constants, exceptions and row schemas for stage B (stages/INTERFACES.md, "Stage B -> B/output/").

Every threshold below is quoted from PREREG.md (freeze b946087); none is a tuning choice.
"""
from typing import Literal

from pydantic import BaseModel, ConfigDict

# ---- registered constants (PREREG §Measured variables, §Other planned analysis S15) ----------
PRIMARY_P1 = 1e-4
PRIMARY_P2 = 1e-4
PRIMARY_P12 = 5e-6
S15A_P12 = 1e-6
S15B_P12 = 1e-5
WINDOW_PRIMARY = 500_000          # ±500 kb, primary colocalization window
WINDOW_WIDE = 1_000_000           # ±1 Mb, extraction window and S15c
MIN_SHARED = 50                   # coloc not run below this many shared variants
PP_H4_THRESHOLD = 0.80            # evidence state
PALINDROME_MAF_LOW = 0.42         # palindromic variants with MAF 0.42-0.58 are dropped
PALINDROME_MAF_HIGH = 0.58
LOW_COVERAGE_FRACTION = 0.50      # low coverage: < 50% of pQTL-window variants retained
PROXY_R2 = 0.80                   # LD proxy threshold, 1000 Genomes EUR
SPLICE_SQTL_PP_H4 = 0.80          # splicing flag: sQTL PP.H4 >= 0.80 in any tissue
SPLICE_EQTL_PP_H4 = 0.50          # and eQTL PP.H4 < 0.50 in all three
SPLICE_TISSUES = 3                # liver, whole blood, the gene's highest-median-TPM GTEx v8 tissue

Direction = Literal["decrease", "increase", "ambiguous"]
InstrumentSource = Literal["ukbppp", "decode", "interval"]
PLATFORM = {"ukbppp": "Olink", "decode": "SomaScan", "interval": "SomaScan"}
TABLE12_SOURCES = ("ukbppp", "decode")   # descriptive table 12: agreement where a gene has both
OutcomeSource = Literal["gwas_catalog", "finngen", "opengwas"]
EvidenceState = Literal["supportive", "contradictory", "inconclusive"]
NotRunReason = Literal["regional_file_unavailable", "outcome_file_unavailable",
                       "fewer_than_50_shared", "direction_ambiguous", ""]
PlatformConcordance = Literal["concordant", "discordant", "untested"]
Build = Literal["GRCh37", "GRCh38"]

# Canonical regional table, one row per variant, alleles upper case, effect on `ea`.
VARIANT_COLUMNS = ["rsid", "chrom", "pos", "ea", "oa", "eaf", "beta", "se", "p", "n"]


# ---- exceptions -------------------------------------------------------------------------------
class StageBError(Exception):
    """Base class for every stage B failure."""


class InputContractError(StageBError):
    """An input file does not satisfy the contract stage B relies on."""


class AmbiguousInstrumentError(StageBError):
    """An instrument cannot be resolved to exactly one regional file."""


class RetrievalError(StageBError):
    """A regional file could not be retrieved after retries (the plan's 'unavailable')."""


class LDReferenceError(StageBError):
    """The 1000 Genomes EUR reference cannot supply LD for a region."""


class StaleCheckpointError(StageBError):
    """A unit checkpoint directory was written under another run fingerprint (stage_b/checkpoint.py)."""


class ColocBackendError(StageBError):
    """The R coloc backend failed or returned a malformed result."""


class _Row(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


# ---- inputs -----------------------------------------------------------------------------------
class HypothesisInput(_Row):
    """The only hypotheses.csv columns stage B reads. `eligible` and every other field stay
    unread, so no phase-derived value enters stage B. `ukbppp_assay_ids` / `decode_assay_ids`
    name the gene's assays in those two sources (empty where it has none, or, for deCODE, where
    the selected outcome GWAS is Icelandic); they feed only the cross-source agreement of
    descriptive table 12."""

    hypothesis_id: str
    gene_symbol: str
    gene_ensembl: str
    direction: Direction
    instrument_source: InstrumentSource
    instrument_assay_id: str
    platform: Literal["Olink", "SomaScan"]
    outcome_accession: str
    outcome_source: OutcomeSource
    outcome_n_case: int
    outcome_n_control: int
    ukbppp_assay_ids: str = ""
    decode_assay_ids: str = ""


HYPOTHESIS_INPUT_COLUMNS = list(HypothesisInput.model_fields)


class Sentinel(_Row):
    """A source's reported cis-pQTL sentinel for one assay, in the source's own build."""

    source: InstrumentSource
    assay_id: str
    rsid: str
    chrom: str
    pos: int
    build: Build
    neg_log10_p: float
    locator: str = ""   # UKB-PPP protein id (GENE:UNIPROT:OID:v1), deCODE SeqId, INTERVAL target full name


class OutcomeSpec(_Row):
    accession: str
    source: OutcomeSource
    n_case: int
    n_control: int
    risk_coded: bool


class InstrumentUnit(_Row):
    """One unit of Modal work: an instrument and every outcome GWAS paired with it."""

    unit_key: str
    source: InstrumentSource
    assay_id: str
    gene_symbol: str
    gene_ensembl: str
    platform: Literal["Olink", "SomaScan"]
    sentinel: Sentinel
    pqtl_locator: str             # Synapse entity id (UKB-PPP), URL (deCODE), OpenGWAS id (INTERVAL)
    decode_smp_url: str = ""      # S16 only
    outcomes: tuple[OutcomeSpec, ...]


# ---- outputs ----------------------------------------------------------------------------------
class EvidenceRow(_Row):
    """evidence.csv, column order of INTERFACES.md. Empty strings stand for missing values."""

    hypothesis_id: str
    coloc_run: bool
    not_run_reason: NotRunReason
    pp_h0: float | str
    pp_h1: float | str
    pp_h2: float | str
    pp_h3: float | str
    pp_h4: float | str
    n_shared: int | str
    frac_pqtl_retained: float | str
    frac_outcome_retained: float | str
    sentinel_or_proxy_retained: bool | str
    low_coverage: bool | str
    lead_variant: str
    genetic_direction: Literal[1, -1, 0]
    evidence_state: EvidenceState
    S: Literal[0, 1]
    E: float
    protein_altering: bool | str
    platform_concordant: PlatformConcordance
    splicing_candidate: bool | str
    s15a_pp_h4: float | str
    s15b_pp_h4: float | str
    s15c_pp_h4: float | str
    s15d_pp_h4: float | str
    s15e_pp_h4: float | str
    s15f_pp_h4: float | str
    s15g_pp_h4: float | str
    s15f_low_coverage_excluded: bool | str
    s16_evidence_state: EvidenceState | Literal[""]
    s17_sentinel_p: float | str
    evidence_state_ukbppp: EvidenceState | Literal[""]
    evidence_state_decode: EvidenceState | Literal[""]


EVIDENCE_COLUMNS = list(EvidenceRow.model_fields)


class RegionalManifestRow(_Row):
    source: str
    protein_or_study: str
    window: str
    variants: int
    sha256: str
    status: Literal["ok", "unavailable"]
    detail: str
    retrieved_utc: str


REGIONAL_MANIFEST_COLUMNS = list(RegionalManifestRow.model_fields)
