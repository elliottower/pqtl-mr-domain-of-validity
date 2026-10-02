"""Constants, exceptions and row schemas for stage B (stages/INTERFACES.md, "Stage B -> B/output/").

Every threshold below is quoted from PREREG.md (freeze b946087); none is a tuning choice.
"""
from typing import Literal

from pydantic import BaseModel, ConfigDict, model_validator

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
# Whole files the collect phase downloads once to the stage B volume (stage_b/collect.py).
CollectSource = Literal["decode", "decode_smp", "ukbppp", "ukbppp_rsid_map", "gwas_catalog"]

# Canonical regional table, one row per variant, alleles upper case, effect on `ea`.
VARIANT_COLUMNS = ["rsid", "chrom", "pos", "ea", "oa", "eaf", "beta", "se", "p", "n"]


# ---- exceptions -------------------------------------------------------------------------------
class StageBError(Exception):
    """Base class for every stage B failure."""


class InputContractError(StageBError):
    """An input file does not satisfy the contract stage B relies on."""


class AmbiguousInstrumentError(StageBError):
    """An instrument cannot be resolved to exactly one regional file."""


class SourceAbsent(StageBError):
    """Definitive absence at the source, the only condition recorded as the plan's 'unavailable':
    HTTP 404 or 410, a file its source's listing does not name, an accession without a harmonised
    file, or, for Ensembl only, its HTTP 400 unknown-identifier message naming a requested rsID
    (an Ensembl 404 is a malformed request, not an absence)."""


class RetryableSourceError(StageBError):
    """A fault that says nothing about whether the source holds the file: authentication or
    authorization (401/403), rate limit (429), server error (5xx), timeout, connection reset, or a
    truncated or corrupt gzip stream. Never recorded as unavailable, however often it repeats: the
    file or unit stays unfinished and a later call retries it. `kind` names the class for the
    status report."""

    def __init__(self, kind: str, message: str):
        super().__init__(kind, message)
        self.kind = kind

    def __str__(self) -> str:
        return f"[{self.kind}] {self.args[1]}"


class CollectError(StageBError):
    """The collect record of a file is missing, or a file on the volume is not the file its
    record (or the pinned listing) describes. Not 'unavailable': the unit stops."""


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


class SourceFile(_Row):
    """A source file as the source's own listing names it. Never a URL: a link or token may be
    re-issued, the file's name, size and ETag stay."""

    name: str
    size: int | None = None
    etag: str = ""


def _not_a_url(value: str, what: str) -> None:
    if any(mark in value for mark in ("://", "/", "?")):
        raise ValueError(f"{what} must be a stable identifier, not a URL or path")


class InstrumentUnit(_Row):
    """One unit of Modal work: an instrument and every outcome GWAS paired with it. The record
    holds identities only (it is hashed into the unit fingerprint): no URL, link or token."""

    unit_key: str
    source: InstrumentSource
    assay_id: str
    gene_symbol: str
    gene_ensembl: str
    platform: Literal["Olink", "SomaScan"]
    sentinel: Sentinel
    pqtl_locator: str                        # UKB-PPP OID, deCODE file name, INTERVAL OpenGWAS id
    pqtl_listing: SourceFile | None = None   # deCODE: the file's record in the pinned folder listing
    smp_listing: SourceFile | None = None    # deCODE SMP-normalized file, S16 only
    outcomes: tuple[OutcomeSpec, ...]

    @model_validator(mode="after")
    def identities_only(self) -> "InstrumentUnit":
        _not_a_url(self.pqtl_locator, "pqtl_locator")
        for listing in (self.pqtl_listing, self.smp_listing):
            if listing is not None:
                _not_a_url(listing.name, "a listed file name")
        return self


class CollectTask(_Row):
    """One whole file for the collect phase: `key` is the SeqId (deCODE), the OID (UKB-PPP), the
    chromosome (UKB-PPP rsID map) or the accession (GWAS Catalog); name, size and ETag where the
    pinned listing gives them."""

    source: CollectSource
    key: str
    name: str = ""
    size: int | None = None
    etag: str = ""


class CollectRecord(_Row):
    """The sidecar of one collect task. `collected`: the file is at `path` (relative to the stage B
    root on the volume) with these bytes and sha256. `absent`: the source does not hold it.
    `remote_indexed`: a GWAS Catalog file with a tabix index, queried by region and not downloaded.
    `source_url` carries no query string and no token."""

    status: Literal["collected", "absent", "remote_indexed"]
    source: CollectSource
    key: str
    name: str = ""
    path: str = ""
    bytes: int | None = None
    sha256: str = ""
    md5: str = ""
    etag: str = ""
    last_modified: str = ""
    source_url: str = ""
    detail: str = ""
    utc: str = ""


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
COLLECTED_FILES_COLUMNS = ["source", "key", "status", "name", "bytes", "sha256", "md5", "etag", "last_modified",
                           "source_url", "detail", "utc"]
