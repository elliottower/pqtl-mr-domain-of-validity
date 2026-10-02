"""Row schemas of every stage A output file (stages/INTERFACES.md, "Stage A -> A/output/")."""
from datetime import date
from typing import Literal

from pydantic import BaseModel, ConfigDict

Direction = Literal["decrease", "increase", "ambiguous"]
MechanismClass = Literal["aligned", "blocking", "other"]
InstrumentSource = Literal["ukbppp", "decode", "interval"]
OutcomeSource = Literal["gwas_catalog", "finngen", "opengwas"]
Overlap = Literal["no", "unknown", "yes"]
S19Arm = Literal["neutralizing_biologic", "small_molecule_blocker", ""]


class _Row(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class HypothesisRow(_Row):
    """One therapeutic hypothesis. Column order is the order of INTERFACES.md.

    `in_s1`, `in_s4` and `in_s10` mark passing the set's stage A rules before the held-out
    filter: S1 is in_s1 & heldout, the pooled set S2 is in_s1, S4 is in_s4 & heldout, S10 is
    in_s10 & heldout. `in_s12` is formed on the held-out S1 only. S20 is S1 plus the rows with
    strict_instrument false. S8, S9 and S21 re-form hypotheses (in_s8 / in_s9 / in_s21, below).

    `in_s8` / `in_s9` / `in_s21` mark passing every S1 rule, before the held-out filter, with
    hypotheses re-formed under the S8 / S9 class table or the S21 restored-conflicts rule (S8 is
    in_s8 & heldout). `conflicted_row_restored` is true exactly for variant s21 rows. A re-formed hypothesis
    identical to a primary row is that row; any other is a row of its own with `variant` s8, s9
    or s21 (id hashed with the variant) and in_s1 false. `ukbppp_assay_ids` / `decode_assay_ids` are the gene's assays
    in each source's list (the list the row was written under), for the cross-source
    evidence-state agreement of descriptive table 12; `decode_assay_ids` is empty where the
    selected outcome GWAS includes Icelandic participants (never paired with deCODE).
    `s19_arm` (S19) is neutralizing_biologic, small_molecule_blocker or empty, from the same
    molecule type, action type and localization as the class (mechanism.s19_arm), set only when
    every program of the hypothesis is in that arm."""

    hypothesis_id: str
    gene_symbol: str
    gene_ensembl: str
    indication_id: str
    indication_name: str
    direction: Direction
    mechanism_class: MechanismClass
    class_reason: str
    drug_program_ids: str
    eligible: bool
    heldout: bool
    specificity_rank: int
    instrument_source: InstrumentSource
    instrument_assay_id: str
    platform: Literal["Olink", "SomaScan"]
    outcome_accession: str
    outcome_source: OutcomeSource
    outcome_tier: Literal[1, 2, 3]
    overlap: Overlap
    outcome_n_case: int
    outcome_n_control: int
    outcome_neff: float
    blood_secreted_hpa: bool
    secreted_uniprot: bool
    neuro_psych: bool
    neuro_only: bool
    psych_only: bool
    oncology: bool
    pilot_indication: bool
    karim_launched: bool
    single_protein_row: bool
    conflicted_row_restored: bool
    strict_instrument: bool
    in_s1: bool
    in_s4: bool
    in_s10: bool
    in_s12: bool
    subtype_restricted: bool
    phenotype_broader: bool
    sample_overlap: bool
    overlap_unknown: bool
    pre_pqtl_publication_date: date
    variant: Literal["primary", "s8", "s9", "s21"]
    in_s8: bool
    in_s9: bool
    in_s21: bool
    ukbppp_assay_ids: str
    decode_assay_ids: str
    s19_arm: S19Arm


class OutcomeTraitCodingRow(_Row):
    """outcome_trait_coding.tsv: one row per selected outcome accession. Stage B reads the first
    two columns (PREREG §Measured variables, "Genetic effect direction")."""

    outcome_accession: str
    risk_coded: bool
    outcome_source: OutcomeSource
    trait: str
    n_case: int
    n_control: int
    coding_basis: str


class FunnelRow(_Row):
    step: str
    remaining: int
    excluded: int
    reason: str


class OutcomeSelectionRow(_Row):
    """Descriptive table 2: per indication, every candidate by tier, with overlap content."""

    indication_id: str
    indication_name: str
    outcome_accession: str
    outcome_source: OutcomeSource
    trait: str
    tier: Literal[1, 2, 3]
    n_case: int
    n_control: int
    neff: float
    includes_ukb: Overlap
    includes_iceland: Overlap
    includes_interval: Overlap
    selected: bool


class MechanismCrosstabRow(_Row):
    """Descriptive table 4: molecule type x action type x localization x class."""

    molecule_type: str
    action_type: str
    neutralizing_binding_agent: bool
    blood_secreted_hpa: bool
    row_class: MechanismClass
    row_direction: Direction
    row_reason: str
    program_target_class: MechanismClass
    program_target_direction: Direction
    conflicted: bool
    n_mechanism_rows: int
    n_program_targets: int
    n_hypotheses_s1: int


HYPOTHESIS_COLUMNS = list(HypothesisRow.model_fields)
FUNNEL_COLUMNS = list(FunnelRow.model_fields)
TRAIT_CODING_COLUMNS = list(OutcomeTraitCodingRow.model_fields)
OUTCOME_SELECTION_COLUMNS = list(OutcomeSelectionRow.model_fields)
CROSSTAB_COLUMNS = list(MechanismCrosstabRow.model_fields)
