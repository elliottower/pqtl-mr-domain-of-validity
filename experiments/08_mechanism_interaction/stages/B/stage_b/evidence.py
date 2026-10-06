"""Evidence rules of PREREG §Measured variables, as pure functions.

Genetic effect direction, evidence state (S), evidence score (E), coverage, not-run reasons and
the protein_altering, platform_concordant and splicing_candidate flags. Every threshold comes
from schemas.py, where each is quoted from the plan.
"""
import math
from collections.abc import Iterable, Mapping
from typing import Literal

import numpy as np
import pandas as pd

from stage_b.coloc_backend import AbfResult
from stage_b.schemas import (LOW_COVERAGE_FRACTION, MIN_SHARED, PP_H4_THRESHOLD, SPLICE_EQTL_PP_H4,
                             SPLICE_SQTL_PP_H4, SPLICE_TISSUES, Direction, EvidenceState, InputContractError, NotRunReason,
                             PlatformConcordance)

# PREREG §Flags: "a missense, in-frame, stop, frameshift or splice consequence in the target gene".
PROTEIN_ALTERING_EXACT = {"missense_variant", "inframe_insertion", "inframe_deletion", "stop_gained",
                          "stop_lost", "frameshift_variant"}
PROTEIN_ALTERING_SUBSTRING = "splice"


# ---- direction and state ------------------------------------------------------------------------

def genetic_direction(abf: AbfResult, h: pd.DataFrame, risk_coded: bool) -> tuple[str, Literal[1, -1, 0]]:
    """(lead variant, direction). The lead is the variant with the largest SNP.PP.H4 that
    survived harmonization (every variant coloc saw did). Direction is the sign of
    beta_outcome / beta_protein there: +1 risk-increasing protein, -1 protective, 0 when the
    ratio is zero or the outcome accession is not coded so that positive beta means higher risk."""
    lead = abf.lead_variant()
    row = h.loc[h["rsid"] == lead]
    if len(row) != 1:
        raise InputContractError(f"lead variant {lead} is not exactly once in the harmonized table")
    ratio = float(row["beta_o"].iloc[0]) / float(row["beta_p"].iloc[0]) if float(row["beta_p"].iloc[0]) != 0 else 0.0
    if not risk_coded or ratio == 0 or not np.isfinite(ratio):
        return lead, 0
    return lead, 1 if ratio > 0 else -1


def direction_match(intervention: Direction, genetic: int) -> Literal[1, -1, 0]:
    """+1 when the genetic effect points the way the drug acts: a decrease requires a
    risk-increasing protein, an increase a protective one. -1 when opposite, 0 when either
    direction is ambiguous."""
    if intervention == "ambiguous" or genetic == 0:
        return 0
    want = 1 if intervention == "decrease" else -1
    return 1 if genetic == want else -1


def evidence_state(coloc_run: bool, pp_h4: float | None, intervention: Direction, genetic: int,
                   threshold: float = PP_H4_THRESHOLD) -> EvidenceState:
    if not coloc_run or pp_h4 is None or pp_h4 < threshold:
        return "inconclusive"
    d = direction_match(intervention, genetic)
    if d == 0:
        return "inconclusive"
    return "supportive" if d == 1 else "contradictory"


def evidence_score(coloc_run: bool, pp_h4: float | None, intervention: Direction, genetic: int) -> float:
    """E = d x PP.H4; d = 0 when colocalization was not run or a direction is ambiguous."""
    if not coloc_run or pp_h4 is None:
        return 0.0
    return float(direction_match(intervention, genetic) * pp_h4)


def sentinel_outcome_p(outcome: pd.DataFrame, rsid: str) -> float | None:
    """S17 (the pilot classifier): the outcome-association p-value at the instrument's sentinel,
    matched by rsID in the outcome region, direction ignored. None when the sentinel is absent
    from the outcome file or listed with more than one p-value."""
    p = pd.to_numeric(outcome.loc[outcome["rsid"] == rsid, "p"], errors="coerce").astype(float)
    vals = set(p[np.isfinite(p)])
    return float(next(iter(vals))) if len(vals) == 1 else None


# ---- not run and coverage -----------------------------------------------------------------------

def not_run_reason(pqtl_available: bool, outcome_available: bool, n_shared: int | None) -> NotRunReason:
    """PREREG §Colocalization: not run if either regional file is unavailable or fewer than 50
    shared variants. The pQTL file is checked first."""
    if not pqtl_available:
        return "regional_file_unavailable"
    if not outcome_available:
        return "outcome_file_unavailable"
    if n_shared is None or n_shared < MIN_SHARED:
        return "fewer_than_50_shared"
    return ""


def n_variants(df: pd.DataFrame) -> int:
    """Distinct variants (position and alleles), so a variant listed under two rsIDs counts once."""
    return int(len(df[["chrom", "pos", "ea", "oa"]].drop_duplicates()))


def coverage(n_pqtl_window: int, n_outcome_window: int, shared_rsids: Iterable[str], sentinel_rsid: str,
             sentinel_proxies: set[str] | None) -> dict:
    """PREREG §Coverage. Low coverage: under 50% of pQTL-window variants retained, or the
    sentinel and its r^2 >= 0.8 proxies all absent from the shared set. `sentinel_proxies` is None
    where the LD panel has no record in the window (pipeline.py, LDUnavailable): the sentinel itself
    retained is still retained, but its absence leaves `sentinel_or_proxy_retained` missing (None),
    and so `low_coverage` wherever the retained fraction alone does not already make it true."""
    shared = set(shared_rsids)
    n_shared = len(shared)
    frac_p = n_shared / n_pqtl_window if n_pqtl_window else 0.0
    frac_o = n_shared / n_outcome_window if n_outcome_window else 0.0
    if sentinel_rsid in shared:
        retained: bool | None = True
    elif sentinel_proxies is None:
        retained = None
    else:
        retained = bool(shared & sentinel_proxies)
    low = True if frac_p < LOW_COVERAGE_FRACTION else (None if retained is None else not retained)
    return {"n_shared": n_shared, "frac_pqtl_retained": frac_p, "frac_outcome_retained": frac_o,
            "sentinel_or_proxy_retained": retained, "low_coverage": low}


# ---- flags --------------------------------------------------------------------------------------

def vep_protein_altering(vep_records: list[dict], gene_ensembl: str) -> bool:
    """True if any VEP transcript consequence in the target gene is protein altering."""
    gene = gene_ensembl.split(".")[0]
    for rec in vep_records:
        for tc in rec.get("transcript_consequences", []) or []:
            if str(tc.get("gene_id", "")).split(".")[0] != gene:
                continue
            for term in tc.get("consequence_terms", []) or []:
                if term in PROTEIN_ALTERING_EXACT or PROTEIN_ALTERING_SUBSTRING in term:
                    return True
    return False


def st29_pav(row: Mapping | None, platform: Literal["Olink", "SomaScan"]) -> bool | None:
    """Eldjarn ST29 `PAV olink` / `PAV soma` for the instrument's platform; None if the protein
    is absent from ST29."""
    if row is None:
        return None
    return str(row["PAV olink" if platform == "Olink" else "PAV soma"]).strip().upper() == "Y"


def protein_altering(vep_hit: bool | None, pav: bool | None) -> bool | Literal[""]:
    """Lead variant or a proxy has a protein-altering consequence in the target gene, or ST29 PAV
    = Y. Missing ("") only when VEP could not be evaluated and ST29 does not say Y."""
    if vep_hit is True or pav is True:
        return True
    if vep_hit is None:
        return ""
    return False


def platform_concordance(row: Mapping | None) -> PlatformConcordance:
    """Eldjarn ST29 `cis pQTL on both and high correlation (> 0.5)`; proteins absent from ST29 or
    tested on one platform are untested, not discordant."""
    if row is None:
        return "untested"
    n_tested = row.get("N platforms tested")
    if n_tested is None or str(n_tested).strip() in ("", "NA") or int(float(n_tested)) < 2:
        return "untested"
    return "concordant" if str(row["cis pQTL on both and high correlation (> 0.5)"]).strip().upper() == "Y" else "discordant"


def lead_sqtl_event(sqtl: pd.DataFrame, gene_ensembl: str) -> str | None:
    """PREREG §Splicing flag: the leafcutter intron with the smallest sQTL p-value in the window
    among the gene's events, chosen by the sQTL statistics alone (ties: lower event ID)."""
    gene = gene_ensembl.split(".")[0]
    g = sqtl[sqtl["gene_id"].astype(str).str.split(".").str[0] == gene]
    g = g[np.isfinite(g["p"].astype(float))]
    if g.empty:
        return None
    best = g.sort_values(["p", "molecular_trait_id"], kind="stable").iloc[0]
    return str(best["molecular_trait_id"])


def splicing_candidate(sqtl_pp_h4: Mapping[str, float | None], eqtl_pp_h4: Mapping[str, float | None],
                       query_failed: bool) -> bool | Literal[""]:
    """Set when PP.H4(pQTL, lead sQTL) >= 0.80 in any of the three tissues and PP.H4(pQTL, eQTL)
    < 0.50 in all three. "Below 0.50 in all three" needs an observed PP.H4 in all three: the flag
    is missing ("") unless `eqtl_pp_h4` holds three finite values, so a tissue where the eQTL
    colocalization could not be computed (no eQTL variants, or fewer than 50 shared) makes the
    flag missing, and the frozen consequence for S11 applies (run on the instruments where the
    flag was computed, missing count reported). Missing also when an eQTL Catalogue query failed.
    A tissue whose sQTL colocalization could not be computed counts as not >= 0.80."""
    if query_failed:
        return ""
    eqtl = [v for v in eqtl_pp_h4.values() if v is not None and math.isfinite(v)]
    if len(eqtl_pp_h4) != SPLICE_TISSUES or len(eqtl) != SPLICE_TISSUES:
        return ""
    any_sqtl = any(v is not None and math.isfinite(v) and v >= SPLICE_SQTL_PP_H4 for v in sqtl_pp_h4.values())
    return any_sqtl and all(v < SPLICE_EQTL_PP_H4 for v in eqtl)
