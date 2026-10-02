"""Harmonization (PREREG §Harmonization, §Genome build).

Variants are matched by rsID, never by position, and no liftover is applied. Outcome effects
are aligned to the pQTL effect allele. Palindromic variants (A/T, C/G) with minor allele
frequency 0.42-0.58 are dropped; other palindromic variants are aligned by frequency. A variant
without an rsID in either file is dropped.
"""
from typing import Literal

import numpy as np
import pandas as pd

from stage_b.schemas import PALINDROME_MAF_HIGH, PALINDROME_MAF_LOW

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
MULTIALLELIC_OTHER = "!"   # deCODE: effect allele tested against all other alleles

HARMONIZED_COLUMNS = ["rsid", "chrom", "pos", "ea", "oa", "eaf_p", "beta_p", "se_p", "p_p", "n_p",
                      "eaf_o", "beta_o", "se_o", "p_o", "n_o"]

DropReason = Literal["allele_mismatch", "palindromic_ambiguous_maf", "palindromic_no_frequency"]


def is_palindromic(a1: str, a2: str) -> bool:
    return len(a1) == 1 and len(a2) == 1 and COMPLEMENT.get(a1) == a2


def minor_allele_frequency(eaf: float) -> float:
    return min(eaf, 1.0 - eaf)


def _complement(allele: str) -> str | None:
    if not allele or any(b not in COMPLEMENT for b in allele):
        return None
    return "".join(COMPLEMENT[b] for b in allele)


def _ambiguous_maf(eaf: float) -> bool:
    m = minor_allele_frequency(eaf)
    return PALINDROME_MAF_LOW <= m <= PALINDROME_MAF_HIGH


def alignment_sign(ea_p: str, oa_p: str, eaf_p: float, ea_o: str, oa_o: str, eaf_o: float) -> int | DropReason:
    """+1 if the outcome effect is already on the pQTL effect allele, -1 if it must be flipped,
    or the reason the variant is dropped."""
    if oa_p == MULTIALLELIC_OTHER:
        if ea_o == ea_p:
            oa_p = oa_o
        elif oa_o == ea_p:
            oa_p = ea_o
        else:
            return "allele_mismatch"
    if is_palindromic(ea_p, oa_p):
        if {ea_o, oa_o} != {ea_p, oa_p}:
            return "allele_mismatch"
        if np.isnan(eaf_p) or np.isnan(eaf_o):
            return "palindromic_no_frequency"
        if _ambiguous_maf(eaf_p) or _ambiguous_maf(eaf_o):
            return "palindromic_ambiguous_maf"
        base = 1 if ea_o == ea_p else -1
        f_same_letter = eaf_o if ea_o == ea_p else 1.0 - eaf_o
        strand_flipped = (f_same_letter - 0.5) * (eaf_p - 0.5) < 0
        return -base if strand_flipped else base
    if ea_o == ea_p and oa_o == oa_p:
        return 1
    if ea_o == oa_p and oa_o == ea_p:
        return -1
    c_ea, c_oa = _complement(ea_o), _complement(oa_o)
    if len(ea_p) == 1 and len(oa_p) == 1 and c_ea is not None and c_oa is not None:
        if c_ea == ea_p and c_oa == oa_p:
            return 1
        if c_ea == oa_p and c_oa == ea_p:
            return -1
    return "allele_mismatch"


def harmonize(pqtl: pd.DataFrame, outcome: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    """Join two canonical tables by rsID and align the outcome to the pQTL effect allele.
    Returns the harmonized table (HARMONIZED_COLUMNS; position and alleles from the pQTL file)
    and a count per drop reason. An rsID that still maps to more than one aligned pair is
    dropped entirely."""
    counts = {"pqtl_no_rsid": int((pqtl["rsid"] == "").sum()),
              "outcome_no_rsid": int((outcome["rsid"] == "").sum()),
              "allele_mismatch": 0, "palindromic_ambiguous_maf": 0, "palindromic_no_frequency": 0,
              "missing_effect": 0, "duplicate_rsid": 0}
    p = pqtl[pqtl["rsid"] != ""]
    o = outcome[outcome["rsid"] != ""]
    m = p.merge(o, on="rsid", suffixes=("_p", "_o"))
    rows = []
    for r in m.itertuples(index=False):
        if not (np.isfinite(r.beta_p) and np.isfinite(r.se_p) and np.isfinite(r.beta_o) and np.isfinite(r.se_o)
                and r.se_p > 0 and r.se_o > 0):
            counts["missing_effect"] += 1
            continue
        s = alignment_sign(r.ea_p, r.oa_p, r.eaf_p, r.ea_o, r.oa_o, r.eaf_o)
        if isinstance(s, str):
            counts[s] += 1
            continue
        rows.append({"rsid": r.rsid, "chrom": r.chrom_p, "pos": r.pos_p, "ea": r.ea_p, "oa": r.oa_p,
                     "eaf_p": r.eaf_p, "beta_p": r.beta_p, "se_p": r.se_p, "p_p": r.p_p, "n_p": r.n_p,
                     "eaf_o": r.eaf_o if s == 1 else 1.0 - r.eaf_o, "beta_o": s * r.beta_o, "se_o": r.se_o,
                     "p_o": r.p_o, "n_o": r.n_o})
    h = pd.DataFrame(rows, columns=HARMONIZED_COLUMNS)
    dup = h["rsid"].duplicated(keep=False)
    counts["duplicate_rsid"] = int(h.loc[dup, "rsid"].nunique())
    h = h[~dup].sort_values("pos", kind="stable").reset_index(drop=True)
    return h, counts
