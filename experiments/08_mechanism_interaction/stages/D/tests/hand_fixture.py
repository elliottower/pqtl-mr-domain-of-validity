"""SYNTHETIC, hand-built stage A/B/C outputs with set memberships worked out by hand (below), for
tests whose oracle must not be stage D's own code. Twelve hypotheses on four genes.

row  gene  blood  ind  class     arm                     programs  rank  heldout in_s1  status_24 (status_12)
r01  G1    yes    I1   aligned   neutralizing_biologic   P1        2     yes     yes    advanced
r02  G1    yes    I2   aligned   ""  (oligonucleotide)   P2        1     yes     yes    no_observed_advancement
r03  G1    yes    I3   blocking  small_molecule_blocker  P3        0     yes     yes    no_phase
r04  G1    yes    I4   aligned   neutralizing_biologic   P1;P4     3     yes     yes    no_observed_advancement
r05  G2    no     I1   blocking  ""                      P5        2     yes     yes    advanced
r06  G2    no     I2   blocking  ""                      P5        1     yes     yes    active (no_observed_advancement)
r07  G2    no     I3   other     ""                      P6        5     yes     yes    no_observed_advancement
r08  G3    yes    I1   blocking  small_molecule_blocker  P7        1     yes     yes    no_observed_advancement
r09  G3    yes    I2   aligned   small_molecule_blocker  P8        1     yes     yes    advanced
r10  G3    yes    I4   aligned   neutralizing_biologic   P9        1     no      yes    advanced
r11  G4    yes    I1   aligned   neutralizing_biologic   P10       1     yes     no     advanced   (in_s4 only)
r12  G4    yes    I2   blocking  small_molecule_blocker  P11       1     yes     yes    business_only

Flags: r01 karim_launched; r06 pilot_indication; r07 psych_only (neuro_psych).

Stage B flags (false unless listed; "missing" is an empty cell):
  protein_altering    r02 true;  r04 missing;  r10 missing (r10 is a pilot key, not in S1)
  low_coverage        r08 true;  r05 missing;  r11 missing (r11 is in S4 only, not in S1)
  splicing_candidate  r09 true;  r06 missing

Hand-worked sets:
  S1  = in_s1 & heldout                      = r01-r09, r12
  S2  = in_s1                                = S1 + r10
  S4  = in_s4 & heldout                      = S1 + r11
  S12 (per gene by rank, indication, lowest program; drop a row whose program is already kept):
        G1: r03 (P3) keep, r02 (P2) keep, r01 (P1) keep, r04 (P1, P4) drop
        G2: r06 (P5) keep, r05 (P5) drop, r07 (P6) keep
        G3: r08 (I1, P7) keep, r09 (I2, P8) keep;  G4: r12 keep
      = r01, r02, r03, r06, r07, r08, r09, r12
  S13 = S1 without karim_launched and pilot_indication = S1 - r01 - r06
  S19 = S1, blood-secreted, arm consistent with class: r01, r03, r04, r08, r12
        (r02 aligned oligonucleotide and r09 aligned-with-blocker-arm excluded; G2 not blood-secreted)
        analysed: r01 (1), r04 (0), r08 (0); excluded: r03 no_phase, r12 business_only
  S22 = S1 without psych_only                = S1 - r07
  Flag restrictions keep a row only when its flag is known to be false; a missing flag is
  excluded and counted, and only S1 rows are counted:
  S5   = S1 with protein_altering false       = S1 - r02 (true) - r04 (missing)
       = r01, r03, r05, r06, r07, r08, r09, r12;   missing excluded: 1 (r04; r10 is outside S1)
  S11  = S1 with splicing_candidate false     = S1 - r09 (true) - r06 (missing)
       = r01, r02, r03, r04, r05, r07, r08, r12;   missing excluded: 1 (r06)
  S15f = S1 with low_coverage false           = S1 - r08 (true) - r05 (missing)
       = r01, r02, r03, r04, r06, r07, r09, r12;   missing excluded: 1 (r05; r11 is outside S1)
  S1 exclusion report (status_24 of rows with no outcome): no_phase 1, active 1, business_only 1
  S14a (12-month): r06 is mature at 12 months, so the report is no_phase 1, business_only 1

Evidence (direction decrease except r07, ambiguous; aligned/blocking always decrease):
  r01 PP.H4 0.78, +1, S15a PP.H4 0.85   primary inconclusive; S15d (0.75) supportive; S15a supportive
  r02 PP.H4 0.95, +1                    supportive at 0.75, 0.80 and 0.90
  r08 PP.H4 0.85, -1                    contradictory at 0.75 and 0.80; inconclusive at 0.90 (S15e)
  others PP.H4 0.10, +1                 inconclusive everywhere
"""
import pandas as pd

ROWS = [  # id, gene, blood, ind, class, arm, programs, rank, heldout, in_s1, status_24, status_12
    ("r01", "G1", True, "I1", "aligned", "neutralizing_biologic", "P1", 2, True, True, "advanced", "advanced"),
    ("r02", "G1", True, "I2", "aligned", "", "P2", 1, True, True, "no_observed_advancement", "no_observed_advancement"),
    ("r03", "G1", True, "I3", "blocking", "small_molecule_blocker", "P3", 0, True, True, "no_phase", "no_phase"),
    ("r04", "G1", True, "I4", "aligned", "neutralizing_biologic", "P1;P4", 3, True, True, "no_observed_advancement",
     "no_observed_advancement"),
    ("r05", "G2", False, "I1", "blocking", "", "P5", 2, True, True, "advanced", "advanced"),
    ("r06", "G2", False, "I2", "blocking", "", "P5", 1, True, True, "active", "no_observed_advancement"),
    ("r07", "G2", False, "I3", "other", "", "P6", 5, True, True, "no_observed_advancement", "no_observed_advancement"),
    ("r08", "G3", True, "I1", "blocking", "small_molecule_blocker", "P7", 1, True, True, "no_observed_advancement",
     "no_observed_advancement"),
    ("r09", "G3", True, "I2", "aligned", "small_molecule_blocker", "P8", 1, True, True, "advanced", "advanced"),
    ("r10", "G3", True, "I4", "aligned", "neutralizing_biologic", "P9", 1, False, True, "advanced", "advanced"),
    ("r11", "G4", True, "I1", "aligned", "neutralizing_biologic", "P10", 1, True, False, "advanced", "advanced"),
    ("r12", "G4", True, "I2", "blocking", "small_molecule_blocker", "P11", 1, True, True, "business_only",
     "business_only"),
]
EVIDENCE = {"r01": (0.78, 1, 0.85), "r02": (0.95, 1, 0.95), "r08": (0.85, -1, 0.85)}   # pp_h4, direction, s15a
S12 = {"r01", "r02", "r03", "r06", "r07", "r08", "r09", "r12"}

EXPECTED = {
    "S1": {"r01", "r02", "r03", "r04", "r05", "r06", "r07", "r08", "r09", "r12"},
    "S2": {"r01", "r02", "r03", "r04", "r05", "r06", "r07", "r08", "r09", "r10", "r12"},
    "S4": {"r01", "r02", "r03", "r04", "r05", "r06", "r07", "r08", "r09", "r11", "r12"},
    "S12": S12,
    "S13": {"r02", "r03", "r04", "r05", "r07", "r08", "r09", "r12"},
    "S19": {"r01", "r03", "r04", "r08", "r12"},
    "S22": {"r01", "r02", "r03", "r04", "r05", "r06", "r08", "r09", "r12"},
    "S5": {"r01", "r03", "r05", "r06", "r07", "r08", "r09", "r12"},
    "S11": {"r01", "r02", "r03", "r04", "r05", "r07", "r08", "r12"},
    "S15f": {"r01", "r02", "r03", "r04", "r06", "r07", "r09", "r12"},
}
# stage B flags: row -> value, "" for a missing flag; every other row is False
PROTEIN_ALTERING = {"r02": True, "r04": "", "r10": ""}
LOW_COVERAGE = {"r08": True, "r05": "", "r11": ""}
SPLICING_CANDIDATE = {"r09": True, "r06": ""}
FLAG_MISSING_EXCLUDED = {"S5": 1, "S11": 1, "S15f": 1}      # among S1 rows: r04; r06; r05
FLAG_MISSING_ANYWHERE = {"S5": 2, "S11": 1, "S15f": 2}      # all twelve rows: r04, r10; r06; r05, r11
S19_OUTCOMES = {"r01": 1.0, "r04": 0.0, "r08": 0.0}
S1_EXCLUDED = {"no_phase": 1, "active": 1, "business_only": 1}
S14A_EXCLUDED = {"no_phase": 1, "business_only": 1}
# (row, set) -> evidence state, worked by hand from the table above
STATES = {("r01", "S1"): "inconclusive", ("r01", "S15d"): "supportive", ("r01", "S15e"): "inconclusive",
          ("r01", "S15a"): "supportive", ("r02", "S1"): "supportive", ("r02", "S15d"): "supportive",
          ("r02", "S15e"): "supportive", ("r08", "S1"): "contradictory", ("r08", "S15d"): "contradictory",
          ("r08", "S15e"): "inconclusive", ("r05", "S15d"): "inconclusive"}


def _state(pp: float, gd: int, direction: str) -> str:
    if direction == "ambiguous" or pp < 0.80:
        return "inconclusive"
    return "supportive" if gd == 1 else "contradictory"


def hand_tables() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    hyp, ev, out = [], [], []
    for (hid, gene, blood, ind, cls, arm, progs, rank, heldout, in_s1, st24, st12) in ROWS:
        direction = "ambiguous" if cls == "other" else "decrease"
        hyp.append({
            "hypothesis_id": hid, "gene_symbol": gene, "gene_ensembl": f"ENSG_{gene}", "indication_id": ind,
            "indication_name": f"indication {ind}", "direction": direction, "mechanism_class": cls,
            "class_reason": "hand", "drug_program_ids": progs, "eligible": True, "heldout": heldout,
            "specificity_rank": rank, "instrument_source": "ukbppp", "instrument_assay_id": f"A_{gene}",
            "platform": "Olink", "outcome_accession": f"GCST_{ind}", "outcome_source": "gwas_catalog",
            "outcome_tier": 1, "overlap": "no", "outcome_n_case": 5000, "outcome_n_control": 50000,
            "outcome_neff": 4 / (1 / 5000 + 1 / 50000) * (1 + rank / 10), "blood_secreted_hpa": blood,
            "secreted_uniprot": blood, "neuro_psych": hid == "r07", "neuro_only": False, "psych_only": hid == "r07",
            "oncology": False, "pilot_indication": hid == "r06", "karim_launched": hid == "r01",
            "single_protein_row": True, "conflicted_row_restored": False, "strict_instrument": True,
            "in_s1": in_s1, "in_s4": in_s1 or hid == "r11", "in_s10": in_s1, "in_s12": hid in S12,
            "subtype_restricted": False, "phenotype_broader": False, "sample_overlap": False,
            "overlap_unknown": hid == "r11", "pre_pqtl_publication_date": "2021-06-01", "variant": "primary",
            "in_s8": in_s1, "in_s9": in_s1, "in_s21": in_s1, "s19_arm": arm})
        pp, gd, s15a = EVIDENCE.get(hid, (0.10, 1, 0.10))
        state = _state(pp, gd, direction)
        e = 0.0 if direction == "ambiguous" else gd * pp
        ev.append({
            "hypothesis_id": hid, "coloc_run": True, "not_run_reason": "", "pp_h0": 0.0, "pp_h1": 0.0, "pp_h2": 0.0,
            "pp_h3": round(1 - pp, 6), "pp_h4": pp, "n_shared": 500, "frac_pqtl_retained": 0.9,
            "frac_outcome_retained": 0.9, "sentinel_or_proxy_retained": True,
            "low_coverage": LOW_COVERAGE.get(hid, False),
            "lead_variant": f"rs{hid}", "genetic_direction": gd, "evidence_state": state,
            "S": int(state == "supportive"), "E": e, "protein_altering": PROTEIN_ALTERING.get(hid, False),
            "platform_concordant": "concordant", "splicing_candidate": SPLICING_CANDIDATE.get(hid, False),
            "s15a_pp_h4": s15a, **{f"s15{x}_pp_h4": pp for x in "bcdefg"},
            "s15f_low_coverage_excluded": False, "s16_evidence_state": "", "s17_sentinel_p": 0.5,
            "evidence_state_ukbppp": state, "evidence_state_decode": ""})
        adv = {"advanced": "1", "no_observed_advancement": "0"}
        out.append({
            "hypothesis_id": hid, "status_24": st24, "advanced_24": adv.get(st24, ""), "status_12": st12,
            "status_36": st24, "status_48": st24, "advanced_12": adv.get(st12, ""), "advanced_36": adv.get(st24, ""),
            "advanced_48": adv.get(st24, ""),
            "stop_code": {"no_observed_advancement": "efficacy", "business_only": "business"}.get(st24, ""),
            "why_stopped_texts": "", "efficacy_coded": st24 == "no_observed_advancement",
            "phase2_start_date": "2015-01-01", "phase3_start_date": "", "time_to_phase3_days": 3000,
            "phase3_event": 0, "phase3_success": "", "approved": False, "chembl_max_phase_for_ind": 2.0,
            "earliest_phase2_start": "2015-01-01", "last_phase2_end_date": "2018-01-01"})
    return pd.DataFrame(hyp), pd.DataFrame(ev), pd.DataFrame(out)
