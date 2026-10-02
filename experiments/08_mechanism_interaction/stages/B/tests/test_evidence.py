import pandas as pd
import pytest

from helpers import table
from stage_b.coloc_backend import AbfResult
from stage_b.evidence import (coverage, direction_match, evidence_score, evidence_state, genetic_direction,
                              lead_sqtl_event, not_run_reason, platform_concordance, protein_altering,
                              sentinel_outcome_p, splicing_candidate, st29_pav, vep_protein_altering)


def abf(pp_h4: float, snp_pp_h4: dict[str, float]) -> AbfResult:
    rest = (1 - pp_h4) / 4
    return AbfResult(pp=(rest, rest, rest, rest, pp_h4), nsnps=len(snp_pp_h4), snp_pp_h4=snp_pp_h4)


def harmonized(rows: list[tuple[str, float, float]]) -> pd.DataFrame:
    return pd.DataFrame([{"rsid": r, "beta_p": bp, "beta_o": bo} for r, bp, bo in rows])


@pytest.mark.parametrize("intervention, genetic, expected", [
    ("decrease", 1, 1), ("decrease", -1, -1), ("increase", -1, 1), ("increase", 1, -1),
    ("ambiguous", 1, 0), ("ambiguous", -1, 0), ("decrease", 0, 0), ("increase", 0, 0),
])
def test_direction_match(intervention, genetic, expected):
    assert direction_match(intervention, genetic) == expected


@pytest.mark.parametrize("run, pp, intervention, genetic, state", [
    (True, 0.80, "decrease", 1, "supportive"),
    (True, 0.95, "increase", -1, "supportive"),
    (True, 0.80, "decrease", -1, "contradictory"),
    (True, 0.99, "increase", 1, "contradictory"),
    (True, 0.7999, "decrease", 1, "inconclusive"),
    (True, 0.99, "ambiguous", 1, "inconclusive"),
    (True, 0.99, "decrease", 0, "inconclusive"),
    (False, None, "decrease", 1, "inconclusive"),
])
def test_evidence_state_table(run, pp, intervention, genetic, state):
    assert evidence_state(run, pp, intervention, genetic) == state


@pytest.mark.parametrize("run, pp, intervention, genetic, e", [
    (True, 0.6, "decrease", 1, 0.6), (True, 0.6, "decrease", -1, -0.6), (True, 0.6, "ambiguous", 1, 0.0),
    (True, 0.6, "increase", 0, 0.0), (False, None, "decrease", 1, 0.0),
])
def test_evidence_score_is_direction_times_pp_h4(run, pp, intervention, genetic, e):
    assert evidence_score(run, pp, intervention, genetic) == pytest.approx(e)


def test_decrease_is_supportive_only_for_a_risk_increasing_protein():
    h = harmonized([("rs1", 0.5, 0.2), ("rs2", 0.5, -0.2)])
    lead_up, gd_up = genetic_direction(abf(0.9, {"rs1": 0.9, "rs2": 0.1}), h, risk_coded=True)
    lead_dn, gd_dn = genetic_direction(abf(0.9, {"rs1": 0.1, "rs2": 0.9}), h, risk_coded=True)
    assert (lead_up, gd_up) == ("rs1", 1)
    assert (lead_dn, gd_dn) == ("rs2", -1)
    assert evidence_state(True, 0.9, "decrease", gd_up) == "supportive"
    assert evidence_state(True, 0.9, "decrease", gd_dn) == "contradictory"
    assert evidence_state(True, 0.9, "increase", gd_dn) == "supportive"


def test_genetic_direction_uses_the_ratio_so_a_negative_pqtl_beta_reverses_it():
    h = harmonized([("rs1", -0.5, 0.2)])
    assert genetic_direction(abf(0.9, {"rs1": 1.0}), h, risk_coded=True) == ("rs1", -1)


def test_accession_not_risk_coded_makes_direction_ambiguous():
    h = harmonized([("rs1", 0.5, 0.2)])
    lead, gd = genetic_direction(abf(0.99, {"rs1": 1.0}), h, risk_coded=False)
    assert (lead, gd) == ("rs1", 0)
    assert evidence_state(True, 0.99, "decrease", gd) == "inconclusive"


def test_lead_variant_ties_break_on_rsid():
    assert abf(0.9, {"rs9": 0.5, "rs10": 0.5}).lead_variant() == "rs10"


@pytest.mark.parametrize("pqtl, outcome, shared, reason", [
    (False, True, 500, "regional_file_unavailable"),
    (False, False, None, "regional_file_unavailable"),
    (True, False, None, "outcome_file_unavailable"),
    (True, True, 49, "fewer_than_50_shared"),
    (True, True, 0, "fewer_than_50_shared"),
    (True, True, 50, ""),
])
def test_not_run_reasons(pqtl, outcome, shared, reason):
    assert not_run_reason(pqtl, outcome, shared) == reason


def test_coverage_fractions_and_low_coverage_boundaries():
    shared = [f"rs{i}" for i in range(50)]
    c = coverage(100, 200, shared, "rs0", set())
    assert (c["n_shared"], c["frac_pqtl_retained"], c["frac_outcome_retained"]) == (50, 0.5, 0.25)
    assert c["sentinel_or_proxy_retained"] and not c["low_coverage"]
    assert coverage(101, 200, shared, "rs0", set())["low_coverage"]


def test_sentinel_absent_but_proxy_present_is_retained():
    shared = [f"rs{i}" for i in range(80)]
    assert not coverage(100, 100, shared, "rsX", {"rsX", "rs7"})["low_coverage"]
    c = coverage(100, 100, shared, "rsX", {"rsX", "rsY"})
    assert not c["sentinel_or_proxy_retained"] and c["low_coverage"]


def test_protein_altering_combines_vep_and_st29():
    recs = [{"transcript_consequences": [{"gene_id": "ENSG1", "consequence_terms": ["intron_variant"]},
                                         {"gene_id": "ENSG2", "consequence_terms": ["missense_variant"]}]}]
    assert not vep_protein_altering(recs, "ENSG1")
    assert vep_protein_altering(recs, "ENSG2.4")
    splice = [{"transcript_consequences": [{"gene_id": "ENSG1", "consequence_terms": ["splice_region_variant"]}]}]
    assert vep_protein_altering(splice, "ENSG1")
    row = {"PAV olink": "Y", "PAV soma": "N"}
    assert st29_pav(row, "Olink") is True and st29_pav(row, "SomaScan") is False and st29_pav(None, "Olink") is None
    assert protein_altering(False, True) is True
    assert protein_altering(True, None) is True
    assert protein_altering(False, False) is False
    assert protein_altering(None, False) == ""


@pytest.mark.parametrize("row, expected", [
    (None, "untested"),
    ({"N platforms tested": 1, "cis pQTL on both and high correlation (> 0.5)": "N"}, "untested"),
    ({"N platforms tested": 2, "cis pQTL on both and high correlation (> 0.5)": "Y"}, "concordant"),
    ({"N platforms tested": 2, "cis pQTL on both and high correlation (> 0.5)": "N"}, "discordant"),
])
def test_platform_concordance(row, expected):
    assert platform_concordance(row) == expected


def test_lead_sqtl_event_is_chosen_by_sqtl_p_within_the_gene_only():
    sq = pd.DataFrame({"gene_id": ["ENSG1", "ENSG1", "ENSG2", "ENSG1"],
                       "molecular_trait_id": ["i1", "i2", "i3", "i2"], "p": [1e-5, 1e-9, 1e-30, 1e-3]})
    assert lead_sqtl_event(sq, "ENSG1") == "i2"
    assert lead_sqtl_event(sq, "ENSG9") is None


T3 = ("liver", "blood", "top")


def tissues(liver, blood, top) -> dict:
    return dict(zip(T3, (liver, blood, top)))


@pytest.mark.parametrize("sq, eq, failed, expected", [
    (tissues(0.85, 0.1, None), tissues(0.2, 0.49, 0.0), False, True),
    (tissues(0.1, None, 0.95), tissues(0.2, 0.49, 0.3), False, True),        # sQTL in the top tissue only
    (tissues(0.85, 0.1, 0.1), tissues(0.2, 0.50, 0.1), False, False),        # 0.50 is not below 0.50
    (tissues(0.79, None, None), tissues(0.2, 0.1, 0.1), False, False),       # no sQTL at 0.80
    (tissues(None, None, None), tissues(0.2, 0.1, 0.1), False, False),       # no sQTL computed anywhere
    # an eQTL colocalization that could not be computed is not an observed PP.H4 below 0.50
    (tissues(0.9, None, None), tissues(None, None, None), False, ""),
    (tissues(0.9, 0.9, 0.9), tissues(0.1, None, 0.1), False, ""),
    (tissues(0.9, 0.9, 0.9), tissues(0.1, float("nan"), 0.1), False, ""),
    (tissues(0.1, 0.1, 0.1), tissues(0.1, None, 0.9), False, ""),            # missing whatever the others say
    ({"liver": 0.9, "blood": 0.9}, {"liver": 0.1, "blood": 0.1}, False, ""),  # only two tissues
    (tissues(0.9, 0.9, 0.9), tissues(0.1, 0.1, 0.1), True, ""),              # eQTL Catalogue query failed
])
def test_splicing_candidate(sq, eq, failed, expected):
    assert splicing_candidate(sq, eq, failed) == expected


def test_sentinel_outcome_p_matches_the_sentinel_rsid_only():
    t = table([{"rsid": "rs1", "p": 0.2}, {"rsid": "rs2", "p": 0.03}, {"rsid": "rs3", "p": 0.5},
               {"rsid": "rs3", "p": 0.01}, {"rsid": "rs4", "p": float("nan")}, {"rsid": "rs5", "p": 0.04},
               {"rsid": "rs5", "p": 0.04}])
    assert sentinel_outcome_p(t, "rs2") == pytest.approx(0.03)
    assert sentinel_outcome_p(t, "rs5") == pytest.approx(0.04)   # the same record listed twice
    assert sentinel_outcome_p(t, "rs3") is None                   # two p-values: not one sentinel record
    assert sentinel_outcome_p(t, "rs4") is None
    assert sentinel_outcome_p(t, "rs9") is None
