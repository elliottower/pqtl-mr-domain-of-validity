import numpy as np
import pytest

from helpers import table
from stage_b.harmonize import alignment_sign, harmonize

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


@pytest.mark.parametrize("ea_o, oa_o, expected", [
    ("A", "G", 1),      # same coding
    ("G", "A", -1),     # alleles swapped
    ("T", "C", 1),      # other strand
    ("C", "T", -1),     # other strand, swapped
    ("A", "C", "allele_mismatch"),
    ("AT", "G", "allele_mismatch"),
])
def test_alignment_sign_non_palindromic(ea_o, oa_o, expected):
    assert alignment_sign("A", "G", 0.3, ea_o, oa_o, 0.3) == expected


@pytest.mark.parametrize("eaf_p, eaf_o", [(0.42, 0.2), (0.58, 0.2), (0.5, 0.2), (0.2, 0.45), (0.2, 0.56)])
def test_palindromic_with_maf_in_042_058_is_dropped(eaf_p, eaf_o):
    assert alignment_sign("A", "T", eaf_p, "A", "T", eaf_o) == "palindromic_ambiguous_maf"


def test_palindromic_just_outside_the_band_is_kept():
    assert alignment_sign("A", "T", 0.41, "A", "T", 0.41) == 1
    assert alignment_sign("A", "T", 0.59, "A", "T", 0.59) == 1


@pytest.mark.parametrize("ea_o, oa_o, eaf_o, expected", [
    ("A", "T", 0.2, 1),    # same letter, same frequency side: same allele
    ("A", "T", 0.8, -1),   # same letter, opposite side: other strand, the letter is the other allele
    ("T", "A", 0.8, -1),   # swapped letters, frequency says the swap is real
    ("T", "A", 0.2, 1),    # swapped letters on the other strand: identical allele
])
def test_palindromic_outside_band_aligned_by_frequency(ea_o, oa_o, eaf_o, expected):
    assert alignment_sign("A", "T", 0.2, ea_o, oa_o, eaf_o) == expected


def test_palindromic_without_frequency_is_dropped():
    assert alignment_sign("C", "G", 0.1, "C", "G", float("nan")) == "palindromic_no_frequency"


def test_decode_multiallelic_other_allele_aligns_on_effect_allele():
    assert alignment_sign("A", "!", 0.1, "A", "C", 0.1) == 1
    assert alignment_sign("A", "!", 0.1, "C", "A", 0.9) == -1
    assert alignment_sign("A", "!", 0.1, "C", "G", 0.9) == "allele_mismatch"


def test_harmonize_flips_beta_and_frequency_and_drops_what_the_plan_drops():
    p = table([{"rsid": "rs1", "pos": 1, "ea": "A", "oa": "G", "eaf": 0.3, "beta": 0.5},
               {"rsid": "rs2", "pos": 2, "ea": "A", "oa": "T", "eaf": 0.45, "beta": 0.5},
               {"rsid": "", "pos": 3},
               {"rsid": "rs4", "pos": 4, "ea": "C", "oa": "T", "eaf": 0.1, "beta": -0.2},
               {"rsid": "rs5", "pos": 5, "ea": "A", "oa": "G"}])
    o = table([{"rsid": "rs1", "ea": "G", "oa": "A", "eaf": 0.7, "beta": 0.3},
               {"rsid": "rs2", "ea": "A", "oa": "T", "eaf": 0.45, "beta": 0.3},
               {"rsid": "rs4", "ea": "G", "oa": "T", "beta": 0.1},
               {"rsid": "rs5", "ea": "A", "oa": "G", "beta": 0.1},
               {"rsid": "rs5", "ea": "A", "oa": "G", "beta": 0.2},
               {"rsid": "rs9", "ea": "A", "oa": "G"}])
    h, counts = harmonize(p, o)
    assert list(h["rsid"]) == ["rs1"]
    assert h.loc[0, "beta_o"] == pytest.approx(-0.3)
    assert h.loc[0, "eaf_o"] == pytest.approx(0.3)
    assert counts["palindromic_ambiguous_maf"] == 1
    assert counts["allele_mismatch"] == 1
    assert counts["pqtl_no_rsid"] == 1
    assert counts["duplicate_rsid"] == 1


def test_harmonize_recovers_the_true_outcome_sign_under_any_reported_coding():
    """Many variants: the outcome file reports each on a random strand and allele order; after
    harmonization every outcome beta must equal the true effect on the pQTL effect allele."""
    rng = np.random.default_rng()
    n = 4000
    bases = np.array(list("ACGT"))
    rows_p, rows_o, truth = [], [], {}
    for i in range(n):
        a1 = rng.choice(bases)
        a2 = rng.choice([b for b in bases if b not in (a1, COMP[a1])])   # never palindromic
        eaf = rng.uniform(0.05, 0.95)
        b_true = rng.normal()
        rs = f"rs{i}"
        truth[rs] = b_true
        rows_p.append({"rsid": rs, "pos": i, "ea": a1, "oa": a2, "eaf": eaf, "beta": rng.normal()})
        ea_o, oa_o, b_o, f_o = a1, a2, b_true, eaf
        if rng.random() < 0.5:
            ea_o, oa_o, b_o, f_o = oa_o, ea_o, -b_o, 1 - f_o
        if rng.random() < 0.5:
            ea_o, oa_o = COMP[ea_o], COMP[oa_o]
        rows_o.append({"rsid": rs, "pos": i, "ea": ea_o, "oa": oa_o, "eaf": f_o, "beta": b_o})
    h, counts = harmonize(table(rows_p), table(rows_o))
    assert len(h) == n
    got = dict(zip(h["rsid"], h["beta_o"]))
    assert all(got[k] == pytest.approx(v) for k, v in truth.items())
    assert h["eaf_o"].to_numpy() == pytest.approx(h["eaf_p"].to_numpy())


def test_rows_without_an_se_are_kept_only_where_the_outcome_mode_uses_no_se():
    p = table([{"rsid": f"rs{i}", "pos": i, "ea": "A", "oa": "G", "eaf": 0.3} for i in range(5)])
    o = table([{"rsid": "rs0", "ea": "A", "oa": "G", "se": float("nan")},                 # no SE: kept in pvalue_coloc only
               {"rsid": "rs1", "ea": "G", "oa": "A", "se": float("nan"), "eaf": 0.8},     # flipped, frequency aligned
               {"rsid": "rs2", "ea": "A", "oa": "G", "se": float("nan"), "p": 0.0},       # p = 0: never in the p-value form
               {"rsid": "rs3", "ea": "A", "oa": "G", "se": float("nan"), "eaf": float("nan")},   # no frequency for the MAF
               {"rsid": "rs4", "ea": "A", "oa": "G", "beta": float("nan")}])               # no direction
    native, counts_native = harmonize(p, o)
    assert list(native["rsid"]) == [] and counts_native["missing_effect"] == 5
    pvalue, counts_pvalue = harmonize(p, o, outcome_se=False)
    assert list(pvalue["rsid"]) == ["rs0", "rs1"] and counts_pvalue["missing_effect"] == 3
    assert pvalue["se_o"].isna().all() and list(pvalue["beta_o"]) == [0.1, -0.1] and list(pvalue["eaf_o"]) == pytest.approx([0.3, 0.2])
