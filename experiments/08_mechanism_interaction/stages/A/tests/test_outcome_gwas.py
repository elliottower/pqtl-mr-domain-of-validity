import random

import pandas as pd
import pytest

from stage_a.outcome_gwas import combine_candidates, effective_n, finngen_candidates, gwas_catalog_candidates, \
    indication_tiers, opengwas_candidates, select_outcome
from tests.synth import disease_row, gcat_study


def _gcat(*studies) -> pd.DataFrame:
    return gwas_catalog_candidates([(pd.DataFrame([s for s, _ in studies]), pd.DataFrame([a for _, a in studies]),
                                     "published")])


def test_gwas_catalog_overlap_from_listed_cohorts():
    c = _gcat(gcat_study("G1", "MONDO_1", "EstBB|FinnGen", 1000, 9000),
              gcat_study("G2", "MONDO_1", "UKB", 1000, 9000),
              gcat_study("G3", "MONDO_1", "", 1000, 9000),
              gcat_study("G4", "MONDO_1", "", 1000, 9000, country="Iceland"),
              gcat_study("G5", "MONDO_1", "INTERVAL", 1000, 9000)).set_index("accession")
    assert tuple(c.loc["G1", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("no", "no", "no")
    assert tuple(c.loc["G2", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("yes", "no", "no")
    assert tuple(c.loc["G3", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("unknown",) * 3
    assert tuple(c.loc["G4", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("unknown", "yes", "unknown")
    assert c.loc["G5", "includes_interval"] == "yes"
    assert c.loc["G1", "mapped_ids"] == ["MONDO_1"]


def test_opengwas_is_never_independent_unless_finngen():
    gcat = _gcat(gcat_study("GCST9", "MONDO_9", "EstBB", 1000, 9000))
    fgn = finngen_candidates(pd.DataFrame([{"studyId": "FINNGEN_R12_K1", "diseaseIds": ["MONDO_2"],
                                            "traitFromSourceMappedIds": None, "nCases": 100, "nControls": 900,
                                            "traitFromSource": "k", "hasSumstats": True}]))
    info = {"ieu-a-7": {"trait": "x", "ncase": 5, "ncontrol": 5, "population": "European", "consortium": "EstBB"},
            "ukb-b-1": {"trait": "x", "ncase": 5, "ncontrol": 5, "population": "European", "consortium": ""},
            "finn-b-K1": {"trait": "x", "ncase": 5, "ncontrol": 5, "population": "European", "consortium": "FinnGen"},
            "ebi-a-GCST9": {"trait": "x", "ncase": 5, "ncontrol": 5, "population": "European", "consortium": ""},
            "prot-a-1": {"trait": "x"}}
    o = opengwas_candidates(info, gcat, fgn, {}).set_index("accession")
    assert "prot-a-1" not in o.index
    assert tuple(o.loc["ieu-a-7", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("unknown",) * 3
    assert o.loc["ukb-b-1", "includes_ukb"] == "yes"
    assert tuple(o.loc["finn-b-K1", ["includes_ukb", "includes_iceland", "includes_interval"]]) == ("no",) * 3
    assert o.loc["finn-b-K1", "mapped_ids"] == ["MONDO_2"]
    assert o.loc["ebi-a-GCST9", "mapped_ids"] == ["MONDO_9"]


def test_opengwas_label_match_is_exact_name_or_synonym():
    labels = {"heart failure": {"MONDO_HF"}, "cardiac failure": {"MONDO_HF"}}
    info = {"ieu-a-1": {"trait": "Cardiac failure", "ncase": 5, "ncontrol": 5, "population": "European"},
            "ieu-a-2": {"trait": "heart failure or related", "ncase": 5, "ncontrol": 5, "population": "European"}}
    empty = pd.DataFrame(columns=["accession", "mapped_ids"])
    o = opengwas_candidates(info, empty, empty, labels).set_index("accession")
    assert o.loc["ieu-a-1", "mapped_ids"] == ["MONDO_HF"]
    assert o.loc["ieu-a-2", "mapped_ids"] == []


def test_candidate_qualification_uses_effective_n_floor():
    assert effective_n(1000, 1000) == pytest.approx(2000.0)
    c = combine_candidates(_gcat(gcat_study("A", "M", "x", 1000, 1000), gcat_study("B", "M", "x", 999, 1000),
                                 gcat_study("C", "M", "x", 9000, 9000, ancestry="East Asian")),
                           finngen_candidates(pd.DataFrame(columns=["studyId"])),
                           pd.DataFrame(columns=["accession"]))
    assert dict(zip(c["accession"], c["qualifies"])) == {"A": True, "B": False, "C": False}


def test_tiers_exact_obsolete_descendant_ancestor():
    dis = pd.DataFrame([disease_row("MONDO_P", "p", desc=["MONDO_C"], anc=["MONDO_G"], obs=["MONDO:OLD"]),
                        disease_row("MONDO_C", "c", anc=["MONDO_P", "MONDO_G"])])
    c = combine_candidates(_gcat(gcat_study("EX", "MONDO_P", "x", 9000, 9000),
                                 gcat_study("OB", "MONDO_OLD", "x", 9000, 9000),
                                 gcat_study("DE", "MONDO_C", "x", 9000, 9000),
                                 gcat_study("AN", "MONDO_G", "x", 9000, 9000),
                                 gcat_study("UN", "MONDO_Z", "x", 9000, 9000)),
                           finngen_candidates(pd.DataFrame(columns=["studyId"])), pd.DataFrame(columns=["accession"]))
    acc = dict(enumerate(c["accession"]))
    t = indication_tiers({"MONDO_P"}, c, dis)["MONDO_P"]
    assert {k: sorted(acc[i] for i in v) for k, v in t.items()} == {1: ["EX", "OB"], 2: ["DE"], 3: ["AN"]}


def _cands(*rows) -> pd.DataFrame:
    cols = ["source", "accession", "neff", "includes_ukb", "includes_iceland", "includes_interval"]
    return pd.DataFrame([dict(zip(cols, r)) for r in rows])


def test_selection_order_tier_then_independence_then_n_then_source_then_accession():
    c = _cands(("gwas_catalog", "T2", 9e9, "no", "no", "no"),        # 0 tier 2
               ("gwas_catalog", "UNK", 9e6, "unknown", "no", "no"),  # 1 tier 1 unknown, larger N
               ("opengwas", "B", 1e4, "no", "no", "no"),             # 2 tier 1 no
               ("finngen", "Z", 1e4, "no", "no", "no"),              # 3 tier 1 no, same N, FinnGen before OpenGWAS
               ("finngen", "A", 1e4, "no", "no", "no"))              # 4 lower accession wins the tie
    sel = select_outcome(["ukbppp"], {1: [1, 2, 3, 4], 2: [0]}, c)
    assert (sel.candidate, sel.tier, sel.overlap, sel.rule) == (4, 1, "no", "first_source_tier1_no")
    c2 = c.assign(neff=[9e9, 9e6, 1e4, 2e4, 1e4])
    assert select_outcome(["ukbppp"], {1: [1, 2, 3, 4], 2: [0]}, c2).candidate == 3


def test_first_source_with_tier1_no_is_used_else_first_source_fallback():
    c = _cands(("gwas_catalog", "UKBSTUDY", 1e6, "yes", "no", "unknown"),
               ("gwas_catalog", "OTHER", 1e4, "unknown", "unknown", "unknown"))
    sel = select_outcome(["ukbppp", "decode"], {1: [0]}, c)
    assert (sel.source, sel.candidate, sel.overlap, sel.rule) == ("decode", 0, "no", "first_source_tier1_no")
    sel = select_outcome(["ukbppp", "interval"], {1: [0, 1]}, c)
    assert (sel.source, sel.candidate, sel.tier, sel.overlap, sel.rule) == ("ukbppp", 1, 1, "unknown",
                                                                             "first_source_fallback")


def test_decode_instrument_skips_icelandic_candidate_even_at_better_tier():
    c = _cands(("gwas_catalog", "ICE", 1e6, "no", "yes", "no"), ("finngen", "FG", 1e4, "no", "no", "no"))
    sel = select_outcome(["decode"], {1: [0], 2: [1]}, c)
    assert (sel.source, sel.candidate, sel.tier, sel.overlap) == ("decode", 1, 2, "no")
    none = select_outcome(["decode"], {1: [0]}, c)
    assert none.candidate is None and none.source == "decode"


def test_decode_instrument_never_pairs_with_icelandic_outcome_randomized():
    sources_all = ["ukbppp", "decode", "interval"]
    flags = ["no", "unknown", "yes"]
    for _ in range(3000):
        n = random.randint(1, 6)
        c = _cands(*[(random.choice(["gwas_catalog", "finngen", "opengwas"]), f"A{i}", random.choice([1e4, 5e4, 1e5]),
                      random.choice(flags), random.choice(flags), random.choice(flags)) for i in range(n)])
        pools: dict[int, list[int]] = {1: [], 2: [], 3: []}
        for i in range(n):
            pools[random.choice([1, 2, 3])].append(i)
        srcs = [s for s in sources_all if random.random() < 0.6] or ["decode"]
        sel = select_outcome(srcs, pools, c)
        if sel.source == "decode" and sel.candidate is not None:
            assert c.at[sel.candidate, "includes_iceland"] != "yes"
        if sel.candidate is not None:
            assert sel.overlap == c.at[sel.candidate, {"ukbppp": "includes_ukb", "decode": "includes_iceland",
                                                         "interval": "includes_interval"}[sel.source]]
