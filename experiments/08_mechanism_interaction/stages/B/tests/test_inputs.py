import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from test_pipeline import fp

from stage_b.collect import collect_tasks
from stage_b.ld import aligned_ld, parse_genotypes, proxies
from stage_b.parsers import (decode_to_canonical, filter_decode, filter_decode_annotation, filter_decode_excluded,
                             filter_finngen, filter_gwas_catalog, filter_ukbppp, normalize_chrom,
                             parse_ukbppp_rsid_map, split_lines, ukbppp_to_canonical)
from stage_b.schemas import (HYPOTHESIS_INPUT_COLUMNS, AmbiguousInstrumentError, HypothesisInput,
                             InputContractError, LDReferenceError, Sentinel, SourceFile)
from stage_b.schemas import InputContractError as _InputContractError
from stage_b.units import load_hypotheses as _load_hypotheses
from stage_b.sentinels import (decode_sentinels, interval_opengwas_id, interval_sentinels, select_assay,
                               ukbppp_sentinels)
from stage_b.units import build_units, load_trait_coding, unit_key


def rows(text: str, sep: str | None = "\t"):
    it = split_lines(text.strip().splitlines(), sep)
    return next(it), it


def test_normalize_chrom():
    assert [normalize_chrom(c) for c in ("chr1", "1", "23", "chrX", "x")] == ["1", "1", "X", "X", "X"]


def test_ukbppp_window_filter_rsid_join_and_log10p():
    header, it = rows("""
CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA
1 1000 1:900:A:G:imp:v1 A G 0.2 1 30000 ADD 0.1 0.01 100 20 NA
1 2000001 1:1990000:C:T:imp:v1 C T 0.2 1 30000 ADD 0.1 0.01 100 3 NA
2 1000 2:900:A:G:imp:v1 A G 0.2 1 30000 ADD 0.1 0.01 100 3 NA
1 900000 1:800000:C:T:imp:v1 C T 0.4 1 30000 ADD -0.1 0.01 100 2 NA
""", None)
    window = filter_ukbppp(it, header, "1", 500_000, 500_000)
    assert [w["ID"] for w in window] == ["1:900:A:G:imp:v1", "1:800000:C:T:imp:v1"]
    mh, mit = rows("ID\tREF\tALT\trsid\tPOS19\tPOS38\n1:900:A:G:imp:v1\tA\tG\trs11\t900\t1000")
    df = ukbppp_to_canonical(window, parse_ukbppp_rsid_map(mit, mh, {w["ID"] for w in window}))
    assert list(df["rsid"]) == ["rs11", ""]
    assert (df.loc[0, "ea"], df.loc[0, "oa"], df.loc[0, "eaf"]) == ("G", "A", 0.2)
    assert df.loc[0, "p"] == pytest.approx(1e-20)


def test_decode_uses_annotation_frequency_corrected_alleles_and_drops_excluded():
    header, it = rows("""
Chrom\tPos\tName\trsids\teffectAllele\totherAllele\tBeta\tPval\tmin_log10_pval\tSE\tN\tImpMAF
chr1\t100\tchr1:100:A:G\trs1\tA\tG\t0.2\t1e-5\t5\t0.01\t35000\t0.4
chr1\t200\tchr1:200:C:T\trs2\tC\tC\t0.1\t1e-3\t3\t0.01\t35000\t0.1
chr1\t300\tchr1:300:G:T\trs3\tG\tT\t0.1\t1e-3\t3\t0.01\t35000\t0.1
chr1\t400\tchr1:400:A:C\trs4,rs44\tA\tC\t0.1\t1e-3\t3\t0.01\t35000\t0.1
chr1\t500\tchr1:500:A:C\trs5\tA\tC\t0.1\t1e-3\t3\t0.01\t35000\t0.1
""")
    window = filter_decode(it, header, "1", 300, 1000)
    names = {w["Name"] for w in window}
    ah, ait = rows("""
Chrom\tPos\tName\trsids\teffectAllele\totherAllele\teffectAlleleFreq
chr1\t100\tchr1:100:A:G\trs1\tA\tG\t0.6
chr1\t200\tchr1:200:C:T\trs2\tC\t!\t0.05
chr1\t400\tchr1:400:A:C\trs4,rs44\tA\tC\t0.3
chr1\t500\tchr1:500:A:C\trs5\tC\tA\t0.3
""")
    eh, eit = rows("Chrom\tPos\tName\trsids\teffectAllele\totherAllele\nchr1\t300\tchr1:300:G:T\trs3\tG\tT")
    df, counts = decode_to_canonical(window, filter_decode_annotation(ait, ah, names), filter_decode_excluded(eit, eh, names))
    assert list(df["rsid"]) == ["rs1", "rs2", "rs4", "rs44"]
    assert df.set_index("rsid").loc["rs1", "eaf"] == 0.6          # annotation, not ImpMAF (0.4)
    assert df.set_index("rsid").loc["rs2", "oa"] == "!"           # multiallelic correction
    assert counts == {"excluded": 1, "not_in_annotation": 0, "annotation_allele_mismatch": 1}


def test_finngen_effect_allele_is_alt_and_rsids_split():
    header, it = rows("""
#chrom\tpos\tref\talt\trsids\tnearest_genes\tpval\tmlogp\tbeta\tsebeta\taf_alt\taf_alt_cases\taf_alt_controls
1\t100\tA\tG\trs1,rs9\tX\t0.01\t2\t0.2\t0.05\t0.3\t0.3\t0.3
1\t9000000\tA\tG\trs2\tX\t0.01\t2\t0.2\t0.05\t0.3\t0.3\t0.3
""")
    df = filter_finngen(it, header, "1", 100, 1000)
    assert list(df["rsid"]) == ["rs1", "rs9"]
    assert (df.loc[0, "ea"], df.loc[0, "oa"], df.loc[0, "eaf"]) == ("G", "A", 0.3)


def test_gwas_catalog_beta_from_odds_ratio_and_contract_errors():
    header, it = rows("""
hm_variant_id\thm_rsid\thm_chrom\thm_pos\thm_other_allele\thm_effect_allele\thm_beta\thm_odds_ratio\thm_effect_allele_frequency\tp_value\tstandard_error
x\trs1\t1\t100\tA\tG\tNA\t2.0\t0.2\t1e-4\t0.1
x\trs2\t1\t200\tA\tG\t0.5\tNA\t0.2\t1e-4\t0.1
""")
    df = filter_gwas_catalog(it, header, "1", 100, 1000)
    assert df["beta"].tolist() == pytest.approx([math.log(2.0), 0.5])
    with pytest.raises(InputContractError):
        filter_gwas_catalog(iter([]), ["hm_chrom", "hm_pos", "hm_effect_allele", "hm_other_allele", "standard_error",
                                       "hm_rsid"], "1", 1, 1)


def test_sentinel_rules_per_source():
    st9 = pd.DataFrame({"UKBPPP ProteinID": ["A:P1:OID1:v1"] * 3 + ["B:P2:OID2:v1"],
                        "rsID": ["rs1", "rs2", "rs3", "rs4"], "CHROM": [1, 1, 1, 2], "GENPOS (hg38)": [10, 20, 30, 40],
                        "log10(p) (discovery)": [50.0, 80.0, 200.0, 12.0], "cis/trans": ["cis", "cis", "trans", "cis"]})
    u = {s.assay_id: s for s in ukbppp_sentinels(st9)}
    assert u["OID1"].rsid == "rs2" and u["OID2"].rsid == "rs4" and u["OID1"].build == "GRCh38"
    st02 = pd.DataFrame({"SeqId": ["1_1", "1_1", "1_1"], "variant": ["rs5", "rs6", "rs7"], "chr\n(var.)": ["chr1"] * 3,
                         "pos\n(var.)": [1, 2, 3], "cis/\ntrans": ["cis", "cis", "trans"],
                         "Rank\n(cond.\nsign.)": [2, 1, 1], "-Log10(P)\n(adj.)": [300.0, 100.0, 500.0]})
    assert decode_sentinels(st02)[0].rsid == "rs6"
    st4 = pd.DataFrame({"SOMAmer ID": ["G.1.2.3", "G.1.2.3"], "Target fullname": ["Prot G"] * 2,
                        "Sentinel variant*": ["rs8", "rs9"], "Chr": [1, 1], "Pos": [5, 6],
                        "cis/ trans": ["cis", "cis"], "meta_p": [1e-10, 1e-40]})
    s = interval_sentinels(st4)[0]
    assert (s.rsid, s.build) == ("rs9", "GRCh37")
    sents = {("decode", "a"): Sentinel(source="decode", assay_id="a", rsid="rs1", chrom="1", pos=1, build="GRCh38", neg_log10_p=10),
             ("decode", "b"): Sentinel(source="decode", assay_id="b", rsid="rs2", chrom="1", pos=1, build="GRCh38", neg_log10_p=90)}
    assert select_assay(["a", "b"], sents, "decode").assay_id == "b"
    with pytest.raises(AmbiguousInstrumentError):
        select_assay(["c"], sents, "decode")


def test_interval_opengwas_id_requires_a_unique_trait_match():
    info = [{"id": "prot-a-1", "trait": "Prot G"}, {"id": "prot-a-2", "trait": "Prot H"},
            {"id": "prot-a-3", "trait": "Prot H"}, {"id": "ebi-a-1", "trait": "Prot G"}]
    assert interval_opengwas_id("Prot G", info) == "prot-a-1"
    with pytest.raises(AmbiguousInstrumentError):
        interval_opengwas_id("Prot H", info)


def hyp(hid: str, assay: str, acc: str, gene: str = "ENSG1") -> HypothesisInput:
    return HypothesisInput(hypothesis_id=hid, gene_symbol="G", gene_ensembl=gene, direction="decrease",
                           instrument_source="decode", instrument_assay_id=assay, platform="SomaScan",
                           outcome_accession=acc, outcome_source="finngen", outcome_n_case=100, outcome_n_control=900)


def test_build_units_groups_by_instrument_and_reports_unresolved():
    sents = {("decode", "a"): Sentinel(source="decode", assay_id="a", rsid="rs1", chrom="1", pos=1, build="GRCh38", neg_log10_p=10)}
    hyps = [hyp("h1", "a", "F1"), hyp("h2", "a", "F1"), hyp("h3", "a", "F2"), hyp("h4", "zz", "F1")]
    listed = {"a": SourceFile(name="1_1_G_G.txt.gz", size=10, etag="e" * 32)}
    units, hyp_unit, unresolved, source_units = build_units(hyps, sents, lambda s: "1_1_G_G.txt.gz",
                                                            {"F1": True, "F2": False}, listed, {})
    assert (units[0].pqtl_locator, units[0].pqtl_listing, units[0].smp_listing) == ("1_1_G_G.txt.gz", listed["a"], None)
    assert len(units) == 1 and [o.accession for o in units[0].outcomes] == ["F1", "F2"]      # one gene, one assay: one unit
    assert (units[0].unit_key, units[0].gene_ensembl) == ("decode__a__ENSG1", "ENSG1")
    assert units[0].outcomes[1].risk_coded is False
    assert hyp_unit == {h: "decode__a__ENSG1" for h in ("h1", "h2", "h3")} and set(unresolved) == {"h4"}
    assert source_units == {"h1": {"decode": "decode__a__ENSG1"}, "h2": {"decode": "decode__a__ENSG1"},
                            "h3": {"decode": "decode__a__ENSG1"}, "h4": {"decode": ""}}
    with pytest.raises(InputContractError):
        build_units(hyps, sents, lambda s: "u", {"F1": True}, {}, {})


def test_an_assay_serving_two_genes_gives_one_unit_per_gene_sharing_the_sentinel_and_the_file():
    sents = {("decode", "a"): Sentinel(source="decode", assay_id="a", rsid="rs1", chrom="1", pos=1, build="GRCh38", neg_log10_p=10)}
    listed = {"a": SourceFile(name="1_1_G_G.txt.gz", size=10, etag="e" * 32)}
    smp = {"a": SourceFile(name="1_1_G_G.txt.gz", size=9, etag="f" * 32)}
    hyps = [hyp("h1", "a", "F1"), hyp("h2", "a", "F1", gene="ENSG2"), hyp("h3", "a", "F2", gene="ENSG2"), hyp("h4", "a", "F3")]
    units, hyp_unit, unresolved, source_units = build_units(hyps, sents, lambda s: listed[s.assay_id].name,
                                                            {"F1": True, "F2": True, "F3": True}, listed, smp)
    by_key = {u.unit_key: u for u in units}
    assert list(by_key) == ["decode__a__ENSG1", "decode__a__ENSG2"]
    one, two = by_key["decode__a__ENSG1"], by_key["decode__a__ENSG2"]
    assert (one.gene_ensembl, [o.accession for o in one.outcomes]) == ("ENSG1", ["F1", "F3"])
    assert (two.gene_ensembl, [o.accession for o in two.outcomes]) == ("ENSG2", ["F1", "F2"])
    for field in ("source", "assay_id", "platform", "sentinel", "pqtl_locator", "pqtl_listing", "smp_listing"):
        assert getattr(one, field) == getattr(two, field), field
    assert (one.pqtl_listing, one.smp_listing) == (listed["a"], smp["a"])
    # each hypothesis goes to the unit of its own gene
    assert hyp_unit == {"h1": "decode__a__ENSG1", "h2": "decode__a__ENSG2", "h3": "decode__a__ENSG2", "h4": "decode__a__ENSG1"}
    assert source_units == {h: {"decode": k} for h, k in hyp_unit.items()} and unresolved == {}
    # the shared file is one collect task (and one for its SMP-normalized release), whatever the number of genes
    assert [(t.source, t.key, t.name, t.size) for t in collect_tasks(units)] == [("decode", "a", "1_1_G_G.txt.gz", 10),
                                                                              ("decode_smp", "a", "1_1_G_G.txt.gz", 9)]
    assert collect_tasks(units) == collect_tasks([one]) == collect_tasks([two])
    assert fp(one) != fp(two)                                        # two checkpoints, never resumed as each other


def test_a_seqid_on_two_genes_lists_gives_each_gene_its_own_unit_in_the_cross_source_pairing():
    sents = {("decode", "a"): Sentinel(source="decode", assay_id="a", rsid="rs1", chrom="4", pos=1, build="GRCh38", neg_log10_p=90),
             ("decode", "b"): Sentinel(source="decode", assay_id="b", rsid="rs3", chrom="4", pos=9, build="GRCh38", neg_log10_p=10),
             ("ukbppp", "OID1"): Sentinel(source="ukbppp", assay_id="OID1", rsid="rs2", chrom="4", pos=5, build="GRCh38",
                                          neg_log10_p=40)}
    selected = hyp("c1", "a", "F1", gene="ENSG_C").model_copy(update={"decode_assay_ids": "a"})
    paired = hyp("b1", "OID1", "F2", gene="ENSG_B").model_copy(update={
        "instrument_source": "ukbppp", "platform": "Olink", "ukbppp_assay_ids": "OID1", "decode_assay_ids": "a;b"})
    units, hyp_unit, unresolved, source_units = build_units([selected, paired], sents, lambda s: f"loc_{s.assay_id}",
                                                            {"F1": True, "F2": True}, {}, {})
    by_key = {u.unit_key: u for u in units}
    assert set(by_key) == {"decode__a__ENSG_C", "decode__a__ENSG_B", "ukbppp__OID1__ENSG_B"}
    assert [o.accession for o in by_key["decode__a__ENSG_C"].outcomes] == ["F1"]
    assert [o.accession for o in by_key["decode__a__ENSG_B"].outcomes] == ["F2"]
    assert by_key["decode__a__ENSG_B"].sentinel == by_key["decode__a__ENSG_C"].sentinel == sents[("decode", "a")]
    assert hyp_unit == {"c1": "decode__a__ENSG_C", "b1": "ukbppp__OID1__ENSG_B"} and unresolved == {}
    assert source_units == {"c1": {"decode": "decode__a__ENSG_C"},
                            "b1": {"ukbppp": "ukbppp__OID1__ENSG_B", "decode": "decode__a__ENSG_B"}}


def test_a_unit_key_is_a_file_name_and_never_stands_for_two_instruments():
    assert unit_key("interval", "LAMA1.LAMB1.LAMC1.2728.62.2", "ENSG00000101680") == (
        "interval__LAMA1.LAMB1.LAMC1.2728.62.2__ENSG00000101680")
    assert unit_key("decode", "15525_294", "ENSG00000196616") == "decode__15525_294__ENSG00000196616"
    assert unit_key("interval", "a/b:c d", "ENS G/1") == "interval__a_b_c_d__ENS_G_1"
    assert re.fullmatch(r"[A-Za-z0-9_.-]+", unit_key("interval", "../x\\y?*", "é"))
    sents = {("decode", a): Sentinel(source="decode", assay_id=a, rsid="rs1", chrom="1", pos=1, build="GRCh38", neg_log10_p=10)
             for a in ("a:b", "a_b")}
    with pytest.raises(InputContractError, match="names two instruments"):       # two assay ids that give one file name
        build_units([hyp("h1", "a:b", "F1"), hyp("h2", "a_b", "F1")], sents, lambda s: "u", {"F1": True}, {}, {})


def test_build_units_pairs_the_outcome_with_the_other_sources_instrument():
    sents = {("decode", "a"): Sentinel(source="decode", assay_id="a", rsid="rs1", chrom="1", pos=1, build="GRCh38", neg_log10_p=10),
             ("ukbppp", "OID1"): Sentinel(source="ukbppp", assay_id="OID1", rsid="rs2", chrom="1", pos=1, build="GRCh38",
                                          neg_log10_p=40)}
    both = hyp("h1", "a", "F1").model_copy(update={"ukbppp_assay_ids": "OID1", "decode_assay_ids": "a"})
    missing_ukb = hyp("h2", "a", "F2").model_copy(update={"ukbppp_assay_ids": "OID9"})
    decode_only = hyp("h3", "a", "F3")
    units, hyp_unit, unresolved, source_units = build_units(
        [both, missing_ukb, decode_only], sents, lambda s: f"loc:{s.assay_id}", {"F1": True, "F2": True, "F3": True}, {}, {})
    by_key = {u.unit_key: u for u in units}
    decode, ukb = "decode__a__ENSG1", "ukbppp__OID1__ENSG1"
    assert set(by_key) == {decode, ukb}
    assert (by_key[ukb].platform, by_key[decode].platform) == ("Olink", "SomaScan")
    assert [o.accession for o in by_key[ukb].outcomes] == ["F1"]
    assert [o.accession for o in by_key[decode].outcomes] == ["F1", "F2", "F3"]
    assert hyp_unit == {"h1": decode, "h2": decode, "h3": decode} and unresolved == {}
    assert source_units == {"h1": {"ukbppp": ukb, "decode": decode}, "h2": {"ukbppp": "", "decode": decode}, "h3": {"decode": decode}}


# ---- LD ---------------------------------------------------------------------------------------

def vcf_lines(genos: dict[tuple[str, str, str], list[str]]) -> list[str]:
    return [f"1\t{i}\t{rs}\t{ref}\t{alt}\t" + "\t".join(g) for i, ((rs, ref, alt), g) in enumerate(genos.items())]


def test_proxies_and_signed_ld():
    lead = ["0|0", "0|1", "1|1", "0|1", "0|0", "1|0"]
    same = list(lead)
    flipped = [{"0|0": "1|1", "1|1": "0|0"}.get(g, g) for g in lead]     # perfectly anti-correlated
    other = ["0|1", "0|0", "0|1", "1|1", "0|1", "0|0"]
    meta, dos = parse_genotypes(vcf_lines({("rs1", "A", "G"): lead, ("rs2", "C", "T"): same,
                                           ("rs3", "A", "C"): flipped, ("rs4", "G", "T"): other,
                                           ("rs5", "A", "G,T"): lead}))
    assert list(meta["rsid"]) == ["rs1", "rs2", "rs3", "rs4"]           # multiallelic record dropped
    assert proxies(meta, dos, "rs1") == {"rs1", "rs2", "rs3"}
    assert proxies(meta, dos, "rs404") == set()
    h = pd.DataFrame({"rsid": ["rs1", "rs3", "rs9"], "ea": ["G", "A", "A"], "oa": ["A", "C", "G"]})
    sub, ld = aligned_ld(h, meta, dos)
    assert list(sub["rsid"]) == ["rs1", "rs3"]
    assert ld[0, 1] == pytest.approx(1.0)      # rs3 effect allele A is REF, so its sign flips back


def test_ld_panel_without_rsids_cannot_be_matched():
    with pytest.raises(LDReferenceError):
        parse_genotypes(vcf_lines({(".", "A", "G"): ["0|1", "1|1"], (".", "C", "T"): ["0|1", "0|0"]}))
    assert np.isfinite(parse_genotypes(vcf_lines({("rs1", "A", "G"): ["0|1", ".|."]}))[1]).all()


# ---- the seam with stage A (INTERFACES.md is the contract both stages are written against) -----

INTERFACES = Path(__file__).resolve().parents[2] / "INTERFACES.md"


def stage_a_columns() -> list[str]:
    section = INTERFACES.read_text().split("## Stage A", 1)[1].split("`outcome_trait_coding.tsv`, one row", 1)[0]
    cols = []
    for line in section.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if line.startswith("| ") and cells[0] not in ("column",) and not set(cells[0]) <= {"-"}:
            cols += [c.strip() for c in cells[0].split(",")]
    return cols


def test_every_column_stage_b_reads_is_one_stage_a_writes():
    assert set(HYPOTHESIS_INPUT_COLUMNS) <= set(stage_a_columns())
    assert {"ukbppp_assay_ids", "decode_assay_ids", "outcome_n_case"} <= set(HYPOTHESIS_INPUT_COLUMNS)


def test_trait_coding_file_in_stage_a_layout_loads(tmp_path):
    p = tmp_path / "outcome_trait_coding.tsv"
    pd.DataFrame({"outcome_accession": ["GCST1", "FINNGEN_R12_X"], "risk_coded": [True, False],
                  "outcome_source": ["gwas_catalog", "finngen"], "trait": ["t", "u"], "n_case": [10, 20],
                  "n_control": [100, 200], "coding_basis": ["a", "b"]}).to_csv(p, sep="\t", index=False)
    assert load_trait_coding(p) == {"GCST1": True, "FINNGEN_R12_X": False}


def test_load_hypotheses_refuses_a_repeated_hypothesis_id(tmp_path):
    row = {c: "x" for c in HYPOTHESIS_INPUT_COLUMNS}
    path = tmp_path / "hypotheses.csv"
    pd.DataFrame([row, row]).to_csv(path, index=False)
    with pytest.raises(_InputContractError, match="repeats hypothesis ids"):
        _load_hypotheses(path)
