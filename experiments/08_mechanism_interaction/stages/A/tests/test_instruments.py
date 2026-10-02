import re
from pathlib import Path

import pandas as pd

from stage_a import loaders
from stage_a.instruments import EPIGRAPHDB_NO_ASSAY, UNIPROT_ACC, InstrumentLists, SourceList, decode_lists, \
    decode_strict_diagnostic, interval_lists, split_ids, ukbppp_lists
from stage_a.loaders import load_decode, load_epigraphdb_genes
from tests.synth import World, write_world


def _olink(*rows) -> pd.DataFrame:
    return pd.DataFrame([dict(zip(["Assay", "OlinkID", "UniProt", "UniProt2", "HGNC.symbol", "ensembl_id"], r))
                         for r in rows])


def test_ukbppp_strict_requires_one_gene_and_one_accession():
    olink = _olink(
        ("A", "OID1", "P11111", "P11111", "A", "ENSG1"),
        ("B_C", "OID2", "P22222_P33333", "P22222", "B", "ENSG2"),
        ("B_C", "OID2", "P22222_P33333", "P33333", "C", "ENSG3"),
        ("D", "OID4", "P44444_P55555", "P44444", "D", "ENSG4"),
        ("E", "OID5", "NTproE", "P66666-3", "E", "ENSG5"),
        ("F", "OID6", "P77777-2", "P77777", "F", "ENSG6"),
        ("NOTLISTED", "OID7", "P88888", "P88888", "NOTLISTED", "ENSG7"),
    )
    strict, inclusive = ukbppp_lists(["A", "B_C", "D", "E", "F"], olink)
    assert set(strict.by_symbol) == {"A", "E", "F"}
    assert strict.assays("A", "ENSG1") == {"OID1"}
    assert strict.has("x", "ENSG5")
    assert {"A", "B", "C", "D", "E", "F"} <= set(inclusive.by_symbol)
    assert "NOTLISTED" not in inclusive.by_symbol
    assert inclusive.assays("C", "") == {"OID2"}


def test_ukbppp_same_assay_on_two_panels_stays_one_gene():
    olink = _olink(("A", "OID1", "P11111", "P11111", "A", "ENSG1"), ("A", "OID9", "P11111", "P11111", "A", "ENSG1"))
    strict, _ = ukbppp_lists(["A"], olink)
    assert strict.assays("A", "ENSG1") == {"OID1", "OID9"}


def test_decode_strict_excludes_multi_gene_and_not_included_aptamers():
    st01 = pd.DataFrame([
        {"SeqId": "1-1", "Gene": "A", "UniProt": "P1", "Ensembl": "ENSG1", "Included": "Yes"},
        {"SeqId": "2-2", "Gene": "B, C", "UniProt": "P2, P3", "Ensembl": "ENSG2, ENSG3", "Included": "Yes"},
        {"SeqId": "3-3", "Gene": "D", "UniProt": "P4", "Ensembl": "ENSG4", "Included": "No"},
        {"SeqId": "4-4", "Gene": "E", "UniProt": "P5 P6", "Ensembl": "ENSG5", "Included": "Yes"},
        {"SeqId": "5-5", "Gene": "F", "UniProt": "P7", "Ensembl": "ENSG6", "Included": "Yes"},
    ])
    st02 = pd.DataFrame([
        {"gene": "A", "SeqId": "1-1", "cis_trans": "cis"}, {"gene": "B.C", "SeqId": "2-2", "cis_trans": "cis"},
        {"gene": "D", "SeqId": "3-3", "cis_trans": "cis"}, {"gene": "E", "SeqId": "4-4", "cis_trans": "cis"},
        {"gene": "F", "SeqId": "5-5", "cis_trans": "trans"},
    ])
    strict, inclusive = decode_lists(st01, st02)
    assert set(strict.by_symbol) == {"A"}
    assert strict.assays("A", "ENSG1") == {"1-1"}
    assert not strict.has("F", "ENSG6")
    assert {"A", "B", "C", "D", "E"} <= set(inclusive.by_symbol)
    assert "F" not in inclusive.by_symbol


def test_interval_strict_single_accession_single_symbol():
    st4 = pd.DataFrame([
        {"somamer_id": "S1", "UniProt": "P1", "cis_trans": "cis"},
        {"somamer_id": "S2", "UniProt": "P2 P3", "cis_trans": "cis"},
        {"somamer_id": "S3", "UniProt": "P4", "cis_trans": "cis"},
        {"somamer_id": "S5", "UniProt": "P5", "cis_trans": "trans"},
    ])
    uni2sym = {"P1": {"A"}, "P2": {"B"}, "P3": {"C"}, "P4": {"D1", "D2"}, "P5": {"E"}}
    strict, inclusive = interval_lists(st4, uni2sym, {"Z"})
    assert set(strict.by_symbol) == {"A"}
    assert set(inclusive.by_symbol) == {"A", "B", "C", "D1", "D2", "Z"}
    assert inclusive.assays("Z", "") == {EPIGRAPHDB_NO_ASSAY}


def test_sources_follow_discovery_order():
    lists = InstrumentLists(ukbppp=SourceList(by_symbol={"G": frozenset({"o"})}),
                            decode=SourceList(by_ensembl={"E": frozenset({"s"})}),
                            interval=SourceList(by_symbol={"G": frozenset({"i"})}))
    assert lists.sources_for("G", "E") == ["ukbppp", "decode", "interval"]
    assert lists.sources_for("H", "E") == ["decode"]
    assert lists.sources_for("H", "X") == []


# ---- deCODE: the registered rule (count_v4.py decode_lists), its port and the diagnostic ------------

DECODE_ST01 = [
    # SeqId, Gene, UniProt, Ensembl, Included
    ("1-1", "A", "P1", "ENSG1", "Yes"),        # ST02 names another gene, B, for this SeqId
    ("2-2", "C", "P2", None, "Yes"),           # blank Ensembl ID
    ("3-3", "D", "P3", "ENSG3", "Yes"),        # the two tables agree
    ("4-4", "E", "P4", "ENSG4", "Yes"),        # ST02 names two genes, E and F
    ("5-5", "G, H", "P5, P6", "ENSG5, ENSG6", "Yes"),   # two ST01 genes: not strict
    ("6-6", "K", "P7", "ENSG7", "No"),         # not included in the analysis: no ST01 row
    ("7-7", "L", "P8 P9", "ENSG8", "Yes"),     # two UniProt accessions: not strict
    ("8-8", "M", "P10", "ENSG9 ENSG10", "Yes"),  # two Ensembl IDs: not strict
    ("9-9", "N", "P11", "ENSG11", "Yes"),      # only a trans pQTL
]
DECODE_ST02 = [
    # gene (prot.), SeqId, cis/trans
    ("B", "1-1", "cis"), ("C", "2-2", "cis"), ("D", "3-3", "cis"), (None, "3-3", "cis"), ("E.F", "4-4", "cis"),
    ("G.H", "5-5", "cis"), ("K", "6-6", "cis"), ("L", "7-7", "cis"), ("M", "8-8", "cis"), ("N", "9-9", "trans"),
    ("D", "3-3", "trans"),
]


def decode_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    return (pd.DataFrame(DECODE_ST01, columns=["SeqId", "Gene", "UniProt", "Ensembl", "Included"]),
            pd.DataFrame(DECODE_ST02, columns=["gene", "SeqId", "cis_trans"]))


def count_v4_decode_lists(st01_rows, st02_rows) -> dict:
    """feasibility/v4_round3/count_v4.py `decode_lists`, lines 232-257, transcribed with the
    worksheet rows passed in as tuples (SeqId, Gene, UniProt, Ensembl, Included) and (gene, SeqId,
    cis/trans). The registered rule: the plan defines the lists "as implemented in" that script."""
    seq = {}
    for s, gene, uni, ens, inc in st01_rows:
        if s is None or str(inc).strip() != "Yes":
            continue
        seq[str(s)] = {"genes": set(split_ids(gene)), "unis": set(split_ids(uni)),
                       "ens": {e for e in split_ids(ens) if e.startswith("ENSG")}}
    out = {k: {"sym": set(), "ens": set()} for k in ("strict", "inclusive")}
    for g, s, ct in st02_rows:
        if ct != "cis":
            continue
        s = str(s)
        info = seq.get(s, {"genes": set(), "unis": set(), "ens": set()})
        sym2 = {x for x in re.split(r"[ ,.;|]+", str(g)) if x and x != "NA"}
        out["inclusive"]["sym"] |= sym2 | info["genes"]
        out["inclusive"]["ens"] |= info["ens"]
        if len(info["genes"]) == 1 and len(info["unis"]) == 1 and len(info["ens"]) <= 1:
            out["strict"]["sym"] |= sym2 | info["genes"]
            out["strict"]["ens"] |= info["ens"]
    return out


def test_decode_strict_list_is_the_registered_rule_on_a_hand_worked_fixture():
    strict, inclusive = decode_lists(*decode_frames())
    # 1-1: one ST01 Gene (A), so the SeqId is strict, and it lists A and the ST02 gene B.
    # 2-2: a blank Ensembl ID is "at most one", so C is strict, by symbol only.
    # 4-4: one ST01 Gene (E); the ST02 cell "E.F" adds F.
    # 3-3: its second cis row has an empty gene cell, which count_v4.py line 248 reads as the symbol "None".
    assert strict.by_symbol == {"A": {"1-1"}, "B": {"1-1"}, "C": {"2-2"}, "D": {"3-3"}, "None": {"3-3"}, "E": {"4-4"},
                                "F": {"4-4"}}
    assert strict.by_ensembl == {"ENSG1": {"1-1"}, "ENSG3": {"3-3"}, "ENSG4": {"4-4"}}
    assert strict.has("C", "") and strict.has("B", "") and not strict.has("G", "ENSG5")
    assert not strict.has("K", "ENSG7") and not strict.has("L", "ENSG8") and not strict.has("M", "ENSG9")
    assert not strict.has("N", "ENSG11")
    assert set(inclusive.by_symbol) == {"A", "B", "C", "D", "None", "E", "F", "G", "H", "K", "L", "M"}
    assert set(inclusive.by_ensembl) == {"ENSG1", "ENSG3", "ENSG4", "ENSG5", "ENSG6", "ENSG8", "ENSG9", "ENSG10"}


def test_decode_port_gives_the_lists_of_count_v4_on_the_fixture_and_on_the_synthetic_world():
    w = World()
    world = ([tuple(r[c] for c in ("SeqId", "Gene", "UniProt", "Ensembl", "Included")) for r in w.decode_st01],
             [tuple(r[c] for c in ("gene", "SeqId", "cis_trans")) for r in w.decode_st02])
    for st01_rows, st02_rows in ((DECODE_ST01, DECODE_ST02), world):
        ref = count_v4_decode_lists(st01_rows, st02_rows)
        strict, inclusive = decode_lists(pd.DataFrame(st01_rows, columns=["SeqId", "Gene", "UniProt", "Ensembl", "Included"]),
                                         pd.DataFrame(st02_rows, columns=["gene", "SeqId", "cis_trans"]))
        for name, got in (("strict", strict), ("inclusive", inclusive)):
            assert ref[name]["sym"] == set(got.by_symbol), name
            assert ref[name]["ens"] == set(got.by_ensembl), name
    assert "None" in count_v4_decode_lists(DECODE_ST01, DECODE_ST02)["strict"]["sym"]     # line 248, kept by the port


def test_decode_diagnostic_counts_disagreeing_seqids_and_the_symbols_they_add():
    st01, st02 = decode_frames()
    d = decode_strict_diagnostic(st01, st02)
    # strict SeqIds: 1-1, 2-2, 3-3, 4-4. ST02 symbols differ from the ST01 Gene on 1-1 ({B} vs {A}), 3-3
    # ({D, None} vs {D}: the empty cell) and 4-4 ({E, F} vs {E}). Strict symbols A, B, C, D, None, E, F; B, None
    # and F are the ST01 Gene of no strict SeqId.
    assert d.model_dump() == {"strict_seqids": 4, "strict_seqids_st02_differs_from_st01": 3, "strict_symbols": 7,
                              "strict_symbols_only_through_disagreement": 3}
    strict, _ = decode_lists(st01, st02)
    assert d.strict_symbols == len(strict.by_symbol)
    assert d.strict_seqids == len(set().union(*strict.by_symbol.values()))


def test_decode_diagnostic_does_not_count_a_symbol_that_another_strict_seqid_names_in_st01():
    st01 = pd.DataFrame([("1-1", "A", "P1", "ENSG1", "Yes"), ("2-2", "B", "P2", "ENSG2", "Yes")],
                        columns=["SeqId", "Gene", "UniProt", "Ensembl", "Included"])
    st02 = pd.DataFrame([("B", "1-1", "cis"), ("B", "2-2", "cis")], columns=["gene", "SeqId", "cis_trans"])
    d = decode_strict_diagnostic(st01, st02)
    assert (d.strict_seqids_st02_differs_from_st01, d.strict_symbols_only_through_disagreement) == (1, 0)
    empty = decode_strict_diagnostic(st01, st02.iloc[0:0])
    assert set(empty.model_dump().values()) == {0}


# ---- UKB-PPP and INTERVAL: count_v4.py's handling, kept where it looks like a quirk ------------------

def count_v4_ukbppp_strict_assays(raw, olink: pd.DataFrame) -> set[str]:
    """count_v4.py `ukbppp_lists`, lines 198-209, with the list and the Olink map passed in."""
    olink = olink[olink["Assay"].isin(raw)]
    strict = set()
    for assay, g in olink.groupby("Assay"):
        syms = set(g["HGNC.symbol"].dropna())
        accs = {a.split("-")[0] for u in g["UniProt"].dropna() for a in u.split("_")
                if UNIPROT_ACC.match(a.split("-")[0])}
        if not accs:
            accs = {a.split("-")[0] for a in g["UniProt2"].dropna() if UNIPROT_ACC.match(a.split("-")[0])}
        if "_" not in assay and len(syms) == 1 and len(accs) == 1:
            strict.add(assay)
    return strict


def test_ukbppp_uniprot2_is_not_split_on_underscores_as_in_count_v4():
    olink = _olink(
        ("A", "OID1", "NTproA", "P11111_notanaccession", "A", "ENSG1"),   # one accession only if split: not strict
        ("B", "OID2", "NTproB", "P22222_P33333", "B", "ENSG2"),           # two if split, none unsplit: not strict
        ("C", "OID3", "NTproC", "P44444-2", "C", "ENSG3"),                # an accession with an isoform suffix: strict
        ("D", "OID4", "P55555_notanaccession", "ignored", "D", "ENSG4"),  # `UniProt` itself is split (lines 202-203)
        ("E", "OID5", "NTproE", None, "E", "ENSG5"),                      # no accession in either column
    )
    cis = ["A", "B", "C", "D", "E"]
    strict, inclusive = ukbppp_lists(cis, olink)
    assert set(strict.by_symbol) == {"C", "D"} == count_v4_ukbppp_strict_assays(cis, olink)
    assert set(strict.by_ensembl) == {"ENSG3", "ENSG4"}
    assert {"A", "B", "C", "D", "E"} <= set(inclusive.by_symbol)


def test_ukbppp_strict_assays_are_those_of_count_v4_on_the_synthetic_world():
    w = World()
    olink = pd.DataFrame(w.olink)
    strict, _ = ukbppp_lists(list(w.ukbppp_cis), olink)
    ref = count_v4_ukbppp_strict_assays(list(w.ukbppp_cis), olink)
    assert ref == {"GA", "GD", "GG"} and set(strict.by_symbol) == ref


def test_interval_keeps_an_empty_epigraphdb_token_and_reads_a_falsy_uniprot_cell_as_no_accession():
    st4 = pd.DataFrame([{"somamer_id": "S1", "UniProt": "P1", "cis_trans": "cis"},
                        {"somamer_id": "S2", "UniProt": None, "cis_trans": "cis"},
                        {"somamer_id": "S3", "UniProt": 0, "cis_trans": "cis"},
                        {"somamer_id": "S4", "UniProt": "", "cis_trans": "cis"}], dtype=object)
    strict, inclusive = interval_lists(st4, {"P1": {"A"}, "0": {"ZERO"}, "None": {"NONE"}}, {"Z", ""})
    assert strict.by_symbol == {"A": {"S1"}}
    # count_v4.py line 293 splits "A;;B" into "A", "" and "B"; the empty token stays a list entry
    assert inclusive.by_symbol == {"A": {"S1"}, "Z": {EPIGRAPHDB_NO_ASSAY}, "": {EPIGRAPHDB_NO_ASSAY}}
    assert not inclusive.has("GENE", "ENSG1")


def test_epigraphdb_loader_keeps_every_token_of_the_split(tmp_path):
    path = tmp_path / "epi.csv"
    pd.DataFrame({"gene_symbol": ["A;;B", "C", None], "mr_pval": ["x", "y", "z"]}).to_csv(path, index=False)
    assert load_epigraphdb_genes(path) == {"A", "", "B", "C"}


def test_count_v4_sheets_are_read_over_the_recorded_dimensions_and_cells_stay_as_given(tmp_path):
    path = write_world(World(decode_st02=[{"gene": None, "SeqId": "1000-1", "cis_trans": "cis"},
                                          {"gene": "GF", "SeqId": "2000-2", "cis_trans": "cis"}]), tmp_path / "in")
    st01, st02 = load_decode(tmp_path / "in" / "decode.xlsx")
    assert path.is_file() and st02["gene"].tolist() == [None, "GF"] and st02.dtypes.eq(object).all()
    assert st01.dtypes.eq(object).all() and st01["SeqId"].tolist() == ["1000-1", "2000-2", "3000-3"]
    strict, _ = decode_lists(st01, st02)
    assert set(strict.by_symbol) == {"GB", "None", "GF"}
    src = (Path(loaders.__file__)).read_text()
    assert "def _sheet(path: Path, sheet: str, header_row: int, columns: dict[str, str],\n           reset_dimensions: bool = False)" in src
    for fn, reset in (("load_decode", False), ("load_sun2018_st4", False), ("load_karim_launched", True)):
        body = src.split(f"def {fn}(")[1].split("\ndef ")[0]
        assert ("reset_dimensions=True" in body) is reset, fn
