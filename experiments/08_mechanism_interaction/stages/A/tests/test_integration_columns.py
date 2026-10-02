from pathlib import Path

import pandas as pd
import pytest

from stage_a.build import run_stage_a
from stage_a.mechanism import S8_BLOCKING, label_row, label_row_s8
from stage_a.outcome_gwas import MixedTraitCoding, outcome_trait_coding, trait_coding
from stage_a.schemas import HYPOTHESIS_COLUMNS, TRAIT_CODING_COLUMNS
from tests.synth import World, gcat_study, moa, world_inputs


INTERFACES = Path(__file__).resolve().parents[2] / "INTERFACES.md"


def table_columns(section: str) -> list[str]:
    cols = []
    for line in section.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if line.startswith("| ") and cells[0] not in ("column",) and not set(cells[0]) <= {"-"}:
            cols += [c.strip() for c in cells[0].split(",")]
    return cols


def test_schemas_match_interfaces_md():
    a = INTERFACES.read_text().split("## Stage A", 1)[1].split("## Stage B", 1)[0]
    hyp_part, trait_part = a.split("`outcome_trait_coding.tsv`, one row", 1)
    assert table_columns(hyp_part) == HYPOTHESIS_COLUMNS
    assert table_columns(trait_part) == TRAIT_CODING_COLUMNS


def _one(h: pd.DataFrame, gene: str, ind: str, cls: str) -> pd.Series:
    rows = h[(h["gene_symbol"] == gene) & (h["indication_id"] == ind) & (h["mechanism_class"] == cls)
             & (h["variant"] == "primary")]
    assert len(rows) == 1, rows
    return rows.iloc[0]


@pytest.fixture(scope="module")
def base():
    return run_stage_a(world_inputs(World()))


# ---- S8 / S9 row labels -----------------------------------------------------------------------

@pytest.mark.parametrize("action,text", [("INHIBITOR", "x inhibitor"), ("ANTAGONIST", "x antagonist"),
                                         ("NEGATIVE ALLOSTERIC MODULATOR", "x"), ("BINDING AGENT", "neutralising x")])
def test_s8_moves_antibody_antagonists_of_non_secreted_targets_only(action, text):
    moved = label_row_s8("Antibody", action, text, blood_secreted=False)
    assert (moved.mechanism_class, moved.direction, moved.reason) == ("blocking", "decrease", S8_BLOCKING)
    assert label_row_s8("Antibody", action, text, blood_secreted=True) == label_row("Antibody", action, text, True)
    for molecule in ("Protein", "Enzyme"):
        assert label_row_s8(molecule, action, text, False).mechanism_class == "other"


@pytest.mark.parametrize("molecule,action,text", [("Antibody", "AGONIST", "x agonist"),
                                                  ("Antibody", "BINDING AGENT", "x binder"),
                                                  ("Antibody drug conjugate", "INHIBITOR", "x inhibitor")])
def test_s8_leaves_other_antibody_rows_alone(molecule, action, text):
    assert label_row_s8(molecule, action, text, False) == label_row(molecule, action, text, False)


# ---- S8 / S9 re-formed hypotheses ----------------------------------------------------------------

def _s(h: pd.DataFrame, flag: str) -> pd.DataFrame:
    return h[h[flag] & h["heldout"]]


def test_s8_moves_the_antibody_into_a_new_blocking_hypothesis(base):
    h = base.hypotheses
    s8 = _s(h, "in_s8")
    keys = set(zip(s8["gene_symbol"], s8["indication_id"], s8["mechanism_class"], s8["direction"]))
    assert ("GB", "MONDO_1001", "other", "decrease") not in keys   # AB11 left other
    gb = s8[(s8["gene_symbol"] == "GB") & (s8["indication_id"] == "MONDO_1001") & (s8["mechanism_class"] == "blocking")]
    assert len(gb) == 1 and gb.iloc[0]["drug_program_ids"] == "AB11;SM2" and gb.iloc[0]["variant"] == "s8"
    # the primary rows it replaces are not in S8; unchanged S1 rows are shared, not copied
    assert not bool(_one(h, "GB", "MONDO_1001", "other")["in_s8"])
    assert not bool(_one(h, "GB", "MONDO_1001", "blocking")["in_s8"])   # SM2 alone is not the S8 hypothesis
    ga = _one(h, "GA", "MONDO_1001", "aligned")
    assert (ga["variant"], ga["in_s8"], ga["in_s9"]) == ("primary", True, True)


def test_s8_merge_gives_one_hypothesis_with_the_union_of_programs():
    """An S8-moved antibody program and a small-molecule blocker, same gene / indication /
    direction: one S8 hypothesis, programs the union, never two."""
    w = World()
    w.molecules.append({"id": "SM12", "drugType": "Small molecule", "parentId": None})
    w.mechanisms.append(moa("SM12", "INHIBITOR", "GB inhibitor", "GB"))
    w.clinical.append(("SM12", "MONDO_1001", "PHASE_2"))
    h = run_stage_a(world_inputs(w)).hypotheses
    prim = _s(h, "in_s1")
    prim_gb = prim[(prim["gene_symbol"] == "GB") & (prim["indication_id"] == "MONDO_1001")]
    assert set(zip(prim_gb["mechanism_class"], prim_gb["drug_program_ids"])) == {("blocking", "SM12;SM2"), ("other", "AB11")}
    s8 = _s(h, "in_s8")
    gb = s8[(s8["gene_symbol"] == "GB") & (s8["indication_id"] == "MONDO_1001") & (s8["direction"] == "decrease")]
    assert len(gb) == 1
    assert (gb.iloc[0]["mechanism_class"], gb.iloc[0]["drug_program_ids"]) == ("blocking", "AB11;SM12;SM2")
    assert h["hypothesis_id"].is_unique
    key_counts = s8.groupby(["gene_symbol", "indication_id", "direction", "mechanism_class"]).size()
    assert (key_counts == 1).all()


def test_s9_reforms_with_uniprot_localization(base):
    h = base.hypotheses
    s9 = _s(h, "in_s9")
    # GA HPA and UniProt secreted: its aligned row is shared. GB not UniProt-secreted: AB11 stays other.
    assert bool(_one(h, "GA", "MONDO_1001", "aligned")["in_s9"])
    assert bool(_one(h, "GB", "MONDO_1001", "other")["in_s9"])
    h2 = run_stage_a(world_inputs(World(uniprot_secreted=("GA", "GB")))).hypotheses
    s9b = _s(h2, "in_s9")
    gb = s9b[(s9b["gene_symbol"] == "GB") & (s9b["indication_id"] == "MONDO_1001")]
    assert set(zip(gb["mechanism_class"], gb["drug_program_ids"], gb["variant"])) == {
        ("aligned", "AB11", "s9"), ("blocking", "SM2", "primary")}
    assert len(s9) == len(set(zip(s9["gene_symbol"], s9["indication_id"], s9["direction"], s9["mechanism_class"])))


def test_conflicting_rows_rule_applies_within_s8():
    """AB11 gains an agonist row. Under S8 its inhibitor row is blocking/decrease and its agonist
    row other/increase: two informative rows disagree, so AB11-GB is other/ambiguous in S8."""
    w = World()
    w.mechanisms.append(moa("AB11", "AGONIST", "GB agonist", "GB"))
    h = run_stage_a(world_inputs(w)).hypotheses
    s8 = _s(h, "in_s8")
    ab = s8[s8["drug_program_ids"].str.split(";").map(lambda x: "AB11" in x)]
    assert set(zip(ab["mechanism_class"], ab["direction"])) == {("other", "ambiguous")}


def test_variant_only_rows_are_outside_every_primary_set(base):
    v = base.hypotheses[base.hypotheses["variant"] != "primary"]
    assert len(v) and not (v["in_s1"] | v["in_s4"] | v["in_s10"] | v["in_s12"]).any()
    assert v["strict_instrument"].all()
    for x in ("s8", "s9", "s21"):
        assert ((v["variant"] == x) == v[f"in_{x}"]).all()
    assert ((base.hypotheses["variant"] == "s21") == base.hypotheses["conflicted_row_restored"]).all()


def test_funnel_counts_variant_rows(base):
    f = base.funnel.set_index("step")
    h = base.hypotheses
    for v in ("s8", "s9", "s21"):
        assert f.loc[f"{v.upper()}:rows_present_only_under_{v}", "remaining"] == int((h["variant"] == v).sum())
        assert f.loc[f"{v.upper()}:rows_shared_with_primary", "remaining"] == int(
            (h[f"in_{v}"] & (h["variant"] == "primary")).sum())


# ---- per-source assays (descriptive table 12) ---------------------------------------------------

def test_source_assay_columns(base):
    h = base.hypotheses
    gb = _one(h, "GB", "MONDO_1001", "blocking")
    assert (gb["instrument_source"], gb["ukbppp_assay_ids"], gb["decode_assay_ids"]) == ("decode", "", "1000-1")
    ga = _one(h, "GA", "MONDO_1001", "aligned")
    assert (ga["ukbppp_assay_ids"], ga["decode_assay_ids"]) == ("OID001", "")  # GA's SomaScan assay is multi-gene


def test_decode_assays_blank_where_the_outcome_is_icelandic():
    w = World()
    w.decode_st01.append({"SeqId": "4000-4", "Gene": "GA", "UniProt": "P00001", "Ensembl": "ENSG00000000001",
                          "Included": "Yes"})
    w.decode_st02.append({"gene": "GA", "SeqId": "4000-4", "cis_trans": "cis"})
    w.gwas_catalog.append(gcat_study("GCSTICE", "MONDO_3000", "deCODE", 200000, 200000, "Iceland"))
    h = run_stage_a(world_inputs(w)).hypotheses
    ice = _one(h, "GA", "MONDO_3000", "aligned")
    assert (ice["instrument_source"], ice["outcome_accession"], ice["decode_assay_ids"]) == ("ukbppp", "GCSTICE", "")
    other = _one(h, "GA", "MONDO_1001", "aligned")
    assert (other["ukbppp_assay_ids"], other["decode_assay_ids"]) == ("OID001", "4000-4")


# ---- outcome trait coding ----------------------------------------------------------------------

def test_trait_coding_rule():
    assert trait_coding("5,000", 50000)[0] is True
    assert trait_coding(12.0, 3.0)[0] is True
    for ncase, ncontrol in ((0, 100), (100, 0), (None, 100), ("NA", 100)):
        assert trait_coding(ncase, ncontrol)[0] is False


def test_trait_coding_table_covers_every_selected_accession_once(base):
    t = base.outcome_trait_coding
    assert t["outcome_accession"].is_unique
    assert set(t["outcome_accession"]) == set(base.hypotheses["outcome_accession"])
    assert t["risk_coded"].all()


def test_trait_coding_refuses_one_accession_coded_both_ways():
    cands = pd.DataFrame([{"source": "gwas_catalog", "accession": "X", "trait": "t", "ncase": 10, "ncontrol": 10},
                          {"source": "opengwas", "accession": "X", "trait": "t", "ncase": 0, "ncontrol": 10}])
    with pytest.raises(MixedTraitCoding):
        outcome_trait_coding({("gwas_catalog", "X"): 0, ("opengwas", "X"): 1}, cands)


def test_s21_restored_hypothesis_with_a_primary_key_but_other_programs_is_kept():
    """GB / MONDO_2000: SM13 is a plain small-molecule inhibitor (primary blocking/decrease);
    SM3 has conflicting rows (primary other/ambiguous). Restoring SM3's rows puts its inhibitor
    row on the same key as SM13's hypothesis, with programs SM13 + SM3. That S21 hypothesis must
    be written, not dropped because the key already exists."""
    w = World()
    w.molecules.append({"id": "SM13", "drugType": "Small molecule", "parentId": None})
    w.mechanisms.append(moa("SM13", "INHIBITOR", "GB inhibitor", "GB"))
    w.clinical.append(("SM13", "MONDO_2000", "PHASE_2"))
    h = run_stage_a(world_inputs(w)).hypotheses
    key = (h["gene_symbol"] == "GB") & (h["indication_id"] == "MONDO_2000") & (h["mechanism_class"] == "blocking") \
        & (h["direction"] == "decrease")
    rows = h[key]
    assert set(zip(rows["variant"], rows["drug_program_ids"])) == {("primary", "SM13"), ("s21", "SM13;SM3")}
    prim, s21 = rows[rows["variant"] == "primary"].iloc[0], rows[rows["variant"] == "s21"].iloc[0]
    assert (prim["in_s1"], prim["in_s21"]) == (True, False)
    assert (s21["in_s1"], s21["in_s21"], s21["conflicted_row_restored"]) == (False, True, True)
    assert s21["hypothesis_id"] != prim["hypothesis_id"]
    s = h[h["in_s21"] & h["heldout"]]
    assert not s.duplicated(["gene_symbol", "indication_id", "direction", "mechanism_class"]).any()
    # the other restored row (SM3 disrupting agent, other/decrease) is also an S21 hypothesis
    assert ((s["gene_symbol"] == "GB") & (s["indication_id"] == "MONDO_2000") & (s["mechanism_class"] == "other")
            & (s["drug_program_ids"] == "SM3")).sum() == 1
