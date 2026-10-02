import numpy as np
import pandas as pd
import pytest

from stage_d.evidence import evidence_score, evidence_state
from stage_d.join import JoinError, derive, join_stages, load_stage_tables
from stage_d.sets import SetError, form_sets, s12_select


def test_evidence_state_truth_table():
    pp = [0.9, 0.9, 0.9, 0.9, 0.79, 0.9, np.nan, 0.8]
    gd = [1, -1, -1, 1, 1, 0, 1, 1]
    idir = ["decrease", "decrease", "increase", "increase", "decrease", "decrease", "decrease", "ambiguous"]
    assert list(evidence_state(pp, gd, idir)) == ["supportive", "contradictory", "supportive", "contradictory",
                                                  "inconclusive", "inconclusive", "inconclusive", "inconclusive"]
    assert list(evidence_score(pp, gd, idir)) == pytest.approx([0.9, -0.9, 0.9, -0.9, 0.79, 0.0, 0.0, 0.0])
    assert list(evidence_state([0.8, 0.75], [1, 1], ["decrease"] * 2)) == ["supportive", "inconclusive"]


def _raw_joined(root):
    return join_stages(*load_stage_tables(root))


def test_join_rejects_id_mismatch(sealed_stages):
    hyp, ev, out = load_stage_tables(sealed_stages[0])
    with pytest.raises(JoinError, match="A ids missing"):
        join_stages(hyp, ev.iloc[1:], out)


def test_derive_rejects_state_inconsistent_with_pp_h4(sealed_stages):
    raw = _raw_joined(sealed_stages[0])
    i = raw.index[raw["evidence_state"] == "supportive"][0]
    raw.loc[i, "pp_h4"] = 0.5
    with pytest.raises(JoinError, match="evidence_state"):
        derive(raw)


def test_derive_rejects_a_selected_source_state_that_is_not_the_primary_state(sealed_stages):
    raw = _raw_joined(sealed_stages[0])
    assert "evidence_state_ukbppp" in raw.attrs["present_columns"]["B"]
    derive(raw)
    i = raw.index[(raw["instrument_source"] == "ukbppp") & (raw["evidence_state"] == "supportive")][0]
    raw.loc[i, "evidence_state_ukbppp"] = "inconclusive"
    with pytest.raises(JoinError, match="evidence_state_ukbppp disagrees"):
        derive(raw)


def test_selected_source_state_check_skipped_when_column_absent(sealed_stages):
    raw = _raw_joined(sealed_stages[0])
    raw.attrs["present_columns"]["B"] = [c for c in raw.attrs["present_columns"]["B"] if c != "evidence_state_decode"]
    i = raw.index[raw["instrument_source"] == "decode"][0]
    raw.loc[i, "evidence_state_decode"] = None
    derive(raw)


def test_pre_pqtl_publication_is_strictly_before_publication(sealed_stages):
    raw = _raw_joined(sealed_stages[0])
    idx = raw.index[:4]
    raw.loc[idx, "pre_pqtl_publication_date"] = pd.to_datetime(["2021-06-01"] * 4)
    raw.loc[idx, "earliest_phase2_start"] = pd.to_datetime(["2021-05-31", "2021-06-01", "2021-06-02", None])
    df, info = derive(raw)
    assert list(df.loc[idx, "pre_pqtl_publication"].astype(object)) == [True, False, False, pd.NA]
    sets = form_sets(df)
    s6 = set(sets["S6"].frame["hypothesis_id"])
    s1 = sets["S1"].frame
    assert s6 == set(s1.loc[(s1["pre_pqtl_publication"] == True).fillna(False), "hypothesis_id"])  # noqa: E712
    assert sets["S6"].info.notes["date_missing_excluded"] == int(s1["pre_pqtl_publication"].isna().sum())
    assert info.pre_pqtl_publication_true == int((df["pre_pqtl_publication"] == True).sum())  # noqa: E712


def test_covariate_standardization_over_heldout_s1(joined):
    df, info = joined
    ref = df["in_s1"] & df["heldout"]
    assert df.loc[ref, "z_log10_neff"].mean() == pytest.approx(0, abs=1e-12)
    assert df.loc[ref, "z_log10_neff"].std(ddof=1) == pytest.approx(1)
    # three fixture rows have no effective N; variant-only copies of them carry it missing too
    copies = int((df["outcome_neff"].isna() & (df["variant"] != "primary")).sum())
    assert info.imputation.neff_missing_all == 3 + copies


def test_primary_and_pooled_sets(joined, analysis_sets):
    df, _ = joined
    s1, s2 = analysis_sets["S1"].frame, analysis_sets["S2"].frame
    assert set(s1["hypothesis_id"]) == set(df.loc[df["in_s1"] & df["heldout"], "hypothesis_id"])
    assert set(s2["hypothesis_id"]) == set(df.loc[df["in_s1"], "hypothesis_id"])
    assert (~s2["heldout"]).sum() > 0
    assert set(analysis_sets["S4"].frame["hypothesis_id"]) == set(df.loc[df["in_s4"] & df["heldout"], "hypothesis_id"])
    an = analysis_sets["S1"].analysed()
    assert an["y"].notna().all()
    assert set(an["status_24"]) <= {"advanced", "no_observed_advancement"}
    assert analysis_sets["S1"].info.n_analysed == len(an)
    assert sum(analysis_sets["S1"].info.excluded_outcome.values()) == len(s1) - len(an)


def test_s7_keeps_advanced_and_efficacy_coded_stops_only(analysis_sets):
    s7 = analysis_sets["S7"].frame
    kept = s7.loc[s7["y"].notna()]
    assert (kept.loc[kept["y"] == 0, "efficacy_coded"]).all()
    assert set(kept.loc[kept["y"] == 1, "hypothesis_id"]) == set(s7.loc[s7["advanced_24"] == 1, "hypothesis_id"])
    dropped = s7.loc[(s7["advanced_24"] == 0) & ~s7["efficacy_coded"]]
    assert len(dropped) > 0 and dropped["y"].isna().all()


def test_s13_s18_s19_s22_restrictions(analysis_sets):
    s1 = analysis_sets["S1"].frame
    s13 = analysis_sets["S13"].frame
    assert not s13["karim_launched"].any() and not s13["pilot_indication"].any()
    assert len(s13) == int((~s1["karim_launched"] & ~s1["pilot_indication"]).sum())
    assert analysis_sets["S18"].frame["single_protein_row"].all()
    s19 = analysis_sets["S19"].frame
    assert s19["blood_secreted_hpa"].all() and len(s19) > 0
    assert set(zip(s19["cls"], s19["s19_arm"])) <= {("aligned", "neutralizing_biologic"),
                                                    ("blocking", "small_molecule_blocker")}
    s22 = analysis_sets["S22"]
    assert s22.info.models == ("h4",) and not s22.frame["psych_only"].any()
    assert len(s22.frame) == int((~s1["psych_only"]).sum())


def test_s14_uses_maturation_outcomes(analysis_sets):
    for sid, m in (("S14a", 12), ("S14b", 36), ("S14c", 48)):
        f = analysis_sets[sid].frame
        assert np.array_equal(f["y"].to_numpy(), f[f"advanced_{m}"].astype(float).to_numpy(), equal_nan=True)


def test_s15_threshold_sets_move_supportive_counts_monotonically(analysis_sets):
    # hand-worked S15 states are in test_hand_sets; here only the ordering across thresholds
    s1 = analysis_sets["S1"].frame
    assert analysis_sets["S15e"].frame["S"].sum() <= s1["S"].sum() <= analysis_sets["S15d"].frame["S"].sum()
    assert not (analysis_sets["S15f"].frame["low_coverage"] == True).any()  # noqa: E712


def test_s16_replaces_state_only_where_available(analysis_sets):
    s1, s16 = analysis_sets["S1"].frame, analysis_sets["S16"].frame
    has = s1["s16_evidence_state"].notna()
    assert (s16.loc[has, "state"] == s1.loc[has, "s16_evidence_state"]).all()
    assert (s16.loc[~has, "state"] == s1.loc[~has, "evidence_state"]).all()


def test_s17_s20_s21(joined, analysis_sets):
    df, _ = joined
    s17 = analysis_sets["S17"].frame
    assert ((s17["s17_sentinel_p"] < 0.05) == (s17["S"] == 1)).all()
    s20 = set(analysis_sets["S20"].frame["hypothesis_id"])
    assert s20 == set(df.loc[df["heldout"] & (df["in_s1"] | ~df["strict_instrument"]), "hypothesis_id"])
    s21 = set(analysis_sets["S21"].frame["hypothesis_id"])
    assert s21 == set(df.loc[df["heldout"] & df["in_s21"], "hypothesis_id"])
    assert ((df["variant"] == "s21") == df["conflicted_row_restored"]).all()


def test_s8_s9_s17_not_formed_without_their_columns(joined):
    df, _ = joined
    df = df.copy()
    df.attrs["present_columns"] = {k: [c for c in v if c not in ("in_s8", "in_s9", "in_s21", "s17_sentinel_p")]
                                   for k, v in df.attrs["present_columns"].items()}
    sets = form_sets(df)
    for sid in ("S8", "S9", "S21", "S17"):
        assert not sets[sid].info.formed and "no column" in sets[sid].info.reason
        with pytest.raises(SetError):
            sets[sid].analysed()


def test_s8_s9_are_the_reformed_rows(joined, analysis_sets):
    df, _ = joined
    for sid, col, v in (("S8", "in_s8", "s8"), ("S9", "in_s9", "s9"), ("S21", "in_s21", "s21")):
        f = analysis_sets[sid].frame
        assert set(f["hypothesis_id"]) == set(df.loc[df[col] & df["heldout"], "hypothesis_id"])
        assert (f["cls"] == f["mechanism_class"]).all()
        assert (f["variant"] == v).any() and set(f["variant"]) <= {"primary", v}
        assert analysis_sets[sid].info.notes["variant_only_rows"] == int((f["variant"] == v).sum())
        # one row per hypothesis key within the set
        assert not f.duplicated(["gene_ensembl", "indication_id", "direction", "cls"]).any()
    s1 = set(analysis_sets["S1"].frame["hypothesis_id"])
    assert not (set(df.loc[df["variant"] != "primary", "hypothesis_id"]) & s1)


def test_s16_is_not_formed_when_no_smp_normalized_state_was_written(joined):
    df, _ = joined
    df = df.copy()
    assert form_sets(df)["S16"].info.formed and form_sets(df)["S16"].info.notes["replaced"] > 0
    df["s16_evidence_state"] = None
    s16 = form_sets(df)["S16"]
    assert not s16.info.formed and s16.info.reason.startswith("no deCODE SMP-normalized statistics were retrieved")
    assert s16.info.reason.endswith("S16 is not run")
    with pytest.raises(SetError):
        s16.analysed()
    one = df.copy()
    one.loc[one.index[(one["in_s1"] & one["heldout"]).to_numpy()][0], "s16_evidence_state"] = "supportive"
    assert form_sets(one)["S16"].info.formed and form_sets(one)["S16"].info.notes == {"replaced": 1}


def test_s11_not_formed_when_flag_missing_everywhere(joined):
    df, _ = joined
    df = df.copy()
    df["splicing_candidate"] = pd.array([pd.NA] * len(df), dtype="boolean")
    assert not form_sets(df)["S11"].info.formed


def test_s11_excludes_flagged_and_missing(joined, analysis_sets):
    s11 = analysis_sets["S11"]
    assert (s11.frame["splicing_candidate"] == False).all()  # noqa: E712
    s1 = analysis_sets["S1"].frame
    assert s11.info.notes["flag_missing_excluded"] == int(s1["splicing_candidate"].isna().sum())


@pytest.mark.parametrize("set_id,flag", [("S5", "protein_altering"), ("S15f", "low_coverage")])
def test_s5_and_s15f_keep_known_false_flags_and_count_the_missing(joined, set_id, flag):
    df, _ = joined
    df = df.copy()
    s1_rows = df.index[(df["in_s1"] & df["heldout"]).to_numpy()]
    false_rows = [i for i in s1_rows if df.at[i, flag] == False][:7]  # noqa: E712
    assert len(false_rows) == 7
    before = form_sets(df)[set_id]
    assert before.info.notes == {"flag_missing_excluded": 0}
    df[flag] = df[flag].astype("boolean")
    df.loc[false_rows, flag] = pd.NA
    after = form_sets(df)[set_id]
    assert after.info.notes == {"flag_missing_excluded": 7}
    assert set(before.frame["hypothesis_id"]) - set(after.frame["hypothesis_id"]) == set(df.loc[false_rows, "hypothesis_id"])
    assert (after.frame[flag] == False).all()  # noqa: E712
    s1 = form_sets(df)["S1"].frame
    assert len(after.frame) == int((s1[flag] == False).sum())  # noqa: E712
    assert len(s1) - len(after.frame) == 7 + int((s1[flag] == True).sum())  # noqa: E712


def test_s12_hand_example():
    f = pd.DataFrame({
        "gene_ensembl": ["g1", "g1", "g1", "g1", "g2"],
        "specificity_rank": [5, 1, 1, 0, 3],
        "indication_id": ["I9", "I2", "I1", "I7", "I1"],
        "programs": [["p3"], ["p1", "p2"], ["p2"], ["p1"], ["p1"]],
    })
    # order in g1: rank 0 (I7: p1) kept; rank 1 I1 (p2) kept; rank 1 I2 (p1,p2) dropped; rank 5 (p3) kept.
    # g2 is a separate gene, so p1 is not yet used there.
    assert list(s12_select(f)) == [True, False, True, True, True]


def test_s12_disagreement_with_stage_a_raises(joined):
    df, _ = joined
    df = df.copy()
    i = df.index[(df["in_s1"] & df["heldout"] & df["in_s12"])][0]
    df.loc[i, "in_s12"] = False
    with pytest.raises(SetError, match="S12"):
        form_sets(df)
