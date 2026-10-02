"""Set membership, S19, no_phase and S15 against the hand-worked expectations of hand_fixture.py,
through the cross-stage path: written as sealed A/B/C outputs, seals logged in a copy of the frozen
PREREG.md, admitted by the run guard, then loaded, joined and formed into sets."""
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from prereg.log import append
from v8_run_guard import RunNotAuthorized

from stage_d.descriptive import outcome_category, table1_funnel
from stage_d.guard import StageInputError, logged_seals, sha256_file, verify_inputs
from stage_d.join import derive, join_stages, load_stage_tables
from stage_d.pipeline import prepare
from stage_d.schemas import OutcomeRow
from stage_d.sets import SET_ORDER, SetError, form_sets
from stage_d.synthetic import SYNTHETIC_SCRIPT_SHA256, synthetic_fingerprint, write_stage_outputs
from tests.hand_fixture import (EXPECTED, FLAG_MISSING_ANYWHERE, FLAG_MISSING_EXCLUDED, S1_EXCLUDED, S14A_EXCLUDED, S19_OUTCOMES, STATES,
                                 hand_tables)

PREREG = Path(os.environ.get("V8_PREREG", Path(__file__).resolve().parents[3] / "PREREG.md"))
TOKEN = "staged-token-0001"


def log(p: Path, event: str) -> None:
    append(p, "2026-10-09", event, "results not opened")


@pytest.fixture
def hand_stages(tmp_path):
    root = tmp_path / "stages"
    expected = write_stage_outputs(root, *hand_tables())
    prereg = tmp_path / "PREREG.md"
    shutil.copyfile(PREREG, prereg)
    return root, expected, prereg


def seal(prereg: Path, expected: dict[str, str]) -> None:
    for s in "ABC":
        log(prereg, f"SEAL stage={s} manifest_sha256={expected[s]}")
    log(prereg, f"RUN_START stage=D token={TOKEN}")


@pytest.fixture
def hand_sets(hand_stages):
    root, expected, prereg = hand_stages
    seal(prereg, expected)
    seals = logged_seals(prereg, TOKEN, root)
    assert seals == expected
    verify_inputs(root, seals)
    df, _ = derive(join_stages(*load_stage_tables(root)))
    return df, form_sets(df)


def ids(aset) -> set[str]:
    return set(aset.frame["hypothesis_id"])


# ---- run guard on the cross-stage path -------------------------------------------------------------

@pytest.mark.parametrize("events", [
    (),
    ("stage D prepared; SEAL stage=A pending",),
    ("SEAL stage=A manifest_sha256=<A>", "SEAL stage=B manifest_sha256=<B>", "SEAL stage=C manifest_sha256=<C>",
     "stage D started token=" + TOKEN),
    ("SEAL stage=A manifest_sha256=<A>", "SEAL stage=B manifest_sha256=<B>", f"RUN_START stage=D token={TOKEN}"),
    (f"RUN_START stage=D token={TOKEN}", "SEAL stage=A manifest_sha256=<A>", "SEAL stage=B manifest_sha256=<B>",
     "SEAL stage=C manifest_sha256=<C>"),
])
def test_stage_d_refuses_without_logged_seals_before_an_exact_run_start(hand_stages, events):
    root, expected, prereg = hand_stages
    for e in events:
        log(prereg, e.replace("<A>", expected["A"]).replace("<B>", expected["B"]).replace("<C>", expected["C"]))
    with pytest.raises(RunNotAuthorized):
        logged_seals(prereg, TOKEN, root)


def test_stage_d_refuses_a_manifest_edited_after_its_seal_was_logged(hand_stages):
    root, expected, prereg = hand_stages
    seal(prereg, expected)
    m = root / "C" / "output" / "MANIFEST.tsv"
    m.write_text(m.read_text().replace(SYNTHETIC_SCRIPT_SHA256, "0" * 64, 1))
    assert sha256_file(m) != expected["C"]
    with pytest.raises(RunNotAuthorized, match="stage C"):
        logged_seals(prereg, TOKEN, root)


def test_prepare_runs_on_the_logged_seals(hand_stages, tmp_path):
    root, expected, prereg = hand_stages
    seal(prereg, expected)
    seals = logged_seals(prereg, TOKEN, root)
    plan = prepare(root, seals, tmp_path / "work", synthetic_fingerprint(seals))
    assert {k: v["manifest_sha256"] for k, v in plan.seals.items()} == expected


def test_a_caller_supplied_hash_that_is_not_the_manifest_still_refuses(hand_stages):
    root, expected, _ = hand_stages
    with pytest.raises(StageInputError, match="logged"):
        verify_inputs(root, {**expected, "A": "0" * 64})


# ---- memberships -----------------------------------------------------------------------------------

@pytest.mark.parametrize("set_id", sorted(EXPECTED))
def test_hand_worked_membership(hand_sets, set_id):
    _, sets = hand_sets
    assert ids(sets[set_id]) == EXPECTED[set_id]


def test_every_registered_set_is_present_and_hand_checked_sets_have_their_members(hand_sets):
    _, sets = hand_sets
    assert set(sets) == set(SET_ORDER) and len(SET_ORDER) == 30
    assert {"S1", "S12", "S19"} <= set(EXPECTED)
    for set_id, members in EXPECTED.items():
        assert ids(sets[set_id]) == members, set_id


def test_s19_is_neutralizing_biologics_vs_small_molecule_blockers_of_blood_secreted_targets(hand_sets):
    _, sets = hand_sets
    s19 = sets["S19"]
    an = s19.analysed()
    assert dict(zip(an["hypothesis_id"], an["y"])) == S19_OUTCOMES
    assert s19.info.excluded_outcome == {"no_phase": 1, "business_only": 1}
    assert s19.info.notes == {"blood_secreted_aligned_not_neutralizing_biologic": 2,
                              "blood_secreted_blocking_not_small_molecule": 0}


def test_s12_hand_worked_selection_agrees_with_stage_a_and_a_wrong_flag_raises(hand_sets):
    df, sets = hand_sets
    assert ids(sets["S12"]) == EXPECTED["S12"]
    bad = df.copy()
    bad.loc[bad["hypothesis_id"] == "r04", "in_s12"] = True
    with pytest.raises(SetError, match="S12"):
        form_sets(bad)


# ---- flag restrictions: S5, S11, S15f ----------------------------------------------------------------

@pytest.mark.parametrize("set_id,flag", [("S5", "protein_altering"), ("S11", "splicing_candidate"), ("S15f", "low_coverage")])
def test_a_flag_restriction_keeps_only_known_false_rows_and_counts_the_missing(hand_sets, set_id, flag):
    df, sets = hand_sets
    s = sets[set_id]
    assert ids(s) == EXPECTED[set_id]
    assert s.frame[flag].notna().all() and not s.frame[flag].astype(bool).any()
    assert s.info.notes == {"flag_missing_excluded": FLAG_MISSING_EXCLUDED[set_id]}
    s1 = sets["S1"].frame
    dropped = s1.loc[~s1["hypothesis_id"].isin(ids(s))]
    assert int(dropped[flag].isna().sum()) == FLAG_MISSING_EXCLUDED[set_id]
    assert int((dropped[flag] == True).sum()) == len(dropped) - FLAG_MISSING_EXCLUDED[set_id]  # noqa: E712
    assert int(df[flag].isna().sum()) == FLAG_MISSING_ANYWHERE[set_id]   # a missing flag outside S1 is not counted


# ---- no_phase -----------------------------------------------------------------------------------------

def test_no_phase_is_a_missing_outcome_counted_in_the_exclusion_report(hand_sets):
    df, sets = hand_sets
    r03 = df.loc[df["hypothesis_id"] == "r03"].iloc[0]
    assert r03["status_24"] == "no_phase" and np.isnan(sets["S1"].frame.set_index("hypothesis_id").loc["r03", "y"])
    assert sets["S1"].info.excluded_outcome == S1_EXCLUDED
    assert sets["S1"].info.n_analysed == len(EXPECTED["S1"]) - sum(S1_EXCLUDED.values())
    assert sets["S14a"].info.excluded_outcome == S14A_EXCLUDED
    funnel = {r["step"]: r for r in table1_funnel(pd.DataFrame(
        {"step": [], "remaining": [], "excluded": [], "reason": []}), sets["S1"].frame)}
    assert funnel["primary outcome: no_phase"]["excluded"] == 1
    assert outcome_category(sets["S1"].frame.set_index("hypothesis_id").loc[["r03"]]).tolist() == ["no_phase"]


@pytest.mark.parametrize("advanced", ["", "0"])
def test_outcome_schema_accepts_no_phase_only_with_an_empty_outcome(advanced):
    rec = {"hypothesis_id": "h", "status_24": "no_phase", "advanced_24": advanced, "status_12": "no_phase",
           "status_36": "no_phase", "status_48": "no_phase", "advanced_12": "", "advanced_36": "", "advanced_48": "",
           "efficacy_coded": "false", "approved": "false"}
    if advanced:
        with pytest.raises(ValueError):
            OutcomeRow.model_validate(rec)
    else:
        assert OutcomeRow.model_validate(rec).status_24 == "no_phase"


# ---- S15 against hand-worked states -------------------------------------------------------------------

@pytest.mark.parametrize("key", sorted(STATES))
def test_s15_states_are_the_hand_worked_ones(hand_sets, key):
    _, sets = hand_sets
    hid, set_id = key
    f = sets[set_id].frame.set_index("hypothesis_id")
    assert f.loc[hid, "state"] == STATES[key]
    assert f.loc[hid, "S"] == float(STATES[key] == "supportive")
