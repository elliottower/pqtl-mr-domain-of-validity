import numpy as np
import pandas as pd
import pytest

from susie_time_limit import (TimingError, classify, kaplan_meier_completed, pair_row, spawns_by_unit, summarize)

LOG = [{"what": "collect", "utc": "2026-10-05T00:00:00+00:00"},
       {"what": "spawn", "utc": "2026-10-08T03:31:57+00:00", "calls": [{"unit_key": "u1", "call_id": "a"},
                                                                    {"unit_key": "u2", "call_id": "b"}]},
       {"what": "spawn", "utc": "2026-10-08T18:18:16+00:00", "calls": [{"unit_key": "u2", "call_id": "c"}]}]
TIMEOUT = "susie did not complete: Rscript ran past 3600 s on tasks ['s15g']"
STOPPED = "susie did not complete: stopped in R: non-finite LD"


def test_classify_reads_the_basis_and_the_cause_in_the_note():
    assert classify("coloc.susie", "") == classify("abf_susie_criterion_not_met", "") == "completed"
    assert classify("abf_susie_did_not_complete", TIMEOUT) == "censored"
    assert classify("abf_susie_did_not_complete", STOPPED) == "stopped_in_r"
    with pytest.raises(TimingError):
        classify("abf_susie_not_run", "")
    with pytest.raises(TimingError):
        classify("abf_susie_did_not_complete", "something else")


def test_durations_are_whole_minutes_and_the_tier_is_the_latest_spawn_before_the_start():
    spawns = spawns_by_unit(LOG)
    r = pair_row("u1", "X", "coloc.susie", "", "2026-10-08T05:00:59+00:00", "2026-10-08T05:12:01+00:00", spawns)
    assert (r["minutes"], r["kind"], r["tier"]) == (12.0, "completed", "2cpu")
    late = pair_row("u2", "Y", "abf_susie_did_not_complete", TIMEOUT, "2026-10-08T19:00:00+00:00",
                    "2026-10-08T20:00:30+00:00", spawns)
    assert (late["minutes"], late["kind"], late["tier"]) == (60.0, "censored", "8cpu")


def test_a_spawn_inside_the_interval_or_more_than_70_minutes_makes_the_pair_interrupted():
    spawns = spawns_by_unit(LOG)
    across = pair_row("u2", "Z", "coloc.susie", "", "2026-10-08T18:00:00+00:00", "2026-10-08T18:40:00+00:00", spawns)
    assert (across["kind"], across["interrupted_because"]) == ("interrupted", "spawn inside the interval")
    long = pair_row("u1", "W", "coloc.susie", "", "2026-10-08T05:00:00+00:00", "2026-10-08T06:11:00+00:00", spawns)
    assert (long["kind"], long["interrupted_because"]) == ("interrupted", "over 70 minutes")
    edge = pair_row("u1", "V", "coloc.susie", "", "2026-10-08T05:00:00+00:00", "2026-10-08T06:10:00+00:00", spawns)
    assert edge["kind"] == "completed"


def test_a_pair_needs_a_spawn_before_it_and_an_end_after_its_start():
    spawns = spawns_by_unit(LOG)
    with pytest.raises(TimingError, match="no spawn"):
        pair_row("u1", "X", "coloc.susie", "", "2026-10-07T00:00:00+00:00", "2026-10-07T00:05:00+00:00", spawns)
    with pytest.raises(TimingError, match="before its outcome region"):
        pair_row("u1", "X", "coloc.susie", "", "2026-10-08T05:10:00+00:00", "2026-10-08T05:00:00+00:00", spawns)


def test_kaplan_meier_matches_a_hand_worked_example():
    # completions at 2, 2, 5; censored at 3 and 60. S: t=2 -> 1 - 2/5 = 0.6; t=5 -> 0.6 * (1 - 1/2) = 0.3.
    minutes = np.array([2, 2, 3, 5, 60], dtype=float)
    completed = np.array([1, 1, 0, 1, 0], dtype=bool)
    got = kaplan_meier_completed(minutes, completed, np.array([0, 1, 2, 4, 5, 60]))
    assert got == pytest.approx([0.0, 0.0, 0.4, 0.4, 0.7, 0.7])


def test_kaplan_meier_without_censoring_is_the_empirical_share():
    rng = np.random.default_rng()
    minutes = rng.integers(0, 61, size=5000).astype(float)
    grid = np.arange(61)
    got = kaplan_meier_completed(minutes, np.ones(len(minutes), dtype=bool), grid)
    assert got == pytest.approx([(minutes <= g).mean() for g in grid])


def test_summary_counts_each_kind_per_tier_and_leaves_interrupted_and_stopped_runs_off_the_curve():
    pairs = pd.DataFrame({"tier": ["2cpu"] * 5 + ["8cpu"] * 2,
                          "kind": ["completed", "completed", "censored", "interrupted", "stopped_in_r", "completed", "censored"],
                          "minutes": [10.0, 50.0, 60.0, 5.0, 1.0, 46.0, 60.0]})
    s = summarize(pairs)
    assert {k: s["2cpu"][k] for k in ("attempted", "completed", "censored", "stopped_in_r", "interrupted")} == {
        "attempted": 5, "completed": 2, "censored": 1, "stopped_in_r": 1, "interrupted": 1}
    assert s["2cpu"]["completed_in_minutes_45_to_60"] == 1 and s["8cpu"]["completed_in_minutes_45_to_60"] == 1
    curve = dict(zip(s["2cpu"]["curve"]["minute"], s["2cpu"]["curve"]["share_completed"]))
    assert curve[9] == 0 and curve[10] == pytest.approx(1 / 3) and curve[60] == pytest.approx(2 / 3)


def test_an_event_and_a_censoring_at_the_same_minute_count_the_event_first():
    # one completion and one censoring at 5, one completion at 10: S(5) = 1 - 1/3, S(10) = 2/3 * (1 - 1/1) = 0
    got = kaplan_meier_completed(np.array([5.0, 5.0, 10.0]), np.array([True, False, True]), np.array([4, 5, 10]))
    assert got == pytest.approx([0.0, 1 / 3, 1.0])


def test_the_hardware_tier_changes_exactly_at_the_8_cpu_respawn():
    log = [{"what": "spawn", "utc": "2026-10-08T17:30:44+00:00", "calls": [{"unit_key": "a", "call_id": "x"}]},
           {"what": "spawn", "utc": "2026-10-08T17:30:45+00:00", "calls": [{"unit_key": "b", "call_id": "y"}]}]
    spawns = spawns_by_unit(log)
    rows = [pair_row(u, "T", "coloc.susie", "", "2026-10-08T17:31:00+00:00", "2026-10-08T17:40:00+00:00", spawns)
            for u in ("a", "b")]
    assert [r["tier"] for r in rows] == ["2cpu", "8cpu"]


def test_a_spawn_at_the_start_is_the_runs_own_and_a_spawn_at_the_end_interrupts_it():
    log = [{"what": "spawn", "utc": "2026-10-08T05:00:00+00:00", "calls": [{"unit_key": "a", "call_id": "x"}]},
           {"what": "spawn", "utc": "2026-10-08T07:20:00+00:00", "calls": [{"unit_key": "b", "call_id": "y"}]},
           {"what": "spawn", "utc": "2026-10-08T06:00:00+00:00", "calls": [{"unit_key": "b", "call_id": "z"}]}]
    spawns = spawns_by_unit(log)
    at_start = pair_row("a", "T", "coloc.susie", "", "2026-10-08T05:00:00+00:00", "2026-10-08T05:20:00+00:00", spawns)
    assert at_start["kind"] == "completed"                                        # the interval is (start, end]
    at_end = pair_row("b", "T", "coloc.susie", "", "2026-10-08T07:00:00+00:00", "2026-10-08T07:20:00+00:00", spawns)
    assert (at_end["kind"], at_end["interrupted_because"]) == ("interrupted", "spawn inside the interval")
