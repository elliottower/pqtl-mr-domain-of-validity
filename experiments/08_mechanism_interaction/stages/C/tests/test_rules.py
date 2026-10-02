from datetime import date

import pytest

from models import ACTIVE_STATUSES, OUTCOME_FREEZE_DATE, StageCError, Status, StopCode, Trial
from rules import classify_why_stopped, code_hypothesis, combine_stop_codes, months_before, post_freeze_counts

F = OUTCOME_FREEZE_DATE
_n = iter(range(10_000_000, 99_999_999))


def trial(phases=("PHASE2",), status="COMPLETED", why="", start=None, end=None, submitted=None,
          updated=None) -> Trial:
    return Trial(nct_id=f"NCT{next(_n):08d}", phases=tuple(phases), overall_status=status, why_stopped=why,
                 start_date=start, completion_date=end, first_submit_date=submitted, last_update_submit_date=updated)



# ---- dates -------------------------------------------------------------------------------

def test_outcome_freeze_date_is_the_freeze_commit_date():
    assert F == date(2026, 9, 30)


def test_months_before_clamps_the_day():
    assert months_before(date(2026, 3, 31), 1) == date(2026, 2, 28)
    assert months_before(F, 24) == date(2024, 9, 30)
    assert months_before(F, 48) == date(2022, 9, 30)


# ---- primary status branches ---------------------------------------------------------------

@pytest.mark.parametrize("stage", ["PHASE_3", "PREAPPROVAL", "APPROVAL"])
def test_advanced_when_any_drug_reached_phase3_or_later_even_with_active_trials(stage):
    o = code_hypothesis(["PHASE_2", stage], [trial(status="RECRUITING")], F)
    assert set(o.statuses.values()) == {Status.advanced}
    assert set(o.advanced.values()) == {1}
    assert o.stop_code == StopCode.none


def test_phase_2_3_is_not_advanced():
    o = code_hypothesis(["PHASE_2_3"], [trial(end=date(2015, 1, 1))], F)
    assert o.statuses[24] == Status.no_observed_advancement
    assert o.advanced[24] == 0


@pytest.mark.parametrize("status", sorted(ACTIVE_STATUSES))
def test_any_ongoing_trial_of_any_phase_makes_the_hypothesis_active(status):
    old_p2 = trial(end=date(2010, 1, 1))
    o = code_hypothesis(["PHASE_2"], [old_p2, trial(phases=("PHASE1",), status=status)], F)
    assert set(o.statuses.values()) == {Status.active}
    assert set(o.advanced.values()) == {None}


@pytest.mark.parametrize("status", ["UNKNOWN", "COMPLETED", "TERMINATED", "WITHDRAWN", "SUSPENDED"])
def test_non_active_statuses_do_not_block_the_zero(status):
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2010, 1, 1)),
                                      trial(phases=("PHASE1",), status=status, why="efficacy")], F)
    assert o.statuses[24] == Status.no_observed_advancement


def test_maturation_windows_differ_by_last_phase2_end():
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2021, 1, 1)), trial(end=date(2023, 6, 1))], F)
    assert o.statuses == {12: Status.no_observed_advancement, 24: Status.no_observed_advancement,
                          36: Status.no_observed_advancement, 48: Status.active}
    assert o.advanced == {12: 0, 24: 0, 36: 0, 48: None}


@pytest.mark.parametrize("end,expected", [(date(2024, 9, 30), Status.no_observed_advancement),
                                          (date(2024, 10, 1), Status.active)])
def test_24_month_boundary_is_inclusive(end, expected):
    assert code_hypothesis(["PHASE_2"], [trial(end=end)], F).statuses[24] == expected


def test_terminated_phase2_end_counts_as_termination_date():
    o = code_hypothesis(["PHASE_2"], [trial(status="TERMINATED", why="slow accrual", end=date(2019, 5, 1))], F)
    assert o.statuses[24] == Status.no_observed_advancement
    assert o.stop_code == StopCode.unclassified


def test_phase_1_2_and_2_3_registrations_count_as_phase2_ends():
    for phases in (("PHASE1", "PHASE2"), ("PHASE2", "PHASE3")):
        o = code_hypothesis(["PHASE_2"], [trial(phases=phases, end=date(2018, 1, 1))], F)
        assert o.statuses[24] == Status.no_observed_advancement


# ---- undated (pre-stage amendments table) --------------------------------------------------

def test_undated_when_no_trials_are_linked():
    o = code_hypothesis(["PHASE_2"], [], F)
    assert set(o.statuses.values()) == {Status.undated}
    assert set(o.advanced.values()) == {None}


def test_undated_when_phase2_trials_have_no_end_date():
    o = code_hypothesis(["PHASE_2"], [trial(end=None), trial(phases=("PHASE1",), end=date(2010, 1, 1))], F)
    assert o.statuses[24] == Status.undated


def test_undated_when_the_only_dated_phase2_end_is_withdrawn_or_suspended():
    o = code_hypothesis(["PHASE_2"], [trial(status="WITHDRAWN", why="no participants", end=date(2012, 1, 1))], F)
    assert o.statuses[24] == Status.undated


def test_undated_does_not_apply_to_advanced_hypotheses():
    assert code_hypothesis(["PHASE_3"], [], F).statuses[24] == Status.advanced


# ---- business-only exclusion -------------------------------------------------------------------

def test_business_only_when_every_stop_is_business():
    trials = [trial(status="TERMINATED", why="Sponsor decision", end=date(2015, 1, 1)),
              trial(status="WITHDRAWN", why="Lack of FUNDING"),
              trial(end=date(2014, 1, 1))]
    o = code_hypothesis(["PHASE_2"], trials, F)
    assert set(o.statuses.values()) == {Status.business_only}
    assert set(o.advanced.values()) == {None}
    assert o.stop_code == StopCode.business
    assert o.why_stopped_texts == "||".join(t.why_stopped for t in sorted(trials[:2], key=lambda t: t.nct_id))


def test_one_non_business_stop_lifts_the_business_exclusion():
    o = code_hypothesis(["PHASE_2"], [trial(status="TERMINATED", why="strategic reasons", end=date(2015, 1, 1)),
                                      trial(status="TERMINATED", why="futility", end=date(2016, 1, 1))], F)
    assert o.statuses[24] == Status.no_observed_advancement
    assert o.stop_code == StopCode.efficacy


def test_business_stop_on_an_active_hypothesis_is_labelled_active():
    o = code_hypothesis(["PHASE_2"], [trial(status="TERMINATED", why="business", end=date(2015, 1, 1)),
                                      trial(status="RECRUITING")], F)
    assert o.statuses[24] == Status.active


# ---- why_stopped keyword table and precedence ---------------------------------------------------

@pytest.mark.parametrize("text,code", [
    ("Business decision", StopCode.business), ("STRATEGIC reasons", StopCode.business),
    ("commercial", StopCode.business), ("portfolio prioritization", StopCode.business),
    ("sponsor decision", StopCode.business), ("funding withdrawn", StopCode.business),
    ("financial", StopCode.business),
    ("lack of efficacy", StopCode.efficacy), ("Futility", StopCode.efficacy), ("deemed futile", StopCode.efficacy),
    ("lack of benefit", StopCode.efficacy), ("did not meet primary endpoint", StopCode.efficacy),
    ("safety", StopCode.safety), ("Adverse events", StopCode.safety), ("hepatotoxicity", StopCode.safety),
    ("poor tolerability", StopCode.safety),
    ("PI left the institution", StopCode.unclassified), ("", StopCode.unclassified),
])
def test_keyword_table(text, code):
    assert classify_why_stopped(text) == code


@pytest.mark.parametrize("text,code", [
    ("sponsor decision after safety review", StopCode.safety),
    ("business decision following futility analysis", StopCode.efficacy),
    ("safety and efficacy concerns", StopCode.efficacy),
])
def test_one_text_matching_several_codes_takes_the_precedence(text, code):
    assert classify_why_stopped(text) == code


@pytest.mark.parametrize("codes,expected", [
    ([StopCode.business, StopCode.unclassified], StopCode.unclassified),
    ([StopCode.business, StopCode.safety], StopCode.safety),
    ([StopCode.safety, StopCode.efficacy, StopCode.business], StopCode.efficacy),
    ([StopCode.unclassified, StopCode.safety], StopCode.safety),
    ([StopCode.business], StopCode.business),
])
def test_precedence_efficacy_safety_unclassified_business(codes, expected):
    assert combine_stop_codes(codes) == expected


def test_combine_needs_a_code():
    with pytest.raises(StageCError):
        combine_stop_codes([])


def test_hypothesis_stop_code_uses_precedence_across_trials():
    o = code_hypothesis(["PHASE_2"], [trial(status="TERMINATED", why="safety", end=date(2015, 1, 1)),
                                      trial(status="SUSPENDED", why="administrative"),
                                      trial(status="TERMINATED", why="business", end=date(2016, 1, 1))], F)
    assert o.stop_code == StopCode.safety


# ---- completed, no successor --------------------------------------------------------------------

def test_completed_without_stops_is_completed_no_successor_not_efficacy():
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2016, 1, 1))], F)
    assert o.stop_code == StopCode.completed_no_successor
    assert o.why_stopped_texts == ""


def test_completed_trial_why_stopped_text_is_never_read():
    o = code_hypothesis(["PHASE_2"], [trial(status="COMPLETED", why="futility", end=date(2016, 1, 1))], F)
    assert o.stop_code == StopCode.completed_no_successor


def test_stops_take_priority_over_completed_no_successor():
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2016, 1, 1)),
                                      trial(status="TERMINATED", why="toxicity", end=date(2017, 1, 1))], F)
    assert o.stop_code == StopCode.safety


def test_no_trials_gives_no_stop_code():
    assert code_hypothesis(["PHASE_2"], [], F).stop_code == StopCode.none


# ---- trial records dated after the freeze -------------------------------------------------------

def test_trial_registered_after_freeze_is_ignored():
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2016, 1, 1)),
                                      trial(status="RECRUITING", submitted=date(2026, 10, 2))], F)
    assert o.statuses[24] == Status.no_observed_advancement


def test_post_freeze_counts_are_a_diagnostic_and_change_no_label():
    trials = [trial(end=date(2016, 1, 1), submitted=date(2015, 1, 1), updated=date(2026, 10, 1)),
              trial(end=date(2017, 1, 1), submitted=date(2016, 1, 1), updated=date(2026, 9, 30)),
              trial(status="RECRUITING", submitted=date(2026, 10, 1), updated=date(2026, 10, 1)),
              trial(end=date(2018, 1, 1))]
    assert post_freeze_counts(trials, F) == (1, 2)
    stripped = [t.model_copy(update={"last_update_submit_date": None}) for t in trials]
    assert code_hypothesis(["PHASE_2"], trials, F) == code_hypothesis(["PHASE_2"], stripped, F)


def test_end_after_freeze_means_the_trial_was_running_at_the_freeze():
    o = code_hypothesis(["PHASE_2"], [trial(end=date(2016, 1, 1)),
                                      trial(status="TERMINATED", why="futility", end=date(2026, 11, 1))], F)
    assert o.statuses[24] == Status.active
    assert o.stop_code == StopCode.completed_no_successor


# ---- no retrievable phase ---------------------------------------------------------------------

def test_no_phase_rows_is_excluded_and_counted():
    o = code_hypothesis([], [trial(end=date(2016, 1, 1))], F)
    assert set(o.statuses.values()) == {Status.no_phase}
    assert set(o.advanced.values()) == {None}


# ---- time to Phase III -------------------------------------------------------------------------

def test_time_to_phase3_event():
    o = code_hypothesis(["PHASE_3"], [trial(start=date(2015, 3, 1)), trial(start=date(2016, 1, 1)),
                                      trial(phases=("PHASE3",), start=date(2018, 3, 1))], F)
    assert (o.phase2_start, o.phase3_start, o.phase3_event) == (date(2015, 3, 1), date(2018, 3, 1), 1)
    assert o.time_to_phase3_days == (date(2018, 3, 1) - date(2015, 3, 1)).days


def test_time_to_phase3_censored_at_freeze():
    o = code_hypothesis(["PHASE_2"], [trial(start=date(2020, 1, 1), end=date(2021, 1, 1))], F)
    assert (o.phase3_start, o.phase3_event) == (None, 0)
    assert o.time_to_phase3_days == (F - date(2020, 1, 1)).days


def test_phase3_start_after_freeze_is_censored():
    o = code_hypothesis(["PHASE_2"], [trial(start=date(2020, 1, 1)),
                                      trial(phases=("PHASE3",), status="NOT_YET_RECRUITING", start=date(2026, 12, 1))], F)
    assert o.phase3_event == 0
    assert o.time_to_phase3_days == (F - date(2020, 1, 1)).days


def test_phase3_before_the_phase2_origin_is_not_the_event():
    o = code_hypothesis(["PHASE_3"], [trial(phases=("PHASE3",), start=date(2010, 1, 1)),
                                      trial(start=date(2012, 1, 1))], F)
    assert o.phase3_event == 0


def test_phase2_3_registration_is_an_origin_not_an_event():
    o = code_hypothesis(["PHASE_2_3"], [trial(phases=("PHASE2", "PHASE3"), start=date(2014, 1, 1))], F)
    assert (o.phase2_start, o.phase3_event) == (date(2014, 1, 1), 0)


def test_no_dated_phase2_start_excludes_from_time_to_phase3():
    o = code_hypothesis(["PHASE_3"], [trial(start=None), trial(phases=("PHASE3",), start=date(2018, 1, 1)),
                                      trial(status="WITHDRAWN", start=date(2011, 1, 1))], F)
    assert (o.phase2_start, o.phase3_start, o.time_to_phase3_days, o.phase3_event) == (None, None, None, None)


# ---- Phase III success and approval ----------------------------------------------------------------

def test_phase3_success_approved():
    o = code_hypothesis(["APPROVAL", "PHASE_2"], [], F)
    assert (o.approved, o.phase3_success) == (True, 1)


def test_phase3_success_zero_when_mature_and_inactive():
    o = code_hypothesis(["PHASE_3"], [trial(phases=("PHASE3",), status="TERMINATED", why="futility",
                                            end=date(2020, 1, 1))], F)
    assert (o.approved, o.phase3_success) == (False, 0)


@pytest.mark.parametrize("p3", [trial(phases=("PHASE3",), status="RECRUITING"),
                                trial(phases=("PHASE3",), end=date(2025, 1, 1))])
def test_phase3_success_pending(p3):
    assert code_hypothesis(["PHASE_3"], [trial(phases=("PHASE3",), end=date(2015, 1, 1)), p3], F).phase3_success is None


def test_phase3_success_pending_when_phase3_end_is_undated():
    assert code_hypothesis(["PHASE_3"], [trial(phases=("PHASE3",), end=None)], F).phase3_success is None


def test_phase3_success_pending_without_phase3_trial_on_record():
    assert code_hypothesis(["PREAPPROVAL"], [], F).phase3_success is None


def test_phase3_success_empty_when_not_advanced():
    assert code_hypothesis(["PHASE_2"], [trial(phases=("PHASE3",), end=date(2015, 1, 1))], F).phase3_success is None


# ---- last Phase II completion or termination (descriptive table 7) --------------------------

def test_last_phase2_end_is_the_latest_ended_phase2_at_the_freeze():
    trials = [
        trial(phases=("PHASE2",), status="COMPLETED", end=date(2019, 5, 1)),
        trial(phases=("PHASE1", "PHASE2"), status="TERMINATED", why="futility", end=date(2020, 2, 1)),
        trial(phases=("PHASE2",), status="COMPLETED", end=date(2027, 1, 1)),        # ends after the freeze: ongoing
        trial(phases=("PHASE3",), status="COMPLETED", end=date(2022, 6, 1)),        # not Phase II
        trial(phases=("PHASE2",), status="WITHDRAWN", end=date(2021, 1, 1)),        # withdrawn is not an end
        trial(phases=("PHASE2",), status="COMPLETED", end=date(2021, 3, 1),
              submitted=date(2026, 10, 2)),                                         # registered after the freeze
    ]
    o = code_hypothesis(["PHASE_2"], trials, F)
    assert o.last_phase2_end == date(2020, 2, 1)
    assert o.statuses[24] == Status.active  # the Phase II ending after the freeze is ongoing


def test_last_phase2_end_is_the_date_the_maturation_rule_reads():
    for end, mature_24 in ((date(2024, 9, 30), True), (date(2024, 10, 1), False)):
        o = code_hypothesis(["PHASE_2"], [trial(end=end)], F)
        assert o.last_phase2_end == end
        assert (o.statuses[24] == Status.no_observed_advancement) is mature_24


def test_last_phase2_end_kept_for_advanced_and_no_phase_hypotheses():
    t = [trial(end=date(2016, 4, 1))]
    assert code_hypothesis(["PHASE_3"], t, F).last_phase2_end == date(2016, 4, 1)
    assert code_hypothesis([], t, F).last_phase2_end == date(2016, 4, 1)
    assert code_hypothesis(["PHASE_2"], [], F).last_phase2_end is None
