"""The rounding-aware CI-versus-p rule (stage_b/outcome_files.py, CI_P_RULE_ROUNDING_AWARE), beside the
frozen rule it is reported with, on hand-built rows and on seeded random rows against a brute-force
oracle."""
import math
from decimal import Decimal
from statistics import NormalDist

import numpy as np
import pytest
from test_outcome_files import SSF_OR_CI_HEADER, exact_or_ci_row, or_ci_row

from stage_b.outcome_files import (AGREE, CI_P_LOG10_TOLERANCE, CI_P_RULE_FROZEN, CI_P_RULE_ROUNDING_AWARE, CI_Z, DISAGREE,
                                   LOWER_LIMIT_AT_ZERO, NOT_CHECKABLE, ROUNDING_UNINFORMATIVE, CiPCheck, SsfReader,
                                   ci_p_rounding_aware_verdict, ci_p_rule_diagnostic, latent_wald_extremes, log10_two_sided_p, printed_digits,
                                   wald_ratio_extremes)

Z99 = NormalDist().inv_cdf(0.995)
NEG_LOG10_OR_CI_HEADER = [c if c != "p_value" else "neg_log_10_p_value" for c in SSF_OR_CI_HEADER]
SEED = 20261005
# OR 1.0 with limits 1.0 and 1.1: the limits' intervals meet at 1.05, inside OR's [0.95, 1.05], which holds 1
UNINFORMATIVE = ("1.0", "1.0", "1.1")


def wald_p(beta: float, se: float) -> float:
    return math.erfc(abs(beta / se) / math.sqrt(2))


def printed_row(beta: float, se: float, decimals: int = 3, p: str | None = None, z: float = CI_Z) -> list[str]:
    """OR = exp(beta) and limits exp(beta -/+ z se), each printed with `decimals` decimals; p the
    Wald p of (beta, se) at full precision unless given."""
    def fmt(x: float) -> str:
        return f"{x:.{decimals}f}"
    return or_ci_row(odds=fmt(math.exp(beta)), lower=fmt(math.exp(beta - z * se)), upper=fmt(math.exp(beta + z * se)),
                     p=repr(wald_p(beta, se)) if p is None else p)


def small_effect_rows(n: int, z: float = CI_Z, p_format=repr) -> list[list[str]]:
    """Small effects (|ln OR| of a few 1e-3) with SE 5e-4 to 2e-3, OR and limits printed with 3 decimals,
    p the exact Wald p of the unrounded values."""
    rng = np.random.default_rng(SEED)
    rows = []
    for i, (b, s) in enumerate(zip(rng.normal(0, 0.003, n), rng.uniform(0.0005, 0.002, n))):
        row = printed_row(b, s, z=z, p=p_format(wald_p(b, s)))
        row[1], row[8] = str(1000 + i), f"rs{i}"
        rows.append(row)
    return rows


def uninformative_row(p: str = "0.5") -> list[str]:
    return or_ci_row(*UNINFORMATIVE, p=p)


def mismatch_rows(n: int) -> list[list[str]]:
    """99% limits read as 95%, printed to 4 decimals: a genuine semantic mismatch under either rule."""
    return [printed_row(0.03 * (5 + i % 4), 0.03, decimals=4, z=Z99) for i in range(n)]


def run(rows, rule: str = CI_P_RULE_ROUNDING_AWARE, header=SSF_OR_CI_HEADER) -> dict:
    reader = SsfReader(header, 0, "GWAS Catalog X", "or_ci_derived_se")
    check = CiPCheck(reader, rule)
    for r in rows:
        got = reader.row(r)
        assert not isinstance(got, str), got
        check.add(got, r)
    return check.result()


def verdict(odds: str, lower: str, upper: str, p: str, neg_log10: bool = False) -> str:
    return ci_p_rounding_aware_verdict((odds, lower, upper), p, neg_log10, CI_Z)


def counts(result: dict) -> tuple[int, int, int, bool]:
    return result["rows_checked"], result["rows_agree"], result["rows_not_checkable"], result["passed"]


def z_of_log10_p(log10_p: float) -> float:
    """The |z| whose two-sided p has log10 `log10_p`, by bisection on log10_two_sided_p."""
    lo, hi = 0.0, 1e4
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if log10_two_sided_p(mid) > log10_p else (lo, mid)
    return (lo + hi) / 2


def exact_texts(beta: float, se: float) -> tuple[str, str, str]:
    return repr(math.exp(beta)), repr(math.exp(beta - CI_Z * se)), repr(math.exp(beta + CI_Z * se))


# ---- the rule on single rows -----------------------------------------------------------------------------------

def test_a_small_effect_row_that_fails_the_frozen_rule_only_for_the_rounding_of_its_limits_passes():
    # beta 0.002, SE 5e-4: z = 4, p = 6.3e-5; printed 1.002 (1.001, 1.003) gives SE 5.09e-4, p 8.7e-5, 0.14 log10 off
    row = printed_row(0.002, 0.0005)
    assert (row[4], row[10], row[9]) == ("1.002", "1.001", "1.003")
    got = run([row])
    assert counts(got) == (1, 1, 0, True)
    assert counts(got["frozen_rule"]) == (1, 0, 0, False)


def test_a_row_whose_rounding_allows_any_p_is_rounding_uninformative_and_never_an_agreement():
    ext = latent_wald_extremes(UNINFORMATIVE)
    assert ext.zero_width_reachable and ext.contains_zero_effect and (ext.low, ext.high) == (0.0, math.inf)
    for p in ("1e-300", "0.01", "0.5", "0.999"):
        assert verdict(*UNINFORMATIVE, p) == ROUNDING_UNINFORMATIVE
    got = run([uninformative_row(p) for p in ("1e-300", "0.01", "0.5", "0.999")])
    assert counts(got) == (0, 0, 4, False)
    assert got["rows_rounding_uninformative"] == 4 and got["rows_lower_limit_rounding_reaches_zero"] == 0


@pytest.mark.parametrize("texts,why", [
    (("1.01", "1.00", "1.01"), "zero width reachable, OR's interval [1.005, 1.015] excludes 1"),
    (("1.00", "0.90", "1.11"), "OR's interval holds 1, the limits' intervals do not meet"),
    (("1.00", "0.9", "1.0"), "OR's interval holds 1, the limits' intervals meet at 0.95, outside OR's interval"),
])
def test_a_row_is_uninformative_only_when_both_halves_of_the_full_span_hold(texts, why):
    ext = latent_wald_extremes(texts)
    assert not (ext.zero_width_reachable and ext.contains_zero_effect), why
    assert verdict(*texts, "0.5") != ROUNDING_UNINFORMATIVE


def test_uninformative_rows_do_not_rescue_a_file_with_a_few_genuine_mismatches():
    rows = [uninformative_row() for _ in range(1000)] + mismatch_rows(5)
    got = run(rows)
    assert counts(got) == (5, 0, 1000, False) and got["rows_rounding_uninformative"] == 1000
    # counted as agreements, the 1000 rows would have carried the file: 1000 of 1005 is 99.5%
    assert 100 * 1000 >= 95 * 1005


def test_a_file_whose_rows_are_all_uninformative_or_not_checkable_fails_closed():
    got = run([uninformative_row() for _ in range(50)] + [or_ci_row(p="#NA"), or_ci_row(p="1")])
    assert counts(got) == (0, 0, 52, False) and got["rows_rounding_uninformative"] == 50


@pytest.mark.parametrize("bad,passed", [(5, True), (6, False)])
def test_the_threshold_is_95_percent_of_the_informative_rows_exactly(bad, passed):
    good = [exact_or_ci_row(0.1 + 0.001 * i, 0.05) for i in range(100 - bad)]
    got = run(good + mismatch_rows(bad) + [uninformative_row() for _ in range(50)])
    assert counts(got) == (100, 100 - bad, 50, passed) and got["rows_rounding_uninformative"] == 50


def test_corner_combinations_no_single_latent_row_can_take_are_not_used():
    # OR 1.10 in [1.095, 1.105]; L 1.0 in [0.95, 1.05]; U 1.1 in [1.05, 1.15]. The corners L_hi = U_lo = 1.05
    # give a zero width, but OR >= 1.095 cannot lie inside [L, U] = [1.05, 1.05]: the narrowest interval
    # that holds OR is [1.05, OR], so |ln OR| / (ln U - ln L) <= ln 1.095 / ln(1.095 / 1.05).
    texts = ("1.10", "1.0", "1.1")
    ext = latent_wald_extremes(texts)
    assert ext.high == pytest.approx(math.log(1.095) / math.log(1.095 / 1.05), rel=1e-12)
    assert ext.low == pytest.approx(math.log(1.095) / math.log(1.15 / 0.95), rel=1e-12)
    p_floor = log10_two_sided_p(ext.high * 2 * CI_Z)                     # about 2.3e-17
    assert verdict(*texts, f"{10 ** p_floor:.2e}") == AGREE
    assert verdict(*texts, "1e-30") == DISAGREE                          # the corners would have allowed it
    # OR 1.00 in [0.995, 1.005] holds 1, and the limits' intervals meet at 0.95: the corners say any p;
    # no latent row puts OR at 0.95, and the implied p is at least 0.67
    texts = ("1.00", "0.9", "1.0")
    ext = latent_wald_extremes(texts)
    assert ext.contains_zero_effect and not ext.zero_width_reachable
    assert ext.high == pytest.approx(math.log(1 / 0.995) / math.log(0.995 / 0.95), rel=1e-12)
    assert (verdict(*texts, "0.9"), verdict(*texts, "0.01")) == (AGREE, DISAGREE)


def test_a_row_whose_or_lies_outside_its_interval_beyond_rounding_has_no_latent_row_and_disagrees():
    assert latent_wald_extremes(("1.4", "1.6", "1.9")) is None                # OR_hi 1.45 < L_lo 1.55
    assert latent_wald_extremes(("2.1", "1.6", "1.9")) is None                # OR_lo 2.05 > U_hi 1.95
    assert verdict("1.4", "1.6", "1.9", "0.01") == verdict("2.1", "1.6", "1.9", "0.01") == DISAGREE
    # OR_hi = L_lo = 1.55: the one feasible OR is 1.55, on its interval's lower limit, whatever the doubles round to
    touching = latent_wald_extremes(("1.5", "1.6", "1.9"))
    assert touching.low == pytest.approx(math.log(1.55) / math.log(1.95 / 1.55), rel=1e-9)
    assert touching.high == pytest.approx(math.log(1.55) / math.log(1.85 / 1.55), rel=1e-9)


@pytest.mark.parametrize("shift,expected", [(0.099, AGREE), (-0.099, AGREE), (0.101, DISAGREE), (-0.101, DISAGREE)])
def test_the_intervals_meet_up_to_the_0_1_log10_widening_and_no_further(shift, expected):
    beta, se = 0.3, 0.05                              # printed at full precision: rounding intervals of ~1e-16
    texts = exact_texts(beta, se)
    implied = log10_two_sided_p(beta / ((math.log(float(texts[2])) - math.log(float(texts[1]))) / (2 * CI_Z)))
    assert verdict(*texts, repr(10 ** (implied + shift))) == expected


@pytest.mark.parametrize("odds,lower,upper,p,expected", [
    ("1.00", "0.90", "1.11", "0.9", AGREE),           # OR's interval [0.995, 1.005] holds 1: implied p reaches 1
    ("1.00", "0.90", "1.11", "0.5", DISAGREE),        # but no lower than 0.92
    ("1.01", "1.00", "1.01", "1e-300", AGREE),        # L_hi = U_lo = 1.005 = OR_lo: SE 0 at OR 1.005, implied p reaches 0
    ("1.01", "1.00", "1.01", "0.9", DISAGREE),        # and no higher than 0.33
    ("1.2", "1.1", "1.3", "1", NOT_CHECKABLE), ("1.2", "1.1", "1.3", "1.0", NOT_CHECKABLE),
    ("1.2", "1.1", "1.3", "0", NOT_CHECKABLE), ("1.2", "1.1", "1.3", "0.0", NOT_CHECKABLE),
    ("1.2", "1.1", "1.3", "0.99999", DISAGREE),       # checkable just below 1: the implied p is at most 0.03
    ("1.3", "1.3", "1.3", "0.01", NOT_CHECKABLE),     # L not below U
])
def test_boundary_rows(odds, lower, upper, p, expected):
    assert verdict(odds, lower, upper, p) == expected


@pytest.mark.parametrize("decimals", [2, 3, 4, 6, 8, 12])
def test_p_just_above_0_and_just_below_1_at_each_printed_precision(decimals):
    se = 0.05
    for p_true in (10.0 ** -decimals, 3 * 10.0 ** -decimals, 1 - 10.0 ** -decimals):
        texts = exact_texts(NormalDist().inv_cdf(1 - p_true / 2) * se, se)
        printed = f"{p_true:.{decimals}f}"
        assert Decimal(printed) == Decimal(printed).quantize(Decimal(1).scaleb(-decimals)) and 0 < Decimal(printed) < 1
        assert verdict(*texts, printed) == AGREE, printed
        if p_true < 0.5:
            assert verdict(*texts, f"{10 * p_true:.{decimals}f}") == DISAGREE
        else:
            assert verdict(*texts, f"{0.5:.{decimals}f}") == DISAGREE
        assert verdict(*texts, f"{1:.{decimals}f}") == verdict(*texts, f"{0:.{decimals}f}") == NOT_CHECKABLE
        neg = f"{-math.log10(p_true):.{decimals}f}"
        if Decimal(neg) > 0:
            assert verdict(*texts, neg, neg_log10=True) == AGREE, neg


@pytest.mark.parametrize("text,expected", [
    ("1.002", (1002, -3)), ("0.9981", (9981, -4)), ("1e-05", (1, -5)), ("1E-05", (1, -5)), ("1.0e-05", (10, -6)),
    ("1.2E-300", (12, -301)), (" 2.50e-07 ", (250, -9)), ("3", (3, 0)), ("0.000", (0, -3)), ("+0.5", (5, -1)),
    ("+1.25E+02", (125, 0)), ("1.4e-5", (14, -6)), ("1.4E-0005", (14, -6)),
])
def test_printed_digits_read_the_last_printed_digit_of_plain_and_scientific_forms(text, expected):
    got = printed_digits(text)
    assert got is not None and got[1:] == expected


@pytest.mark.parametrize("text", ["#NA", "NA", "", "x", "inf", "nan", "-0.1", "-0", "-0.0", "-0e-5", "Infinity"])
def test_printed_digits_give_nothing_for_a_missing_negative_signed_zero_or_non_finite_field(text):
    assert printed_digits(text) is None


def test_signs_exponent_case_and_long_mantissas_do_not_change_a_verdict():
    se = 0.05
    beta = NormalDist().inv_cdf(1 - 1.4e-5 / 2) * se
    texts = exact_texts(beta, se)
    for p in ("1.4e-05", "1.4E-05", "1.4e-5", "1.4E-0005", "+1.4e-05", "0.000014"):
        assert verdict(*texts, p) == AGREE, p
    plus = tuple("+" + t for t in texts)
    assert verdict(*plus, "+1.4e-05") == AGREE and verdict(*plus, "1.4e-04") == DISAGREE
    for p in ("-0", "-0.0", "+0", "+0.0", "0e-5"):
        assert verdict(*texts, p) == NOT_CHECKABLE, p
        assert verdict(*texts, p, neg_log10=True) == NOT_CHECKABLE, p
    # 60-digit mantissas: OR and limits to 40 significant digits, p to 60
    long_texts = tuple(f"{Decimal(t):.40f}" for t in texts)
    long_p = f"{Decimal(10) ** Decimal(log10_two_sided_p(beta / se)):.60e}"
    assert len(long_p) > 60 and verdict(*long_texts, long_p) == AGREE
    assert verdict(*long_texts, f"{Decimal(long_p) * 10:.60e}") == DISAGREE
    assert printed_digits("1." + "3" * 400)[1] == int("1" + "3" * 400)
    huge = printed_digits("1." + "3" * 5000)                   # past Python's int-from-str digit limit
    assert huge[2] == -5000 and huge[1] % 1000 == 333 and huge[1].bit_length() > 16000
    assert verdict(*texts, "1." + "3" * 5000) == NOT_CHECKABLE and verdict(*texts, "1." + "3" * 5000, neg_log10=True) == DISAGREE
    assert verdict("1." + "2" * 400, "1." + "1" * 400, "1." + "3" * 400, "0." + "0" * 300 + "1") in (AGREE, DISAGREE)


def test_a_very_large_minus_log10_p_is_checked_in_log_space():
    neg = 123456.789
    se = 0.01
    beta = z_of_log10_p(-neg) * se                                   # z about 754
    texts = exact_texts(beta, se)
    assert verdict(*texts, f"{neg:.3f}", neg_log10=True) == AGREE
    assert verdict(*texts, f"{neg + 1:.3f}", neg_log10=True) == DISAGREE
    for huge in ("1e400", "1E+999999", "9" * 500):                   # -log10 p beyond the doubles: p is 0
        assert verdict(*texts, huge, neg_log10=True) == DISAGREE
    assert verdict(*texts, "1e-99999") == DISAGREE
    assert verdict("1.01", "1.00", "1.01", "1e400", neg_log10=True) == AGREE     # SE 0 reachable at OR 1.005: z unbounded


def test_a_scientific_p_carries_the_rounding_of_its_last_mantissa_digit():
    se = 0.05
    texts = exact_texts(NormalDist().inv_cdf(1 - 1.4e-5 / 2) * se, se)  # implied p 1.4e-5 at full precision
    assert verdict(*texts, "1e-05") == AGREE                  # [0.5e-5, 1.5e-5] holds 1.4e-5
    assert verdict(*texts, "1E-05") == AGREE
    assert verdict(*texts, "1.0e-05") == DISAGREE             # [0.95e-5, 1.05e-5]: 0.125 log10 below
    assert verdict(*texts, "1.3E-05") == AGREE


def test_an_underflow_sized_p_is_checked_in_log_space():
    beta, se = 0.4, 0.01                                      # z = 40: p about 3.6e-350, below the doubles
    log10_p = log10_two_sided_p(beta / se)
    exponent = math.floor(log10_p)
    p_text = f"{10 ** (log10_p - exponent):.4f}e{exponent}"
    row = exact_or_ci_row(beta, se, p=0.0)
    row[7] = p_text
    got = run([row])
    assert counts(got) == (1, 1, 0, True)
    assert counts(got["frozen_rule"]) == (0, 0, 1, False)     # read as p = 0 by the reader: the frozen rule cannot check it
    row[7] = "3.6e-300"
    assert counts(run([row])) == (1, 0, 0, False)
    neg = row[:7] + [f"{-log10_p:.4f}"] + row[8:]
    assert counts(run([neg], header=NEG_LOG10_OR_CI_HEADER)) == (1, 1, 0, True)
    neg[7] = "300.0"
    assert counts(run([neg], header=NEG_LOG10_OR_CI_HEADER)) == (1, 0, 0, False)


def test_rows_without_a_p_in_0_1_are_not_checkable_under_either_rule_and_no_checkable_row_fails_the_file():
    got = run([or_ci_row(p="#NA"), or_ci_row(p="0"), or_ci_row(p="1"), or_ci_row(p="NA")])
    assert counts(got) == (0, 0, 4, False) and got["rows_lower_limit_rounding_reaches_zero"] == 0
    assert counts(got["frozen_rule"]) == (0, 0, 4, False)
    assert verdict("1.2", "1.1", "1.3", "#NA") == NOT_CHECKABLE
    assert verdict("1.2", "0", "1.3", "0.01") == LOWER_LIMIT_AT_ZERO          # the reader rejects such a row before this
    assert verdict("1.2", "0.000", "1.3", "0.01") == LOWER_LIMIT_AT_ZERO
    assert verdict("0", "0.1", "1.3", "0.01") == NOT_CHECKABLE


# ---- the file-level result -------------------------------------------------------------------------------------

def test_every_row_printed_from_true_values_agrees_or_is_uninformative_and_the_frozen_rule_fails_the_file():
    for p_format in (repr, lambda p: f"{p:.3g}", lambda p: f"{p:.1e}"):
        got = run(small_effect_rows(2000, p_format=p_format))
        frozen = got["frozen_rule"]
        assert got["rows_checked"] == got["rows_agree"] >= 1500
        # a p that prints as 1 (three significant digits of 0.9996) is not checkable under either rule
        assert got["rows_checked"] + got["rows_rounding_uninformative"] == frozen["rows_checked"]
        assert got["rows_not_checkable"] == frozen["rows_not_checkable"] + got["rows_rounding_uninformative"]
        assert got["passed"] is True and frozen["passed"] is False


@pytest.mark.parametrize("mismatch", ["ci_99_read_as_95", "p_of_another_test"])
def test_a_genuine_semantic_mismatch_fails_both_rules(mismatch):
    rng = np.random.default_rng(SEED)
    rows = []
    for zt, s, sign in zip(rng.uniform(3, 8, 300), rng.uniform(0.02, 0.05, 300), rng.choice([-1, 1], 300)):
        beta = sign * zt * s
        if mismatch == "ci_99_read_as_95":           # 99% limits; p the Wald p of (beta, se)
            rows.append(printed_row(beta, s, decimals=4, z=Z99))
        else:                                         # 95% limits; p three times the Wald p, 0.48 log10 off
            rows.append(printed_row(beta, s, decimals=4, p=repr(3 * wald_p(beta, s))))
    got = run(rows)
    assert counts(got) == (300, 0, 0, False)
    assert counts(got["frozen_rule"]) == (300, 0, 0, False)


def test_the_rounding_aware_result_reports_counts_and_pass_fail_only_with_the_frozen_rule_beside_it():
    rows = [exact_or_ci_row(0.1 + 0.001 * i, 0.05) for i in range(50)]
    got = run(rows)
    assert set(got) == {"rule", "rule_version", "ci_level", "z", "log10_tolerance", "min_agreement_percent", "rows_checked",
                        "rows_agree", "rows_disagree", "rows_not_checkable", "rows_lower_limit_rounding_reaches_zero",
                        "rows_rounding_uninformative", "passed", "frozen_rule", "row_crosstab", "diagnostic"}
    assert (got["rule"], got["rule_version"]) == ("or_ci_vs_p_rounding_aware", CI_P_RULE_ROUNDING_AWARE)
    assert got["frozen_rule"] == run(rows, CI_P_RULE_FROZEN)
    assert all(isinstance(v, int) for v in got["row_crosstab"].values()) and sum(got["row_crosstab"].values()) == 50
    with pytest.raises(ValueError, match="CI-versus-p rule is one of"):
        CiPCheck(SsfReader(SSF_OR_CI_HEADER, 0, "X", "or_ci_derived_se"), "loose")


def test_the_per_file_diagnostic_cross_tabulates_each_row_by_its_two_verdicts():
    small = small_effect_rows(400)
    frozen, alone = run(small, CI_P_RULE_FROZEN), run(small)
    # rows printed from latent Wald rows: no rounding-aware disagreement, so every frozen disagreement is
    # either accepted or made uninformative
    assert alone["diagnostic"]["disagree_under_both"] == 0 and alone["diagnostic"]["frozen_agree_rounding_aware_disagree"] == 0
    rescued = alone["diagnostic"]["frozen_disagree_rounding_aware_agree"]
    assert rescued + alone["row_crosstab"]["disagree__rounding_uninformative"] == frozen["rows_disagree"] and rescued > 0
    got = run(small + mismatch_rows(30) + [uninformative_row() for _ in range(11)] + [or_ci_row(p="#NA")])
    assert got["diagnostic"] == {"frozen_disagree_rounding_aware_agree": rescued, "disagree_under_both": 30,
                                 "rounding_uninformative": 11 + alone["rows_rounding_uninformative"],
                                 "frozen_agree_rounding_aware_disagree": 0}
    assert sum(got["row_crosstab"].values()) == 400 + 30 + 11 + 1
    assert got["row_crosstab"]["not_checkable__not_checkable"] == 1 + alone["row_crosstab"]["not_checkable__not_checkable"]


def test_the_diagnostic_reads_its_four_counts_from_the_cross_tabulation():
    crosstab = {"disagree__agree": 5, "disagree__disagree": 3, "agree__rounding_uninformative": 2,
                "disagree__rounding_uninformative": 4, "not_checkable__rounding_uninformative": 1, "agree__disagree": 6,
                "agree__agree": 100, "not_checkable__not_checkable": 9}
    assert ci_p_rule_diagnostic(crosstab) == {"frozen_disagree_rounding_aware_agree": 5, "disagree_under_both": 3,
                                              "rounding_uninformative": 7, "frozen_agree_rounding_aware_disagree": 6}


# ---- properties over seeded random rows ------------------------------------------------------------------------

def _box(text: str) -> tuple[float, float]:
    _, d, e = printed_digits(text)
    return math.log((d - 0.5) * 10.0 ** e), math.log((d + 0.5) * 10.0 ** e)


def oracle_extremes(texts: tuple[str, str, str], fill: int = 14) -> tuple[float, float] | None:
    """Brute force: |o| / (u - l) over a grid of each rounding box in log space with l <= o <= u. The
    grid is one shared set of values (every box's endpoints, 0, and `fill` points across each box) so
    that coordinates can coincide; the true extremes have coordinates among those values. None when no
    grid point is feasible."""
    o_box, l_box, u_box = (_box(t) for t in texts)
    shared = np.unique(np.concatenate([np.linspace(*b, fill) for b in (o_box, l_box, u_box)] + [[0.0]]))
    grid = [shared[(shared >= b[0]) & (shared <= b[1])] for b in (o_box, l_box, u_box)]
    o, lw, up = np.meshgrid(*grid, indexing="ij")
    feasible = (lw <= o) & (o <= up)
    if not feasible.any():
        return None
    width, size = (up - lw)[feasible], np.abs(o)[feasible]
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(width > 0, size / width, np.where(size > 0, np.inf, 0.0))
    return float(ratio.min()), float(ratio.max())


def random_printed_row(rng) -> tuple[tuple[str, str, str], float, float]:
    """OR and its 95% limits of a latent Wald row, printed to 1-5 decimals; (texts, beta, se)."""
    beta = rng.normal(0, 0.05) if rng.random() < 0.7 else rng.normal(0, 0.5)
    se = math.exp(rng.uniform(math.log(2e-4), math.log(0.3)))
    d = int(rng.integers(1, 6))
    vals = (math.exp(beta), math.exp(beta - CI_Z * se), math.exp(beta + CI_Z * se))
    return tuple(f"{v:.{d}f}" for v in vals), beta, se


def _considered(texts) -> bool:
    parsed = [printed_digits(t) for t in texts]
    return parsed[1][1] > 0 and parsed[0][0] > 0 and parsed[1][0] < parsed[2][0]


def test_the_exact_extremes_equal_a_brute_force_grid_oracle_over_random_rounded_rows():
    rng = np.random.default_rng(SEED)
    checked = infinite = uninformative = 0
    for _ in range(400):
        texts, _, _ = random_printed_row(rng)
        if not _considered(texts):
            continue
        exact, oracle = latent_wald_extremes(texts), oracle_extremes(texts)
        assert oracle is not None and exact is not None, texts
        checked += 1
        assert exact.low == pytest.approx(oracle[0], rel=1e-9, abs=1e-15), texts
        if math.isinf(exact.high):
            infinite += 1
            assert math.isinf(oracle[1]), texts
        else:
            assert exact.high == pytest.approx(oracle[1], rel=1e-9), texts
        uninformative += exact.zero_width_reachable and exact.contains_zero_effect
    assert checked > 300 and infinite > 10 and uninformative > 5


def near_limit_row(rng) -> tuple[str, str, str]:
    """OR, L and U drawn apart and printed to 1-3 decimals, the OR near one of its limits, where the
    constraint L <= OR <= U binds."""
    d = int(rng.integers(1, 4))
    lower = rng.uniform(0.3, 3.0)
    upper = lower + rng.uniform(0.0, 0.3)
    odds = (lower if rng.random() < 0.5 else upper) + rng.normal(0, 0.5 * 10.0 ** -d)
    return tuple(f"{v:.{d}f}" for v in (odds, lower, upper))


def test_the_exact_extremes_lie_inside_the_independent_corner_extremes():
    rng = np.random.default_rng(SEED + 1)
    tighter = {"latent": 0, "near_limit": 0}
    for i in range(4000):
        kind = "latent" if i % 2 else "near_limit"
        texts = random_printed_row(rng)[0] if kind == "latent" else near_limit_row(rng)
        if not _considered(texts):
            continue
        exact = latent_wald_extremes(texts)
        if exact is None:
            continue
        (o0, o1), (l0, l1), (u0, u1) = (_box(t) for t in texts)
        corner_low = (0.0 if o0 <= 0 <= o1 else min(abs(o0), abs(o1))) / (u1 - l0)
        corner_high = max(abs(o0), abs(o1)) / (u0 - l1) if u0 > l1 else math.inf
        # 1e-9: the oracle's logs of (D -/+ 1/2) 10^e and the rule's ln(D -/+ 1/2) + e ln 10 differ in the last bits,
        # which a width of a few 1e-4 magnifies
        assert exact.low >= corner_low * (1 - 1e-9) and exact.high <= corner_high * (1 + 1e-9), texts
        tighter[kind] += exact.high < corner_high * (1 - 1e-6) or exact.low > corner_low * (1 + 1e-6)
    assert tighter["near_limit"] > 100


def test_random_infeasible_and_feasible_boxes_agree_with_the_oracle_on_feasibility():
    rng = np.random.default_rng(SEED + 2)
    infeasible = 0
    for _ in range(600):
        d = int(rng.integers(1, 4))
        lo = rng.uniform(0.2, 3.0)
        vals = sorted([lo, lo + rng.uniform(0.0, 0.5)])
        odds = rng.uniform(0.1, 4.0)
        texts = (f"{odds:.{d}f}", f"{vals[0]:.{d}f}", f"{vals[1]:.{d}f}")
        if not _considered(texts):
            continue
        exact, oracle = latent_wald_extremes(texts), oracle_extremes(texts)
        assert (exact is None) == (oracle is None), texts
        if exact is None:
            infeasible += 1
            assert verdict(*texts, "0.01") == DISAGREE
        else:
            assert exact.low == pytest.approx(oracle[0], rel=1e-9, abs=1e-15)
            assert math.isinf(exact.high) == math.isinf(oracle[1])
    assert infeasible > 50


def test_a_row_printed_from_a_latent_wald_row_never_disagrees():
    rng = np.random.default_rng(SEED + 3)
    seen = {AGREE: 0, ROUNDING_UNINFORMATIVE: 0}
    for _ in range(3000):
        texts, beta, se = random_printed_row(rng)
        if not _considered(texts):
            continue
        log10_p = log10_two_sided_p(beta / se)
        if log10_p >= -1e-12:
            continue
        sig = int(rng.integers(1, 6))
        mantissa_exp = math.floor(log10_p)
        p_text = f"{10 ** (log10_p - mantissa_exp):.{sig - 1}f}e{mantissa_exp}"
        if Decimal(p_text) >= 1:
            continue
        got = verdict(*texts, p_text)
        assert got in seen, (texts, p_text, beta, se)
        seen[got] += 1
        neg_text = f"{-log10_p:.{sig}f}"
        if Decimal(neg_text) > 0:
            assert verdict(*texts, neg_text, neg_log10=True) in seen, (texts, neg_text)
    assert seen[AGREE] > 2000 and seen[ROUNDING_UNINFORMATIVE] > 10


def test_the_verdict_matches_the_one_the_oracle_extremes_give_over_random_rows_and_random_p():
    rng = np.random.default_rng(SEED + 4)
    seen = set()
    for _ in range(400):
        texts, beta, se = random_printed_row(rng)
        if not _considered(texts):
            continue
        p_text = f"{10 ** (log10_two_sided_p(beta / se) + rng.uniform(-1.5, 0.3)):.2e}"
        if not 0 < Decimal(p_text) < 1:
            continue
        low, high = oracle_extremes(texts)
        _, d, e = printed_digits(p_text)
        published = (math.log10(d - 0.5) + e, math.log10(d + 0.5) + e)
        if math.isinf(high) and low == 0.0:
            expected = ROUNDING_UNINFORMATIVE
        else:
            implied = (log10_two_sided_p(high * 2 * CI_Z), log10_two_sided_p(low * 2 * CI_Z))
            gaps = (published[0] - implied[1] - CI_P_LOG10_TOLERANCE, implied[0] - CI_P_LOG10_TOLERANCE - published[1])
            if min(abs(g) for g in gaps) < 1e-9:      # on the boundary to within the oracle's float error
                continue
            expected = AGREE if max(gaps) <= 0 else DISAGREE
        got = verdict(*texts, p_text)
        assert got == expected, (texts, p_text)
        seen.add(got)
    assert seen == {AGREE, DISAGREE, ROUNDING_UNINFORMATIVE}


def test_more_printed_precision_never_widens_the_implied_p_interval():
    """Wherever every rounding interval of the finer print lies inside the coarser one, the finer
    print's range of |ln OR| / SE lies inside the coarser one's, and so does its implied-p interval."""
    rng = np.random.default_rng(SEED + 5)
    nested = 0
    for _ in range(1500):
        beta, se = rng.normal(0, 0.05), math.exp(rng.uniform(math.log(5e-4), math.log(0.1)))
        vals = (math.exp(beta), math.exp(beta - CI_Z * se), math.exp(beta + CI_Z * se))
        coarse_d = int(rng.integers(1, 5))
        fine_d = coarse_d + int(rng.integers(1, 4))
        coarse, fine = (tuple(f"{v:.{d}f}" for v in vals) for d in (coarse_d, fine_d))
        if not (_considered(coarse) and _considered(fine)):
            continue
        if not all(c[0] <= f[0] and f[1] <= c[1] for c, f in zip(map(_box, coarse), map(_box, fine))):
            continue
        nested += 1
        wide, narrow = latent_wald_extremes(coarse), latent_wald_extremes(fine)
        assert narrow.low >= wide.low * (1 - 1e-12) and narrow.high <= wide.high * (1 + 1e-12), (coarse, fine)
        p_true = log10_two_sided_p(beta / se)
        if p_true < -1e-9:
            p_text = repr(10 ** p_true)
            if verdict(*fine, p_text) == AGREE:
                assert verdict(*coarse, p_text) in (AGREE, ROUNDING_UNINFORMATIVE), (coarse, fine, p_text)
    assert nested > 300


def test_the_ratio_extremes_on_hand_boxes():
    # no constraint binds: o in [0.1, 0.2], l in [-0.5, -0.4], u in [0.6, 0.7]
    got = wald_ratio_extremes((0.1, 0.2), (-0.5, -0.4), (0.6, 0.7))
    assert (got.low, got.high) == (pytest.approx(0.1 / 1.2), pytest.approx(0.2 / 1.0))
    # o forced to l: o in [0.3, 0.5], l in [0.45, 0.6], u in [0.7, 0.9]; feasible o in [0.45, 0.5]
    got = wald_ratio_extremes((0.3, 0.5), (0.45, 0.6), (0.7, 0.9))
    assert got.low == pytest.approx(0.45 / 0.45) and got.high == pytest.approx(0.5 / 0.2)
    assert wald_ratio_extremes((0.3, 0.4), (0.45, 0.6), (0.7, 0.9)) is None
    # a zero width at o = 0 only: the supremum is the ratio near it, not infinite
    got = wald_ratio_extremes((-0.1, 0.1), (-0.2, 0.0), (0.0, 0.2))
    assert not got.zero_width_reachable and got.contains_zero_effect and math.isfinite(got.high)
