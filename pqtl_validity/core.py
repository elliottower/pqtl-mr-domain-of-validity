"""Core validity test: classify MR results by mechanism and score informativeness."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from scipy.stats import binomtest

from pqtl_validity.taxonomy import Mechanism, classify_gene


@dataclass(frozen=True)
class ValidityResult:
    """Result of classifying one drug-target pair."""

    gene: str
    mechanism: Optional[Mechanism]
    mr_p: Optional[float]
    mr_significant: bool
    prediction: str  # "INFORMATIVE" or "UNINFORMATIVE"
    outcome: Optional[str]  # "SUCCESS" or "FAILURE" if drug approval known

    def __repr__(self) -> str:
        mech = self.mechanism.value if self.mechanism else "unknown"
        return (
            f"ValidityResult(gene={self.gene!r}, mechanism={mech!r}, "
            f"prediction={self.prediction!r})"
        )

    @property
    def mr_informative_expected(self) -> Optional[bool]:
        if self.mechanism is None:
            return None
        return self.mechanism.mr_informative


def classify_pair(
    gene: str,
    mr_p: Optional[float] = None,
    alpha: float = 0.05,
) -> ValidityResult:
    """Classify a drug-target pair and predict MR informativeness.

    Parameters
    ----------
    gene : str
        Gene symbol of the drug target (e.g. "ACE", "TNF", "EGFR").
    mr_p : float, optional
        P-value from cis-pQTL MR analysis, if available.
    alpha : float
        Significance threshold for the MR test.

    Returns
    -------
    ValidityResult
        Classification with prediction about MR informativeness.
    """
    mechanism = classify_gene(gene)
    mr_significant = mr_p is not None and mr_p < alpha

    if mechanism is None:
        prediction = "UNKNOWN"
    elif mechanism.mr_informative:
        prediction = "INFORMATIVE"
    else:
        prediction = "UNINFORMATIVE"

    return ValidityResult(
        gene=gene,
        mechanism=mechanism,
        mr_p=mr_p,
        mr_significant=mr_significant,
        prediction=prediction,
        outcome=None,
    )


@dataclass(frozen=True)
class ScoringResult:
    """Balanced accuracy evaluation across a set of pairs."""

    n: int
    n_success: int
    n_failure: int
    tp: int
    tn: int
    fp: int
    fn: int
    sensitivity: float
    specificity: float
    balanced_accuracy: float
    ba_ci_low: float
    ba_ci_high: float
    ci_level: float
    n_dropped: int = 0
    dropped_unclassified: tuple[str, ...] = ()
    dropped_mixed: tuple[str, ...] = ()
    dropped_no_outcome: tuple[str, ...] = ()

    def __repr__(self) -> str:
        dropped = f", dropped={self.n_dropped}" if self.n_dropped else ""
        return (
            f"ScoringResult(BA={self.balanced_accuracy:.3f} "
            f"[{self.ba_ci_low:.3f}, {self.ba_ci_high:.3f}], "
            f"n={self.n}{dropped})"
        )


def _coerce_mechanism(value) -> Optional[Mechanism]:
    """Accept a Mechanism, a taxonomy value, or a coarse stratum label."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    if isinstance(value, Mechanism):
        return value

    key = str(value).strip().lower()
    try:
        return Mechanism(key)
    except ValueError:
        pass
    # Enum member names, e.g. "A1_SOLUBLE_LIGAND" or just "A1".
    for m in Mechanism:
        if key == m.name.lower() or key == m.name.split("_")[0].lower():
            return m
    # Coarse three-class labels used in the published data tables. These carry
    # the informative/uninformative distinction without the subcategory, which
    # is all scoring needs.
    coarse = {
        "abundance_modulating": Mechanism.A3_INDIRECT_ABUNDANCE,
        "abundance": Mechanism.A3_INDIRECT_ABUNDANCE,
        "activity_blocking": Mechanism.B1_ENZYME_INHIBITOR,
        "activity": Mechanism.B1_ENZYME_INHIBITOR,
        "mixed": Mechanism.MIXED,
    }
    if key in coarse:
        return coarse[key]
    raise ValueError(
        f"Unrecognized mechanism {value!r}. Pass a Mechanism, a taxonomy value "
        f"such as 'A1_soluble_ligand', or one of {sorted(coarse)}."
    )


def score_pairs(
    genes: Sequence[str],
    outcomes: Sequence[str],
    mr_ps: Optional[Sequence[float]] = None,
    mechanisms: Optional[Sequence] = None,
    alpha: float = 0.05,
    n_bootstrap: int = 10_000,
    ci_level: float = 0.90,
    seed: Optional[int] = None,
    rule: str = "mr_significance",
    on_unclassified: str = "error",
) -> ScoringResult:
    """Score a set of drug-target pairs.

    Two decision rules are available.

    ``rule="mr_significance"`` (default) is the screen the paper evaluates:
    predict SUCCESS when MR is significant (p < alpha), whatever the mechanism.
    This is the standard pQTL-MR screen under test. Applied to all pairs it
    reproduces the composite result; applied to a mechanism stratum it
    reproduces that stratum's result. Mechanism assignments are not required,
    because stratification is performed by passing a subset.

    ``rule="mechanism_gated"`` additionally requires the mechanism to be
    abundance-mediated (A1/A2/A3), predicting FAILURE for activity-blocking
    mechanisms (B1/B2) on the grounds that the null is expected there. This
    rule is not evaluated in the paper and requires a mechanism per pair.

    Parameters
    ----------
    genes : sequence of str
        Gene symbols for each pair.
    outcomes : sequence of str
        "SUCCESS" or "FAILURE" for each pair (drug approval status).
    mr_ps : sequence of float, optional
        MR p-values. If None, classification uses mechanism alone.
    mechanisms : sequence, optional
        Per-pair mechanism assignments, as ``Mechanism`` members or strings.
        Used only by ``rule="mechanism_gated"``. When omitted, mechanisms are
        looked up from the frozen taxonomy by gene symbol, which covers fewer
        genes than the curated data tables do; pass the curated column to
        avoid a reduced denominator.
    alpha : float
        MR significance threshold.
    n_bootstrap : int
        Bootstrap iterations for BA confidence interval.
    ci_level : float
        Confidence interval level (default 0.90 matching paper).
    seed : int, optional
        Random seed for bootstrap.
    on_unclassified : {"error", "skip"}
        Behaviour when a pair has no mechanism assignment. ``"error"`` (the
        default) raises and names the offending genes. ``"skip"`` excludes
        them, recording which in ``ScoringResult.dropped_unclassified``.

    Returns
    -------
    ScoringResult
        Balanced accuracy with bootstrap confidence interval. Any pair excluded
        from scoring is counted in ``n_dropped`` and listed by reason, so a
        reduced denominator is always visible in the result.

    Raises
    ------
    ValueError
        If input lengths disagree, if ``on_unclassified="error"`` and some gene
        has no mechanism assignment, or if no pair survives filtering.
    """
    if len(genes) != len(outcomes):
        raise ValueError("genes and outcomes must have same length")
    if mr_ps is not None and len(mr_ps) != len(genes):
        raise ValueError("mr_ps must have same length as genes")
    if mechanisms is not None and len(mechanisms) != len(genes):
        raise ValueError("mechanisms must have same length as genes")
    if on_unclassified not in ("error", "skip"):
        raise ValueError("on_unclassified must be 'error' or 'skip'")
    if rule not in ("mr_significance", "mechanism_gated"):
        raise ValueError("rule must be 'mr_significance' or 'mechanism_gated'")
    if rule == "mr_significance" and mr_ps is None:
        raise ValueError("rule='mr_significance' requires mr_ps")

    rng = np.random.default_rng(seed)
    gated = rule == "mechanism_gated"

    if gated:
        resolved = [
            _coerce_mechanism(mechanisms[i]) if mechanisms is not None
            else classify_gene(gene)
            for i, gene in enumerate(genes)
        ]
        unclassified = tuple(
            sorted({str(g) for g, m in zip(genes, resolved) if m is None})
        )
        if unclassified and on_unclassified == "error":
            raise ValueError(
                f"{len(unclassified)} gene(s) have no mechanism assignment: "
                f"{', '.join(unclassified)}. Pass `mechanisms=` with the curated "
                f"assignments, or `on_unclassified='skip'` to exclude them, which "
                f"reduces the denominator."
            )
    else:
        resolved = [None] * len(genes)
        unclassified = ()

    predictions = []
    dropped_mixed, dropped_no_outcome = [], []
    for i, gene in enumerate(genes):
        if outcomes[i] not in ("SUCCESS", "FAILURE"):
            # e.g. EXCLUDED or PENDING adjudications.
            dropped_no_outcome.append(str(gene))
            continue

        if not gated:
            predictions.append(
                ("SUCCESS" if mr_ps[i] < alpha else "FAILURE", outcomes[i])
            )
            continue

        mech = resolved[i]
        if mech is None:
            continue
        if mech == Mechanism.MIXED:
            # mr_informative is undefined for MIXED, so the pair is unscorable
            # under this rule.
            dropped_mixed.append(str(gene))
            continue

        if mech.mr_informative:
            significant = mr_ps is None or mr_ps[i] < alpha
            predictions.append(("SUCCESS" if significant else "FAILURE", outcomes[i]))
        else:
            predictions.append(("FAILURE", outcomes[i]))

    if not predictions:
        raise ValueError("No pairs could be scored after filtering.")

    tp = sum(1 for p, o in predictions if p == "SUCCESS" and o == "SUCCESS")
    tn = sum(1 for p, o in predictions if p == "FAILURE" and o == "FAILURE")
    fp = sum(1 for p, o in predictions if p == "SUCCESS" and o == "FAILURE")
    fn = sum(1 for p, o in predictions if p == "FAILURE" and o == "SUCCESS")

    n_s = tp + fn
    n_f = tn + fp
    sens = tp / n_s if n_s > 0 else 0.0
    spec = tn / n_f if n_f > 0 else 0.0
    ba = (sens + spec) / 2.0

    # Bootstrap CI for balanced accuracy
    ba_boots = np.empty(n_bootstrap)
    pred_arr = np.array(predictions, dtype=object)
    for b in range(n_bootstrap):
        idx = rng.integers(0, len(predictions), size=len(predictions))
        sample = pred_arr[idx]
        b_tp = sum(1 for p, o in sample if p == "SUCCESS" and o == "SUCCESS")
        b_fn = sum(1 for p, o in sample if p == "FAILURE" and o == "SUCCESS")
        b_tn = sum(1 for p, o in sample if p == "FAILURE" and o == "FAILURE")
        b_fp = sum(1 for p, o in sample if p == "SUCCESS" and o == "FAILURE")
        b_ns = b_tp + b_fn
        b_nf = b_tn + b_fp
        b_sens = b_tp / b_ns if b_ns > 0 else 0.0
        b_spec = b_tn / b_nf if b_nf > 0 else 0.0
        ba_boots[b] = (b_sens + b_spec) / 2.0

    tail = (1.0 - ci_level) / 2.0
    ci_low = float(np.percentile(ba_boots, tail * 100))
    ci_high = float(np.percentile(ba_boots, (1.0 - tail) * 100))

    return ScoringResult(
        n=len(predictions),
        n_success=n_s,
        n_failure=n_f,
        tp=tp,
        tn=tn,
        fp=fp,
        fn=fn,
        sensitivity=sens,
        specificity=spec,
        balanced_accuracy=ba,
        ba_ci_low=ci_low,
        ba_ci_high=ci_high,
        ci_level=ci_level,
        n_dropped=len(genes) - len(predictions),
        dropped_unclassified=unclassified,
        dropped_mixed=tuple(dropped_mixed),
        dropped_no_outcome=tuple(dropped_no_outcome),
    )
