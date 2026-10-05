"""Analysis sets S1-S22 (PREREG §Study design, Sets; §Other planned analysis).

Each set is a pre-outcome frame (membership fixed before outcomes matter) with four analysis
columns: `cls` (mechanism class used by H1), `state` (evidence state), `S` (1 if supportive),
`E` (continuous score) and `y` (outcome; NaN where the set's outcome rule excludes the row).
The analysed rows of a set are those with a non-missing `y`.

S1 is `in_s1 & heldout`. Assumptions about stage A flags where the plan leaves the encoding
open are stated beside each set and in the stage D report.

Flag restrictions keep only rows whose flag is known to be false: S5 keeps `protein_altering ==
False`, S11 `splicing_candidate == False`, S15f `low_coverage == False`. A row with the flag
missing is excluded from that set and counted in its `flag_missing_excluded` note.

S16 replaces the evidence state by `s16_evidence_state` where stage B wrote one (a deCODE instrument
colocalized with its SMP-normalized statistics). When no S1 hypothesis has one, because no
SMP-normalized file was retrieved, S16 is not formed and is reported as not run.

DEVIATION_SETS are not registered sets. S_frozen_ci_rule is the sensitivity analysis of the
post-freeze amendment of stage B's CI-versus-p validation rule (PRESTAGE_NOTES_DRAFT.md, stage B):
the held-out S1 rows, with every hypothesis whose outcome file fails the frozen rule
(`outcome_file_frozen_ci_p_pass` false) made inconclusive (S = 0, E = 0), as it would have been had
the frozen rule excluded the file. Every other row keeps its state. It is fitted with H1 and H4 at
the primary prior and reported apart from table 10. A file that fails the amended rule and passes
the frozen one is unavailable in both, because stage B ran no colocalization on it. Not formed when
stage B wrote no such column, or when no S1 hypothesis has a failing outcome file (the set would
equal S1).
"""
import numpy as np
import pandas as pd
from pydantic import BaseModel

from stage_d.evidence import PRIMARY_PP_H4_THRESHOLD, evidence_state


class SetError(ValueError):
    pass


class SetInfo(BaseModel):
    set_id: str
    definition: str
    models: tuple[str, ...]
    formed: bool
    reason: str = ""
    n_pre_outcome: int = 0
    n_analysed: int = 0
    excluded_outcome: dict[str, int] = {}
    notes: dict[str, int] = {}


class AnalysisSet(BaseModel):
    model_config = {"arbitrary_types_allowed": True}
    info: SetInfo
    frame: pd.DataFrame | None = None

    def analysed(self) -> pd.DataFrame:
        if self.frame is None:
            raise SetError(f"{self.info.set_id} was not formed: {self.info.reason}")
        return self.frame.loc[self.frame["y"].notna()].copy()


S15_VARIANTS = {
    "S15a": ("s15a_pp_h4", PRIMARY_PP_H4_THRESHOLD, "coloc p12 = 1e-6"),
    "S15b": ("s15b_pp_h4", PRIMARY_PP_H4_THRESHOLD, "coloc p12 = 1e-5"),
    "S15c": ("s15c_pp_h4", PRIMARY_PP_H4_THRESHOLD, "coloc +/-1 Mb window"),
    "S15d": ("pp_h4", 0.75, "PP.H4 threshold 0.75"),
    "S15e": ("pp_h4", 0.90, "PP.H4 threshold 0.90"),
    "S15g": ("s15g_pp_h4", PRIMARY_PP_H4_THRESHOLD, "coloc.susie where either trait has >= 2 credible sets"),
}
MATURATION = {"S14a": 12, "S14b": 36, "S14c": 48}


def s12_select(frame: pd.DataFrame) -> pd.Series:
    """One hypothesis per drug-gene program: within each gene, order by specificity_rank, then
    indication id, then lowest program id; keep a hypothesis only if none of its programs is in a
    hypothesis already kept for that gene."""
    keep = pd.Series(False, index=frame.index)
    order = frame.assign(_min_prog=frame["programs"].map(min)).sort_values(
        ["gene_ensembl", "specificity_rank", "indication_id", "_min_prog"], kind="mergesort")
    for _, g in order.groupby("gene_ensembl", sort=False):
        used: set[str] = set()
        for idx, progs in zip(g.index, g["programs"]):
            if not used.intersection(progs):
                keep[idx] = True
                used.update(progs)
    return keep


def _make(df: pd.DataFrame, mask, set_id: str, definition: str, models=("h1",), *, cls=None,
          state=None, y=None, status="status_24", notes=None) -> AnalysisSet:
    """`status` names the outcome status column whose values label the rows `y` leaves missing
    (the set's exclusion report: active, business_only, undated, no_phase, ...)."""
    frame = df.loc[np.asarray(mask, dtype=bool)].copy()
    frame["cls"] = frame["mechanism_class"] if cls is None else pd.Series(cls, index=df.index).loc[frame.index]
    st = frame["evidence_state"] if state is None else pd.Series(state, index=df.index).loc[frame.index]
    frame["state"] = st.astype(object)
    frame["S"] = (frame["state"] == "supportive").astype(float)
    frame["E"] = frame["E"].astype(float)
    frame["C"] = frame["neuro_psych"].astype(float)
    y_series = frame["advanced_24"] if y is None else pd.Series(y, index=df.index).loc[frame.index]
    frame["y"] = pd.to_numeric(y_series, errors="raise").astype(float)
    excluded = frame.loc[frame["y"].isna(), status].value_counts().to_dict()
    info = SetInfo(set_id=set_id, definition=definition, models=tuple(models), formed=True,
                   n_pre_outcome=len(frame), n_analysed=int(frame["y"].notna().sum()),
                   excluded_outcome={str(k): int(v) for k, v in excluded.items()}, notes=notes or {})
    return AnalysisSet(info=info, frame=frame)


def _not_formed(set_id, definition, models, reason) -> AnalysisSet:
    return AnalysisSet(info=SetInfo(set_id=set_id, definition=definition, models=tuple(models),
                                    formed=False, reason=reason))


def form_sets(df: pd.DataFrame) -> dict[str, AnalysisSet]:
    present = {c for cols in df.attrs.get("present_columns", {}).values() for c in cols}
    base = (df["in_s1"] & df["heldout"]).to_numpy()
    sets: dict[str, AnalysisSet] = {}

    sets["S1"] = _make(df, base, "S1", "held-out; advancement outcome; 24-month maturation")
    # S2: in_s1 regardless of heldout. INTERFACES.md (stage A): in_s1 marks passing every S1 rule
    # before the held-out filter, so S1 = in_s1 & heldout and the pooled set S2 = in_s1.
    pooled = df["in_s1"].to_numpy()
    sets["S2"] = _make(df, pooled, "S2", "pooled: held-out plus pilot keys under V8 rules",
                       notes={"pilot_rows": int((pooled & ~df["heldout"].to_numpy()).sum())})
    sets["S3"] = _make(df, base & (df["platform_concordant"] == "concordant").fillna(False).to_numpy(),
                       "S3", "S1 restricted to platform_concordant proteins")
    sets["S4"] = _make(df, (df["in_s4"] & df["heldout"]).to_numpy(), "S4", "S1 plus overlap_unknown hypotheses")
    pa = df["protein_altering"]
    sets["S5"] = _make(df, base & (pa == False).fillna(False).to_numpy(), "S5",  # noqa: E712
                       "S1 without protein_altering (instruments with the flag computed)",
                       notes={"flag_missing_excluded": int((base & pa.isna().to_numpy()).sum())})
    pre = df["pre_pqtl_publication"]
    sets["S6"] = _make(df, base & (pre == True).fillna(False).to_numpy(), "S6",  # noqa: E712
                       "S1 restricted to pre_pqtl_publication",
                       notes={"date_missing_excluded": int((base & pre.isna().to_numpy()).sum())})
    adv = df["advanced_24"].astype(float)
    y7 = adv.where((adv == 1.0) | df["efficacy_coded"].to_numpy(dtype=bool))
    sets["S7"] = _make(df, base, "S7", "S1 with no-observed-advancement restricted to efficacy-coded stops", y=y7)
    # S8 / S9 / S21: stage A re-forms hypotheses under each set's class table (S21: restored
    # conflicting rows); in_s8 / in_s9 / in_s21 mark the re-formed set (shared primary rows plus
    # variant-only rows) before the held-out filter.
    for sid, col, text in (("S8", "in_s8",
                            "S1 with antibody antagonists of non-blood-secreted targets moved from other to function-blocking"),
                           ("S9", "in_s9", "S1 with UniProt Secreted in place of HPA blood-secreted"),
                           ("S21", "in_s21", "S1 with conflicted program-target rows restored")):
        if col in present:
            mask = (df[col].fillna(False) & df["heldout"]).to_numpy(dtype=bool)
            sets[sid] = _make(df, mask, sid, text,
                              notes={"variant_only_rows": int((mask & (df["variant"] != "primary").to_numpy()).sum())})
        else:
            sets[sid] = _not_formed(sid, text, ("h1",), f"stage A output has no column {col}")
    sets["S10"] = _make(df, (df["in_s10"] & df["heldout"]).to_numpy(), "S10",
                        "S1 plus tier-2 and tier-3 outcome GWAS")
    spl = df["splicing_candidate"]
    computed = spl.notna().to_numpy()
    if not (base & computed).any():
        sets["S11"] = _not_formed("S11", "S1 without splicing_candidate", ("h1",),
                                  "splicing_candidate missing for every S1 instrument")
    else:
        sets["S11"] = _make(df, base & computed & ~(spl == True).fillna(False).to_numpy(), "S11",  # noqa: E712
                            "S1 without splicing_candidate (instruments with the flag computed)",
                            notes={"flag_missing_excluded": int((base & ~computed).sum())})
    s12 = pd.Series(False, index=df.index)
    s12.loc[df.index[base]] = s12_select(df.loc[base])
    a_s12 = (df["in_s12"] & base).to_numpy()
    if not np.array_equal(s12.to_numpy(), a_s12):
        raise SetError(f"S12: D selection and stage A in_s12 disagree on {int((s12.to_numpy() != a_s12).sum())} S1 rows")
    sets["S12"] = _make(df, s12.to_numpy(), "S12", "one hypothesis per drug-gene program")
    sets["S13"] = _make(df, base & ~df["karim_launched"].to_numpy() & ~df["pilot_indication"].to_numpy(), "S13",
                        "S1 without karim_launched and without pilot_indication")
    for sid, months in MATURATION.items():
        sets[sid] = _make(df, base, sid, f"S1 with {months}-month maturation", y=df[f"advanced_{months}"],
                          status=f"status_{months}")
    for sid, (col, thr, text) in S15_VARIANTS.items():
        st = evidence_state(df[col], df["genetic_direction"], df["direction"], threshold=thr)
        sets[sid] = _make(df, base, sid, text, state=st)
    low = df["low_coverage"]
    sets["S15f"] = _make(df, base & (low == False).fillna(False).to_numpy(), "S15f",  # noqa: E712
                         "low_coverage excluded (hypotheses with coverage computed)",
                         notes={"flag_missing_excluded": int((base & low.isna().to_numpy()).sum())})
    has_smp = base & df["s16_evidence_state"].notna().to_numpy()
    if not has_smp.any():
        sets["S16"] = _not_formed("S16", "deCODE SMP-normalized statistics where available", ("h1",),
                                  "no deCODE SMP-normalized statistics were retrieved (s16_evidence_state is empty "
                                  "for every S1 hypothesis); S16 is not run")
    else:
        s16 = df["s16_evidence_state"].fillna(df["evidence_state"])
        sets["S16"] = _make(df, base, "S16", "deCODE SMP-normalized statistics where available", state=s16,
                            notes={"replaced": int(has_smp.sum())})
    if "s17_sentinel_p" in present:
        p = df["s17_sentinel_p"]
        st17 = np.where(p.fillna(np.inf).to_numpy() < 0.05, "supportive", "inconclusive")
        sets["S17"] = _make(df, base, "S17", "evidence = outcome-association p < 0.05 at the sentinel, direction ignored",
                            state=st17, notes={"p_missing": int((base & p.isna().to_numpy()).sum())})
    else:
        sets["S17"] = _not_formed("S17", "evidence = outcome-association p < 0.05 at the sentinel", ("h1",),
                                  "stage B output has no column s17_sentinel_p")
    sets["S18"] = _make(df, base & df["single_protein_row"].to_numpy(), "S18",
                        "S1 restricted to single-protein mechanism rows")
    # S19: neutralizing biologics vs small-molecule blockers of blood-secreted proteins. The arm is
    # stage A's s19_arm; a row enters only when its arm matches its class.
    arm = df["s19_arm"].astype(object)
    arm_ok = (((df["mechanism_class"] == "aligned") & (arm == "neutralizing_biologic"))
              | ((df["mechanism_class"] == "blocking") & (arm == "small_molecule_blocker"))).to_numpy(dtype=bool)
    blood = base & df["blood_secreted_hpa"].to_numpy(dtype=bool)
    sets["S19"] = _make(df, blood & arm_ok, "S19",
                        "H1 on blood-secreted targets: neutralizing biologics vs small-molecule blockers",
                        notes={"blood_secreted_aligned_not_neutralizing_biologic": int(
                                   (blood & (df["mechanism_class"] == "aligned").to_numpy() & ~arm_ok).sum()),
                               "blood_secreted_blocking_not_small_molecule": int(
                                   (blood & (df["mechanism_class"] == "blocking").to_numpy() & ~arm_ok).sum())})
    # S20: A writes non-strict rows only where they pass every other S1 rule.
    sets["S20"] = _make(df, df["heldout"].to_numpy() & (df["in_s1"] | ~df["strict_instrument"]).to_numpy(), "S20",
                        "S1 with source-native inclusive instrument lists")
    sets["S22"] = _make(df, base & ~df["psych_only"].to_numpy(), "S22",
                        "H4 with psychiatric-only indications removed", models=("h4",))
    sets[FROZEN_CI_RULE_SET] = frozen_ci_rule_set(df, base, present)
    return sets


FROZEN_CI_RULE_SET = "S_frozen_ci_rule"
FROZEN_CI_RULE_COLUMN = "outcome_file_frozen_ci_p_pass"
FROZEN_CI_RULE_DEFINITION = ("deviation sensitivity: S1 with every hypothesis whose outcome file fails the frozen "
                             "CI-versus-p rule made inconclusive")


def frozen_ci_rule_set(df: pd.DataFrame, base: np.ndarray, present: set[str]) -> AnalysisSet:
    """S_frozen_ci_rule (module docstring): the S1 rows, the state of each row whose outcome file fails
    the frozen CI-versus-p rule set to inconclusive and its E to 0."""
    models = ("h1", "h4")
    if FROZEN_CI_RULE_COLUMN not in present:
        return _not_formed(FROZEN_CI_RULE_SET, FROZEN_CI_RULE_DEFINITION, models,
                           f"stage B output has no column {FROZEN_CI_RULE_COLUMN}")
    flag = df[FROZEN_CI_RULE_COLUMN]
    fails = base & (flag == False).fillna(False).to_numpy(dtype=bool)  # noqa: E712
    if not fails.any():
        return _not_formed(FROZEN_CI_RULE_SET, FROZEN_CI_RULE_DEFINITION, models,
                           "no held-out S1 hypothesis has an outcome file that fails the frozen CI-versus-p rule; "
                           "the set equals S1")
    state = df["evidence_state"].astype(object).where(~fails, "inconclusive")
    changed = fails & (df["evidence_state"] != "inconclusive").to_numpy(dtype=bool)
    notes = {"outcome_file_fails_frozen_rule": int(fails.sum()),
             "outcome_file_passes_frozen_rule": int((base & (flag == True).fillna(False).to_numpy(dtype=bool)).sum()),  # noqa: E712
             "no_ci_p_check": int((base & flag.isna().to_numpy()).sum()),
             "made_inconclusive": int(changed.sum()),
             **{f"made_inconclusive_from_{s}": int((fails & (df["evidence_state"] == s).to_numpy(dtype=bool)).sum())
                for s in ("supportive", "contradictory")}}
    return _make(df.assign(E=df["E"].where(~fails, 0.0)), base, FROZEN_CI_RULE_SET, FROZEN_CI_RULE_DEFINITION, models,
                 state=state, notes=notes)


SET_ORDER = ("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13",
             "S14a", "S14b", "S14c", "S15a", "S15b", "S15c", "S15d", "S15e", "S15f", "S15g",
             "S16", "S17", "S18", "S19", "S20", "S21", "S22")
# Sets that are not registered: the sensitivity analysis of a logged deviation (module docstring).
DEVIATION_SETS = (FROZEN_CI_RULE_SET,)
