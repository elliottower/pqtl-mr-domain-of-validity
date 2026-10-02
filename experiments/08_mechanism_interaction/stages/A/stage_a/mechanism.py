"""Mechanism class, intervention direction and conflicting mechanism rows
(PREREG §Measured variables, "Mechanism class", "Intervention direction";
§Study design, "Conflicting mechanism rows", "Drug programs").

Direction tokens follow `experiments/PREREG_direction_concordance.md` §3.2-§3.4, adopted
verbatim by the plan. Ported from `count_v4.py` (`token_direction`, `is_neutralizing`,
`effective_action`, `mechanism_class`, `direction_of`, `resolve_drug_targets`).

As in count_v4.py, the class table and the BINDING AGENT test read the action type exactly as Open
Targets gives it (lines 144-162); only the direction rule normalizes its token (line 134, §3.2).
"""
import pandas as pd
from pydantic import BaseModel, ConfigDict

from stage_a.schemas import Direction, MechanismClass, S19Arm

# Mechanism class table, first matching row wins.
BIO_TYPES = frozenset({"Antibody", "Protein", "Enzyme"})
BIO_ACTIONS = frozenset({"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR"})
OLIGO_TYPE = "Oligonucleotide"
OLIGO_ACTIONS = frozenset({"ANTISENSE INHIBITOR", "RNAI INHIBITOR"})
DEGRADER = "DEGRADER"
SMALL_MOLECULE = "Small molecule"
SM_ACTIONS = frozenset({"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR",
                        "INVERSE AGONIST"})
BINDING_AGENT = "BINDING AGENT"
NEUTRALIZING_TOKENS = ("inhibitor", "antagonist", "blocker", "neutraliz", "neutralis", "sequestr")

ROW1 = "row1_biologic_inhibitor_blood_secreted"
ROW2 = "row2_oligonucleotide_blood_secreted"
ROW3 = "row3_degrader_blood_secreted"
ROW4 = "row4_small_molecule_blocker"
OTHER_BIOLOGIC_NOT_BLOOD = "other:biologic_target_not_blood_secreted"
OTHER_OLIGO_NOT_BLOOD = "other:oligonucleotide_target_not_blood_secreted"
OTHER_DEGRADER_NOT_BLOOD = "other:degrader_target_not_blood_secreted"
OTHER_BINDING_AGENT = "other:binding_agent_not_neutralizing"
OTHER_NOT_IN_TABLE = "other:molecule_or_action_type_not_in_table"
OTHER_CONFLICT = "other:conflicting_mechanism_rows"
# S8 (PREREG §Other planned analysis): "antibody antagonists of targets that are not
# blood-secreted moved from other to function-blocking". Read as the rows row 1 sends to other
# for localization alone (molecule type Antibody; INHIBITOR, ANTAGONIST, NEGATIVE ALLOSTERIC
# MODULATOR or a neutralizing BINDING AGENT); Protein and Enzyme rows stay other.
S8_MOLECULE_TYPE = "Antibody"
S8_BLOCKING = "s8:antibody_antagonist_target_not_blood_secreted"

# S19 (PREREG §Other planned analysis): "neutralizing biologics vs small-molecule blockers of
# blood-secreted proteins". A neutralizing biologic is a row-1 row (Antibody, Protein or Enzyme;
# INHIBITOR, ANTAGONIST, NEGATIVE ALLOSTERIC MODULATOR or a neutralizing BINDING AGENT; target
# blood-secreted). A small-molecule blocker is a row-4 row whose target is blood-secreted. Row 2
# (oligonucleotides) and row 3 (degraders) are abundance-aligned but in neither arm.
S19_NEUTRALIZING = "neutralizing_biologic"
S19_SMALL_MOLECULE = "small_molecule_blocker"

# PREREG_direction_concordance.md §3.3.
TIER1_OVERRIDES = {"INVERSE AGONIST": "BLOCKING", "NEGATIVE ALLOSTERIC MODULATOR": "BLOCKING",
                   "POSITIVE ALLOSTERIC MODULATOR": "ACTIVATING", "PARTIAL AGONIST": "ACTIVATING"}
BLOCKING_SUBSTR = ("INHIBITOR", "ANTAGONIST", "BLOCKER", "DEGRADER", "DISRUPT", "SUPPRESSOR")
ACTIVATING_SUBSTR = ("AGONIST", "ACTIVATOR", "OPENER", "STABILISER", "STABILIZER", "RELEASING AGENT")


class RowLabel(BaseModel):
    model_config = ConfigDict(frozen=True)

    mechanism_class: MechanismClass
    direction: Direction
    reason: str
    neutralizing_binding_agent: bool


def normalize_token(token: str | None) -> str:
    """§3.2 step 2: strip, uppercase, collapse internal whitespace."""
    return " ".join(str(token if token is not None else "").upper().split())


def token_direction(token: str | None) -> str:
    """§3.3: BLOCKING, ACTIVATING or AMBIGUOUS, tiers evaluated in order."""
    t = normalize_token(token)
    if not t:
        return "AMBIGUOUS"
    if t in TIER1_OVERRIDES:
        return TIER1_OVERRIDES[t]
    if any(s in t for s in BLOCKING_SUBSTR):
        return "BLOCKING"
    if any(s in t for s in ACTIVATING_SUBSTR):
        return "ACTIVATING"
    return "AMBIGUOUS"


def is_neutralizing(action_type: str | None, moa_text: str | None) -> bool:
    """count_v4.py line 145: the action type compared as given, the mechanism text lower-cased."""
    return action_type == BINDING_AGENT and any(k in str(moa_text).lower() for k in NEUTRALIZING_TOKENS)


def effective_action(action_type: str | None, moa_text: str | None) -> str | None:
    """A neutralizing BINDING AGENT is treated as INHIBITOR for class and direction; any other
    action type is returned as given (count_v4.py line 149)."""
    return "INHIBITOR" if is_neutralizing(action_type, moa_text) else action_type


def direction_of(action_eff: str | None) -> Direction:
    """§3.4 applied to one mechanism row: a single token (count_v4.py lines 165-167)."""
    d = token_direction(action_eff)
    if d == "BLOCKING":
        return "decrease"
    if d == "ACTIVATING":
        return "increase"
    return "ambiguous"


def mechanism_class(molecule_type: str, action_eff: str | None, blood_secreted: bool) -> tuple[MechanismClass, str]:
    """(class, reason): the table row that fired, or why the row is other."""
    if molecule_type in BIO_TYPES and action_eff in BIO_ACTIONS:
        return ("aligned", ROW1) if blood_secreted else ("other", OTHER_BIOLOGIC_NOT_BLOOD)
    if molecule_type == OLIGO_TYPE and action_eff in OLIGO_ACTIONS:
        return ("aligned", ROW2) if blood_secreted else ("other", OTHER_OLIGO_NOT_BLOOD)
    if action_eff == DEGRADER and blood_secreted:
        return "aligned", ROW3
    if molecule_type == SMALL_MOLECULE and action_eff in SM_ACTIONS:
        return "blocking", ROW4
    if action_eff == DEGRADER:
        return "other", OTHER_DEGRADER_NOT_BLOOD
    if action_eff == BINDING_AGENT:
        return "other", OTHER_BINDING_AGENT
    return "other", OTHER_NOT_IN_TABLE


def label_row(molecule_type: str, action_type: str | None, moa_text: str | None,
              blood_secreted: bool) -> RowLabel:
    eff = effective_action(action_type, moa_text)
    cls, reason = mechanism_class(molecule_type, eff, blood_secreted)
    return RowLabel(mechanism_class=cls, direction=direction_of(eff), reason=reason,
                    neutralizing_binding_agent=is_neutralizing(action_type, moa_text))


def label_row_s8(molecule_type: str, action_type: str | None, moa_text: str | None,
                 blood_secreted: bool) -> RowLabel:
    """The primary label, except that an antibody row which is other only because its target is
    not blood-secreted is function-blocking (S8). Direction is unchanged."""
    base = label_row(molecule_type, action_type, moa_text, blood_secreted)
    if base.reason == OTHER_BIOLOGIC_NOT_BLOOD and molecule_type == S8_MOLECULE_TYPE:
        return base.model_copy(update={"mechanism_class": "blocking", "reason": S8_BLOCKING})
    return base


def s19_arm(reason: str, blood_secreted: bool) -> S19Arm:
    """The S19 arm of a program-target (or of a restored S21 row) from the class-table rows that
    decided its class (`;`-joined reasons): every reason row 1 -> neutralizing biologic; every
    reason row 4 and the target blood-secreted -> small-molecule blocker; otherwise neither."""
    reasons = set(str(reason).split(";"))
    if reasons == {ROW1}:
        return S19_NEUTRALIZING
    if reasons == {ROW4} and blood_secreted:
        return S19_SMALL_MOLECULE
    return ""


def hypothesis_s19_arm(program_arms: dict[str, set[str]]) -> S19Arm:
    """A hypothesis is in an arm only if every one of its programs is in that arm (a program whose
    rows fall in different arms, or in none, is in no arm)."""
    arms = {next(iter(a)) if len(a) == 1 else "" for a in program_arms.values()}
    if len(arms) == 1 and "" not in arms:
        return arms.pop()
    return ""


def resolve_program_targets(rows: pd.DataFrame) -> pd.DataFrame:
    """One class and direction per (program, gene).

    `rows` carries program, gene, row_class, row_direction, row_reason, one row per distinct
    mechanism row (exact duplicates already collapsed). Informative rows are those whose
    direction is not ambiguous. If the informative rows agree in class and direction, that
    class and direction are used; if they conflict, the program-target is other/ambiguous. With
    no informative row, every row is other/ambiguous (no aligned or blocking row has an
    ambiguous direction), and the program-target is other/ambiguous.
    """
    out = []
    for (program, gene), grp in rows.groupby(["program", "gene"], sort=True):
        informative = grp[grp["row_direction"] != "ambiguous"]
        pairs = set(zip(informative["row_class"], informative["row_direction"]))
        if len(pairs) == 1:
            cls, direction = next(iter(pairs))
            reason = ";".join(sorted(set(informative["row_reason"])))
        elif len(pairs) == 0:
            cls, direction = "other", "ambiguous"
            reason = ";".join(sorted(set(grp["row_reason"])))
        else:
            cls, direction, reason = "other", "ambiguous", OTHER_CONFLICT
        out.append({"program": program, "gene": gene, "pt_class": cls, "pt_direction": direction,
                    "pt_reason": reason, "conflict": len(pairs) > 1})
    return pd.DataFrame(out, columns=["program", "gene", "pt_class", "pt_direction", "pt_reason", "conflict"])
