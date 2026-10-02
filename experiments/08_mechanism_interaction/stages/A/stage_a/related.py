"""Related-indication rule and S12 (PREREG §Study design, "Related indications",
"Specificity rank"; §Other planned analysis, S12). Ported from `count_v4.py`
(`related_indication_rule`, `s12_one_per_program`) with `specificity_rank` ordering and
lower-ID tie-breaking.
"""
import pandas as pd

from stage_a.indications import specificity_order_key


def related_indication_rule(hyp: pd.DataFrame, descendants: dict[str, set[str]],
                            rank: dict[str, int]) -> pd.Series:
    """True = retained. `hyp` carries gene, direction, mechanism_class, indication_id and
    programs (frozenset of Phase II+ canonical programs). Within one gene, direction and class,
    terms are processed most specific first; X is dropped when every program of X is a program
    of some already-retained hypothesis on an ontology descendant of X."""
    keep = pd.Series(True, index=hyp.index)
    for _, g in hyp.groupby(["gene", "direction", "mechanism_class"], sort=True):
        if len(g) < 2:
            continue
        order = sorted(g.index, key=lambda i: specificity_order_key(g.at[i, "indication_id"], rank))
        retained: list = []
        for i in order:
            desc = descendants.get(g.at[i, "indication_id"], set())
            covered: set = set()
            for j in retained:
                if g.at[j, "indication_id"] in desc:
                    covered |= g.at[j, "programs"]
            programs = g.at[i, "programs"]
            if programs and programs <= covered:
                keep[i] = False
            else:
                retained.append(i)
    return keep


def s12_one_per_program(s1: pd.DataFrame, rank: dict[str, int]) -> pd.Series:
    """True = kept. Within each gene, S1 hypotheses ordered by specificity rank, then indication
    ID, then lowest program ID; a hypothesis is kept only if none of its programs is in a
    hypothesis already kept."""
    keep = pd.Series(False, index=s1.index)
    for _, g in s1.groupby("gene", sort=True):
        order = sorted(g.index, key=lambda i: (*specificity_order_key(g.at[i, "indication_id"], rank),
                                               min(g.at[i, "programs"])))
        used: set = set()
        for i in order:
            if not (g.at[i, "programs"] & used):
                keep[i] = True
                used |= g.at[i, "programs"]
    return keep
