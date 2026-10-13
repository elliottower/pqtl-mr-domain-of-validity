"""Q4 (Open Targets GWAS credible-set support G and G_pre beside S) and Q5 (S_0.5 in place of S)."""
from dataclasses import replace

import numpy as np
import pandas as pd

from followup.common import FollowupError, wilson
from followup.data import check_ensembl
from stage_d.designs import Design
from stage_d.evidence import _match_opposite

L2G_MIN = 0.5
EVIDENCE_COLUMNS = ["studyLocusId", "targetId", "diseaseId", "resourceScore", "score", "publicationDate"]
PP_BINS = (0.0, 0.1, 0.5, 0.8)


def check_evidence(ev: pd.DataFrame) -> dict:
    """Every row: resourceScore present and equal to score (the documented L2G score)."""
    if ev["resourceScore"].isna().any():
        raise FollowupError(f"{int(ev['resourceScore'].isna().sum())} credible-set rows have no resourceScore")
    if not np.array_equal(ev["resourceScore"].to_numpy(), ev["score"].to_numpy()):
        raise FollowupError(f"{int((ev['resourceScore'] != ev['score']).sum())} credible-set rows have resourceScore != score")
    return {"rows": int(len(ev)), "resourceScore_min": float(ev["resourceScore"].min()),
            "resourceScore_max": float(ev["resourceScore"].max())}


def closure(disease: pd.DataFrame) -> dict[str, set[str]]:
    """D(i) = {i} plus its descendants in the 26.09 disease table."""
    return {i: {i, *(d or [])} for i, d in zip(disease["id"], disease["descendants"])}


def support(df: pd.DataFrame, ev: pd.DataFrame, closures: dict[str, set[str]]) -> pd.DataFrame:
    """G, G_loci and G_pre per hypothesis (pseudocode in PREREG §Measured variables)."""
    check_ensembl(df["gene_ensembl"], "V8 gene_ensembl")
    q = ev.loc[ev["resourceScore"] >= L2G_MIN].drop_duplicates(["studyLocusId", "targetId", "diseaseId", "publicationDate"])
    check_ensembl(q["targetId"], "Open Targets targetId")
    q = q.assign(pub=pd.to_datetime(q["publicationDate"], errors="raise", format="%Y-%m-%d"))
    by_target = {t: g for t, g in q.groupby("targetId")}
    missing_ind = sorted(set(df["indication_id"]) - set(closures))
    if missing_ind:
        raise FollowupError(f"{len(missing_ind)} indications not in the disease table (first {missing_ind[0]})")
    rows = []
    for hid, gene, ind, start in zip(df["hypothesis_id"], df["gene_ensembl"], df["indication_id"],
                                     df["earliest_phase2_start"]):
        g = by_target.get(gene)
        hits = g.loc[g["diseaseId"].isin(closures[ind])] if g is not None else None
        loci = set() if hits is None else set(hits["studyLocusId"])
        if pd.isna(start):
            pre = pd.NA
        else:
            pre = int(hits is not None and bool((hits["pub"].notna() & (hits["pub"] < start)).any()))
        rows.append({"hypothesis_id": hid, "G": int(bool(loci)), "G_loci": len(loci), "G_pre": pre})
    return pd.DataFrame(rows).astype({"G_pre": "Int64"})


def augment(d: Design, new_id: str, extra: dict[str, np.ndarray]) -> Design:
    """The design with extra fixed-effect columns appended (focal and role columns unchanged)."""
    X = np.column_stack([d.X, *extra.values()])
    return replace(d, design_id=new_id, X=X, columns=d.columns + tuple(extra), focal=np.concatenate([d.focal, np.zeros(len(extra))]))


def g_columns(d: Design, g_by_id: pd.Series, with_interaction: bool) -> dict[str, np.ndarray]:
    g = g_by_id.loc[d.hypothesis_ids].to_numpy(dtype=float)
    if np.isnan(g).any():
        raise FollowupError(f"{d.design_id}: G missing for {int(np.isnan(g).sum())} analysed rows")
    cols = {"G": g}
    if with_interaction:
        cols["GxA"] = g * d.X[:, d.stratum_col]
    return cols


def cell_table(frame: pd.DataFrame, g_col: str) -> list[dict]:
    """Advancement rate with its Wilson interval in every class x S x G cell of the analysed rows."""
    out = []
    a = frame.loc[frame["y"].notna() & frame[g_col].notna()]
    for (cls, s, g), grp in a.groupby(["cls", "S", g_col]):
        k, n = int(grp["y"].sum()), int(len(grp))
        out.append({"class": cls, "S": int(s), g_col: int(g), "n": n, "advanced": k, "rate": k / n,
                    "wilson95": wilson(k, n)})
    return out


def pp_bin_table(frame: pd.DataFrame) -> list[dict]:
    """Advancement by PP.H4 bin x direction state, analysed rows with a completed colocalization."""
    a = frame.loc[frame["y"].notna() & frame["coloc_run"].astype(bool)]
    match, opposite = _match_opposite(a["genetic_direction"], a["direction"])
    direction = np.where(match, "concordant", np.where(opposite, "discordant", "ambiguous"))
    edges = list(PP_BINS) + [np.inf]
    labels = ["[0,0.1)", "[0.1,0.5)", "[0.5,0.8)", "[0.8,1]"]
    bins = pd.cut(a["pp_h4"].astype(float), edges, right=False, labels=labels)
    out = []
    for (b, dstate), grp in a.assign(bin=bins, dir=direction).groupby(["bin", "dir"], observed=False):
        k, n = int(grp["y"].sum()), int(len(grp))
        out.append({"pp_h4_bin": str(b), "direction": dstate, "n": n, "advanced": k,
                    "rate": k / n if n else None, "wilson95": wilson(k, n)})
    return out
