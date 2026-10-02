"""Design matrices for every registered model (PREREG §Statistical models).

Model kinds
  h1       advance ~ S + A + S:A + covariates; aligned and function-blocking rows; focal S:A (> 0)
  h2       advance ~ S + covariates; all rows; focal S (> 0)
  h4       advance ~ S + C + S:C + covariates; all rows; focal S:C (< 0)
  h3       advance ~ supportive + contradictory + covariates (reference inconclusive); all rows;
           focal supportive - contradictory (> 0)
  h1_E, h2_E  h1 / h2 with S replaced by the continuous score E
  within   h1 on the genes with hypotheses in both classes, gene fixed effect in place of the
           gene random intercept
Random part in every kind: gene intercept (except `within`), uncorrelated gene slope on the
evidence regressor (S, E, or the supportive indicator for h3), indication intercept and a
multiple-membership drug-program intercept with weights 1/k.

Covariates: platform (SomaScan = 1), z(log10 N_eff), oncology. A covariate constant within the
analysed rows is dropped and the drop is recorded.
"""
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components

COVARIATES = (("platform_somascan", "platform"), ("z_log10_neff", "z_log10_neff"), ("oncology_f", "oncology"))
KINDS = ("h1", "h2", "h4", "h3", "h1_E", "h2_E", "within")


class DesignError(ValueError):
    pass


@dataclass(frozen=True)
class Design:
    design_id: str
    kind: str
    y: np.ndarray
    X: np.ndarray
    columns: tuple[str, ...]
    focal: np.ndarray
    focal_sign: float
    ev_col: int
    stratum_col: int
    inter_col: int
    slope: np.ndarray
    gene: np.ndarray
    n_gene: int
    indication: np.ndarray
    n_indication: int
    W: sp.csr_matrix
    k: np.ndarray
    component: np.ndarray
    group: np.ndarray
    hypothesis_ids: np.ndarray
    gene_fixed: bool
    dropped_covariates: tuple[str, ...] = field(default=())

    @property
    def n(self) -> int:
        return len(self.y)


def component_labels(genes, program_lists) -> np.ndarray:
    """Connected components of the bipartite gene / drug-program graph, one label per row
    (the clustering of power_v9_fast.component_clusters, built from a frame)."""
    gidx = {g: i for i, g in enumerate(sorted(set(genes)))}
    progs = sorted({p for ps in program_lists for p in ps})
    pidx = {p: len(gidx) + j for j, p in enumerate(progs)}
    rows, cols = [], []
    for g, ps in zip(genes, program_lists):
        for p in ps:
            rows.append(gidx[g])
            cols.append(pidx[p])
    n = len(gidx) + len(pidx)
    adj = sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    _, labels = connected_components(adj, directed=False)
    return np.array([labels[gidx[g]] for g in genes])


def program_matrix(program_lists) -> tuple[sp.csr_matrix, np.ndarray]:
    progs = sorted({p for ps in program_lists for p in ps})
    col = {p: j for j, p in enumerate(progs)}
    rows, cols, vals = [], [], []
    for i, ps in enumerate(program_lists):
        for p in ps:
            rows.append(i)
            cols.append(col[p])
            vals.append(1.0 / len(ps))
    k = np.array([len(ps) for ps in program_lists], dtype=float)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(program_lists), len(progs))), k


def crossover_genes(frame: pd.DataFrame) -> list[str]:
    h1 = frame.loc[frame["cls"].isin(["aligned", "blocking"])]
    by_gene = h1.groupby("gene_ensembl")["cls"].nunique()
    return sorted(by_gene.index[by_gene == 2])


def within_gene_estimable(frame: pd.DataFrame) -> tuple[bool, dict]:
    """S:A with a gene fixed effect is identified only from crossover genes with a supportive
    hypothesis in both classes (OPEN_QUESTIONS item 45)."""
    sub = frame.loc[frame["gene_ensembl"].isin(crossover_genes(frame)) & frame["cls"].isin(["aligned", "blocking"])]
    sup = sub.loc[sub["S"] == 1].groupby("gene_ensembl")["cls"].nunique()
    both = int((sup == 2).sum())
    return both > 0, {"crossover_genes": int(sub["gene_ensembl"].nunique()), "rows": len(sub),
                      "genes_supportive_in_both_classes": both}


def rows_for(kind: str, frame: pd.DataFrame) -> pd.DataFrame:
    if kind in ("h1", "h1_E"):
        return frame.loc[frame["cls"].isin(["aligned", "blocking"])]
    if kind == "within":
        return frame.loc[frame["cls"].isin(["aligned", "blocking"]) & frame["gene_ensembl"].isin(crossover_genes(frame))]
    return frame


def build_design(design_id: str, kind: str, frame: pd.DataFrame, y_col: str = "y") -> Design:
    if kind not in KINDS:
        raise DesignError(f"unknown model kind {kind}")
    rows = rows_for(kind, frame)
    rows = rows.loc[rows[y_col].notna()].reset_index(drop=True)
    if rows.empty:
        raise DesignError(f"{design_id}: no analysed rows")
    ev_name = "E" if kind.endswith("_E") else "S"
    ev = rows[ev_name].to_numpy(dtype=float)
    cols: list[np.ndarray] = []
    names: list[str] = []
    ev_col = stratum_col = inter_col = -1
    if kind == "h3":
        sup = (rows["state"] == "supportive").to_numpy(dtype=float)
        con = (rows["state"] == "contradictory").to_numpy(dtype=float)
        cols += [sup, con]
        names += ["supportive", "contradictory"]
        slope = sup
    else:
        cols.append(ev)
        names.append(ev_name)
        ev_col = 0
        slope = ev
        if kind in ("h1", "h1_E", "within", "h4"):
            strat_name = "C" if kind == "h4" else "A"
            strat = rows["C"].to_numpy(dtype=float) if kind == "h4" else (rows["cls"] == "aligned").to_numpy(dtype=float)
            cols += [strat, ev * strat]
            names += [strat_name, f"{ev_name}x{strat_name}"]
            stratum_col, inter_col = 1, 2
    dropped = []
    for col, label in COVARIATES:
        v = rows[col].to_numpy(dtype=float)
        if np.ptp(v) == 0:
            dropped.append(label)
            continue
        cols.append(v)
        names.append(label)
    gene_codes, gene = np.unique(rows["gene_ensembl"].to_numpy(), return_inverse=True)
    if kind == "within":
        for j in range(1, len(gene_codes)):
            cols.append((gene == j).astype(float))
            names.append(f"gene[{gene_codes[j]}]")
    X = np.column_stack(cols)
    focal = np.zeros(X.shape[1])
    if kind == "h3":
        focal[0], focal[1] = 1.0, -1.0
    elif inter_col >= 0:
        focal[inter_col] = 1.0
    else:
        focal[ev_col] = 1.0
    _, ind = np.unique(rows["indication_id"].to_numpy(), return_inverse=True)
    W, k = program_matrix(list(rows["programs"]))
    group = (rows["cls"].astype(str) + "|" + rows["state"].astype(str)).to_numpy()
    return Design(
        design_id=design_id, kind=kind, y=rows[y_col].to_numpy(dtype=float), X=X, columns=tuple(names),
        focal=focal, focal_sign=-1.0 if kind == "h4" else 1.0, ev_col=ev_col, stratum_col=stratum_col,
        inter_col=inter_col, slope=slope, gene=gene, n_gene=len(gene_codes), indication=ind,
        n_indication=int(ind.max()) + 1, W=W, k=k,
        component=component_labels(list(rows["gene_ensembl"]), list(rows["programs"])),
        group=group, hypothesis_ids=rows["hypothesis_id"].to_numpy(), gene_fixed=kind == "within",
        dropped_covariates=tuple(dropped))


def save_design(d: Design, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path, design_id=d.design_id, kind=d.kind, y=d.y, X=d.X, columns=np.array(d.columns), focal=d.focal,
        focal_sign=d.focal_sign, roles=np.array([d.ev_col, d.stratum_col, d.inter_col]), slope=d.slope,
        gene=d.gene, n_gene=d.n_gene, indication=d.indication, n_indication=d.n_indication,
        W_data=d.W.data, W_indices=d.W.indices, W_indptr=d.W.indptr, W_shape=np.array(d.W.shape), k=d.k,
        component=d.component, group=d.group.astype(str), hypothesis_ids=d.hypothesis_ids.astype(str),
        gene_fixed=d.gene_fixed, dropped=np.array(d.dropped_covariates, dtype=str))


def load_design(path: Path) -> Design:
    z = np.load(path, allow_pickle=False)
    roles = z["roles"]
    return Design(
        design_id=str(z["design_id"]), kind=str(z["kind"]), y=z["y"], X=z["X"], columns=tuple(str(c) for c in z["columns"]),
        focal=z["focal"], focal_sign=float(z["focal_sign"]), ev_col=int(roles[0]), stratum_col=int(roles[1]),
        inter_col=int(roles[2]), slope=z["slope"], gene=z["gene"], n_gene=int(z["n_gene"]),
        indication=z["indication"], n_indication=int(z["n_indication"]),
        W=sp.csr_matrix((z["W_data"], z["W_indices"], z["W_indptr"]), shape=tuple(z["W_shape"])), k=z["k"],
        component=z["component"], group=z["group"], hypothesis_ids=z["hypothesis_ids"],
        gene_fixed=bool(z["gene_fixed"]), dropped_covariates=tuple(str(c) for c in z["dropped"]))
