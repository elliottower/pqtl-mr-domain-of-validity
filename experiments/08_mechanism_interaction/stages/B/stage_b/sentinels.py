"""Instrument sentinels from each source's own published cis-pQTL table (PREREG §Instrument
lists: "The sentinel and regional statistics used in stage B come from the same source").

- UKB-PPP: Sun et al. 2023 ST9 (significant pQTLs, discovery cohort); cis rows of the assay's
  UKBPPP ProteinID; the row with the largest log10(p) (discovery). Position GENPOS (hg38).
- deCODE: Ferkingstad 2021 ST02; cis rows of the SeqId; the sentinel is Rank (cond. sign.) = 1,
  ties by the larger -Log10(P) (adj.). Position chr/pos (var.), hg38.
- INTERVAL: Sun et al. 2018 ST4; cis rows of the SOMAmer; the smallest meta-analysis p.
  Position Chr/Pos, GRCh37.

Where a hypothesis names more than one assay, the assay whose sentinel has the largest
-log10 p is used (select_assay).
"""
import math
import re
from pathlib import Path

import openpyxl
import pandas as pd

from stage_b.parsers import normalize_chrom
from stage_b.schemas import AmbiguousInstrumentError, InputContractError, Sentinel

OID = re.compile(r"(OID\d+)")


def _pick(df: pd.DataFrame, key: str, score: str) -> dict[str, pd.Series]:
    out = {}
    for k, g in df.groupby(key, sort=True):
        g = g.sort_values([score, "rsid"], ascending=[False, True], kind="stable")
        out[str(k)] = g.iloc[0]
    return out


def ukbppp_sentinels(st9: pd.DataFrame) -> list[Sentinel]:
    """`st9` columns: UKBPPP ProteinID, rsID, CHROM, GENPOS (hg38), log10(p) (discovery), cis/trans."""
    need = ["UKBPPP ProteinID", "rsID", "CHROM", "GENPOS (hg38)", "log10(p) (discovery)", "cis/trans"]
    miss = [c for c in need if c not in st9.columns]
    if miss:
        raise InputContractError(f"UKB-PPP ST9 lacks {miss}")
    d = st9[st9["cis/trans"].astype(str).str.strip() == "cis"].copy()
    d["assay"] = d["UKBPPP ProteinID"].astype(str).str.extract(OID, expand=False)
    d["rsid"] = d["rsID"].astype(str)
    d = d.dropna(subset=["assay"])
    return [Sentinel(source="ukbppp", assay_id=a, rsid=r["rsid"], chrom=normalize_chrom(r["CHROM"]),
                     pos=int(r["GENPOS (hg38)"]), build="GRCh38", neg_log10_p=float(r["log10(p) (discovery)"]),
                     locator=str(r["UKBPPP ProteinID"]))
            for a, r in _pick(d, "assay", "log10(p) (discovery)").items()]


def decode_sentinels(st02: pd.DataFrame) -> list[Sentinel]:
    """`st02` columns (header newlines collapsed to spaces): SeqId, variant, chr (var.), pos (var.),
    cis/ trans, Rank (cond. sign.), -Log10(P) (adj.)."""
    cols = {c: " ".join(str(c).split()) for c in st02.columns}
    d = st02.rename(columns=cols)
    need = ["SeqId", "variant", "chr (var.)", "pos (var.)", "cis/ trans", "Rank (cond. sign.)", "-Log10(P) (adj.)"]
    miss = [c for c in need if c not in d.columns]
    if miss:
        raise InputContractError(f"deCODE ST02 lacks {miss}")
    d = d[(d["cis/ trans"].astype(str).str.strip() == "cis") & (pd.to_numeric(d["Rank (cond. sign.)"]) == 1)].copy()
    d["rsid"] = d["variant"].astype(str)
    return [Sentinel(source="decode", assay_id=a, rsid=r["rsid"], chrom=normalize_chrom(r["chr (var.)"]),
                     pos=int(r["pos (var.)"]), build="GRCh38", neg_log10_p=float(r["-Log10(P) (adj.)"]), locator=a)
            for a, r in _pick(d, "SeqId", "-Log10(P) (adj.)").items()]


def interval_sentinels(st4: pd.DataFrame) -> list[Sentinel]:
    """`st4` columns: SOMAmer ID, Target fullname, Sentinel variant*, Chr, Pos, cis/ trans,
    meta_p (the meta-analysis p column, named by the loader)."""
    need = ["SOMAmer ID", "Target fullname", "Sentinel variant*", "Chr", "Pos", "cis/ trans", "meta_p"]
    miss = [c for c in need if c not in st4.columns]
    if miss:
        raise InputContractError(f"INTERVAL ST4 lacks {miss}")
    d = st4[st4["cis/ trans"].astype(str).str.strip() == "cis"].copy()
    d["rsid"] = d["Sentinel variant*"].astype(str)
    d["score"] = [-math.log10(float(p)) if float(p) > 0 else math.inf for p in d["meta_p"]]
    return [Sentinel(source="interval", assay_id=a, rsid=r["rsid"], chrom=normalize_chrom(r["Chr"]),
                     pos=int(r["Pos"]), build="GRCh37", neg_log10_p=float(r["score"]),
                     locator=str(r["Target fullname"]))
            for a, r in _pick(d, "SOMAmer ID", "score").items()]


def select_assay(assay_ids: list[str], sentinels: dict[tuple[str, str], Sentinel], source: str) -> Sentinel:
    """One sentinel for a hypothesis. Several assays: the strongest sentinel; none: raise."""
    found = [sentinels[(source, a)] for a in assay_ids if (source, a) in sentinels]
    if not found:
        raise AmbiguousInstrumentError(f"no {source} sentinel for assays {assay_ids}")
    return min(found, key=lambda s: (-s.neg_log10_p, s.assay_id))


def interval_opengwas_id(target_fullname: str, gwasinfo: list[dict]) -> str:
    """The OpenGWAS prot-a-* dataset whose trait equals the SOMAmer's target full name. Zero or
    several matches raise AmbiguousInstrumentError (the instrument then has no regional file)."""
    hits = sorted({g["id"] for g in gwasinfo
                   if str(g.get("id", "")).startswith("prot-a-") and str(g.get("trait", "")).strip() == target_fullname.strip()})
    if len(hits) != 1:
        raise AmbiguousInstrumentError(f"INTERVAL target '{target_fullname}' matches {len(hits)} prot-a datasets: {hits}")
    return hits[0]


# ---- workbook loaders (pinned supplementary tables; read only the named sheets) ------------------

def _sheet(path: Path, sheet: str, header_row: int) -> pd.DataFrame:
    ws = openpyxl.load_workbook(path, read_only=True)[sheet]
    rows = ws.iter_rows(min_row=header_row, values_only=True)
    header = [str(c) if c is not None else f"_col{i}" for i, c in enumerate(next(rows))]
    return pd.DataFrame([r for r in rows if any(c is not None for c in r)], columns=header)


def load_ukbppp_st9(path: Path) -> pd.DataFrame:
    return _sheet(path, "ST9", 5)


def load_decode_st02(path: Path) -> pd.DataFrame:
    return _sheet(path, "ST02", 3)


def load_interval_st4(path: Path) -> pd.DataFrame:
    df = _sheet(path, "ST4 - pQTL summary", 5)
    cols = list(df.columns)
    if "Meta-analysis" not in cols:
        raise InputContractError("INTERVAL ST4 header lacks 'Meta-analysis'")
    df["meta_p"] = df.iloc[:, cols.index("Meta-analysis") + 2]   # sub-header row: beta, SE, p
    return df[pd.to_numeric(df["meta_p"], errors="coerce").notna()].reset_index(drop=True)
