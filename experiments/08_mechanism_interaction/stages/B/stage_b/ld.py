"""LD from the 1000 Genomes phase 3 EUR panel, GRCh38 (PREREG §Data collection, §Coverage,
§Flags, S15g). Variants are matched to the panel by rsID, like everything else in stage B.

Input is the text of `bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'` over the EUR
samples; fetch.py produces it.
"""
from collections.abc import Iterable

import numpy as np
import pandas as pd

from stage_b.schemas import PROXY_R2, LDReferenceError

MIN_RSID_FRACTION = 0.5   # a panel whose ID column is mostly not rsIDs cannot be matched by rsID


def parse_genotypes(lines: Iterable[str]) -> tuple[pd.DataFrame, np.ndarray]:
    """(meta with rsid/ref/alt, dosage matrix variants x samples). Multiallelic records and
    rsIDs seen more than once are dropped; a missing genotype takes the variant's mean dosage."""
    meta, rows = [], []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or "," in f[4]:
            continue
        g = []
        for gt in f[5:]:
            a = gt.replace("|", "/").split("/")
            g.append(float("nan") if "." in a else float(sum(int(x) for x in a)))
        meta.append({"chrom": f[0], "pos": int(f[1]), "rsid": f[2], "ref": f[3].upper(), "alt": f[4].upper()})
        rows.append(g)
    if not rows:
        raise LDReferenceError("no biallelic records returned for the region")
    m = pd.DataFrame(meta)
    d = np.array(rows, dtype=float)
    rs_frac = float(m["rsid"].str.startswith("rs").mean())
    if rs_frac < MIN_RSID_FRACTION:
        raise LDReferenceError(f"only {rs_frac:.2f} of panel records carry an rsID; rsID matching is impossible")
    keep = m["rsid"].str.startswith("rs") & ~m["rsid"].duplicated(keep=False)
    m, d = m[keep.to_numpy()].reset_index(drop=True), d[keep.to_numpy()]
    col_mean = np.nanmean(d, axis=1, keepdims=True)
    d = np.where(np.isnan(d), col_mean, d)
    return m, d


def _corr_with(d: np.ndarray, i: int) -> np.ndarray:
    x = d - d.mean(axis=1, keepdims=True)
    sd = np.sqrt((x ** 2).sum(axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (x @ x[i]) / (sd * sd[i])
    return np.where(np.isfinite(r), r, 0.0)


def proxies(meta: pd.DataFrame, dosage: np.ndarray, rsid: str, r2_min: float = PROXY_R2) -> set[str]:
    """rsIDs at r^2 >= r2_min with `rsid` in the panel, the variant itself included. Empty when
    the variant is not in the panel or is monomorphic there."""
    hit = np.flatnonzero(meta["rsid"].to_numpy() == rsid)
    if hit.size == 0:
        return set()
    r = _corr_with(dosage, int(hit[0]))
    if np.allclose(dosage[hit[0]], dosage[hit[0]].mean()):
        return set()
    return set(meta.loc[r ** 2 >= r2_min, "rsid"]) | {rsid}


def aligned_ld(h: pd.DataFrame, meta: pd.DataFrame, dosage: np.ndarray) -> tuple[pd.DataFrame, np.ndarray]:
    """Signed LD for the harmonized variants present in the panel, oriented to the pQTL effect
    allele (the allele both datasets' effects refer to). Variants absent, allele-inconsistent or
    monomorphic in the panel are dropped. Returns (subset of h in panel order, r matrix)."""
    pos = {r: i for i, r in enumerate(meta["rsid"])}
    keep, sign, idx = [], [], []
    for k, r in enumerate(h.itertuples(index=False)):
        i = pos.get(r.rsid)
        if i is None:
            continue
        ref, alt = meta.at[i, "ref"], meta.at[i, "alt"]
        if r.ea == alt and r.oa == ref:
            s = 1.0
        elif r.ea == ref and r.oa == alt:
            s = -1.0
        else:
            continue
        if np.allclose(dosage[i], dosage[i].mean()):
            continue
        keep.append(k)
        sign.append(s)
        idx.append(i)
    if len(idx) < 2:
        raise LDReferenceError(f"only {len(idx)} harmonized variants found in the LD panel")
    x = dosage[idx] * np.array(sign)[:, None]
    ld = np.corrcoef(x)
    return h.iloc[keep].reset_index(drop=True), ld
