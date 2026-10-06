"""LD from the 1000 Genomes phase 3 EUR panel, GRCh38 (PREREG §Data collection, §Coverage,
§Flags, S15g). Variants are matched to the panel by rsID, like everything else in stage B.

Input is the text of `bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n'` over the EUR
samples; fetch.py produces it. The registered panel (the 20190312 biallelic SNV and indel release
on GRCh38) carries "." in its ID column, so a record's rsID is not read from the panel: it is the
rsID the UKB-PPP rsID map of the chromosome (Synapse syn51396727, the map the plan registers for
rsID matching) gives for the record's GRCh38 position (the map's POS38) and its unordered pair of
REF and ALT alleles, as for a GWAS-SSF file without an rsID column (outcome_files.map_rsids,
`map_rsid`). No liftover and no strand flip: a record whose pair is the complement of the map's has
no rsID. A record with no such map row, or with rows naming more than one rsID, has no rsID and is
dropped. The share of EUR-polymorphic records (`eur_polymorphic`) that get an rsID must reach
MIN_RSID_FRACTION.
"""
from collections.abc import Iterable, Mapping

import numpy as np
import pandas as pd

from stage_b.outcome_files import map_rsid
from stage_b.schemas import PROXY_R2, LDReferenceError

# The check exists to detect a panel whose records cannot be matched by rsID. It is taken over the records
# polymorphic among the EUR samples only: an EUR-monomorphic record carries no LD information, and most of the
# panel's records without a map rsID are rare variants that are monomorphic in EUR.
MIN_RSID_FRACTION = 0.5


def eur_polymorphic(genotypes: Iterable[str]) -> bool:
    """True when the called alleles of the EUR genotypes hold at least one reference (0) and at least one
    alternate allele; a missing allele (.) counts as neither. Haploid calls count their one allele."""
    alleles = {x for gt in genotypes for x in gt.replace("|", "/").split("/") if x != "."}
    return "0" in alleles and bool(alleles - {"0"})


def parse_genotypes(lines: Iterable[str],
                    rsids: Mapping[tuple[int, frozenset], set[str]]) -> tuple[pd.DataFrame, np.ndarray]:
    """(meta with rsid/ref/alt, dosage matrix variants x samples). `rsids` is the UKB-PPP map at the
    records' GRCh38 positions (outcome_files.map_rsids on GRCh38). Each record takes its rsID from
    it (module docstring); the panel's ID column is not read. Multiallelic records, records without
    an rsID and rsIDs seen on more than one record are dropped; a missing genotype takes the
    variant's mean dosage."""
    meta, rows, polymorphic = [], [], []
    for line in lines:
        f = line.rstrip("\n").split("\t")
        if len(f) < 6 or "," in f[4]:
            continue
        g = []
        for gt in f[5:]:
            a = gt.replace("|", "/").split("/")
            g.append(float("nan") if "." in a else float(sum(int(x) for x in a)))
        pos, ref, alt = int(f[1]), f[3].upper(), f[4].upper()
        meta.append({"chrom": f[0], "pos": pos, "rsid": map_rsid(rsids, pos, ref, alt), "ref": ref, "alt": alt})
        rows.append(g)
        polymorphic.append(eur_polymorphic(f[5:]))
    if not rows:
        raise LDReferenceError("no biallelic records returned for the region")
    m = pd.DataFrame(meta)
    d = np.array(rows, dtype=float)
    poly = np.array(polymorphic)
    if not poly.any():
        raise LDReferenceError("no record of the region is polymorphic among the EUR samples")
    rs_frac = float((m["rsid"].to_numpy()[poly] != "").mean())
    if rs_frac < MIN_RSID_FRACTION:
        raise LDReferenceError(f"only {rs_frac:.2f} of the EUR-polymorphic panel records carry an rsID; "
                               "rsID matching is impossible")
    keep = (m["rsid"] != "") & ~m["rsid"].duplicated(keep=False)
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
