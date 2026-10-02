"""Streaming parsers: source-specific rows -> the canonical regional table (schemas.VARIANT_COLUMNS).

Every parser takes a header and an iterator of already split rows, keeps only rows inside the
window, and never materializes the whole file. Column names are those of each source's own
documentation (planning/perplexity_v8_design/F_data_readmes/NOTES.md and READMEproteomics.txt;
UKB-PPP REGENIE output as described in Sun et al. 2023 ST9 abbreviations; FinnGen R12 and GWAS
Catalog harmonised headers). A header missing a required column raises InputContractError.
"""
import math
import re
from collections.abc import Iterable, Iterator

import numpy as np
import pandas as pd

from stage_b.schemas import VARIANT_COLUMNS, InputContractError

MISSING = {"", "NA", "nan", "NaN", "None", ".", "-"}


def normalize_chrom(value: str | int) -> str:
    """'chr1' / '1' -> '1'; '23' / 'chrX' -> 'X'."""
    c = str(value).strip()
    c = c[3:] if c.lower().startswith("chr") else c
    if re.fullmatch(r"\d+\.0+", c):     # a chromosome column read as float by a spreadsheet loader: 5.0
        c = c.split(".", 1)[0]
    return {"23": "X", "24": "Y", "x": "X", "y": "Y"}.get(c, c)


def in_window(chrom: str, pos: int, center_chrom: str, center_pos: int, half_width: int) -> bool:
    return normalize_chrom(chrom) == normalize_chrom(center_chrom) and abs(int(pos) - center_pos) <= half_width


def split_lines(lines: Iterable[str], sep: str | None) -> Iterator[list[str]]:
    for line in lines:
        line = line.rstrip("\r\n")
        if line and not line.startswith("##"):
            yield line.split(sep) if sep is not None else line.split()


def _index(header: list[str], required: list[str], what: str) -> dict[str, int]:
    idx = {h.lstrip("#"): i for i, h in enumerate(header)}
    missing = [c for c in required if c not in idx]
    if missing:
        raise InputContractError(f"{what}: header lacks {missing}; header was {header}")
    return idx


def _first(idx: dict[str, int], options: list[str], what: str, required: bool = True) -> int | None:
    for o in options:
        if o in idx:
            return idx[o]
    if required:
        raise InputContractError(f"{what}: header has none of {options}")
    return None


def _num(x: str) -> float:
    return float("nan") if x in MISSING else float(x)


def _frame(rows: list[dict]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=VARIANT_COLUMNS)
    for c in ("pos",):
        df[c] = df[c].astype("int64")
    for c in ("eaf", "beta", "se", "p", "n"):
        df[c] = df[c].astype("float64")
    return df


def _p_from_log10(neg_log10_p: float) -> float:
    return float(10.0 ** (-neg_log10_p)) if math.isfinite(neg_log10_p) else float("nan")


# ---- UKB-PPP (REGENIE per-chromosome file inside the per-protein tar) ---------------------------

def parse_ukbppp_rsid_map(rows: Iterable[list[str]], header: list[str], wanted_ids: set[str]) -> dict[str, str]:
    """UKB-PPP rsID map (syn51396727): variant ID -> rsID, for the IDs in `wanted_ids` only."""
    idx = _index(header, ["ID", "rsid"], "UKB-PPP rsID map")
    out = {}
    for r in rows:
        vid, rs = r[idx["ID"]], r[idx["rsid"]]
        if vid in wanted_ids and rs.startswith("rs"):
            out[vid] = rs
    return out


def filter_ukbppp(rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                  half_width: int) -> list[dict]:
    """Window rows of a UKB-PPP REGENIE file as dicts (ID kept for the rsID join). Positions are
    GENPOS, GRCh38; ALLELE1 is the effect allele and A1FREQ its frequency; LOG10P is -log10 p."""
    idx = _index(header, ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"],
                 "UKB-PPP REGENIE file")
    out = []
    for r in rows:
        if not in_window(r[idx["CHROM"]], int(r[idx["GENPOS"]]), chrom, center, half_width):
            continue
        out.append({"ID": r[idx["ID"]], "chrom": normalize_chrom(r[idx["CHROM"]]), "pos": int(r[idx["GENPOS"]]),
                    "ea": r[idx["ALLELE1"]].upper(), "oa": r[idx["ALLELE0"]].upper(),
                    "eaf": _num(r[idx["A1FREQ"]]), "beta": _num(r[idx["BETA"]]), "se": _num(r[idx["SE"]]),
                    "p": _p_from_log10(_num(r[idx["LOG10P"]])), "n": _num(r[idx["N"]])})
    return out


def ukbppp_to_canonical(window_rows: list[dict], rsid_map: dict[str, str]) -> pd.DataFrame:
    """Attach rsIDs from the UKB-PPP map; a variant absent from the map keeps an empty rsID (it
    counts in the pQTL window and is dropped at harmonization)."""
    return _frame([{**{k: v for k, v in r.items() if k != "ID"}, "rsid": rsid_map.get(r["ID"], "")}
                   for r in window_rows])


# ---- deCODE (Ferkingstad 2021 non-normalized per-protein file + annotation + excluded) ----------

def filter_decode(rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                  half_width: int) -> list[dict]:
    idx = _index(header, ["Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "SE", "N"],
                 "deCODE per-protein file")
    out = []
    for r in rows:
        if not in_window(r[idx["Chrom"]], int(r[idx["Pos"]]), chrom, center, half_width):
            continue
        out.append({"Name": r[idx["Name"]], "rsids": r[idx["rsids"]], "chrom": normalize_chrom(r[idx["Chrom"]]),
                    "pos": int(r[idx["Pos"]]), "ea": r[idx["effectAllele"]].upper(),
                    "oa": r[idx["otherAllele"]].upper(), "beta": _num(r[idx["Beta"]]),
                    "se": _num(r[idx["SE"]]), "p": _num(r[idx["Pval"]]), "n": _num(r[idx["N"]])})
    return out


def filter_decode_annotation(rows: Iterable[list[str]], header: list[str], names: set[str]) -> dict[str, tuple[str, str, float]]:
    """assocvariants.annotated.txt.gz: Name -> (effectAllele, otherAllele, effectAlleleFreq)."""
    idx = _index(header, ["Name", "effectAllele", "otherAllele", "effectAlleleFreq"], "deCODE annotation file")
    return {r[idx["Name"]]: (r[idx["effectAllele"]].upper(), r[idx["otherAllele"]].upper(),
                             _num(r[idx["effectAlleleFreq"]]))
            for r in rows if r[idx["Name"]] in names}


def filter_decode_excluded(rows: Iterable[list[str]], header: list[str], names: set[str]) -> set[str]:
    idx = _index(header, ["Name"], "deCODE excluded-variants file")
    return {r[idx["Name"]] for r in rows if r[idx["Name"]] in names}


def decode_to_canonical(window_rows: list[dict], annotation: dict[str, tuple[str, str, float]],
                        excluded: set[str]) -> tuple[pd.DataFrame, dict]:
    """PREREG §deCODE release: effect-allele frequency from the annotation file (not ImpMAF),
    multiallelic rows take its corrected alleles (otherAllele '!'), excluded variants dropped.
    A row whose annotation effect allele differs from the file's is dropped and counted. A row
    with several rsIDs yields one row per rsID; rsID duplicates are resolved at harmonization."""
    counts = {"excluded": 0, "not_in_annotation": 0, "annotation_allele_mismatch": 0}
    out = []
    for r in window_rows:
        if r["Name"] in excluded:
            counts["excluded"] += 1
            continue
        ann = annotation.get(r["Name"])
        if ann is None:
            counts["not_in_annotation"] += 1
            continue
        ea, oa, eaf = ann
        if ea != r["ea"]:
            counts["annotation_allele_mismatch"] += 1
            continue
        rsids = [s for s in r["rsids"].replace(";", ",").split(",") if s.startswith("rs")] or [""]
        for rs in rsids:
            out.append({"rsid": rs, "chrom": r["chrom"], "pos": r["pos"], "ea": ea, "oa": oa, "eaf": eaf,
                        "beta": r["beta"], "se": r["se"], "p": r["p"], "n": r["n"]})
    return _frame(out), counts


# ---- OpenGWAS (INTERVAL prot-a-* and OpenGWAS outcomes; GRCh37) ---------------------------------

def opengwas_to_canonical(records: list[dict]) -> pd.DataFrame:
    """Records of POST /associations: rsid, chr, position, ea, nea, eaf, beta, se, p, n."""
    out = []
    for r in records:
        for k in ("rsid", "chr", "position", "ea", "nea", "beta", "se"):
            if k not in r:
                raise InputContractError(f"OpenGWAS association record lacks '{k}': {r}")
        out.append({"rsid": str(r["rsid"]) if str(r["rsid"]).startswith("rs") else "",
                    "chrom": normalize_chrom(r["chr"]), "pos": int(r["position"]),
                    "ea": str(r["ea"]).upper(), "oa": str(r["nea"]).upper(),
                    "eaf": _num(str(r.get("eaf", ""))), "beta": _num(str(r["beta"])), "se": _num(str(r["se"])),
                    "p": _num(str(r.get("p", ""))), "n": _num(str(r.get("n", "")))})
    return _frame(out)


# ---- GWAS Catalog harmonised summary statistics (GRCh38) ----------------------------------------

def filter_gwas_catalog(rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                        half_width: int) -> pd.DataFrame:
    """Old (hm_*) and GWAS-SSF harmonised layouts. beta is hm_beta/beta, else log of the odds
    ratio; the rsID is hm_rsid/rsid, else variant_id where it is an rsID."""
    idx = {h.lstrip("#"): i for i, h in enumerate(header)}
    w = "GWAS Catalog harmonised file"
    i_chr = _first(idx, ["hm_chrom", "chromosome"], w)
    i_pos = _first(idx, ["hm_pos", "base_pair_location"], w)
    i_ea = _first(idx, ["hm_effect_allele", "effect_allele"], w)
    i_oa = _first(idx, ["hm_other_allele", "other_allele"], w)
    i_beta = _first(idx, ["hm_beta", "beta"], w, required=False)
    i_or = _first(idx, ["hm_odds_ratio", "odds_ratio"], w, required=False)
    if i_beta is None and i_or is None:
        raise InputContractError(f"{w}: neither beta nor odds ratio in header {header}")
    i_se = _first(idx, ["standard_error"], w)
    i_eaf = _first(idx, ["hm_effect_allele_frequency", "effect_allele_frequency"], w, required=False)
    i_p = _first(idx, ["p_value"], w, required=False)
    i_rs = _first(idx, ["hm_rsid", "rsid", "rs_id", "variant_id"], w)
    i_n = _first(idx, ["n"], w, required=False)
    out = []
    for r in rows:
        if r[i_chr] in MISSING or r[i_pos] in MISSING:
            continue
        if not in_window(r[i_chr], int(float(r[i_pos])), chrom, center, half_width):
            continue
        beta = _num(r[i_beta]) if i_beta is not None else float("nan")
        if not math.isfinite(beta) and i_or is not None:
            orv = _num(r[i_or])
            beta = math.log(orv) if math.isfinite(orv) and orv > 0 else float("nan")
        rs = r[i_rs]
        out.append({"rsid": rs if rs.startswith("rs") else "", "chrom": normalize_chrom(r[i_chr]),
                    "pos": int(float(r[i_pos])), "ea": r[i_ea].upper(), "oa": r[i_oa].upper(),
                    "eaf": _num(r[i_eaf]) if i_eaf is not None else float("nan"), "beta": beta,
                    "se": _num(r[i_se]), "p": _num(r[i_p]) if i_p is not None else float("nan"),
                    "n": _num(r[i_n]) if i_n is not None else float("nan")})
    return _frame(out)


# ---- FinnGen R12 (GRCh38) ------------------------------------------------------------------------

def filter_finngen(rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                   half_width: int) -> pd.DataFrame:
    """FinnGen: effect allele is `alt`, frequency `af_alt`; `rsids` may list several."""
    idx = _index(header, ["chrom", "pos", "ref", "alt", "rsids", "pval", "beta", "sebeta", "af_alt"], "FinnGen file")
    out = []
    for r in rows:
        if not in_window(r[idx["chrom"]], int(r[idx["pos"]]), chrom, center, half_width):
            continue
        rsids = [s for s in r[idx["rsids"]].split(",") if s.startswith("rs")] or [""]
        for rs in rsids:
            out.append({"rsid": rs, "chrom": normalize_chrom(r[idx["chrom"]]), "pos": int(r[idx["pos"]]),
                        "ea": r[idx["alt"]].upper(), "oa": r[idx["ref"]].upper(), "eaf": _num(r[idx["af_alt"]]),
                        "beta": _num(r[idx["beta"]]), "se": _num(r[idx["sebeta"]]), "p": _num(r[idx["pval"]]),
                        "n": float("nan")})
    return _frame(out)


# ---- eQTL Catalogue tabix files (GTEx v8 ge and leafcutter; GRCh38) ----------------------------

def filter_eqtl_catalogue(rows: Iterable[list[str]], header: list[str], chrom: str, center: int,
                          half_width: int) -> pd.DataFrame:
    """eQTL Catalogue summary-statistics rows. The ALT allele is always the effect allele (eQTL
    Catalogue FAQ); its frequency is ac/an; N is an/2. molecular_trait_id and gene_id are kept."""
    idx = _index(header, ["molecular_trait_id", "chromosome", "position", "ref", "alt", "pvalue", "beta", "se",
                          "ac", "an", "gene_id", "rsid"], "eQTL Catalogue file")
    out = []
    for r in rows:
        if not in_window(r[idx["chromosome"]], int(r[idx["position"]]), chrom, center, half_width):
            continue
        an, ac = _num(r[idx["an"]]), _num(r[idx["ac"]])
        rs = r[idx["rsid"]]
        out.append({"rsid": rs if rs.startswith("rs") else "", "chrom": normalize_chrom(r[idx["chromosome"]]),
                    "pos": int(r[idx["position"]]), "ea": r[idx["alt"]].upper(), "oa": r[idx["ref"]].upper(),
                    "eaf": ac / an if an > 0 else float("nan"), "beta": _num(r[idx["beta"]]),
                    "se": _num(r[idx["se"]]), "p": _num(r[idx["pvalue"]]),
                    "n": an / 2 if math.isfinite(an) else float("nan"),
                    "molecular_trait_id": r[idx["molecular_trait_id"]], "gene_id": r[idx["gene_id"]]})
    df = pd.DataFrame(out, columns=VARIANT_COLUMNS + ["molecular_trait_id", "gene_id"])
    df["pos"] = df["pos"].astype("int64")
    for c in ("eaf", "beta", "se", "p", "n"):
        df[c] = df[c].astype("float64")
    return df


def restrict_window(df: pd.DataFrame, chrom: str, center: int, half_width: int) -> pd.DataFrame:
    """Subset a canonical table to a narrower window (the ±500 kb primary is a subset of ±1 Mb)."""
    keep = (df["chrom"].map(normalize_chrom) == normalize_chrom(chrom)) & (np.abs(df["pos"] - center) <= half_width)
    return df[keep].reset_index(drop=True)
