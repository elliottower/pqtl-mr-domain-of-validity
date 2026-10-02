"""Outcome-GWAS candidates, tiers, sample overlap and the joint instrument/outcome selection
(PREREG §Measured variables, "Outcome-GWAS candidates and tiers", "Sample overlap",
"Instrument and outcome GWAS, selected together").

Candidate mapping is that of `count_all_indications.py` (`gwas_catalog_candidates`,
`finngen_candidates`, `opengwas_candidates`, `indication_tiers`), as the plan requires, rewritten
to take loaded tables. The selection is `count_v4.py` `select_outcome`. No association statistic
is read: candidates carry metadata only.
"""
import re

import pandas as pd
from pydantic import BaseModel, ConfigDict

from stage_a.schemas import InstrumentSource, Overlap

NEFF_FLOOR = 2000.0
UKB_PATTERN = re.compile(r"\bUKB\b|UKBB|UK ?Biobank|UKBiobank", re.IGNORECASE)
ICELAND_PATTERN = re.compile(r"deCODE|Iceland", re.IGNORECASE)
INTERVAL_PATTERN = re.compile(r"\bINTERVAL\b")
OPENGWAS_SKIP_PREFIX = ("eqtl-a-", "prot-", "met-", "ubm-")

OVERLAP_COL = {"ukbppp": "includes_ukb", "decode": "includes_iceland", "interval": "includes_interval"}
OVERLAP_RANK = {"no": 0, "unknown": 1, "yes": 2}
CAND_SOURCE_RANK = {"gwas_catalog": 0, "finngen": 1, "opengwas": 2}
CANDIDATE_COLUMNS = ["source", "accession", "trait", "mapped_ids", "mapping_method", "ncase", "ncontrol",
                     "neff", "european", "full_sumstats", "case_control", "includes_ukb",
                     "includes_iceland", "includes_interval"]


class Selection(BaseModel):
    """The selected (instrument source, outcome candidate) for one gene-indication pair.
    `candidate` is a row label of the candidate table; None when no source has a candidate."""

    model_config = ConfigDict(frozen=True)

    source: InstrumentSource | None
    candidate: int | None
    tier: int | None
    overlap: Overlap | None
    rule: str


def curie(uri: str) -> str:
    return uri.strip().rstrip("/").split("/")[-1].replace(":", "_")


def norm_label(s) -> str:
    s = str(s).lower().replace("'", "").replace("’", "")
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def effective_n(ncase: float, ncontrol: float) -> float | None:
    if not ncase or not ncontrol or ncase <= 0 or ncontrol <= 0:
        return None
    return 4.0 / (1.0 / ncase + 1.0 / ncontrol)


def to_num(x) -> float | None:
    try:
        v = float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return None
    return v if v == v else None


def parse_cases_controls(text) -> tuple[float, float]:
    ncase = ncontrol = 0.0
    for num, word in re.findall(r"(\d[\d,]*)\s[^\d]*?\b(cases?|controls?)\b", str(text)):
        v = float(num.replace(",", ""))
        if word.startswith("case"):
            ncase += v
        else:
            ncontrol += v
    return ncase, ncontrol


def overlap_flag(pattern: re.Pattern, text: str, cohort_listed: bool) -> Overlap:
    """yes if the metadata names the cohort; no if the source lists its cohorts and none
    matches; otherwise unknown."""
    if pattern.search(text):
        return "yes"
    return "no" if cohort_listed else "unknown"


def label_index(disease: pd.DataFrame) -> dict[str, set[str]]:
    idx: dict[str, set[str]] = {}
    for did, name, syn in zip(disease["id"], disease["name"], disease["exactSynonyms"]):
        for lab in [name, *(list(syn) if syn is not None else [])]:
            idx.setdefault(norm_label(lab), set()).add(did)
    return idx


def gwas_catalog_candidates(frames: list[tuple[pd.DataFrame, pd.DataFrame, str]]) -> pd.DataFrame:
    """`frames`: (studies, ancestries, 'published' | 'unpublished') with upper-cased column names."""
    rows = []
    for st, an, kind in frames:
        uri_col = "MAPPED_TRAIT_URI" if "MAPPED_TRAIT_URI" in st.columns else next(
            c for c in st.columns if "TRAIT" in c and "URI" in c and "BACKGROUND" not in c)
        full_col = next((c for c in st.columns if c.startswith("FULL SUMMARY")), None)
        loc_col = next((c for c in st.columns if "SUMMARY STATS LOCATION" in c), None)
        size_col = next((c for c in st.columns if c.startswith("INITIAL SAMPLE")), None)
        trait_col = next(c for c in st.columns if c in ("DISEASE/TRAIT", "TRAIT"))
        cohort_col = "COHORT" if "COHORT" in st.columns else None
        an_init = an[an["STAGE"].str.lower().eq("initial")] if "STAGE" in an.columns else an
        anc = an_init.groupby("STUDY ACCESSION").agg(
            ancestries=("BROAD ANCESTRAL CATEGORY", lambda s: "|".join(sorted(set(s.dropna())))),
            countries=("COUNTRY OF RECRUITMENT", lambda s: "|".join(sorted(set(s.dropna())))),
            ncase_anc=("NUMBER OF CASES", lambda s: sum(to_num(x) or 0 for x in s)),
            ncontrol_anc=("NUMBER OF CONTROLS", lambda s: sum(to_num(x) or 0 for x in s)),
        )
        for r in st.itertuples(index=False):
            rec = dict(zip(st.columns, r))
            acc = rec.get("STUDY ACCESSION")
            ids = sorted({curie(u) for u in str(rec.get(uri_col) or "").split(",")
                          if u.strip() and u.strip().lower() != "nan"})
            a = anc.loc[acc] if acc in anc.index else None
            ancestries = a["ancestries"] if a is not None else ""
            ncase = a["ncase_anc"] if a is not None else 0
            ncontrol = a["ncontrol_anc"] if a is not None else 0
            if not (ncase and ncontrol):
                ncase, ncontrol = parse_cases_controls(rec.get(size_col, ""))
            full = str(rec.get(full_col, "")).strip().lower() == "yes" if full_col else True
            if kind == "unpublished":
                full = full or bool(str(rec.get(loc_col, "") or "").strip() not in ("", "nan"))
            cohort = str(rec.get(cohort_col, "") or "") if cohort_col else ""
            cohort = "" if cohort.lower() == "nan" else cohort
            text = " ".join([cohort, str(rec.get(size_col, "")), a["countries"] if a is not None else ""])
            rows.append({
                "source": "gwas_catalog", "accession": acc, "trait": rec.get(trait_col),
                "mapped_ids": ids, "mapping_method": "gwas_catalog_mapped_trait_uri",
                "ncase": ncase, "ncontrol": ncontrol, "neff": effective_n(ncase, ncontrol),
                "european": ancestries == "European", "full_sumstats": full,
                "case_control": bool(ncase and ncontrol),
                "includes_ukb": overlap_flag(UKB_PATTERN, text, bool(cohort)),
                "includes_iceland": overlap_flag(ICELAND_PATTERN, text, bool(cohort)),
                "includes_interval": overlap_flag(INTERVAL_PATTERN, text, bool(cohort)),
            })
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def finngen_candidates(study: pd.DataFrame) -> pd.DataFrame:
    """FinnGen studies of the Open Targets study index with their OT disease mappings. FinnGen
    lists its cohort, so every overlap flag is no."""
    fg = study[study["studyId"].str.upper().str.startswith("FINNGEN")]
    rows = []
    for rec in fg.to_dict(orient="records"):
        ids: set[str] = set()
        for c in ("diseaseIds", "traitFromSourceMappedIds"):
            v = rec.get(c)
            if v is not None and not isinstance(v, float):
                ids |= {curie(x) for x in v}
        ncase, ncontrol = to_num(rec.get("nCases")) or 0, to_num(rec.get("nControls")) or 0
        rows.append({
            "source": "finngen", "accession": rec["studyId"], "trait": rec.get("traitFromSource"),
            "mapped_ids": sorted(ids), "mapping_method": "open_targets_study_index_finngen",
            "ncase": ncase, "ncontrol": ncontrol, "neff": effective_n(ncase, ncontrol),
            "european": True, "full_sumstats": bool(rec.get("hasSumstats", True)),
            "case_control": bool(ncase and ncontrol),
            "includes_ukb": "no", "includes_iceland": "no", "includes_interval": "no",
        })
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def opengwas_candidates(info: dict, gcat: pd.DataFrame, fgn: pd.DataFrame,
                        labels: dict[str, set[str]]) -> pd.DataFrame:
    """OpenGWAS `gwasinfo`: ebi-a-* through the GWAS Catalog accession, finn-b-* through the
    FinnGen endpoint, otherwise the ontology field, otherwise an exact label/synonym match.
    A `consortium` field is not a cohort list, so only finn-b-* can be no."""
    gcat_ids = dict(zip(gcat["accession"], gcat["mapped_ids"]))
    fg_code: dict[str, set[str]] = {}
    for acc, ids in zip(fgn["accession"], fgn["mapped_ids"]):
        code = re.sub(r"^FINNGEN_R\d+_", "", acc, flags=re.IGNORECASE)
        fg_code.setdefault(code.upper(), set()).update(ids)
    rows = []
    for gid, x in info.items():
        if gid.startswith(OPENGWAS_SKIP_PREFIX):
            continue
        ids: set[str] = set()
        method = "unmapped"
        if gid.startswith("ebi-a-") and gid[6:] in gcat_ids:
            ids, method = set(gcat_ids[gid[6:]]), "gwas_catalog_accession"
        elif gid.startswith("finn-b-") and gid[7:].upper() in fg_code:
            ids, method = set(fg_code[gid[7:].upper()]), "open_targets_finngen_endpoint_code"
        ont = str(x.get("ontology") or "")
        if not ids and ont not in ("", "NA", "None", "nan"):
            ids, method = {curie(o) for o in re.split(r"[;,| ]+", ont) if o}, "opengwas_ontology_field"
        if not ids and norm_label(x.get("trait", "")) in labels:
            ids, method = set(labels[norm_label(x["trait"])]), "exact_label_or_synonym"
        ncase, ncontrol = to_num(x.get("ncase")) or 0, to_num(x.get("ncontrol")) or 0
        text = " ".join(str(x.get(k, "")) for k in ("consortium", "author", "note"))
        finn = gid.startswith("finn-b-")
        rows.append({
            "source": "opengwas", "accession": gid, "trait": x.get("trait"),
            "mapped_ids": sorted(ids), "mapping_method": method,
            "ncase": ncase, "ncontrol": ncontrol, "neff": effective_n(ncase, ncontrol),
            "european": str(x.get("population")) == "European", "full_sumstats": True,
            "case_control": bool(ncase and ncontrol),
            "includes_ukb": "yes" if (gid.startswith("ukb-") or UKB_PATTERN.search(text)) else (
                "no" if finn else "unknown"),
            "includes_iceland": "yes" if ICELAND_PATTERN.search(text) else ("no" if finn else "unknown"),
            "includes_interval": "yes" if INTERVAL_PATTERN.search(text) else ("no" if finn else "unknown"),
        })
    return pd.DataFrame(rows, columns=CANDIDATE_COLUMNS)


def combine_candidates(gcat: pd.DataFrame, fgn: pd.DataFrame, ogw: pd.DataFrame) -> pd.DataFrame:
    """Candidate table with the `qualifies` flag: case-control, European, full summary
    statistics, effective N >= 2,000."""
    frames = [f for f in (gcat, fgn, ogw) if len(f)]
    cands = (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()).reindex(columns=CANDIDATE_COLUMNS)
    cands["qualifies"] = (cands["case_control"].astype(bool) & cands["european"].astype(bool)
                          & cands["full_sumstats"].astype(bool)
                          & (cands["neff"].astype(float).fillna(0) >= NEFF_FLOOR))
    return cands


def _as_set(x) -> set[str]:
    return set(x) if x is not None else set()


def indication_tiers(indications: set[str], cands: pd.DataFrame, disease: pd.DataFrame) -> dict[str, dict[int, list[int]]]:
    """indication -> {1: [...], 2: [...], 3: [...]} qualifying candidate row labels. Tier 1: a
    mapped ID equals the indication or an obsolete term OT lists for it; tier 2: a descendant;
    tier 3: an ancestor. A candidate in tier 1 is not repeated in tiers 2 or 3."""
    q = cands[cands["qualifies"]]
    by_id: dict[str, list[int]] = {}
    for i, ids in zip(q.index, q["mapped_ids"]):
        for d in ids:
            by_id.setdefault(d, []).append(int(i))
    rel = disease.set_index("id")
    out = {}
    for d in sorted(indications):
        if d in rel.index:
            r = rel.loc[d]
            desc, anc, obs = _as_set(r["descendants"]), _as_set(r["ancestors"]), _as_set(r["obsoleteTerms"])
        else:
            desc = anc = obs = set()
        exact = {d} | {curie(o) for o in obs}
        t1 = sorted({i for x in exact for i in by_id.get(x, [])})
        t2 = sorted({i for x in desc for i in by_id.get(x, [])} - set(t1))
        t3 = sorted({i for x in anc for i in by_id.get(x, [])} - set(t1))
        out[d] = {1: t1, 2: t2, 3: t3}
    return out


def select_outcome(sources: list[str], pools: dict[int, list[int]], cands: pd.DataFrame) -> Selection:
    """Joint selection for one gene-indication pair. Per source (registered order), the first
    candidate by tier, independence (no < unknown < yes), larger effective N, source order GWAS
    Catalog < FinnGen < OpenGWAS, lower accession. A deCODE instrument is never paired with a
    candidate that includes Icelandic participants. The first source whose selection is tier 1
    and no is used; otherwise the first source's selection."""
    per_source = []
    for s in sources:
        col = OVERLAP_COL[s]
        ranked = []
        for tier in (1, 2, 3):
            for i in pools.get(tier, []):
                ov = cands.at[i, col]
                if s == "decode" and cands.at[i, "includes_iceland"] == "yes":
                    continue
                ranked.append((tier, OVERLAP_RANK[ov], -float(cands.at[i, "neff"]),
                               CAND_SOURCE_RANK[cands.at[i, "source"]], str(cands.at[i, "accession"]), i))
        if ranked:
            best = min(ranked)
            per_source.append(Selection(source=s, candidate=best[5], tier=best[0], overlap=cands.at[best[5], col],
                                        rule="per_source"))
        else:
            per_source.append(Selection(source=s, candidate=None, tier=None, overlap=None, rule="per_source"))
    for sel in per_source:
        if sel.tier == 1 and sel.overlap == "no":
            return sel.model_copy(update={"rule": "first_source_tier1_no"})
    if not per_source:
        return Selection(source=None, candidate=None, tier=None, overlap=None, rule="no_instrument_source")
    return per_source[0].model_copy(update={"rule": "first_source_fallback"})


# ---- outcome trait coding (PREREG §Measured variables, "Genetic effect direction") ------------

RISK_CODED_BASIS = "case_control: effect on case status, positive beta = higher risk"
NOT_RISK_CODED_BASIS = "not case-control: no case status to code risk on"


class MixedTraitCoding(ValueError):
    """One outcome accession names candidates whose trait coding differs."""


def trait_coding(ncase, ncontrol) -> tuple[bool, str]:
    """(risk_coded, basis). The plan's outcome traits are case-control disease, coded so that a
    positive beta means higher risk: every source's case-control summary statistics report the
    effect on case status (log odds, or a linear model on the case indicator), so the sign is
    risk-coded whatever the effect unit. A candidate without both case and control counts has no
    case status and is not risk-coded; the qualifying rule never selects one, so a false value
    in the table flags a contract violation upstream."""
    if (to_num(ncase) or 0) > 0 and (to_num(ncontrol) or 0) > 0:
        return True, RISK_CODED_BASIS
    return False, NOT_RISK_CODED_BASIS


def outcome_trait_coding(selected: dict[tuple[str, str], int], cands: pd.DataFrame) -> pd.DataFrame:
    """One row per selected outcome accession, the file stage B reads (`units.load_trait_coding`
    wants outcome_accession, risk_coded as its first two columns). `selected`:
    (outcome source, accession) -> candidate row label."""
    rows, by_acc = [], {}
    for (source, acc), i in sorted(selected.items()):
        c = cands.loc[i]
        coded, basis = trait_coding(c["ncase"], c["ncontrol"])
        if by_acc.setdefault(acc, coded) != coded:
            raise MixedTraitCoding(f"outcome accession {acc} has candidates coded both ways")
        rows.append({"outcome_accession": acc, "risk_coded": coded, "outcome_source": source,
                     "trait": str(c["trait"] if c["trait"] is not None else ""),
                     "n_case": int(round(float(c["ncase"]))), "n_control": int(round(float(c["ncontrol"]))),
                     "coding_basis": basis})
    return pd.DataFrame(rows, columns=["outcome_accession", "risk_coded", "outcome_source", "trait", "n_case",
                                       "n_control", "coding_basis"]).drop_duplicates("outcome_accession")
