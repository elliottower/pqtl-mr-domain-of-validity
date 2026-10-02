"""Held-out keys, foreknowledge flags and publication dates (PREREG §Study design, "Sets";
§Measured variables, "Flags"; §Explanation of foreknowledge, items 4-5).

The pilot label -> Open Targets ID map is `count_v4.py`'s: frozen_candidates_v34.csv disease
names matched exactly (case-insensitive) to OT 26.09 disease names, with the multiple-sclerosis
override of count_coverage.py.
"""
from datetime import date

import pandas as pd

from stage_a.instruments import SOURCE_ORDER

MS_OVERRIDE = {"multiple sclerosis": "EFO_0803536"}  # count_coverage.py OT_ID_OVERRIDES

# Journal online publication dates (Crossref `published-online`) of the three cis-pQTL sources, for
# `pre_pqtl_publication_date`. They are logged below the line of PREREG.md as a pre-stage note
# before stage A runs; `run_stage_a.py` refuses to run while any value is None.
#
# PREREG §Measured variables: `pre_pqtl_publication` holds when the earliest Phase II start of the
# hypothesis precedes the online publication of every source reporting a cis-pQTL for the gene,
# which is to precede the earliest such publication. Stage A writes that earliest date
# (`earliest_publication`, over the sources whose lists hold the gene); the flag itself is computed
# in stage D (stage_d/join.py `derive`: earliest_phase2_start < pre_pqtl_publication_date, strict),
# because the start date comes from stage C.
PUBLICATION_DATES: dict[str, date | None] = {
    "ukbppp": date(2023, 10, 4),    # Sun et al. 2023, Nature 622:329-338, doi 10.1038/s41586-023-06592-6
    "decode": date(2021, 12, 2),    # Ferkingstad et al. 2021, Nat Genet 53:1712-1721, doi 10.1038/s41588-021-00978-w
    "interval": date(2018, 6, 6),   # Sun et al. 2018, Nature 558:73-79, doi 10.1038/s41586-018-0175-2
}


class PilotKeyMappingError(ValueError):
    """A pilot disease label did not map to exactly one Open Targets disease ID."""


class PublicationDatesMissing(ValueError):
    """A cis-pQTL source has no frozen publication date."""


def norm_old(s: str) -> str:
    """Key normalization of count_coverage.py, used for the pilot labels."""
    return str(s).lower().replace("'", "").replace("-", " ").strip()


def pilot_label_map(frozen_candidates: pd.DataFrame, disease: pd.DataFrame) -> dict[str, str]:
    """normalized pilot label -> OT disease ID. `frozen_candidates`: disease, ot_disease_id
    (the latter holds the Open Targets disease *name*, as in count_coverage.py)."""
    name2id: dict[str, list[str]] = {}
    for i, n in zip(disease["id"], disease["name"]):
        name2id.setdefault(str(n).lower(), []).append(i)
    out, problems = {}, {}
    for lab, name in zip(frozen_candidates["disease"], frozen_candidates["ot_disease_id"]):
        key = str(name).lower()
        if key in MS_OVERRIDE:
            out[norm_old(lab)] = MS_OVERRIDE[key]
            continue
        hits = name2id.get(key, [])
        if len(hits) == 1:
            out[norm_old(lab)] = hits[0]
        else:
            problems[lab] = hits
    if problems:
        raise PilotKeyMappingError(f"pilot labels without exactly one OT disease ID: {problems}")
    return out


def pilot_keys(classifications: list[pd.DataFrame], label_map: dict[str, str]) -> set[tuple[str, str]]:
    """(gene, indication ID) keys of classification_v5.csv and classification_v5_1.csv."""
    seen = pd.concat([c[["gene", "disease"]] for c in classifications], ignore_index=True)
    ids = seen["disease"].map(lambda x: label_map.get(norm_old(x)))
    unmapped = sorted(set(seen.loc[ids.isna(), "disease"]))
    if unmapped:
        raise PilotKeyMappingError(f"pilot classification labels absent from the label map: {unmapped}")
    return set(zip(seen["gene"], ids))


def karim_launched_keys(st17: pd.DataFrame, disease: pd.DataFrame) -> set[tuple[str, str]]:
    """(gene, OT disease ID) for every Karim et al. 2026 launched pQTL-supported pair. Karim
    indications are MeSH IDs; an OT disease matches when its dbXRefs list that MeSH ID."""
    mesh2ot: dict[str, set[str]] = {}
    for did, xrefs in zip(disease["id"], disease["dbXRefs"]):
        for x in (list(xrefs) if xrefs is not None else []):
            prefix, _, local = str(x).partition(":")
            if prefix.upper() == "MESH" and local:
                mesh2ot.setdefault(local, set()).add(did)
    return {(g, d) for g, m in zip(st17["gene"], st17["indication_mesh_id"])
            for d in mesh2ot.get(str(m).strip(), set())}


def require_publication_dates(dates: dict[str, date | None]) -> dict[str, date]:
    missing = [s for s in SOURCE_ORDER if dates.get(s) is None]
    if missing:
        raise PublicationDatesMissing(
            f"no frozen publication date for {missing}; set flags.PUBLICATION_DATES by a logged "
            "pre-stage amendment before running stage A")
    return {s: dates[s] for s in SOURCE_ORDER}


def earliest_publication(sources: list[str], dates: dict[str, date]) -> date:
    """Publication date of the earliest source reporting a cis-pQTL for the gene."""
    if not sources:
        raise ValueError("a hypothesis without an instrument source has no publication date")
    return min(dates[s] for s in sources)
