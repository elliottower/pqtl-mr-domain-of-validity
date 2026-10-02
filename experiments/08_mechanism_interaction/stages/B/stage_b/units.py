"""Hypotheses -> units of Modal work. One unit per instrument and gene (source + assay + gene),
carrying every outcome GWAS that gene's hypotheses selected. An assay that serves several genes (a
SOMAmer for a protein complex, a SeqId two genes' lists name) gives one unit per gene: the units
share the assay's sentinel, window and regional file, which the collect phase downloads once
(collect.collect_tasks), and each takes its own gene for the gene-specific steps (the VEP
consequence in the target gene and the splicing flag).

Reads only the hypotheses.csv columns in schemas.HYPOTHESIS_INPUT_COLUMNS.
"""
import re
from collections.abc import Callable, Mapping
from pathlib import Path

import pandas as pd

from stage_b.schemas import (HYPOTHESIS_INPUT_COLUMNS, PLATFORM, TABLE12_SOURCES, AmbiguousInstrumentError,
                             HypothesisInput, InputContractError, InstrumentUnit, OutcomeSpec, Sentinel, SourceFile)
from stage_b.sentinels import select_assay

TRUE = {"true", "1", "yes", "y"}
FALSE = {"false", "0", "no", "n"}


def load_hypotheses(path: Path) -> list[HypothesisInput]:
    header = pd.read_csv(path, nrows=0).columns
    missing = [c for c in HYPOTHESIS_INPUT_COLUMNS if c not in header]
    if missing:
        raise InputContractError(f"{path} lacks columns {missing}")
    df = pd.read_csv(path, usecols=HYPOTHESIS_INPUT_COLUMNS, dtype=str, keep_default_na=False)
    dup = sorted(df.loc[df["hypothesis_id"].duplicated(), "hypothesis_id"].unique())
    if dup:
        raise InputContractError(f"{path} repeats hypothesis ids {dup[:5]} ({len(dup)} in all)")
    return [HypothesisInput(**r) for r in df.to_dict(orient="records")]


def load_trait_coding(path: Path) -> dict[str, bool]:
    """TSV: outcome_accession, risk_coded. PREREG §Genetic effect direction: stage A records each
    selected accession's trait coding; hypotheses.csv (INTERFACES.md) has no column for it, so it
    arrives here as its own file."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if list(df.columns[:2]) != ["outcome_accession", "risk_coded"]:
        raise InputContractError(f"{path}: expected columns outcome_accession, risk_coded")
    out = {}
    for acc, v in zip(df["outcome_accession"], df["risk_coded"]):
        v = v.strip().lower()
        if v not in TRUE | FALSE:
            raise InputContractError(f"{path}: risk_coded '{v}' for {acc} is not a boolean")
        out[acc] = v in TRUE
    return out


def unit_key(source: str, assay_id: str, gene_ensembl: str) -> str:
    """`<source>__<assay>__<gene>`, each part with every character outside [A-Za-z0-9_.-] replaced
    by '_'. It names the unit's checkpoint directory and is never parsed back."""
    return "__".join(re.sub(r"[^A-Za-z0-9_.-]", "_", part) for part in (source, assay_id, gene_ensembl))


def _assays(ids: str) -> list[str]:
    return [a for a in ids.split(";") if a]


def build_units(hyps: list[HypothesisInput], sentinels: dict[tuple[str, str], Sentinel],
                locate: Callable[[Sentinel], str], risk_coded: dict[str, bool],
                decode_files: Mapping[str, SourceFile], smp_files: Mapping[str, SourceFile]
                ) -> tuple[list[InstrumentUnit], dict[str, str], dict[str, str], dict[str, dict[str, str]]]:
    """Returns (units, hypothesis_id -> unit_key, hypothesis_id -> reason it has no unit,
    hypothesis_id -> {source: unit_key} for UKB-PPP and deCODE).

    `locate(sentinel)` gives the pQTL locator (UKB-PPP OID, deCODE file name, INTERVAL OpenGWAS id;
    never a URL) and raises AmbiguousInstrumentError when the source cannot name exactly one
    regional file; such a hypothesis has no regional file and takes the inconclusive state (PREREG
    §Missing data). `decode_files` and `smp_files` map a SeqId to its record in the pinned deCODE
    folder listing (non-normalized, and SMP-normalized for S16); a deCODE unit carries them, so
    its fingerprint names the file by name, size and ETag.

    A unit is one instrument for one gene. The sentinel is selected by source and assay alone
    (sentinels.select_assay), so the units of an assay that serves several genes carry the same
    sentinel, locator and listing and differ in gene and outcomes.

    Descriptive table 12 compares evidence states where a gene has cis-pQTLs in both UKB-PPP and
    deCODE, so each hypothesis's outcome GWAS is also paired with its instrument in each of those
    sources it names assays for (`ukbppp_assay_ids`, `decode_assay_ids`), besides the selected
    one. The last map gives that unit per source; "" where the source names assays but no
    regional file resolves. A source without assays is absent from the map."""
    missing_coding = sorted({h.outcome_accession for h in hyps} - set(risk_coded))
    if missing_coding:
        raise InputContractError(f"no trait coding for outcome accessions {missing_coding[:10]} "
                                 f"({len(missing_coding)} in all)")
    by_unit: dict[str, dict] = {}

    def attach(h: HypothesisInput, source: str, assays: list[str]) -> str:
        s = select_assay(assays, sentinels, source)
        loc = locate(s)
        k = unit_key(s.source, s.assay_id, h.gene_ensembl)
        u = by_unit.setdefault(k, {"unit_key": k, "source": s.source, "assay_id": s.assay_id,
                                   "gene_symbol": h.gene_symbol, "gene_ensembl": h.gene_ensembl,
                                   "platform": PLATFORM[s.source], "sentinel": s, "pqtl_locator": loc,
                                   "pqtl_listing": decode_files.get(s.assay_id) if s.source == "decode" else None,
                                   "smp_listing": smp_files.get(s.assay_id) if s.source == "decode" else None,
                                   "outcomes": {}})
        if (u["source"], u["assay_id"], u["gene_ensembl"]) != (s.source, s.assay_id, h.gene_ensembl):
            raise InputContractError(f"unit key {k} names two instruments: {u['source']} {u['assay_id']} for "
                                     f"{u['gene_ensembl']}, and {s.source} {s.assay_id} for {h.gene_ensembl}")
        spec = OutcomeSpec(accession=h.outcome_accession, source=h.outcome_source, n_case=h.outcome_n_case,
                           n_control=h.outcome_n_control, risk_coded=risk_coded[h.outcome_accession])
        prev = u["outcomes"].setdefault(spec.accession, spec)
        if prev != spec:
            raise InputContractError(f"outcome {spec.accession} has inconsistent records across hypotheses")
        return k

    hyp_unit, unresolved, source_units = {}, {}, {}
    for h in hyps:
        try:
            hyp_unit[h.hypothesis_id] = attach(h, h.instrument_source, _assays(h.instrument_assay_id))
        except AmbiguousInstrumentError as e:
            unresolved[h.hypothesis_id] = f"regional_file_unavailable: {e}"
        per_source = {}
        for src in TABLE12_SOURCES:
            if src == h.instrument_source:
                per_source[src] = hyp_unit.get(h.hypothesis_id, "")
            elif _assays(getattr(h, f"{src}_assay_ids")):
                try:
                    per_source[src] = attach(h, src, _assays(getattr(h, f"{src}_assay_ids")))
                except AmbiguousInstrumentError:
                    per_source[src] = ""
        if per_source:
            source_units[h.hypothesis_id] = per_source
    units = [InstrumentUnit(**{**u, "outcomes": tuple(u["outcomes"][a] for a in sorted(u["outcomes"]))})
             for u in by_unit.values()]
    return sorted(units, key=lambda u: u.unit_key), hyp_unit, unresolved, source_units
