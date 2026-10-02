"""`plan`: hypotheses.csv and each source's instrument metadata -> the units of stage B.

Reads the two sealed stage A files and five pinned tables, and nothing else: no regional file, no
outcome statistic, no network.

    ukbppp_st9         Sun et al. 2023 supplementary workbook, sheet ST9 (UKB-PPP sentinels)
    decode_st02        Ferkingstad et al. 2021 supplementary workbook, sheet ST02 (deCODE sentinels)
    interval_st4       Sun et al. 2018 supplementary workbook, sheet ST4 (INTERVAL sentinels)
    opengwas_gwasinfo  the OpenGWAS listing (INTERVAL prot-a-* dataset ids)
    decode_listing     the listing of the deCODE proteomics folder (file name, size and ETag per SeqId)

Each is a file of the `pqtl-v8-inputs` volume pinned in modal_inputs_manifest.json; `pinned_tables`
hashes the file it is about to read and raises unless it is the pinned one. A sixth table,
`decode_smp_listing` (the folder of the SMP-normalized release, S16), is read when the manifest
pins it and is otherwise absent, in which case no unit carries an SMP file and S16 is empty.

A deCODE instrument is located by its SeqId in the folder listing. A SeqId the listing does not
name, or names more than once, has no regional file (AmbiguousInstrumentError): the hypothesis is
unresolved and takes the inconclusive state, as for any source that cannot name one file. The unit
records the file's name, size and ETag and no address: the link and token of the folder are a
credential of the collect phase (stage_b/fetch.py) and appear in no plan, unit or fingerprint.

`make_plan` returns the text of units.jsonl and the content of unit_plan.json, which records the
sha256 of every file read (`a_outputs`, `tables`) and of units.jsonl.
"""
import hashlib
import json
import re
from collections import Counter
from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import Path

from v8_manifest import sha256_file

from stage_b.collect import collect_tasks
from stage_b.launch import A_FILES
from stage_b.schemas import AmbiguousInstrumentError, InputContractError, Sentinel, SourceFile
from stage_b.sentinels import (decode_sentinels, interval_opengwas_id, interval_sentinels, load_decode_st02,
                               load_interval_st4, load_ukbppp_st9, ukbppp_sentinels)
from stage_b.units import build_units, load_hypotheses, load_trait_coding

PLAN_TABLES = {
    "ukbppp_st9": "inputs/ukbppp/sun2023_MOESM3_ESM.xlsx",
    "decode_st02": "inputs/decode/ferkingstad2021_MOESM4_ESM.xlsx",
    "interval_st4": "feasibility/v2_all_indications/inputs/interval/sun2018_MOESM4_ESM.xlsx",
    "opengwas_gwasinfo": "feasibility/v2_all_indications/inputs/opengwas/gwasinfo_all.json",
    "decode_listing": "inputs/decode/decode_proteomics_folder_listing_2026-09-30.json",
}
OPTIONAL_TABLES = {"decode_smp_listing": "inputs/decode/decode_smp_folder_listing.json"}
DECODE_KEY = re.compile(r"(\d+_\d+)_[^/]+\.txt\.gz")


def pinned_tables(inputs_manifest: Path, experiment_root: Path) -> dict[str, Path]:
    """name -> file under `experiment_root`, each hashed and equal to its pin in the manifest."""
    pins = {f["path"]: f["sha256"] for f in json.loads(inputs_manifest.read_text())["files"]}
    missing = [rel for rel in PLAN_TABLES.values() if rel not in pins]
    if missing:
        raise InputContractError(f"{inputs_manifest} pins no sha256 for {missing}")
    out = {}
    for name, rel in {**PLAN_TABLES, **{n: r for n, r in OPTIONAL_TABLES.items() if r in pins}}.items():
        got = sha256_file(experiment_root / rel)
        if got != pins[rel]:
            raise InputContractError(f"{experiment_root / rel}: sha256 {got} differs from the pin {pins[rel]}")
        out[name] = experiment_root / rel
    return out


def decode_listing_files(listing: Mapping) -> dict[str, SourceFile]:
    """SeqId -> the one file the folder listing names for it. A SeqId with several files is left
    out, so it resolves to no regional file."""
    by_seqid: dict[str, list[SourceFile]] = {}
    for f in listing["files"]:
        m = DECODE_KEY.fullmatch(str(f["Key"]))
        if m is not None:
            by_seqid.setdefault(m.group(1), []).append(
                SourceFile(name=str(f["Key"]), size=int(f["Size"]), etag=str(f.get("ETag", "")).strip('"')))
    return {seqid: files[0] for seqid, files in by_seqid.items() if len(files) == 1}


def make_plan(hypotheses: Path, trait_coding: Path, tables: Mapping[str, Path]) -> tuple[str, dict]:
    """(text of units.jsonl, content of unit_plan.json). `tables` holds the files of PLAN_TABLES
    and, optionally, of OPTIONAL_TABLES, by name."""
    absent = [name for name in PLAN_TABLES if name not in tables]
    if absent:
        raise InputContractError(f"plan needs the tables {absent}")
    hyps = load_hypotheses(hypotheses)
    coding = load_trait_coding(trait_coding)
    sents = (ukbppp_sentinels(load_ukbppp_st9(tables["ukbppp_st9"])) + decode_sentinels(load_decode_st02(tables["decode_st02"]))
             + interval_sentinels(load_interval_st4(tables["interval_st4"])))
    gwasinfo = json.loads(tables["opengwas_gwasinfo"].read_text())
    gwasinfo = list(gwasinfo.values()) if isinstance(gwasinfo, dict) else gwasinfo
    decode_files = decode_listing_files(json.loads(tables["decode_listing"].read_text()))
    smp_files = (decode_listing_files(json.loads(tables["decode_smp_listing"].read_text()))
                 if "decode_smp_listing" in tables else {})

    def locate(s: Sentinel) -> str:
        if s.source == "ukbppp":
            return s.assay_id
        if s.source == "decode":
            if s.assay_id not in decode_files:
                raise AmbiguousInstrumentError(f"the deCODE folder listing names no single file for SeqId {s.assay_id}")
            return decode_files[s.assay_id].name
        return interval_opengwas_id(s.locator, gwasinfo)

    units, hyp_unit, unresolved, source_units = build_units(hyps, {(s.source, s.assay_id): s for s in sents}, locate,
                                                            coding, decode_files, smp_files)
    units_text = "".join(u.model_dump_json() + "\n" for u in units)
    plan = {"generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"), "hypotheses": len(hyps),
            "units": len(units), "units_by_source": dict(sorted(Counter(u.source for u in units).items())),
            "collect_tasks_by_source": dict(sorted(Counter(t.source for t in collect_tasks(units)).items())),
            "a_outputs": dict(zip(A_FILES, (sha256_file(hypotheses), sha256_file(trait_coding)))),
            "tables": {name: sha256_file(path) for name, path in sorted(tables.items())},
            "units_sha256": hashlib.sha256(units_text.encode()).hexdigest(),
            "hypothesis_unit": hyp_unit, "unresolved": unresolved, "hypothesis_source_units": source_units,
            "unit_ids": {u.unit_key: [u.source, u.assay_id] for u in units}}
    return units_text, plan
