"""Unit results -> evidence.csv, regional_manifest.tsv, MANIFEST.tsv (INTERFACES.md, stage B).

Evidence is computed per hypothesis from its instrument unit's result for its selected outcome
accession and its own intervention direction. Eldjarn ST29 supplies the platform-concordance
flag and the platform PAV half of protein_altering. `s17_sentinel_p` is the outcome p-value at
the sentinel (S17); `evidence_state_ukbppp` / `evidence_state_decode` are the primary evidence
rule applied with each source's instrument (descriptive table 12), empty where the hypothesis
names no assay in that source. `s16_evidence_state` is the evidence rule applied to a deCODE
instrument's SMP-normalized statistics (S16) and is empty wherever no such colocalization exists:
for every instrument of another source, and for every hypothesis when no SMP-normalized listing is
pinned.

`collected_files.tsv` lists every record of the collect phase (stage_b/collect.py): each whole file
with its size, sha256 and ETag, and each file recorded absent or queried by region. Its `source_url`
column holds the host of the record's address only (`source_host`), for every source; the record on
the volume keeps the address it has.

The output directory holds tables of hashes and derived summary results only. No regional extract
(a unit directory's `*.tsv.gz` or `ld.npz`, which hold rows of the source files) and no whole file
is written there: `regional_manifest.tsv` and `collected_files.tsv` carry their sha256 values.

`write_outputs` also writes run_info.json (run token, repository commit, code digest, the tool
versions the units ran under) and puts the code digest and the commit in INPUTS.tsv as the
`stage_code` row.
"""
import csv
import json
import urllib.parse
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Literal

import openpyxl
from v8_manifest import (RUN_INFO_NAME, InputRecord, code_record, code_sha256, guard_files, portable_path,
                         relative_files, sha256_file, write_manifest)

from stage_b.checkpoint import FINGERPRINT_NAME
from stage_b.evidence import (evidence_score, evidence_state, platform_concordance, protein_altering,
                              st29_pav)
from stage_b.pipeline import safe
from stage_b.schemas import (COLLECTED_FILES_COLUMNS, EVIDENCE_COLUMNS, REGIONAL_MANIFEST_COLUMNS, TABLE12_SOURCES,
                             CollectRecord, EvidenceRow, EvidenceState, HypothesisInput, InputContractError, InstrumentUnit,
                             RegionalManifestRow, StaleCheckpointError)

ST29_SHEET = "ST29_protein_classification"


def load_st29(path: Path) -> dict[tuple[str, str], dict]:
    """Eldjarn 2023 ST29 only. Keys: ('olink', OlinkID) and ('soma', SeqId)."""
    ws = openpyxl.load_workbook(path, read_only=True)[ST29_SHEET]
    rows = ws.iter_rows(min_row=3, values_only=True)
    header = [str(c) for c in next(rows)]
    out: dict[tuple[str, str], dict] = {}
    for r in rows:
        rec = dict(zip(header, r))
        if rec.get("OlinkID"):
            out.setdefault(("olink", str(rec["OlinkID"])), rec)
        if rec.get("SeqId"):
            out.setdefault(("soma", str(rec["SeqId"])), rec)
    return out


def st29_key(source: str, assay_id: str) -> tuple[str, str] | None:
    """UKB-PPP OID -> OlinkID; deCODE SeqId as is; INTERVAL SOMAmer 'GENE.14151.4.3' -> '14151_4'."""
    if source == "ukbppp":
        return ("olink", assay_id)
    if source == "decode":
        return ("soma", assay_id)
    parts = assay_id.split(".")
    if len(parts) >= 3 and parts[-3].isdigit() and parts[-2].isdigit():
        return ("soma", f"{parts[-3]}_{parts[-2]}")
    return None


def _blank_row(h: HypothesisInput, reason: str, pa, pc, spl, s16: str) -> dict:
    return {"hypothesis_id": h.hypothesis_id, "coloc_run": False, "not_run_reason": reason,
            **{f"pp_h{i}": "" for i in range(5)}, "n_shared": "", "frac_pqtl_retained": "",
            "frac_outcome_retained": "", "sentinel_or_proxy_retained": "", "low_coverage": "",
            "lead_variant": "", "genetic_direction": 0, "evidence_state": "inconclusive", "S": 0, "E": 0.0,
            "protein_altering": pa, "platform_concordant": pc, "splicing_candidate": spl,
            **{f"s15{x}_pp_h4": "" for x in "abcdefg"}, "s15f_low_coverage_excluded": "", "s16_evidence_state": s16,
            "s17_sentinel_p": "", "evidence_state_ukbppp": "", "evidence_state_decode": ""}


def s16_state(h: HypothesisInput, source: str, result: Mapping) -> EvidenceState | Literal[""]:
    """The evidence rule on the SMP-normalized colocalization of a deCODE instrument (S16); ""
    where the unit holds none for the hypothesis's outcome, which is every unit of another source
    and every deCODE unit planned without an SMP-normalized file."""
    s = result["s16"].get(h.outcome_accession) if source == "decode" else None
    if s is None:
        return ""
    return evidence_state(s["coloc_run"], s.get("pp_h4"), h.direction, s.get("genetic_direction", 0))


def evidence_row(h: HypothesisInput, source: str, assay_id: str, result: Mapping | None,
                 st29: Mapping[tuple[str, str], dict]) -> EvidenceRow:
    key = st29_key(source, assay_id) if assay_id else None
    st = st29.get(key) if key else None
    pc = platform_concordance(st)
    pav = st29_pav(st, h.platform)
    if result is None or not result["pqtl_available"]:
        return EvidenceRow(**_blank_row(h, "regional_file_unavailable", protein_altering(None, pav), pc, "", ""))
    spl = result["splicing"]["splicing_candidate"] if result.get("splicing") else ""
    o = result["outcomes"].get(h.outcome_accession)
    if o is None:
        raise InputContractError(f"unit {result['unit_key']} has no record for {h.outcome_accession}")
    if o["not_run_reason"] == "outcome_file_unavailable":
        return EvidenceRow(**_blank_row(h, "outcome_file_unavailable", protein_altering(None, pav), pc, spl, ""))
    cov = o["coverage"]
    base = _blank_row(h, o["not_run_reason"], protein_altering(None, pav), pc, spl, s16_state(h, source, result))
    p17 = o.get("s17_sentinel_p")
    base.update({"s17_sentinel_p": "" if p17 is None else p17, "n_shared": cov["n_shared"], "frac_pqtl_retained": cov["frac_pqtl_retained"],
                 "frac_outcome_retained": cov["frac_outcome_retained"],
                 "sentinel_or_proxy_retained": cov["sentinel_or_proxy_retained"], "low_coverage": cov["low_coverage"],
                 "s15f_low_coverage_excluded": cov["low_coverage"]})
    if not o["coloc_run"]:
        return EvidenceRow(**base)
    pp = o["pp"]
    gd = int(o["genetic_direction"])
    state = evidence_state(True, pp[4], h.direction, gd)
    vep = result["vep"].get(o["lead_variant"], {})
    base.update({"coloc_run": True, **{f"pp_h{i}": pp[i] for i in range(5)}, "lead_variant": o["lead_variant"],
                 "genetic_direction": gd, "evidence_state": state, "S": int(state == "supportive"),
                 "E": evidence_score(True, pp[4], h.direction, gd),
                 "protein_altering": protein_altering(vep.get("hit"), pav),
                 "s15a_pp_h4": o["s15a_pp_h4"], "s15b_pp_h4": o["s15b_pp_h4"],
                 "s15c_pp_h4": "" if o["s15c_pp_h4"] is None else o["s15c_pp_h4"],
                 "s15d_pp_h4": pp[4], "s15e_pp_h4": pp[4], "s15f_pp_h4": pp[4], "s15g_pp_h4": o["s15g_pp_h4"]})
    return EvidenceRow(**base)


def source_state(h: HypothesisInput, k: str, results: Mapping[str, dict]) -> EvidenceState:
    """The primary evidence rule for `h` with the instrument of unit `k` ("" = no resolvable
    regional file, which is inconclusive as in the primary state)."""
    if not k:
        return "inconclusive"
    if k not in results:
        raise InputContractError(f"unit {k} has no result; stage B is not complete")
    r = results[k]
    if not r["pqtl_available"]:
        return "inconclusive"
    o = r["outcomes"].get(h.outcome_accession)
    if o is None:
        raise InputContractError(f"unit {k} has no record for {h.outcome_accession}")
    if not o["coloc_run"]:
        return "inconclusive"
    return evidence_state(True, o["pp"][4], h.direction, int(o["genetic_direction"]))


def build_evidence(hyps: list[HypothesisInput], hyp_unit: Mapping[str, str], unit_ids: Mapping[str, tuple[str, str]],
                   results: Mapping[str, dict], st29: Mapping[tuple[str, str], dict],
                   source_units: Mapping[str, Mapping[str, str]] | None = None) -> list[EvidenceRow]:
    """One row per hypothesis in hypotheses.csv, in its order. A hypothesis without a unit (no
    sentinel or no resolvable regional file) is regional_file_unavailable. `source_units` is the
    fourth value of units.build_units (hypothesis -> {ukbppp|decode: unit_key})."""
    rows = []
    for h in hyps:
        k = hyp_unit.get(h.hypothesis_id)
        if k is None:
            row = evidence_row(h, h.instrument_source, h.instrument_assay_id.split(";")[0], None, st29)
        elif k not in results:
            raise InputContractError(f"unit {k} has no result; stage B is not complete")
        else:
            source, assay = unit_ids[k]
            row = evidence_row(h, source, assay, results[k], st29)
        per = (source_units or {}).get(h.hypothesis_id, {})
        row = row.model_copy(update={f"evidence_state_{src}": source_state(h, per[src], results)
                                     for src in TABLE12_SOURCES if src in per})
        rows.append(EvidenceRow.model_validate(row.model_dump()))
    return rows


# ---- files ---------------------------------------------------------------------------------------


def _fmt(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, float):
        return repr(v)
    return str(v)


def write_table(path: Path, columns: list[str], rows: list[dict], sep: str) -> int:
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter=sep, lineterminator="\n")
        w.writerow(columns)
        for r in rows:
            w.writerow([_fmt(r[c]) for c in columns])
    return len(rows)


def source_host(address: str) -> str:
    """What collected_files.tsv publishes of a collect record's `source_url`: the host of an address
    that has one (`https://download.decode.is/s3/download?...` -> `download.decode.is`), without
    port, path or query; the scheme of an address without a host (`synapse:syn1` -> `synapse`); ""
    for a record without an address."""
    parts = urllib.parse.urlsplit(address)
    return parts.hostname or parts.scheme


def regional_rows(unit_metas: Mapping[str, list[dict]]) -> list[RegionalManifestRow]:
    """unit_key -> list of {source, protein_or_study, window, variants, sha256, status, detail,
    retrieved_utc} as collected from each unit's *.meta.json."""
    return [RegionalManifestRow(**m) for k in sorted(unit_metas) for m in unit_metas[k]]


def write_outputs(out_dir: Path, evidence: list[EvidenceRow], regional: list[RegionalManifestRow],
                  script_paths: list[Path], input_paths: Mapping[str, Path], roots: Sequence[tuple[str, Path]], *,
                  script_root: Path, run_token: str, repo_commit: str, tools: Mapping[str, str],
                  collected: Sequence[CollectRecord] = ()) -> Path:
    """evidence.csv, regional_manifest.tsv and collected_files.tsv first, then run_info.json, INPUTS.tsv (each input,
    path relative to `roots`, and the `stage_code` row) and MANIFEST.tsv in the shared format
    (stages/run_guard/v8_manifest.py). `script_paths` lie under `script_root` (the stage B
    directory) and are hashed by their path relative to it, with the shared guard modules."""
    out_dir.mkdir(parents=True, exist_ok=True)
    script = code_sha256({**relative_files(script_root, script_paths), **guard_files()})
    files = {
        out_dir / "evidence.csv": write_table(out_dir / "evidence.csv", EVIDENCE_COLUMNS,
                                              [r.model_dump() for r in evidence], ","),
        out_dir / "regional_manifest.tsv": write_table(out_dir / "regional_manifest.tsv", REGIONAL_MANIFEST_COLUMNS,
                                                       [r.model_dump() for r in regional], "\t"),
        out_dir / "collected_files.tsv": write_table(
            out_dir / "collected_files.tsv", COLLECTED_FILES_COLUMNS,
            [{**r.model_dump(), "bytes": "" if r.bytes is None else r.bytes, "source_url": source_host(r.source_url)}
             for r in sorted(collected, key=lambda r: (r.source, r.key))], "\t"),
    }
    run_info = out_dir / RUN_INFO_NAME
    run_info.write_text(json.dumps({"stage": "B", "run_token": run_token, "repo_commit": repo_commit,
                                    "code_sha256": script, "tools": dict(sorted(tools.items()))}, indent=1, sort_keys=True))
    inputs = [InputRecord(name=k, path=portable_path(p, roots), sha256=sha256_file(p))
              for k, p in sorted(input_paths.items())]
    inputs.append(code_record("B", script, repo_commit))
    return write_manifest(out_dir, [*files, run_info], script, inputs)


def collect_unit_dir(unit: InstrumentUnit, unit_dir: Path, fingerprint: str, collect_sha256: str) -> tuple[dict, list[dict]]:
    """(result.json, regional-manifest rows) from one unit's checkpoint directory, as copied
    from the volume (`modal volume get pqtl-v8-stage-b /stage_b/units <dir>`). The directory's
    FINGERPRINT.json and its result.json must both carry `fingerprint`
    (checkpoint.unit_fingerprint of this unit under this run) and `collect_sha256`, the collect
    digest that fingerprint was computed with (checkpoint.collect_digest, recomputed by the caller
    from the collect records it holds)."""
    res_path = unit_dir / "result.json"
    if not res_path.exists():
        raise InputContractError(f"unit {unit.unit_key} has no result.json; stage B is not complete")
    result = json.loads(res_path.read_text())
    path = unit_dir / FINGERPRINT_NAME
    bound = json.loads(path.read_text()) if path.exists() else {}
    found = [(d.get("fingerprint"), d.get("collect_sha256")) for d in (bound, result)]
    if found != [(fingerprint, collect_sha256)] * 2:
        raise StaleCheckpointError(f"unit {unit.unit_key}: {unit_dir} carries fingerprints and collect digests {found}, "
                                   f"this run is {(fingerprint, collect_sha256)}")
    names = {"pqtl": (unit.source, unit.assay_id), "pqtl_smp": (f"{unit.source}_smp", unit.assay_id)}
    names.update({f"outcome__{safe(o.accession)}": (o.source, o.accession) for o in unit.outcomes})
    metas = []
    for stem, (src, what) in sorted(names.items()):
        p = unit_dir / f"{stem}.meta.json"
        if p.exists():
            m = json.loads(p.read_text())
            metas.append({"source": src, "protein_or_study": what, "window": m["window"], "variants": m["variants"],
                          "sha256": m["sha256"], "status": m["status"], "detail": m["detail"],
                          "retrieved_utc": m["retrieved_utc"]})
    return result, metas


def load_collect_records(collect_dir: Path) -> list[CollectRecord]:
    """Every record under a copy of the volume's collect directory
    (`modal volume get pqtl-v8-stage-b /stage_b/collect <dir>`), sorted by source and key."""
    records = [CollectRecord.model_validate_json(p.read_text()) for p in sorted(collect_dir.rglob("*.json"))]
    return sorted(records, key=lambda r: (r.source, r.key))
