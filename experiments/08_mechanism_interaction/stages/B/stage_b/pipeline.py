"""One instrument unit, start to finish, with a checkpoint after every step.

The unit is: sentinel positions -> pQTL region (±1 Mb) -> 1000G EUR LD -> each outcome region
-> each colocalization -> VEP for each lead variant -> splicing flag -> S16 (deCODE only) ->
result.json. Every step writes its file to the store and commits before the next starts; a
restarted unit reads what is on disk and resumes at the first missing step. Retrieval is behind
the Fetcher protocol (fetch.VolumeFetcher on Modal, fakes in tests) and R behind ColocBackend.

A step is recorded as unavailable only when the fetcher raises SourceAbsent, a definitive absence
at the source (stage_b/remote.py). Every other retrieval fault (RetryableSourceError: an expired
credential, a rate limit, a server error, a timeout, a broken or corrupt stream) and a missing
collect record (CollectError) propagate out of `process_unit`: nothing is written for the step, the
unit has no result.json, and a later call resumes at that step.

Resuming is allowed only under the fingerprint the directory was written with
(stage_b/checkpoint.py: the full unit record, the source pins, the code digest, the frozen plan
hash and the tool versions of the image). `DirStore.bind` writes FINGERPRINT.json (the fingerprint
and the tool versions it covers) before the first step and raises
StaleCheckpointError on a directory that holds files under another fingerprint or under none.
"""
import gzip
import hashlib
import io
import json
import os
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Protocol

import numpy as np
import pandas as pd

from stage_b.checkpoint import FINGERPRINT_NAME
from stage_b.coloc_backend import (AbfResult, ColocBackend, ColocTask, SusieResult, outcome_dataset,
                                   pqtl_dataset, qtl_dataset)
from stage_b.evidence import (coverage, genetic_direction, lead_sqtl_event, n_variants, not_run_reason,
                              sentinel_outcome_p, splicing_candidate, vep_protein_altering)
from stage_b.harmonize import harmonize
from stage_b.ld import aligned_ld, proxies
from stage_b.parsers import restrict_window
from stage_b.schemas import (PRIMARY_P1, ColocBackendError, PRIMARY_P2, PRIMARY_P12, S15A_P12, S15B_P12, VARIANT_COLUMNS,
                             WINDOW_PRIMARY, WINDOW_WIDE, Build, InstrumentUnit, LDReferenceError,
                             OutcomeSpec, Sentinel, SourceAbsent, StaleCheckpointError)

SOURCE_BUILD: dict[str, Build] = {"ukbppp": "GRCh38", "decode": "GRCh38", "interval": "GRCh37"}
OUTCOME_BUILD: dict[str, Build] = {"gwas_catalog": "GRCh38", "finngen": "GRCh38", "opengwas": "GRCh37"}


class Fetcher(Protocol):
    def positions(self, sentinel: Sentinel) -> dict[str, int | None]: ...
    def pqtl_region(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int, smp: bool = False) -> pd.DataFrame: ...
    def outcome_region(self, spec: OutcomeSpec, chrom: str, center: int, half_width: int) -> pd.DataFrame: ...
    def ld_panel(self, chrom: str, center_grch38: int, half_width: int) -> tuple[pd.DataFrame, np.ndarray]: ...
    def vep(self, rsids: list[str], build: Build) -> list[dict]: ...
    def qtl_regions(self, gene_ensembl: str, chrom: str, center_grch38: int, half_width: int) -> dict[str, dict[str, pd.DataFrame]]: ...


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


class DirStore:
    """Checkpoint directory for one unit. Writes are atomic (temp file + rename) and each is
    followed by `commit` (the Modal volume commit on Modal, a no-op locally)."""

    def __init__(self, root: Path, commit: Callable[[], None] = lambda: None):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.commit = commit

    def has(self, name: str) -> bool:
        return (self.root / name).exists()

    def bind(self, fingerprint: str, tools: Mapping[str, str]) -> None:
        """Bind the directory to `fingerprint`: write FINGERPRINT.json (with `tools`, the tool
        versions the fingerprint covers) when the directory is empty, accept a directory already
        bound to it, refuse any other state."""
        path = self.root / FINGERPRINT_NAME
        if path.exists():
            found = json.loads(path.read_text()).get("fingerprint")
            if found != fingerprint:
                raise StaleCheckpointError(f"{self.root} was written under fingerprint {found}, this run is "
                                           f"{fingerprint}; the checkpoint is not resumed")
            return
        others = sorted(p.name for p in self.root.iterdir() if p.name != f".{FINGERPRINT_NAME}.tmp")
        if others:
            raise StaleCheckpointError(f"{self.root} holds {others[:5]} and no {FINGERPRINT_NAME}; "
                                       "the checkpoint is not resumed")
        self.put_json(FINGERPRINT_NAME, {"fingerprint": fingerprint, "tools": dict(sorted(tools.items()))})

    def _write(self, name: str, data: bytes) -> str:
        tmp = self.root / f".{name}.tmp"
        tmp.write_bytes(data)
        os.replace(tmp, self.root / name)
        self.commit()
        return sha256_bytes(data)

    def put_json(self, name: str, obj) -> str:
        return self._write(name, json.dumps(obj, indent=1, sort_keys=True, default=_jsonable).encode())

    def json(self, name: str):
        return json.loads((self.root / name).read_text())

    def put_table(self, name: str, df: pd.DataFrame) -> str:
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
            gz.write(df.to_csv(sep="\t", index=False).encode())
        return self._write(name, buf.getvalue())

    def table(self, name: str) -> pd.DataFrame:
        df = pd.read_csv(self.root / name, sep="\t", dtype={"rsid": str, "chrom": str, "ea": str, "oa": str},
                         keep_default_na=False, na_values=[""])
        df["rsid"] = df["rsid"].fillna("")
        return df

    def put_npz(self, name: str, **arrays) -> str:
        buf = io.BytesIO()
        np.savez_compressed(buf, **arrays)
        return self._write(name, buf.getvalue())

    def npz(self, name: str) -> dict:
        with np.load(self.root / name, allow_pickle=False) as z:
            return {k: z[k] for k in z.files}


def _jsonable(x):
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.bool_,)):
        return bool(x)
    if isinstance(x, set):
        return sorted(x)
    raise TypeError(f"not JSON serializable: {type(x)}")


def safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def _abf(tid: str, h: pd.DataFrame, d2, p12: float) -> ColocTask:
    return ColocTask(id=tid, method="abf", p1=PRIMARY_P1, p2=PRIMARY_P2, p12=p12, d1=pqtl_dataset(h), d2=d2)


# ---- steps ---------------------------------------------------------------------------------------

def _region_step(store: DirStore, name: str, fetch: Callable[[], pd.DataFrame], window: str) -> dict:
    """Extract (or reload) one regional table; the meta file records status and sha256."""
    meta_name = f"{name}.meta.json"
    if store.has(meta_name):
        return store.json(meta_name)
    try:
        df = fetch()
        sha = store.put_table(f"{name}.tsv.gz", df[VARIANT_COLUMNS + [c for c in df.columns if c not in VARIANT_COLUMNS]])
        meta = {"status": "ok", "variants": int(len(df)), "sha256": sha, "window": window, "detail": ""}
    except SourceAbsent as e:
        meta = {"status": "unavailable", "variants": 0, "sha256": "", "window": window, "detail": str(e)}
    store.put_json(meta_name, {**meta, "retrieved_utc": pd.Timestamp.now(tz="UTC").isoformat()})
    return store.json(meta_name)


def colocalize_pair(pqtl_wide: pd.DataFrame, outcome_wide: pd.DataFrame, sentinel: Sentinel, outcome_center: int,
                    spec: OutcomeSpec, ld_meta: pd.DataFrame, dosage: np.ndarray, sentinel_proxies: set[str],
                    backend: ColocBackend) -> dict:
    """Primary coloc.abf on ±500 kb and the S15a-c, S15g variants for one instrument x outcome."""
    p500 = restrict_window(pqtl_wide, sentinel.chrom, sentinel.pos, WINDOW_PRIMARY)
    o500 = restrict_window(outcome_wide, sentinel.chrom, outcome_center, WINDOW_PRIMARY)
    h, counts = harmonize(p500, o500)
    cov = coverage(n_variants(p500), n_variants(o500), h["rsid"], sentinel.rsid, sentinel_proxies)
    rec = {"accession": spec.accession, "harmonization": counts, "coverage": cov,
           "n_pqtl_window": n_variants(p500), "n_outcome_window": n_variants(o500),
           "s17_sentinel_p": sentinel_outcome_p(outcome_wide, sentinel.rsid)}
    reason = not_run_reason(True, True, len(h))
    if reason:
        return {**rec, "coloc_run": False, "not_run_reason": reason}
    d2 = outcome_dataset(h, spec.n_case, spec.n_control)
    tasks = [_abf("primary", h, d2, PRIMARY_P12), _abf("s15a", h, d2, S15A_P12), _abf("s15b", h, d2, S15B_P12)]
    h1m, _ = harmonize(pqtl_wide, outcome_wide)
    s15c_run = not not_run_reason(True, True, len(h1m))
    if s15c_run:
        tasks.append(_abf("s15c", h1m, outcome_dataset(h1m, spec.n_case, spec.n_control), PRIMARY_P12))
    susie_note = ""
    try:
        hs, ld = aligned_ld(h, ld_meta, dosage)
        tasks.append(ColocTask(id="s15g", method="susie", p1=PRIMARY_P1, p2=PRIMARY_P2, p12=PRIMARY_P12,
                               d1=pqtl_dataset(hs), d2=outcome_dataset(hs, spec.n_case, spec.n_control),
                               LD=ld.tolist()))
    except LDReferenceError as e:
        susie_note = f"susie not run: {e}"
    res = backend.run(tasks)
    prim = res["primary"]
    if not isinstance(prim, AbfResult):
        raise ColocBackendError("primary task did not return a coloc.abf result")
    lead, gd = genetic_direction(prim, h, spec.risk_coded)
    susie = res.get("s15g")
    if isinstance(susie, SusieResult) and max(susie.n_cs1, susie.n_cs2) >= 2 and susie.max_pp_h4 is not None:
        s15g, s15g_method = susie.max_pp_h4, "coloc.susie"
    else:
        s15g, s15g_method = prim.pp[4], "coloc.abf (fewer than two credible sets in both traits, or susie not run)"
    return {**rec, "coloc_run": True, "not_run_reason": "", "pp": list(prim.pp), "nsnps": prim.nsnps,
            "lead_variant": lead, "lead_snp_pp_h4": prim.snp_pp_h4[lead], "genetic_direction": gd,
            "s15a_pp_h4": res["s15a"].pp[4], "s15b_pp_h4": res["s15b"].pp[4],
            "s15c_pp_h4": res["s15c"].pp[4] if s15c_run else None, "s15c_n_shared": int(len(h1m)),
            "s15g_pp_h4": s15g, "s15g_method": s15g_method, "susie_note": susie_note,
            "susie": susie.model_dump() if isinstance(susie, SusieResult) else None}


def splicing_step(unit: InstrumentUnit, p500: pd.DataFrame, qtl: dict[str, dict[str, pd.DataFrame]],
                  backend: ColocBackend) -> dict:
    """Per tissue: coloc.abf (primary priors) of the pQTL with the gene's lead leafcutter intron
    and with the gene's eQTL."""
    gene = unit.gene_ensembl.split(".")[0]
    tasks, per_tissue = [], {}
    for tissue, frames in sorted(qtl.items()):
        info = {"lead_event": None, "sqtl_n_shared": 0, "eqtl_n_shared": 0}
        ev = lead_sqtl_event(frames["leafcutter"], gene)
        info["lead_event"] = ev
        if ev is not None:
            sq = frames["leafcutter"][frames["leafcutter"]["molecular_trait_id"] == ev]
            hs, _ = harmonize(p500, sq[VARIANT_COLUMNS])
            info["sqtl_n_shared"] = int(len(hs))
            if not not_run_reason(True, True, len(hs)):
                tasks.append(_abf(f"sqtl::{tissue}", hs, qtl_dataset(hs), PRIMARY_P12))
        ge = frames["ge"]
        eq = ge[ge["gene_id"].astype(str).str.split(".").str[0] == gene]
        he, _ = harmonize(p500, eq[VARIANT_COLUMNS])
        info["eqtl_n_shared"] = int(len(he))
        if not not_run_reason(True, True, len(he)):
            tasks.append(_abf(f"eqtl::{tissue}", he, qtl_dataset(he), PRIMARY_P12))
        per_tissue[tissue] = info
    res = backend.run(tasks)
    sq_pp = {t: (res[f"sqtl::{t}"].pp[4] if f"sqtl::{t}" in res else None) for t in per_tissue}
    eq_pp = {t: (res[f"eqtl::{t}"].pp[4] if f"eqtl::{t}" in res else None) for t in per_tissue}
    return {"tissues": per_tissue, "sqtl_pp_h4": sq_pp, "eqtl_pp_h4": eq_pp, "query_failed": False, "detail": "",
            "splicing_candidate": splicing_candidate(sq_pp, eq_pp, False)}


def _checked_result(store: DirStore, fingerprint: str) -> dict:
    result = store.json("result.json")
    if result.get("fingerprint") != fingerprint:
        raise StaleCheckpointError(f"{store.root / 'result.json'} carries fingerprint {result.get('fingerprint')}, "
                                   f"this run is {fingerprint}")
    return result


def process_unit(unit: InstrumentUnit, fetcher: Fetcher, backend: ColocBackend, store: DirStore,
                 fingerprint: str, tools: Mapping[str, str]) -> dict:
    """`fingerprint` is checkpoint.unit_fingerprint of this unit under this run, computed with
    `tools` (checkpoint.tool_versions of the container)."""
    store.bind(fingerprint, tools)
    if store.has("result.json"):
        return _checked_result(store, fingerprint)
    s = unit.sentinel
    if not store.has("positions.json"):
        store.put_json("positions.json", fetcher.positions(s))
    pos = store.json("positions.json")
    pos38 = pos.get("GRCh38")

    pmeta = _region_step(store, "pqtl", lambda: fetcher.pqtl_region(unit, s.chrom, s.pos, WINDOW_WIDE),
                         f"{s.chrom}:{s.pos - WINDOW_WIDE}-{s.pos + WINDOW_WIDE} {SOURCE_BUILD[unit.source]}")
    if pmeta["status"] != "ok":
        result = {"unit_key": unit.unit_key, "fingerprint": fingerprint, "pqtl_available": False,
                  "pqtl_detail": pmeta["detail"], "outcomes": {}, "splicing": None, "vep": {}, "s16": {}}
        store.put_json("result.json", result)
        return result
    pqtl = store.table("pqtl.tsv.gz")
    p500 = restrict_window(pqtl, s.chrom, s.pos, WINDOW_PRIMARY)

    if pos38 is None:
        raise LDReferenceError(f"sentinel {s.rsid} has no GRCh38 position; 1000G EUR LD cannot be read")
    if not store.has("ld.npz"):
        m, d = fetcher.ld_panel(s.chrom, pos38, WINDOW_WIDE)
        store.put_npz("ld.npz", rsid=m["rsid"].to_numpy(str), ref=m["ref"].to_numpy(str),
                      alt=m["alt"].to_numpy(str), pos=m["pos"].to_numpy(int), dosage=d)
    z = store.npz("ld.npz")
    ld_meta = pd.DataFrame({"rsid": z["rsid"], "ref": z["ref"], "alt": z["alt"], "pos": z["pos"]})
    dosage = z["dosage"]
    sent_prox = proxies(ld_meta, dosage, s.rsid)

    outcomes = {}
    for spec in unit.outcomes:
        tag = safe(spec.accession)
        center = pos.get(OUTCOME_BUILD[spec.source])
        if center is None:
            ometa = {"status": "unavailable", "detail": f"sentinel not mapped in {OUTCOME_BUILD[spec.source]}"}
        else:
            ometa = _region_step(store, f"outcome__{tag}",
                                 lambda spec=spec, center=center: fetcher.outcome_region(spec, s.chrom, center, WINDOW_WIDE),
                                 f"{s.chrom}:{center - WINDOW_WIDE}-{center + WINDOW_WIDE} {OUTCOME_BUILD[spec.source]}")
        if ometa["status"] != "ok":
            outcomes[spec.accession] = {"accession": spec.accession, "coloc_run": False,
                                        "not_run_reason": "outcome_file_unavailable", "detail": ometa["detail"]}
            continue
        cname = f"coloc__{tag}.json"
        if not store.has(cname):
            store.put_json(cname, colocalize_pair(pqtl, store.table(f"outcome__{tag}.tsv.gz"), s, center, spec,
                                                  ld_meta, dosage, sent_prox, backend))
        outcomes[spec.accession] = store.json(cname)

    vep = {}
    for rec in outcomes.values():
        lead = rec.get("lead_variant")
        if not lead or lead in vep:
            continue
        vname = f"vep__{safe(lead)}.json"
        if not store.has(vname):
            rsids = sorted(proxies(ld_meta, dosage, lead) | {lead})
            try:
                hit = vep_protein_altering(fetcher.vep(rsids, SOURCE_BUILD[unit.source]), unit.gene_ensembl)
                store.put_json(vname, {"rsids": rsids, "hit": hit, "detail": ""})
            except SourceAbsent as e:
                store.put_json(vname, {"rsids": rsids, "hit": None, "detail": str(e)})
        vep[lead] = store.json(vname)

    if not store.has("splicing.json"):
        try:
            qtl = fetcher.qtl_regions(unit.gene_ensembl, s.chrom, pos38, WINDOW_PRIMARY)
            store.put_json("splicing.json", splicing_step(unit, p500, qtl, backend))
        except SourceAbsent as e:
            store.put_json("splicing.json", {"query_failed": True, "detail": str(e), "splicing_candidate": ""})
    splicing = store.json("splicing.json")

    s16 = {}
    if unit.source == "decode" and unit.smp_listing is not None:
        smeta = _region_step(store, "pqtl_smp", lambda: fetcher.pqtl_region(unit, s.chrom, s.pos, WINDOW_WIDE, smp=True),
                             f"{s.chrom}:{s.pos - WINDOW_WIDE}-{s.pos + WINDOW_WIDE} GRCh38")
        for spec in unit.outcomes:
            if smeta["status"] != "ok" or outcomes[spec.accession].get("not_run_reason") == "outcome_file_unavailable":
                continue
            tag = safe(spec.accession)
            name = f"s16__{tag}.json"
            if not store.has(name):
                smp = restrict_window(store.table("pqtl_smp.tsv.gz"), s.chrom, s.pos, WINDOW_PRIMARY)
                o500 = restrict_window(store.table(f"outcome__{tag}.tsv.gz"), s.chrom,
                                       pos[OUTCOME_BUILD[spec.source]], WINDOW_PRIMARY)
                h, _ = harmonize(smp, o500)
                if not_run_reason(True, True, len(h)):
                    store.put_json(name, {"coloc_run": False, "n_shared": int(len(h))})
                else:
                    r = backend.run([_abf("s16", h, outcome_dataset(h, spec.n_case, spec.n_control), PRIMARY_P12)])["s16"]
                    lead, gd = genetic_direction(r, h, spec.risk_coded)
                    store.put_json(name, {"coloc_run": True, "n_shared": int(len(h)), "pp_h4": r.pp[4],
                                          "lead_variant": lead, "genetic_direction": gd})
            s16[spec.accession] = store.json(name)

    result = {"unit_key": unit.unit_key, "fingerprint": fingerprint, "pqtl_available": True, "pqtl_detail": "",
              "positions": pos,
              "sentinel_proxies": sorted(sent_prox), "outcomes": outcomes, "vep": vep, "splicing": splicing,
              "s16": s16, "coloc_session": getattr(backend, "session", {})}
    store.put_json("result.json", result)
    return store.json("result.json")


def reset_unavailable(store: DirStore) -> list[str]:
    """Remove the unit's result.json and every regional step recorded as unavailable (and the
    steps that consumed them), so a rerun retries them. Only a definitive absence is recorded as
    unavailable, so this is for a source shown to hold the file after all; each use is logged
    below the line of PREREG.md."""
    removed = []
    stale = [p for p in store.root.glob("*.meta.json") if json.loads(p.read_text())["status"] != "ok"]
    for meta in stale:
        stem = meta.name[: -len(".meta.json")]
        tag = stem.split("__", 1)[1] if "__" in stem else None
        victims = [meta] + ([store.root / f"coloc__{tag}.json", store.root / f"s16__{tag}.json"] if tag else [])
        if stem in ("pqtl", "pqtl_smp"):
            victims += list(store.root.glob("s16__*.json"))
        for v in victims:
            if v.exists():
                v.unlink()
                removed.append(v.name)
    spl = store.root / "splicing.json"
    if spl.exists() and json.loads(spl.read_text()).get("query_failed"):
        spl.unlink()
        removed.append(spl.name)
    if removed and (store.root / "result.json").exists():
        (store.root / "result.json").unlink()
        removed.append("result.json")
    if removed:
        store.commit()
    return removed
