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
hash, the tool versions of the image and the collect digest, the sha256 over the collect records
of the unit's whole files). `DirStore.bind` writes FINGERPRINT.json (the fingerprint, and the tool
versions and collect digest it covers) before the first step and raises StaleCheckpointError on a
directory that holds files under another fingerprint or under none.

A source that keeps failing is never turned into "unavailable": no step counts failures, and no
code path writes an unavailable record, a result or a collect record after any number of them.

Each outcome is colocalized in the uncertainty mode its collect record names (`Fetcher.outcome_mode`,
schemas.UncertaintyMode; native_se for every source but a GWAS Catalog file the validation gave
another mode). In pvalue_coloc, harmonization keeps rows without an SE, the outcome dataset is
coloc.abf's p-value form (coloc_backend.pvalue_dataset), beta gives the direction only, and S15g is
not run with coloc.susie, which takes beta and varbeta only: its PP.H4 is the primary coloc.abf
value, with the reason in `susie_note`, as where susie is not run for want of LD.

Where the registered LD panel has no record in the unit's window (fetch.VolumeFetcher.ld_panel raises
LDUnavailable: chrX outside the pseudoautosomal regions), the unit is not failed. ld_unavailable.json
records the reason and the unit runs without LD: the sentinel has no known proxies, so
`sentinel_or_proxy_retained` is missing unless the sentinel itself is shared and `low_coverage` is
missing where it would depend on that (evidence.coverage); coloc.susie is not run, so S15g takes the
primary coloc.abf PP.H4 with the reason in `susie_note`, as for a pvalue_coloc outcome (the 2026-10-05
amendment and the plan's S15g); the VEP step asks for the lead variant alone. The primary coloc.abf, S15a-f,
S16 and S17 do not use LD and run as usual. The result carries the reason as `ld_unavailable`.

A sentinel given as several rsIDs is resolved through Ensembl (fetch.resolve_sentinel): the unit
then runs on Ensembl's name and position for it, and a regional row whose rsID is any of its aliases
(the listed ids and Ensembl's names and synonyms) is read as the sentinel (`with_sentinel_name`). A
sentinel Ensembl does not resolve (SentinelUnresolved) gives a result with `pqtl_available` false, the
registered consequence of an unavailable regional file, and the reason in `pqtl_detail`.
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
from stage_b.schemas import (PRIMARY_P1, PRIMARY_P2, PRIMARY_P12, S15A_P12, S15B_P12, VARIANT_COLUMNS,
                             WINDOW_PRIMARY, WINDOW_WIDE, Build, ColocBackendError, InstrumentUnit, LDReferenceError,
                             LDUnavailable, OutcomeSpec, Sentinel, SentinelUnresolved, SourceAbsent, StaleCheckpointError,
                             UncertaintyMode)

SOURCE_BUILD: dict[str, Build] = {"ukbppp": "GRCh38", "decode": "GRCh38", "interval": "GRCh37"}


class Fetcher(Protocol):
    def positions(self, sentinel: Sentinel) -> dict[str, int | None]: ...
    def pqtl_region(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int, smp: bool = False) -> pd.DataFrame: ...
    def outcome_region(self, spec: OutcomeSpec, chrom: str, center: int, half_width: int) -> pd.DataFrame: ...
    def outcome_build(self, spec: OutcomeSpec) -> Build: ...
    def outcome_mode(self, spec: OutcomeSpec) -> UncertaintyMode: ...
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

    def bind(self, fingerprint: str, tools: Mapping[str, str], collect_sha256: str) -> None:
        """Bind the directory to `fingerprint`: write FINGERPRINT.json (with `tools` and
        `collect_sha256`, the tool versions and the collect digest the fingerprint covers) when the
        directory is empty, accept a directory already bound to it, refuse any other state."""
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
        self.put_json(FINGERPRINT_NAME, {"fingerprint": fingerprint, "tools": dict(sorted(tools.items())),
                                         "collect_sha256": collect_sha256})

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


def with_sentinel_name(df: pd.DataFrame, aliases: set[str], name: str) -> pd.DataFrame:
    """`df` with every row whose rsID is one of the sentinel's `aliases` renamed to `name`; rows that are
    then the same variant twice (same rsID, chromosome, position and alleles, as a deCODE row listing two
    of the aliases gives) are kept once. Unchanged when `aliases` is empty."""
    if not aliases:
        return df
    hit = df["rsid"].isin(aliases).to_numpy()
    if not hit.any():
        return df
    out = df.copy()
    out.loc[hit, "rsid"] = name
    repeat = out.duplicated(subset=[c for c in ("rsid", "chrom", "pos", "ea", "oa", "ref", "alt") if c in out.columns])
    return out[~(repeat.to_numpy() & hit)].reset_index(drop=True)


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


PVALUE_COLOC_SUSIE_NOTE = ("susie not run: the outcome file gives no standard error (pvalue_coloc), and coloc.susie "
                           "takes beta and varbeta only")


def colocalize_pair(pqtl_wide: pd.DataFrame, outcome_wide: pd.DataFrame, sentinel: Sentinel, outcome_center: int,
                    spec: OutcomeSpec, ld_meta: pd.DataFrame, dosage: np.ndarray, sentinel_proxies: set[str] | None,
                    backend: ColocBackend, mode: UncertaintyMode, ld_unavailable: str = "") -> dict:
    """Primary coloc.abf on ±500 kb and the S15a-c, S15g variants for one instrument x outcome, the
    outcome read in its uncertainty `mode` (module docstring). With `ld_unavailable` (the reason the
    panel has no record in the window), `sentinel_proxies` is None and coloc.susie is not run: S15g is
    the primary coloc.abf PP.H4, with the reason in `susie_note`."""
    outcome_se = mode != "pvalue_coloc"
    p500 = restrict_window(pqtl_wide, sentinel.chrom, sentinel.pos, WINDOW_PRIMARY)
    o500 = restrict_window(outcome_wide, sentinel.chrom, outcome_center, WINDOW_PRIMARY)
    h, counts = harmonize(p500, o500, outcome_se)
    cov = coverage(n_variants(p500), n_variants(o500), h["rsid"], sentinel.rsid, sentinel_proxies)
    rec = {"accession": spec.accession, "uncertainty_mode": mode, "harmonization": counts, "coverage": cov,
           "n_pqtl_window": n_variants(p500), "n_outcome_window": n_variants(o500),
           "s17_sentinel_p": sentinel_outcome_p(outcome_wide, sentinel.rsid)}
    reason = not_run_reason(True, True, len(h))
    if reason:
        return {**rec, "coloc_run": False, "not_run_reason": reason}
    d2 = outcome_dataset(h, spec.n_case, spec.n_control, mode)
    tasks = [_abf("primary", h, d2, PRIMARY_P12), _abf("s15a", h, d2, S15A_P12), _abf("s15b", h, d2, S15B_P12)]
    h1m, _ = harmonize(pqtl_wide, outcome_wide, outcome_se)
    s15c_run = not not_run_reason(True, True, len(h1m))
    if s15c_run:
        tasks.append(_abf("s15c", h1m, outcome_dataset(h1m, spec.n_case, spec.n_control, mode), PRIMARY_P12))
    susie_note = f"susie not run: {ld_unavailable}" if ld_unavailable else "" if outcome_se else PVALUE_COLOC_SUSIE_NOTE
    if outcome_se and not ld_unavailable:
        try:
            hs, ld = aligned_ld(h, ld_meta, dosage)
            tasks.append(ColocTask(id="s15g", method="susie", p1=PRIMARY_P1, p2=PRIMARY_P2, p12=PRIMARY_P12,
                                   d1=pqtl_dataset(hs), d2=outcome_dataset(hs, spec.n_case, spec.n_control, mode),
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


def _checked_result(store: DirStore, fingerprint: str, collect_sha256: str) -> dict:
    result = store.json("result.json")
    found = (result.get("fingerprint"), result.get("collect_sha256"))
    if found != (fingerprint, collect_sha256):
        raise StaleCheckpointError(f"{store.root / 'result.json'} carries fingerprint and collect digest {found}, "
                                   f"this run is {(fingerprint, collect_sha256)}")
    return result


def process_unit(unit: InstrumentUnit, fetcher: Fetcher, backend: ColocBackend, store: DirStore,
                 fingerprint: str, tools: Mapping[str, str], collect_sha256: str) -> dict:
    """`fingerprint` is checkpoint.unit_fingerprint of this unit under this run, computed with
    `tools` (checkpoint.tool_versions of the container) and `collect_sha256`
    (checkpoint.collect_digest of the unit's collect records)."""
    store.bind(fingerprint, tools, collect_sha256)
    if store.has("result.json"):
        return _checked_result(store, fingerprint, collect_sha256)
    s = unit.sentinel
    if not store.has("positions.json"):
        try:
            store.put_json("positions.json", fetcher.positions(s))
        except SentinelUnresolved as e:
            store.put_json("positions.json", {"GRCh37": None, "GRCh38": None, "unresolved": str(e)})
    pos = store.json("positions.json")
    if pos.get("unresolved"):
        result = {"unit_key": unit.unit_key, "fingerprint": fingerprint, "collect_sha256": collect_sha256,
                  "pqtl_available": False, "pqtl_detail": f"no resolvable sentinel: {pos['unresolved']}",
                  "sentinel_unresolved": pos["unresolved"], "outcomes": {}, "splicing": None, "vep": {}, "s16": {}}
        store.put_json("result.json", result)
        return result
    resolved = pos.get("sentinel") or {}
    aliases = set(resolved.get("aliases", []))
    if resolved:
        s = s.model_copy(update={"rsid": resolved["rsid"], "pos": resolved["pos"]})
    pos38 = pos.get("GRCh38")

    def table(name: str) -> pd.DataFrame:
        return with_sentinel_name(store.table(name), aliases, s.rsid)

    pmeta = _region_step(store, "pqtl", lambda: fetcher.pqtl_region(unit, s.chrom, s.pos, WINDOW_WIDE),
                         f"{s.chrom}:{s.pos - WINDOW_WIDE}-{s.pos + WINDOW_WIDE} {SOURCE_BUILD[unit.source]}")
    if pmeta["status"] != "ok":
        result = {"unit_key": unit.unit_key, "fingerprint": fingerprint, "collect_sha256": collect_sha256,
                  "pqtl_available": False, "pqtl_detail": pmeta["detail"], "outcomes": {}, "splicing": None, "vep": {},
                  "s16": {}}
        store.put_json("result.json", result)
        return result
    pqtl = table("pqtl.tsv.gz")
    p500 = restrict_window(pqtl, s.chrom, s.pos, WINDOW_PRIMARY)

    if pos38 is None:
        raise LDReferenceError(f"sentinel {s.rsid} has no GRCh38 position; 1000G EUR LD cannot be read")
    if not store.has("ld.npz") and not store.has("ld_unavailable.json"):
        try:
            m, d = fetcher.ld_panel(s.chrom, pos38, WINDOW_WIDE)
            store.put_npz("ld.npz", rsid=m["rsid"].to_numpy(str), ref=m["ref"].to_numpy(str),
                          alt=m["alt"].to_numpy(str), pos=m["pos"].to_numpy(int), dosage=d)
        except LDUnavailable as e:
            store.put_json("ld_unavailable.json", {"reason": str(e)})
    ld_unavailable = store.json("ld_unavailable.json")["reason"] if store.has("ld_unavailable.json") else ""
    sent_prox: set[str] | None
    if ld_unavailable:
        ld_meta = pd.DataFrame({"rsid": pd.Series(dtype=str), "ref": pd.Series(dtype=str), "alt": pd.Series(dtype=str),
                                "pos": pd.Series(dtype=int)})
        dosage, sent_prox = np.zeros((0, 0)), None
    else:
        z = store.npz("ld.npz")
        rsid = np.where(np.isin(z["rsid"], sorted(aliases)), s.rsid, z["rsid"]) if aliases else z["rsid"]   # rows kept: dosage aligned
        ld_meta = pd.DataFrame({"rsid": rsid, "ref": z["ref"], "alt": z["alt"], "pos": z["pos"]})
        dosage = z["dosage"]
        sent_prox = proxies(ld_meta, dosage, s.rsid)

    outcomes = {}
    builds = {spec.accession: fetcher.outcome_build(spec) for spec in unit.outcomes}
    modes = {spec.accession: fetcher.outcome_mode(spec) for spec in unit.outcomes}
    for spec in unit.outcomes:
        tag = safe(spec.accession)
        center = pos.get(builds[spec.accession])
        if center is None:
            ometa = {"status": "unavailable", "detail": f"sentinel not mapped in {builds[spec.accession]}"}
        else:
            ometa = _region_step(store, f"outcome__{tag}",
                                 lambda spec=spec, center=center: fetcher.outcome_region(spec, s.chrom, center, WINDOW_WIDE),
                                 f"{s.chrom}:{center - WINDOW_WIDE}-{center + WINDOW_WIDE} {builds[spec.accession]}")
        if ometa["status"] != "ok":
            outcomes[spec.accession] = {"accession": spec.accession, "coloc_run": False,
                                        "not_run_reason": "outcome_file_unavailable", "detail": ometa["detail"]}
            continue
        cname = f"coloc__{tag}.json"
        if not store.has(cname):
            store.put_json(cname, colocalize_pair(pqtl, table(f"outcome__{tag}.tsv.gz"), s, center, spec,
                                                  ld_meta, dosage, sent_prox, backend, modes[spec.accession], ld_unavailable))
        outcomes[spec.accession] = store.json(cname)

    vep = {}
    for rec in outcomes.values():
        lead = rec.get("lead_variant")
        if not lead or lead in vep:
            continue
        vname = f"vep__{safe(lead)}.json"
        if not store.has(vname):
            rsids = [lead] if ld_unavailable else sorted(proxies(ld_meta, dosage, lead) | {lead})
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
                smp = restrict_window(table("pqtl_smp.tsv.gz"), s.chrom, s.pos, WINDOW_PRIMARY)
                o500 = restrict_window(table(f"outcome__{tag}.tsv.gz"), s.chrom,
                                       pos[builds[spec.accession]], WINDOW_PRIMARY)
                mode = modes[spec.accession]
                h, _ = harmonize(smp, o500, mode != "pvalue_coloc")
                if not_run_reason(True, True, len(h)):
                    store.put_json(name, {"coloc_run": False, "n_shared": int(len(h))})
                else:
                    r = backend.run([_abf("s16", h, outcome_dataset(h, spec.n_case, spec.n_control, mode),
                                          PRIMARY_P12)])["s16"]
                    lead, gd = genetic_direction(r, h, spec.risk_coded)
                    store.put_json(name, {"coloc_run": True, "n_shared": int(len(h)), "pp_h4": r.pp[4],
                                          "lead_variant": lead, "genetic_direction": gd})
            s16[spec.accession] = store.json(name)

    result = {"unit_key": unit.unit_key, "fingerprint": fingerprint, "collect_sha256": collect_sha256,
              "pqtl_available": True, "pqtl_detail": "", "positions": pos, "ld_unavailable": ld_unavailable,
              "sentinel_proxies": None if sent_prox is None else sorted(sent_prox), "outcomes": outcomes, "vep": vep, "splicing": splicing,
              "s16": s16, "coloc_session": getattr(backend, "session", {})}
    store.put_json("result.json", result)
    return store.json("result.json")


SUPERSEDED_UNITS = Path("superseded") / "units"
SUPERSEDED_UNIT_ERRORS = Path("superseded") / "errors" / "units"


def bound_fingerprint(unit_dir: Path) -> str | None:
    """The fingerprint FINGERPRINT.json binds `unit_dir` to; None when it has none."""
    path = unit_dir / FINGERPRINT_NAME
    return json.loads(path.read_text()).get("fingerprint") if path.is_file() else None


def supersede_units(root: Path, current: Mapping[str, str], stamp: str, commit: Callable[[], None] = lambda: None) -> dict:
    """Archive the unit directories of `current` (unit key -> the fingerprint this run computes for the
    unit) that `run_unit` would refuse as stale, so a new call starts them afresh. Nothing is deleted.

    Refused, with nothing moved, when any of the units holds result.json: a finished unit is never
    superseded. A directory bound to its current fingerprint is kept (a new call resumes it), and so is
    its error marker. A directory bound to another fingerprint, or holding files and no FINGERPRINT.json,
    is moved to <root>/superseded/units/<stamp>/<unit key>/, and the unit's error marker
    (status.error_marker) to <root>/superseded/errors/units/<stamp>/<unit key>.json; the marker of a
    unit without a directory is archived the same way. Returns the unit keys by what was done."""
    units = root / "units"
    finished = sorted(k for k in current if (units / k / "result.json").exists())
    if finished:
        raise StaleCheckpointError(f"{len(finished)} units hold result.json (first: {finished[0]}); a finished unit is "
                                   "never superseded, and nothing was moved")
    dest_units, dest_errors = root / SUPERSEDED_UNITS / stamp, root / SUPERSEDED_UNIT_ERRORS / stamp
    taken = [p for p in (dest_units, dest_errors) if p.exists()]
    if taken:
        raise StaleCheckpointError(f"{taken[0]} already exists; a superseded directory is never overwritten")
    out: dict[str, list[str]] = {"moved": [], "kept_current": [], "markers_archived": [], "untouched": []}
    for key, fingerprint in sorted(current.items()):
        unit_dir, marker = units / key, root / "errors" / "units" / f"{key}.json"
        if unit_dir.is_dir() and any(unit_dir.iterdir()) and bound_fingerprint(unit_dir) == fingerprint:
            out["kept_current"].append(key)
            continue
        if unit_dir.is_dir():
            dest_units.mkdir(parents=True, exist_ok=True)
            os.replace(unit_dir, dest_units / key)
            out["moved"].append(key)
        if marker.is_file():
            dest_errors.mkdir(parents=True, exist_ok=True)
            os.replace(marker, dest_errors / marker.name)
            out["markers_archived"].append(key)
        if key not in out["moved"] and key not in out["markers_archived"]:
            out["untouched"].append(key)
        commit()
    return out


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
