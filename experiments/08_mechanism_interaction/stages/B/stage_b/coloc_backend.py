"""Colocalization behind one interface: the R `coloc` package called through `Rscript`.

Why an Rscript subprocess and not rpy2: the plan names the R package, so R must run the
arithmetic; a subprocess exchanging JSON couples Python and R only through a file format, so
the Python wheel set and the R build can be pinned independently (rpy2 compiles against the R
headers present at install time and breaks when either side moves). An R failure surfaces as a
non-zero exit or a per-task `ok: false` with R's message, never as a crashed interpreter. One
Rscript call per instrument batches every task of that unit, so start-up cost is paid once.

Tests substitute any object with the same `run` method (the ColocBackend protocol).
"""
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Literal, Protocol

import numpy as np
import pandas as pd
from pydantic import BaseModel, ConfigDict

from stage_b.schemas import ColocBackendError

R_SCRIPT = Path(__file__).resolve().parent / "coloc_run.R"


class ColocDataset(BaseModel):
    model_config = ConfigDict(frozen=True)

    snp: list[str]
    beta: list[float]
    varbeta: list[float]
    N: float
    type: Literal["quant", "cc"]
    MAF: list[float] | None = None
    sdY: float | None = None
    s: float | None = None


class ColocTask(BaseModel):
    model_config = ConfigDict(frozen=True)

    id: str
    method: Literal["abf", "susie"]
    p1: float
    p2: float
    p12: float
    d1: ColocDataset
    d2: ColocDataset
    LD: list[list[float]] | None = None


class AbfResult(BaseModel):
    model_config = ConfigDict(frozen=True)

    pp: tuple[float, float, float, float, float]   # PP.H0 .. PP.H4
    nsnps: int
    snp_pp_h4: dict[str, float]

    def lead_variant(self) -> str:
        """The variant with the largest SNP.PP.H4 (ties: lower rsID, so the choice is stable)."""
        return min(self.snp_pp_h4.items(), key=lambda kv: (-kv[1], kv[0]))[0]


class SusieResult(BaseModel):
    model_config = ConfigDict(frozen=True)

    n_cs1: int
    n_cs2: int
    n_pairs: int
    max_pp_h4: float | None


class ColocBackend(Protocol):
    def run(self, tasks: list[ColocTask]) -> dict[str, AbfResult | SusieResult]: ...


def _maf(eaf: pd.Series) -> list[float] | None:
    m = np.minimum(eaf.to_numpy(float), 1 - eaf.to_numpy(float))
    return [float(x) for x in m] if np.all(np.isfinite(m)) and np.all(m > 0) else None


def pqtl_dataset(h: pd.DataFrame) -> ColocDataset:
    """pQTL side from a harmonized table. Effects are in SD units in all three sources (UKB-PPP
    rank-inverse-normal NPX, deCODE 'Beta (in standard deviations)', INTERVAL inverse-normal),
    so sdY = 1."""
    n = float(np.nanmedian(h["n_p"])) if h["n_p"].notna().any() else float("nan")
    if not np.isfinite(n):
        raise ColocBackendError("pQTL sample size missing for every variant")
    return ColocDataset(snp=list(h["rsid"]), beta=list(h["beta_p"].astype(float)),
                        varbeta=list((h["se_p"].astype(float)) ** 2), N=n, type="quant",
                        MAF=_maf(h["eaf_p"]), sdY=1.0)


def outcome_dataset(h: pd.DataFrame, n_case: int, n_control: int) -> ColocDataset:
    """Case-control outcome: s is the case fraction of the source record's counts."""
    n = n_case + n_control
    if n_case <= 0 or n_control <= 0:
        raise ColocBackendError(f"outcome case/control counts invalid: {n_case}/{n_control}")
    return ColocDataset(snp=list(h["rsid"]), beta=list(h["beta_o"].astype(float)),
                        varbeta=list((h["se_o"].astype(float)) ** 2), N=float(n), type="cc",
                        MAF=_maf(h["eaf_o"]) or _maf(h["eaf_p"]), s=n_case / n)


def qtl_dataset(h: pd.DataFrame) -> ColocDataset:
    """eQTL / sQTL side (eQTL Catalogue): quantitative, scale unknown, so coloc estimates sdY
    from MAF, N and varbeta."""
    n = float(np.nanmedian(h["n_o"])) if h["n_o"].notna().any() else float("nan")
    maf = _maf(h["eaf_o"])
    if not np.isfinite(n) or maf is None:
        raise ColocBackendError("QTL dataset needs N and a finite MAF for every variant")
    return ColocDataset(snp=list(h["rsid"]), beta=list(h["beta_o"].astype(float)),
                        varbeta=list((h["se_o"].astype(float)) ** 2), N=n, type="quant", MAF=maf)


class RscriptColoc:
    """Runs coloc_run.R. `r_libs` prepends a library directory (R_LIBS_USER) when given."""

    def __init__(self, rscript: str = "Rscript", r_libs: str | None = None, timeout_s: int = 3600):
        self.rscript = rscript
        self.r_libs = r_libs
        self.timeout_s = timeout_s
        self.session: dict = {}

    def run(self, tasks: list[ColocTask]) -> dict[str, AbfResult | SusieResult]:
        if not tasks:
            return {}
        env = dict(os.environ)
        if self.r_libs:
            env["R_LIBS_USER"] = self.r_libs
        with tempfile.TemporaryDirectory() as tmp:
            req, resp = Path(tmp) / "request.json", Path(tmp) / "response.json"
            req.write_text(json.dumps({"tasks": [t.model_dump(exclude_none=True) for t in tasks]}))
            proc = subprocess.run([self.rscript, str(R_SCRIPT), str(req), str(resp)], env=env,
                                  capture_output=True, text=True, timeout=self.timeout_s)
            if proc.returncode != 0 or not resp.exists():
                raise ColocBackendError(f"Rscript exited {proc.returncode}: {proc.stderr[-2000:]}")
            payload = json.loads(resp.read_text())
        self.session = payload.get("session", {})
        return parse_response(payload, {t.id for t in tasks})


def parse_response(payload: dict, expected_ids: set[str]) -> dict[str, AbfResult | SusieResult]:
    out: dict[str, AbfResult | SusieResult] = {}
    for t in payload.get("tasks", []):
        if not t.get("ok"):
            raise ColocBackendError(f"coloc task {t.get('id')} failed in R: {t.get('error')}")
        if t["method"] == "abf":
            snp = t["snp"] if isinstance(t["snp"], list) else [t["snp"]]
            h4 = t["snp_pp_h4"] if isinstance(t["snp_pp_h4"], list) else [t["snp_pp_h4"]]
            out[t["id"]] = AbfResult(pp=tuple(float(x) for x in t["pp"]), nsnps=int(t["nsnps"]),
                                     snp_pp_h4=dict(zip(snp, (float(x) for x in h4))))
        else:
            mx = t.get("max_pp_h4")
            out[t["id"]] = SusieResult(n_cs1=int(t["n_cs1"]), n_cs2=int(t["n_cs2"]), n_pairs=int(t["n_pairs"]),
                                       max_pp_h4=None if mx is None else float(mx))
    missing = expected_ids - set(out)
    if missing:
        raise ColocBackendError(f"coloc response lacks tasks {sorted(missing)}")
    return out
