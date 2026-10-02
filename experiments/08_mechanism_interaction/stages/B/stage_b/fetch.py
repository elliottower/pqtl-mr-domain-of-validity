"""Retrieval for stage B, in two halves that share one error classification (stage_b/remote.py).

`RemoteSources` is the collect phase's view of the sources of whole files (stage_b/collect.py):
where each file is, under which name, and how to open it from a byte offset. It is the only code
that holds a deCODE folder token or asks Synapse for a download link, and neither ever leaves it:
`source_url` is the address without query string and with the token replaced by `<token>`.

`VolumeFetcher` is the Fetcher of pipeline.py for the analyze phase. Whole files are read only from
the stage B volume, through their collect records (`collect.read_record`, `collect.collected_file`);
a missing record raises CollectError and a record of `absent` raises SourceAbsent. It queries by
region what is served by region: OpenGWAS associations, tabix on FinnGen, on GWAS Catalog files that
have an index and on the eQTL Catalogue, bcftools on 1000 Genomes, Ensembl and GTEx. Each of those
raises SourceAbsent only on a definitive absence and RetryableSourceError on anything else. For
Ensembl the only definitive absence is its HTTP 400 whose JSON error names a requested rsID as
unknown (`ensembl_json`, `ensembl_unknown_id`); every other Ensembl status, 404 included, raises.

Endpoints are the ones named in PREREG §Data collection. Two are not API calls any more:
- eQTL Catalogue: the REST API is deprecated ("The RESTful API has now been deprecated and is no
  longer available", ebi.ac.uk/eqtl/Data_access, read 2026-09-30); regional queries use tabix on
  the eQTL Catalogue's own per-dataset files, rate limited as that page asks. The path table
  gives them as ftp://ftp.ebi.ac.uk/...; the same host serves the same paths over HTTPS, and both
  the file and its index are read there (`https_path`), so HTTP statuses are classified as for
  every other source. Set EQTL_CATALOGUE_ENABLED = False to apply the frozen consequence instead
  (splicing_candidate missing for every instrument).
- GTEx v8 median TPM (highest-expression tissue) comes from the GTEx Portal API v2.

DECODE_FILE_URL is the address of one file of a deCODE folder link. deCODE mails
`https://download.decode.is/folder/<token>`; its download-form script (main.a632ff70.js, sha256
7ee76c3ae286a588592bce128eed92c17975368913a2fa33d74b1f08867c880f) builds each file's link as
`/s3/download?token=<token>&file=<Key>`, Key being the file's key in the folder listing. The reply
is followed whichever way it hands over the file: a redirect to a signed address, or a JSON body
carrying one (`RemoteFile.open_at`). Size, ETag and range support are those of the final response;
the signed address is used once and is never logged or stored. Whatever the address, a file is
accepted only with the size and ETag of the pinned listing (collect.collect_one).
"""
import csv
import gzip
import io
import json
import re
import subprocess
import tarfile
import tempfile
import time
import urllib.parse
from collections.abc import Callable, Collection, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol

import numpy as np
import pandas as pd
import requests
from pydantic import BaseModel, ConfigDict

from stage_b.collect import collected_file, read_record
from stage_b.ld import parse_genotypes
from stage_b.parsers import (decode_to_canonical, filter_decode, filter_decode_annotation, filter_decode_excluded,
                             filter_eqtl_catalogue, filter_finngen, filter_gwas_catalog, filter_ukbppp,
                             normalize_chrom, opengwas_to_canonical, parse_ukbppp_rsid_map, split_lines,
                             ukbppp_to_canonical)
from stage_b.remote import ABSENT_STATUS, CORRUPT, attempt, http, http_json, status_kind
from stage_b.schemas import (Build, CollectError, CollectTask, InstrumentUnit, LDReferenceError, OutcomeSpec,
                             RetryableSourceError, Sentinel, SourceAbsent)

EQTL_CATALOGUE_ENABLED = True

SYNAPSE_UKBPPP_EUR = "syn51365303"          # European discovery, one tar per protein
SYNAPSE_RSID_MAPS = "syn51396727"           # GRCh38 position -> rsID maps, one file per chromosome
OPENGWAS_API = "https://api.opengwas.io/api"
ENSEMBL_REST = {"GRCh38": "https://rest.ensembl.org", "GRCh37": "https://grch37.rest.ensembl.org"}
GWASCAT_FTP = "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
FINNGEN_R12 = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_{endpoint}.gz"
KG_VCF = ("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/"
          "20190312_biallelic_SNV_and_INDEL/ALL.chr{chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz")
KG_PANEL = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
EQTLCAT_PATHS = ("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/"
                 "tabix/tabix_ftp_paths.tsv")
GTEX_API = "https://gtexportal.org/api/v2"
EQTLCAT_MIN_INTERVAL_S = 2.0   # eQTL Catalogue asks for at most one tabix request every couple of seconds

# GTEx v8 tissueSiteDetailId -> eQTL Catalogue GTEx sample_group (tabix_ftp_paths.tsv). GTEx tissues
# without an eQTL Catalogue dataset map to None.
GTEX_TO_EQTLCAT = {
    "Adipose_Subcutaneous": "adipose_subcutaneous", "Adipose_Visceral_Omentum": "adipose_visceral",
    "Adrenal_Gland": "adrenal_gland", "Artery_Aorta": "artery_aorta", "Artery_Coronary": "artery_coronary",
    "Artery_Tibial": "artery_tibial", "Bladder": None, "Brain_Amygdala": "brain_amygdala",
    "Brain_Anterior_cingulate_cortex_BA24": "brain_anterior_cingulate_cortex",
    "Brain_Caudate_basal_ganglia": "brain_caudate", "Brain_Cerebellar_Hemisphere": "brain_cerebellar_hemisphere",
    "Brain_Cerebellum": "brain_cerebellum", "Brain_Cortex": "brain_cortex",
    "Brain_Frontal_Cortex_BA9": "brain_frontal_cortex", "Brain_Hippocampus": "brain_hippocampus",
    "Brain_Hypothalamus": "brain_hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia": "brain_nucleus_accumbens",
    "Brain_Putamen_basal_ganglia": "brain_putamen", "Brain_Spinal_cord_cervical_c-1": "brain_spinal_cord",
    "Brain_Substantia_nigra": "brain_substantia_nigra", "Breast_Mammary_Tissue": "breast",
    "Cells_Cultured_fibroblasts": "fibroblast", "Cells_EBV-transformed_lymphocytes": "LCL",
    "Cervix_Ectocervix": None, "Cervix_Endocervix": None, "Colon_Sigmoid": "colon_sigmoid",
    "Colon_Transverse": "colon_transverse", "Esophagus_Gastroesophageal_Junction": "esophagus_gej",
    "Esophagus_Mucosa": "esophagus_mucosa", "Esophagus_Muscularis": "esophagus_muscularis",
    "Fallopian_Tube": None, "Heart_Atrial_Appendage": "heart_atrial_appendage",
    "Heart_Left_Ventricle": "heart_left_ventricle", "Kidney_Cortex": "kidney_cortex", "Kidney_Medulla": None,
    "Liver": "liver", "Lung": "lung", "Minor_Salivary_Gland": "minor_salivary_gland", "Muscle_Skeletal": "muscle",
    "Nerve_Tibial": "nerve_tibial", "Ovary": "ovary", "Pancreas": "pancreas", "Pituitary": "pituitary",
    "Prostate": "prostate", "Skin_Not_Sun_Exposed_Suprapubic": "skin_not_sun_exposed",
    "Skin_Sun_Exposed_Lower_leg": "skin_sun_exposed", "Small_Intestine_Terminal_Ileum": "small_intestine",
    "Spleen": "spleen", "Stomach": "stomach", "Testis": "testis", "Thyroid": "thyroid", "Uterus": "uterus",
    "Vagina": "vagina", "Whole_Blood": "blood",
}
FIXED_TISSUES = ("liver", "blood")   # GTEx v8 liver and whole blood


DECODE_FILE_URL = "https://download.decode.is/s3/download?token={token}&file={key}"   # token and key URL-encoded
# Ensembl REST error replies, from its documentation and source (read 2026-10-02):
# - github.com/Ensembl/ensembl-rest/wiki/HTTP-Response-Codes: 400 "Occurs during exceptional circumstances
#   such as the service is unable to find an ID ... the JSON object is an exception hash with the message
#   keyed under error"; 404 "Indicates a badly formatted request. Check your URL". An Ensembl 404 is
#   therefore a fault of the request, not an absence.
# - GET /variation/:species/:id, ensembl-rest lib/EnsEMBL/REST/Model/Variation.pm (fetch_variation):
#   `Catalyst::Exception->throw("$variation_id not found for $species")`, the species as the URL gives it.
# - POST /vep/:species/id, ensembl-rest lib/EnsEMBL/REST/Controller/VEP.pm (get_consequences) returns the
#   first VEP warning as the error when no input id gave a consequence; ensembl-vep
#   modules/Bio/EnsEMBL/VEP/Parser/ID.pm warns "No variant found with ID '$id'". Ids VEP does not know
#   are left out of a 200 reply that holds consequences for other ids.
# The pages rest.ensembl.org/documentation/info/variation_id and .../vep_id_post document parameters
# only. Both messages were also returned verbatim by rest.ensembl.org and grch37.rest.ensembl.org for an
# invented rsID on 2026-10-02: {"error":"rs99999999999999 not found for human"} and
# {"error":"No variant found with ID 'rs99999999999999'"}.
# The variation message ends with the species segment of the request URL, which is always `human` here.
ENSEMBL_UNKNOWN_VARIATION = re.compile(r"(?P<id>rs[0-9]+) not found for human")
ENSEMBL_UNKNOWN_VEP_ID = re.compile(r"No variant found with ID '(?P<id>rs[0-9]+)'")
ENSEMBL_ERROR_BYTES = 4096      # an error reply longer than this is not one of the two messages
EBI_FTP, EBI_HTTPS = "ftp://ftp.ebi.ac.uk/", "https://ftp.ebi.ac.uk/"


class Endpoints(BaseModel):
    """Where each source is reached. The defaults are the registered sources; the dry run and the
    tests point them at a local server."""
    model_config = ConfigDict(frozen=True)

    opengwas_api: str = OPENGWAS_API
    ensembl_rest: dict[str, str] = ENSEMBL_REST
    gwascat_ftp: str = GWASCAT_FTP
    finngen: str = FINNGEN_R12
    kg_vcf: str = KG_VCF
    kg_panel: str = KG_PANEL
    eqtlcat_paths: str = EQTLCAT_PATHS
    gtex_api: str = GTEX_API
    decode_file: str = DECODE_FILE_URL
    decode_smp_file: str = DECODE_FILE_URL


class SynapseLike(Protocol):
    def children(self, parent: str) -> list[dict]:
        """The files of a Synapse folder, each {"name", "id"}."""

    def file(self, entity_id: str) -> dict:
        """{"name", "size", "md5", "url"} of a file entity; `url` is a fresh pre-signed link."""


def _text_lines(raw) -> Iterator[str]:
    return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8")


def _header_and_rows(lines: Iterator[str], sep: str | None, what: str) -> tuple[list[str], Iterator[list[str]]]:
    rows = split_lines(lines, sep)
    try:
        header = next(rows)
    except StopIteration:
        raise RetryableSourceError("corrupt", f"{what}: the file has no header line") from None
    return header, rows


def _local(read: Callable[[], object], what: str):
    """A read of a file on the volume; a truncated or corrupt stream raises, it is not 'unavailable'."""
    try:
        return read()
    except CORRUPT as err:
        raise RetryableSourceError("corrupt", f"{what}: {type(err).__name__}") from None


def https_path(path: str) -> str:
    """An EBI FTP path as the same path on the same host over HTTPS; any other address unchanged."""
    return EBI_HTTPS + path[len(EBI_FTP):] if path.startswith(EBI_FTP) else path


def ensembl_unknown_id(body: bytes, message: re.Pattern, requested: Collection[str]) -> bool:
    """True only when `body` is Ensembl's error JSON, an object whose one key `error` holds a string,
    and that string is, from first to last character, the unknown-identifier message `message` naming
    one of the `requested` rsIDs. A body that does not parse, has another structure, says anything
    more or anything else (request size, syntax, a backend fault, a page not found) or names an rsID
    that was not requested is not an answer about the requested identifier."""
    try:
        doc = json.loads(body)
    except ValueError:
        return False
    if not isinstance(doc, dict) or set(doc) != {"error"} or not isinstance(doc["error"], str):
        return False
    match = message.fullmatch(doc["error"])
    return match is not None and match.group("id") in requested


def ensembl_json(method: str, url: str, what: str, unknown: re.Pattern, requested: Collection[str], **kwargs):
    """The JSON body of an Ensembl REST request for the rsIDs `requested`. HTTP 400 with the
    unknown-identifier message `unknown` naming one of them (`ensembl_unknown_id`) is Ensembl's
    definitive answer and raises SourceAbsent. Every other status that is not 2xx raises
    RetryableSourceError: any other 400 as `protocol`, and 404 or 410 as `protocol` too, because
    Ensembl answers 404 for a badly formed URL and 400, not 404, for an id it does not know."""
    def get():
        with requests.request(method, url, **kwargs) as r:
            if (r.status_code == 400 and len(r.content) <= ENSEMBL_ERROR_BYTES
                    and ensembl_unknown_id(r.content, unknown, requested)):
                raise SourceAbsent(f"{what}: HTTP 400, the id is not known to Ensembl")
            if not 200 <= r.status_code < 300:
                kind = "protocol" if r.status_code in ABSENT_STATUS else status_kind(r.status_code)
                raise RetryableSourceError(kind, f"{what}: HTTP {r.status_code}")
            return r.json()
    return attempt(get, what)


def _link_in(reply) -> str | None:
    """The first http(s) address in a JSON reply, wherever it sits."""
    if isinstance(reply, str):
        return reply if reply.startswith(("https://", "http://")) else None
    for value in reply.values() if isinstance(reply, dict) else reply if isinstance(reply, list) else ():
        found = _link_in(value)
        if found is not None:
            return found
    return None


def url_header(url: str, what: str) -> list[str]:
    """First non-'##' line of a gzipped remote file."""
    def first(lines: Iterator[str]) -> list[str]:
        for line in lines:
            if not line.startswith("##"):
                return line.rstrip("\n").split("\t")
        raise RetryableSourceError("corrupt", f"{what}: the file has no header line")

    def get() -> list[str]:
        with http("GET", url, what, stream=True, timeout=120) as r:
            return first(_text_lines(r.raw))
    return attempt(get, f"header of {what}")


def tabix_rows(url: str, chrom: str, start: int, end: int, workdir: Path, what: str) -> list[list[str]]:
    """Rows of a tabix-indexed remote file in the region; tries the chromosome names in use
    ('1', 'chr1'; for X also '23'). The file is known to exist (its header was read), so a tabix
    failure is a retryable fault."""
    names = [chrom, f"chr{chrom}"] + (["23", "chr23"] if chrom == "X" else [])
    for name in names:
        proc = attempt(lambda n=name: subprocess.run(["tabix", url, f"{n}:{max(start, 1)}-{end}"], cwd=workdir,
                                                     capture_output=True, text=True, timeout=1800, check=True),
                       f"tabix {what}")
        if proc.stdout:
            return [line.split("\t") for line in proc.stdout.splitlines() if line and not line.startswith("#")]
    return []


# ---- collect phase: where the whole files are ------------------------------------------------------

@dataclass(frozen=True)
class RemoteFile:
    """collect.Resolved: one whole file and how to open it. `url` is called at every (re)connection,
    so a link that expires during a download is replaced by a fresh one. With `link_reply`, the
    address may answer with the file (directly or through redirects, which are followed) or with a
    JSON body carrying the file's address, which is then opened; neither address is kept."""
    name: str
    source_url: str
    what: str
    url: Callable[[], str]
    indexed: bool = False
    md5: str = ""
    link_reply: bool = False

    def open_at(self, offset: int) -> requests.Response:
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        r = http("GET", self.url(), self.what, headers=headers, stream=True, timeout=600)
        if self.link_reply and "json" in r.headers.get("Content-Type", "").lower():
            with r:
                link = _link_in(r.json())
            if link is None:
                raise RetryableSourceError("protocol", f"{self.what}: a JSON reply without a download address")
            r = http("GET", link, self.what, headers=headers, stream=True, timeout=600)
        return r


def gwas_catalog_dir(base: str, accession: str) -> str:
    num = int(re.sub(r"\D", "", accession))
    lo = (num - 1) // 1000 * 1000 + 1
    return f"{base}/GCST{lo:06d}-GCST{lo + 999:06d}/{accession}/harmonised/"


class RemoteSources:
    """collect.Sources over the registered sources. `synapse` is called once, on the first UKB-PPP
    task, so a deCODE or GWAS Catalog call needs no Synapse login."""

    def __init__(self, cache_dir: Path, decode_token: str, synapse: Callable[[], SynapseLike],
                 endpoints: Endpoints = Endpoints(), decode_smp_token: str = ""):
        self.cache, self.endpoints = cache_dir, endpoints
        self.tokens = {"decode": decode_token, "decode_smp": decode_smp_token}
        self._make_synapse, self._synapse = synapse, None

    def resolve(self, task: CollectTask) -> RemoteFile:
        if task.source in ("decode", "decode_smp"):
            return self._decode(task)
        if task.source == "ukbppp":
            return self._synapse_file(SYNAPSE_UKBPPP_EUR, lambda name: f"_{task.key}_" in name, f"UKB-PPP {task.key}")
        if task.source == "ukbppp_rsid_map":
            tags = {f"_chr{task.key}_"} | ({"_chr23_"} if task.key == "X" else set())
            return self._synapse_file(SYNAPSE_RSID_MAPS, lambda name: any(t in name for t in tags),
                                      f"UKB-PPP rsID map chr{task.key}")
        return self._gwas_catalog(task.key)

    def _decode(self, task: CollectTask) -> RemoteFile:
        if not task.name:
            raise CollectError(f"{task.source} {task.key}: the task carries no file name from the pinned listing")
        token = self.tokens[task.source]
        if not token:
            raise RetryableSourceError("auth", f"{task.source} {task.key}: no folder token is set")
        template = self.endpoints.decode_file if task.source == "decode" else self.endpoints.decode_smp_file
        key = urllib.parse.quote(task.name, safe="")
        return RemoteFile(name=task.name, source_url=template.format(token="<token>", key=key), link_reply=True,
                          what=f"{task.source} {task.key}",
                          url=lambda: template.format(token=urllib.parse.quote(token, safe=""), key=key))

    def synapse(self) -> SynapseLike:
        if self._synapse is None:
            self._synapse = self._make_synapse()
        return self._synapse

    def _children(self, parent: str) -> list[dict]:
        path = self.cache / f"synapse_{parent}_children.json"
        if not path.is_file():
            kids = [{"name": k["name"], "id": k["id"]} for k in self.synapse().children(parent)]
            self.cache.mkdir(parents=True, exist_ok=True)
            tmp = path.with_name(f".{path.name}.tmp")
            tmp.write_text(json.dumps(kids))
            tmp.replace(path)
        return json.loads(path.read_text())

    def _synapse_file(self, parent: str, match: Callable[[str], bool], what: str) -> RemoteFile:
        hits = [c for c in self._children(parent) if match(c["name"])]
        if len(hits) != 1:
            raise SourceAbsent(f"{what}: the Synapse folder {parent} lists {len(hits)} files for it")
        entity = hits[0]["id"]
        meta = self.synapse().file(entity)
        return RemoteFile(name=meta["name"], source_url=f"synapse:{entity}", what=what, md5=meta.get("md5") or "",
                          url=lambda: self.synapse().file(entity)["url"])

    def _gwas_catalog(self, accession: str) -> RemoteFile:
        base = gwas_catalog_dir(self.endpoints.gwascat_ftp, accession)
        what = f"GWAS Catalog {accession}"
        with http("GET", base, f"{what} listing", timeout=120) as r:
            names = set(re.findall(r'href="([^"]+)"', r.text))
        files = sorted(n for n in names if n.endswith(".h.tsv.gz"))
        if len(files) != 1:
            raise SourceAbsent(f"{what}: {len(files)} harmonised files")
        url = base + files[0]
        return RemoteFile(name=files[0], source_url=url, what=what, url=lambda: url, indexed=f"{files[0]}.tbi" in names)


# ---- analyze phase ---------------------------------------------------------------------------------

class VolumeFetcher:
    """pipeline.Fetcher: whole files from the stage B volume, regional queries from the network."""

    def __init__(self, root: Path, opengwas_token: str, decode_annotation: Path, decode_excluded: Path,
                 endpoints: Endpoints = Endpoints()):
        self.root, self.endpoints = root, endpoints
        self.opengwas_token = opengwas_token
        self.decode_annotation = decode_annotation
        self.decode_excluded = decode_excluded
        self.cache = root / "cache"
        self.cache.mkdir(parents=True, exist_ok=True)
        self._last_eqtlcat = 0.0

    def _collected(self, source: str, key: str) -> Path:
        """The verified file of a collect task; SourceAbsent when collect recorded it absent."""
        record = read_record(self.root, source, key)
        if record.status == "absent":
            raise SourceAbsent(record.detail or f"{source} {key}: absent at the source")
        return collected_file(self.root, record)

    # ---- positions -------------------------------------------------------------------------
    def positions(self, sentinel: Sentinel) -> dict[str, int | None]:
        out: dict[str, int | None] = {"GRCh37": None, "GRCh38": None, sentinel.build: sentinel.pos}
        other: Build = "GRCh37" if sentinel.build == "GRCh38" else "GRCh38"
        if sentinel.rsid.startswith("rs"):
            try:
                doc = ensembl_json("GET", f"{self.endpoints.ensembl_rest[other]}/variation/human/{sentinel.rsid}",
                                   f"Ensembl {other} {sentinel.rsid}", ENSEMBL_UNKNOWN_VARIATION, {sentinel.rsid},
                                   headers={"Content-Type": "application/json"}, timeout=60)
            except SourceAbsent:       # Ensembl reports the rsID as unknown on this build: no position there
                doc = {}
            maps = [m for m in doc.get("mappings", [])
                    if normalize_chrom(m.get("seq_region_name", "")) == sentinel.chrom
                    and m.get("assembly_name") == other]
            out[other] = int(maps[0]["start"]) if len(maps) == 1 else None
        return out

    # ---- pQTL --------------------------------------------------------------------------------
    def _ukbppp(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        what = f"UKB-PPP {unit.assay_id}"
        tar_path = self._collected("ukbppp", unit.assay_id)
        chr_tags = {f"chr{chrom}_"} | ({"chr23_"} if chrom == "X" else set())

        def window() -> list[dict]:
            with tarfile.open(tar_path) as tf:
                members = [m for m in tf.getmembers() if any(t in Path(m.name).name for t in chr_tags)]
                if len(members) != 1:
                    raise SourceAbsent(f"{what}: the tar holds {len(members)} members for chr{chrom}")
                header, rows = _header_and_rows(_text_lines(tf.extractfile(members[0])), None, what)
                return filter_ukbppp(rows, header, chrom, center, half_width)

        def rsids(wanted: set[str]) -> dict[str, str]:
            with gzip.open(self._collected("ukbppp_rsid_map", chrom), "rt") as fh:
                header, rows = _header_and_rows(fh, "\t", f"UKB-PPP rsID map chr{chrom}")
                return parse_ukbppp_rsid_map(rows, header, wanted)

        rows = _local(window, what)
        return ukbppp_to_canonical(rows, _local(lambda: rsids({r["ID"] for r in rows}), f"UKB-PPP rsID map chr{chrom}"))

    def _decode(self, unit: InstrumentUnit, smp: bool, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        source = "decode_smp" if smp else "decode"
        what = f"{source} {unit.assay_id}"
        path = self._collected(source, unit.assay_id)

        def read(p: Path, parse: Callable):
            with gzip.open(p, "rt") as fh:
                header, rows = _header_and_rows(fh, "\t", what)
                return parse(rows, header)

        window = _local(lambda: read(path, lambda rows, header: filter_decode(rows, header, chrom, center, half_width)), what)
        names = {r["Name"] for r in window}
        ann = read(self.decode_annotation, lambda rows, header: filter_decode_annotation(rows, header, names))
        excl = read(self.decode_excluded, lambda rows, header: filter_decode_excluded(rows, header, names))
        df, _counts = decode_to_canonical(window, ann, excl)
        return df

    def _opengwas(self, dataset: str, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        what = f"OpenGWAS {dataset}"
        recs = http_json("POST", f"{self.endpoints.opengwas_api}/associations", what, timeout=600,
                         headers={"Authorization": f"Bearer {self.opengwas_token}"},
                         json={"variant": [f"{chrom}:{max(center - half_width, 1)}-{center + half_width}"],
                               "id": [dataset], "proxies": 0})
        if not isinstance(recs, list):
            raise RetryableSourceError("protocol", f"{what}: the response is not a list of associations")
        return opengwas_to_canonical(recs)

    def pqtl_region(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int, smp: bool = False) -> pd.DataFrame:
        if unit.source == "ukbppp":
            return self._ukbppp(unit, chrom, center, half_width)
        if unit.source == "decode":
            return self._decode(unit, smp, chrom, center, half_width)
        return self._opengwas(unit.pqtl_locator, chrom, center, half_width)

    # ---- outcomes ------------------------------------------------------------------------------
    def outcome_region(self, spec: OutcomeSpec, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        lo, hi = center - half_width, center + half_width
        if spec.source == "opengwas":
            return self._opengwas(spec.accession, chrom, center, half_width)
        if spec.source == "finngen":
            what = f"FinnGen {spec.accession}"
            url = self.endpoints.finngen.format(endpoint=re.sub(r"^FINNGEN_R\d+_", "", spec.accession, flags=re.IGNORECASE))
            header = url_header(url, what)       # a 404 here is the endpoint's absence
            with tempfile.TemporaryDirectory() as tmp:
                return filter_finngen(tabix_rows(url, chrom, lo, hi, Path(tmp), what), header, chrom, center, half_width)
        what = f"GWAS Catalog {spec.accession}"
        record = read_record(self.root, "gwas_catalog", spec.accession)
        if record.status == "absent":
            raise SourceAbsent(record.detail or f"{what}: absent at the source")
        if record.status == "remote_indexed":
            header = url_header(record.source_url, what)
            with tempfile.TemporaryDirectory() as tmp:
                return filter_gwas_catalog(tabix_rows(record.source_url, chrom, lo, hi, Path(tmp), what), header, chrom,
                                           center, half_width)
        path = collected_file(self.root, record)

        def read() -> pd.DataFrame:
            with gzip.open(path, "rt") as fh:
                header, rows = _header_and_rows(fh, "\t", what)
                return filter_gwas_catalog(rows, header, chrom, center, half_width)
        return _local(read, what)

    # ---- 1000 Genomes EUR ------------------------------------------------------------------------
    def _eur_samples(self) -> Path:
        path = self.cache / "1000g_eur_samples.txt"
        if not path.exists():
            def get() -> str:
                with http("GET", self.endpoints.kg_panel, "1000G panel", timeout=120) as r:
                    return r.text
            rows = list(csv.DictReader(io.StringIO(attempt(get, "1000G panel")), delimiter="\t"))
            eur = sorted(r["sample"] for r in rows if r.get("super_pop") == "EUR")
            if len(eur) != 503:
                raise LDReferenceError(f"1000G panel lists {len(eur)} EUR samples, expected 503")
            path.write_text("\n".join(eur) + "\n")
        return path

    def ld_panel(self, chrom: str, center_grch38: int, half_width: int) -> tuple[pd.DataFrame, np.ndarray]:
        samples = self._eur_samples()
        url = self.endpoints.kg_vcf.format(chrom=chrom)
        fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"
        with tempfile.TemporaryDirectory() as tmp:
            for name in (chrom, f"chr{chrom}"):
                region = f"{name}:{max(center_grch38 - half_width, 1)}-{center_grch38 + half_width}"
                cmd = (f"bcftools view -r {region} -S {samples} --force-samples -m2 -M2 -Ou '{url}' | "
                       f"bcftools query -f '{fmt}'")
                proc = attempt(lambda c=cmd: subprocess.run(["bash", "-o", "pipefail", "-c", c], cwd=tmp,
                                                            capture_output=True, text=True, timeout=3600, check=True),
                               f"1000G {region}")
                if proc.stdout:
                    return parse_genotypes(proc.stdout.splitlines())
        raise LDReferenceError(f"1000G EUR returned no records for chr{chrom}:{center_grch38}")

    # ---- VEP ---------------------------------------------------------------------------------------
    def vep(self, rsids: list[str], build: Build) -> list[dict]:
        """VEP records of `rsids`, asked in batches of 200. Ensembl leaves an id it does not know
        out of a reply that holds records for other ids, and answers HTTP 400 naming an unknown id
        only when no id of the batch gave a record; such a batch adds nothing, like the ids left
        out of another batch's reply. SourceAbsent is raised only when every batch was answered
        that way, so the outcome does not depend on where a batch boundary falls."""
        batches = [b for b in ([r for r in rsids[i:i + 200] if r.startswith("rs")] for i in range(0, len(rsids), 200)) if b]
        out, unknown = [], []
        for batch in batches:
            try:
                out += ensembl_json("POST", f"{self.endpoints.ensembl_rest[build]}/vep/human/id", f"VEP {build}",
                                    ENSEMBL_UNKNOWN_VEP_ID, set(batch), json={"ids": batch}, timeout=300,
                                    headers={"Content-Type": "application/json", "Accept": "application/json"})
            except SourceAbsent as err:
                unknown.append(str(err))
        if batches and len(unknown) == len(batches):
            raise SourceAbsent(unknown[0])
        return out

    # ---- eQTL Catalogue (GTEx v8) ------------------------------------------------------------------
    def _eqtlcat_paths(self) -> pd.DataFrame:
        path = self.cache / "eqtlcat_tabix_ftp_paths.tsv"
        if not path.exists():
            def get() -> str:
                with http("GET", self.endpoints.eqtlcat_paths, "eQTL Catalogue path table", timeout=120) as r:
                    return r.text
            path.write_text(attempt(get, "eQTL Catalogue path table"))
        return pd.read_csv(path, sep="\t", dtype=str)

    def _top_tpm_tissue(self, gene_ensembl: str) -> str:
        genes = http_json("GET", f"{self.endpoints.gtex_api}/reference/gene", f"GTEx gene {gene_ensembl}", timeout=60,
                          params={"geneId": gene_ensembl.split(".")[0], "gencodeVersion": "v26",
                                  "genomeBuild": "GRCh38/hg38"})["data"]
        if len(genes) != 1:
            raise SourceAbsent(f"GTEx: {len(genes)} gencode ids for {gene_ensembl}")
        rows = http_json("GET", f"{self.endpoints.gtex_api}/expression/medianGeneExpression",
                         f"GTEx median TPM {gene_ensembl}", timeout=60,
                         params={"gencodeId": genes[0]["gencodeId"], "datasetId": "gtex_v8"})["data"]
        if not rows:
            raise SourceAbsent(f"GTEx: no median TPM for {gene_ensembl}")
        top = max(rows, key=lambda x: (float(x["median"]), x["tissueSiteDetailId"]))["tissueSiteDetailId"]
        if top not in GTEX_TO_EQTLCAT:
            raise SourceAbsent(f"GTEx tissue {top} not in the eQTL Catalogue crosswalk")
        group = GTEX_TO_EQTLCAT[top]
        if group is None:
            raise SourceAbsent(f"highest-median-TPM tissue {top} has no eQTL Catalogue GTEx dataset")
        return group

    def qtl_regions(self, gene_ensembl: str, chrom: str, center_grch38: int, half_width: int) -> dict[str, dict[str, pd.DataFrame]]:
        if not EQTL_CATALOGUE_ENABLED:
            raise SourceAbsent("eQTL Catalogue regional queries disabled (frozen consequence applies)")
        paths = self._eqtlcat_paths()
        gtex = paths[paths["study_label"] == "GTEx"]
        tissues = sorted(set(FIXED_TISSUES) | {self._top_tpm_tissue(gene_ensembl)})
        out: dict[str, dict[str, pd.DataFrame]] = {}
        with tempfile.TemporaryDirectory() as tmp:
            for t in tissues:
                out[t] = {}
                for qm in ("ge", "leafcutter"):
                    row = gtex[(gtex["sample_group"] == t) & (gtex["quant_method"] == qm)]
                    if len(row) != 1:
                        raise SourceAbsent(f"eQTL Catalogue: {len(row)} GTEx {qm} datasets for {t}")
                    url, what = https_path(row["ftp_path"].iloc[0]), f"eQTL Catalogue GTEx {t} {qm}"
                    wait = EQTLCAT_MIN_INTERVAL_S - (time.monotonic() - self._last_eqtlcat)
                    if wait > 0:
                        time.sleep(wait)
                    header = url_header(url, what)
                    rows = tabix_rows(url, chrom, center_grch38 - half_width, center_grch38 + half_width, Path(tmp), what)
                    self._last_eqtlcat = time.monotonic()
                    out[t][qm] = filter_eqtl_catalogue(rows, header, chrom, center_grch38, half_width)
        return out
