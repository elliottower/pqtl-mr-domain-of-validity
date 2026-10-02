"""Remote retrieval for stage B (the Fetcher of pipeline.py). Runs inside the Modal image only:
it needs synapseclient, `tabix` and `bcftools` on PATH, and credentials in the environment.

Everything is streamed and filtered to the window; no full summary-statistics file is kept.
A UKB-PPP tar and an rsID map are downloaded to the container's scratch disk, read once, and
deleted. Retrieval failures after retries raise RetrievalError, which pipeline.py records as the
plan's "unavailable".

Endpoints are the ones named in PREREG §Data collection. Two are not API calls any more:
- eQTL Catalogue: the REST API is deprecated ("The RESTful API has now been deprecated and is no
  longer available", ebi.ac.uk/eqtl/Data_access, read 2026-09-30); regional queries use tabix on
  the eQTL Catalogue's own per-dataset files, rate limited as that page asks. Set
  EQTL_CATALOGUE_ENABLED = False to apply the frozen consequence instead (splicing_candidate
  missing for every instrument).
- GTEx v8 median TPM (highest-expression tissue) comes from the GTEx Portal API v2.
"""
import csv
import gzip
import io
import re
import subprocess
import tarfile
import tempfile
import time
import urllib.request
from collections.abc import Callable, Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import synapseclient
from synapseclient.core.exceptions import SynapseError

from stage_b.ld import parse_genotypes
from stage_b.parsers import (decode_to_canonical, filter_decode, filter_decode_annotation, filter_decode_excluded,
                             filter_eqtl_catalogue, filter_finngen, filter_gwas_catalog, filter_ukbppp,
                             normalize_chrom, opengwas_to_canonical, parse_ukbppp_rsid_map, split_lines,
                             ukbppp_to_canonical)
from stage_b.schemas import Build, InstrumentUnit, LDReferenceError, OutcomeSpec, RetrievalError, Sentinel

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


def retry(fn: Callable, what: str, attempts: int = 4, wait_s: float = 5.0):
    last: Exception | None = None
    for i in range(attempts):
        try:
            return fn()
        except (requests.RequestException, OSError, subprocess.SubprocessError, SynapseError) as e:
            last = e
            time.sleep(wait_s * (2 ** i))
    raise RetrievalError(f"{what}: {last}")


def _text_lines(raw) -> Iterator[str]:
    return io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding="utf-8")


def _header_and_rows(lines: Iterator[str], sep: str | None) -> tuple[list[str], Iterator[list[str]]]:
    rows = split_lines(lines, sep)
    try:
        header = next(rows)
    except StopIteration as e:
        raise RetrievalError("empty file") from e
    return header, rows


def url_header(url: str) -> list[str]:
    """First non-'##' line of a gzipped remote file (http, https or ftp)."""
    def get():
        with urllib.request.urlopen(url, timeout=120) as resp:
            for line in _text_lines(resp):
                if not line.startswith("##"):
                    return line.rstrip("\n").split("\t")
        raise RetrievalError(f"{url}: no header line")
    return retry(get, f"header of {url}")


def tabix_rows(url: str, chrom: str, start: int, end: int, workdir: Path) -> list[list[str]]:
    """Rows of a tabix-indexed remote file in the region; tries the chromosome names in use
    ('1', 'chr1'; for X also '23')."""
    names = [chrom, f"chr{chrom}"] + (["23", "chr23"] if chrom == "X" else [])
    for name in names:
        proc = retry(lambda n=name: subprocess.run(["tabix", url, f"{n}:{max(start, 1)}-{end}"], cwd=workdir,
                                                   capture_output=True, text=True, timeout=1800, check=True),
                     f"tabix {url}")
        if proc.stdout:
            return [line.split("\t") for line in proc.stdout.splitlines() if line and not line.startswith("#")]
    return []


class RemoteFetcher:
    def __init__(self, synapse_token: str, opengwas_token: str, decode_annotation: Path, decode_excluded: Path,
                 cache_dir: Path):
        self.syn = synapseclient.Synapse(silent=True)
        self.syn.login(authToken=synapse_token)
        self.opengwas_token = opengwas_token
        self.decode_annotation = decode_annotation
        self.decode_excluded = decode_excluded
        self.cache = cache_dir
        self.cache.mkdir(parents=True, exist_ok=True)
        self._last_eqtlcat = 0.0

    # ---- positions -------------------------------------------------------------------------
    def positions(self, sentinel: Sentinel) -> dict[str, int | None]:
        out: dict[str, int | None] = {"GRCh37": None, "GRCh38": None, sentinel.build: sentinel.pos}
        other: Build = "GRCh37" if sentinel.build == "GRCh38" else "GRCh38"
        if sentinel.rsid.startswith("rs"):
            def get():
                r = requests.get(f"{ENSEMBL_REST[other]}/variation/human/{sentinel.rsid}",
                                 headers={"Content-Type": "application/json"}, timeout=60)
                r.raise_for_status()
                return r.json()
            maps = [m for m in retry(get, f"Ensembl {other} {sentinel.rsid}").get("mappings", [])
                    if normalize_chrom(m.get("seq_region_name", "")) == sentinel.chrom
                    and m.get("assembly_name") == other]
            out[other] = int(maps[0]["start"]) if len(maps) == 1 else None
        return out

    # ---- pQTL --------------------------------------------------------------------------------
    def _synapse_children(self, parent: str) -> list[dict]:
        path = self.cache / f"{parent}_children.csv"
        if not path.exists():
            kids = retry(lambda: list(self.syn.getChildren(parent, includeTypes=["file"])), f"list {parent}")
            pd.DataFrame([{"name": k["name"], "id": k["id"]} for k in kids]).to_csv(path, index=False)
        return pd.read_csv(path, dtype=str).to_dict(orient="records")

    def _synapse_download(self, entity_id: str, tmp: Path) -> Path:
        ent = retry(lambda: self.syn.get(entity_id, downloadLocation=str(tmp)), f"download {entity_id}")
        return Path(ent.path)

    def _ukbppp(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        oid = unit.assay_id
        tars = [c for c in self._synapse_children(SYNAPSE_UKBPPP_EUR) if f"_{oid}_" in c["name"]]
        if len(tars) != 1:
            raise RetrievalError(f"UKB-PPP: {len(tars)} tars for {oid}")
        chr_tags = {f"chr{chrom}_"} | ({"chr23_"} if chrom == "X" else set())
        with tempfile.TemporaryDirectory() as tmp:
            tar_path = self._synapse_download(tars[0]["id"], Path(tmp))
            with tarfile.open(tar_path) as tf:
                members = [m for m in tf.getmembers() if any(t in Path(m.name).name for t in chr_tags)]
                if len(members) != 1:
                    raise RetrievalError(f"UKB-PPP tar {tars[0]['name']}: {len(members)} members for chr{chrom}")
                header, rows = _header_and_rows(_text_lines(tf.extractfile(members[0])), None)
                window = filter_ukbppp(rows, header, chrom, center, half_width)
            tar_path.unlink()
            maps = [c for c in self._synapse_children(SYNAPSE_RSID_MAPS)
                    if any(re.search(rf"_{t}", c["name"]) for t in chr_tags)]
            if len(maps) != 1:
                raise RetrievalError(f"UKB-PPP rsID map: {len(maps)} files for chr{chrom}")
            map_path = self._synapse_download(maps[0]["id"], Path(tmp))
            with gzip.open(map_path, "rt") as fh:
                header, rows = _header_and_rows(fh, "\t")
                rsids = parse_ukbppp_rsid_map(rows, header, {r["ID"] for r in window})
            map_path.unlink()
        return ukbppp_to_canonical(window, rsids)

    def _decode(self, url: str, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        def get():
            with requests.get(url, stream=True, timeout=600) as r:
                r.raise_for_status()
                header, rows = _header_and_rows(_text_lines(r.raw), "\t")
                return filter_decode(rows, header, chrom, center, half_width)
        window = retry(get, f"deCODE {url.split('?')[0]}")
        names = {r["Name"] for r in window}
        with gzip.open(self.decode_annotation, "rt") as fh:
            header, rows = _header_and_rows(fh, "\t")
            ann = filter_decode_annotation(rows, header, names)
        with gzip.open(self.decode_excluded, "rt") as fh:
            header, rows = _header_and_rows(fh, "\t")
            excl = filter_decode_excluded(rows, header, names)
        df, _counts = decode_to_canonical(window, ann, excl)
        return df

    def _opengwas(self, dataset: str, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        def get():
            r = requests.post(f"{OPENGWAS_API}/associations", timeout=600,
                              headers={"Authorization": f"Bearer {self.opengwas_token}"},
                              json={"variant": [f"{chrom}:{max(center - half_width, 1)}-{center + half_width}"],
                                    "id": [dataset], "proxies": 0})
            r.raise_for_status()
            return r.json()
        recs = retry(get, f"OpenGWAS {dataset}")
        if not isinstance(recs, list):
            raise RetrievalError(f"OpenGWAS {dataset}: unexpected response {str(recs)[:200]}")
        return opengwas_to_canonical(recs)

    def pqtl_region(self, unit: InstrumentUnit, chrom: str, center: int, half_width: int, smp: bool = False) -> pd.DataFrame:
        if unit.source == "ukbppp":
            return self._ukbppp(unit, chrom, center, half_width)
        if unit.source == "decode":
            return self._decode(unit.decode_smp_url if smp else unit.pqtl_locator, chrom, center, half_width)
        return self._opengwas(unit.pqtl_locator, chrom, center, half_width)

    # ---- outcomes ------------------------------------------------------------------------------
    def _gwas_catalog_file(self, acc: str) -> tuple[str, bool]:
        num = int(re.sub(r"\D", "", acc))
        lo = (num - 1) // 1000 * 1000 + 1
        base = f"{GWASCAT_FTP}/GCST{lo:06d}-GCST{lo + 999:06d}/{acc}/harmonised/"

        def get():
            r = requests.get(base, timeout=120)
            r.raise_for_status()
            return r.text
        names = set(re.findall(r'href="([^"]+)"', retry(get, f"GWAS Catalog listing {acc}")))
        files = sorted(n for n in names if n.endswith(".h.tsv.gz"))
        if len(files) != 1:
            raise RetrievalError(f"GWAS Catalog {acc}: {len(files)} harmonised files")
        return base + files[0], f"{files[0]}.tbi" in names

    def outcome_region(self, spec: OutcomeSpec, chrom: str, center: int, half_width: int) -> pd.DataFrame:
        lo, hi = center - half_width, center + half_width
        if spec.source == "opengwas":
            return self._opengwas(spec.accession, chrom, center, half_width)
        if spec.source == "finngen":
            url = FINNGEN_R12.format(endpoint=re.sub(r"^FINNGEN_R\d+_", "", spec.accession, flags=re.IGNORECASE))
            with tempfile.TemporaryDirectory() as tmp:
                return filter_finngen(tabix_rows(url, chrom, lo, hi, Path(tmp)), url_header(url), chrom, center, half_width)
        url, indexed = self._gwas_catalog_file(spec.accession)
        header = url_header(url)
        if indexed:
            with tempfile.TemporaryDirectory() as tmp:
                return filter_gwas_catalog(tabix_rows(url, chrom, lo, hi, Path(tmp)), header, chrom, center, half_width)

        def stream():
            with requests.get(url, stream=True, timeout=1800) as r:
                r.raise_for_status()
                h, rows = _header_and_rows(_text_lines(r.raw), "\t")
                return filter_gwas_catalog(rows, h, chrom, center, half_width)
        return retry(stream, f"GWAS Catalog {spec.accession}")

    # ---- 1000 Genomes EUR ------------------------------------------------------------------------
    def _eur_samples(self) -> Path:
        path = self.cache / "1000g_eur_samples.txt"
        if not path.exists():
            def get():
                r = requests.get(KG_PANEL, timeout=120)
                r.raise_for_status()
                return r.text
            rows = list(csv.DictReader(io.StringIO(retry(get, "1000G panel")), delimiter="\t"))
            eur = sorted(r["sample"] for r in rows if r.get("super_pop") == "EUR")
            if len(eur) != 503:
                raise LDReferenceError(f"1000G panel lists {len(eur)} EUR samples, expected 503")
            path.write_text("\n".join(eur) + "\n")
        return path

    def ld_panel(self, chrom: str, center_grch38: int, half_width: int) -> tuple[pd.DataFrame, np.ndarray]:
        samples = self._eur_samples()
        url = KG_VCF.format(chrom=chrom)
        fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n"
        with tempfile.TemporaryDirectory() as tmp:
            for name in (chrom, f"chr{chrom}"):
                region = f"{name}:{max(center_grch38 - half_width, 1)}-{center_grch38 + half_width}"
                cmd = (f"bcftools view -r {region} -S {samples} --force-samples -m2 -M2 -Ou '{url}' | "
                       f"bcftools query -f '{fmt}'")
                proc = retry(lambda c=cmd: subprocess.run(["bash", "-o", "pipefail", "-c", c], cwd=tmp,
                                                          capture_output=True, text=True, timeout=3600, check=True),
                             f"1000G {region}")
                if proc.stdout:
                    return parse_genotypes(proc.stdout.splitlines())
        raise LDReferenceError(f"1000G EUR returned no records for chr{chrom}:{center_grch38}")

    # ---- VEP ---------------------------------------------------------------------------------------
    def vep(self, rsids: list[str], build: Build) -> list[dict]:
        out = []
        for i in range(0, len(rsids), 200):
            batch = [r for r in rsids[i:i + 200] if r.startswith("rs")]

            def post(b=batch):
                r = requests.post(f"{ENSEMBL_REST[build]}/vep/human/id", json={"ids": b}, timeout=300,
                                  headers={"Content-Type": "application/json", "Accept": "application/json"})
                r.raise_for_status()
                return r.json()
            if batch:
                out += retry(post, f"VEP {build}")
        return out

    # ---- eQTL Catalogue (GTEx v8) ------------------------------------------------------------------
    def _eqtlcat_paths(self) -> pd.DataFrame:
        path = self.cache / "eqtlcat_tabix_ftp_paths.tsv"
        if not path.exists():
            def get():
                r = requests.get(EQTLCAT_PATHS, timeout=120)
                r.raise_for_status()
                return r.text
            path.write_text(retry(get, "eQTL Catalogue path table"))
        return pd.read_csv(path, sep="\t", dtype=str)

    def _top_tpm_tissue(self, gene_ensembl: str) -> str:
        def ref():
            r = requests.get(f"{GTEX_API}/reference/gene", timeout=60,
                             params={"geneId": gene_ensembl.split(".")[0], "gencodeVersion": "v26",
                                     "genomeBuild": "GRCh38/hg38"})
            r.raise_for_status()
            return r.json()["data"]
        genes = retry(ref, f"GTEx gene {gene_ensembl}")
        if len(genes) != 1:
            raise RetrievalError(f"GTEx: {len(genes)} gencode ids for {gene_ensembl}")

        def med():
            r = requests.get(f"{GTEX_API}/expression/medianGeneExpression", timeout=60,
                             params={"gencodeId": genes[0]["gencodeId"], "datasetId": "gtex_v8"})
            r.raise_for_status()
            return r.json()["data"]
        rows = retry(med, f"GTEx median TPM {gene_ensembl}")
        if not rows:
            raise RetrievalError(f"GTEx: no median TPM for {gene_ensembl}")
        top = max(rows, key=lambda x: (float(x["median"]), x["tissueSiteDetailId"]))["tissueSiteDetailId"]
        if top not in GTEX_TO_EQTLCAT:
            raise RetrievalError(f"GTEx tissue {top} not in the eQTL Catalogue crosswalk")
        group = GTEX_TO_EQTLCAT[top]
        if group is None:
            raise RetrievalError(f"highest-median-TPM tissue {top} has no eQTL Catalogue GTEx dataset")
        return group

    def qtl_regions(self, gene_ensembl: str, chrom: str, center_grch38: int, half_width: int) -> dict[str, dict[str, pd.DataFrame]]:
        if not EQTL_CATALOGUE_ENABLED:
            raise RetrievalError("eQTL Catalogue regional queries disabled (frozen consequence applies)")
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
                        raise RetrievalError(f"eQTL Catalogue: {len(row)} GTEx {qm} datasets for {t}")
                    url = row["ftp_path"].iloc[0]
                    wait = EQTLCAT_MIN_INTERVAL_S - (time.monotonic() - self._last_eqtlcat)
                    if wait > 0:
                        time.sleep(wait)
                    header = url_header(url)
                    rows = tabix_rows(url, chrom, center_grch38 - half_width, center_grch38 + half_width, Path(tmp))
                    self._last_eqtlcat = time.monotonic()
                    out[t][qm] = filter_eqtl_catalogue(rows, header, chrom, center_grch38, half_width)
        return out
