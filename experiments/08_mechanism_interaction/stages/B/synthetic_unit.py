"""Synthetic sources with known answers, for the stage B dry run and the end-to-end tests. No study
data and no real source: every file is generated here and served by fake_remote.FakeRemote on
127.0.0.1, in each source's own format, so the real collect and analyze code (stage_b/collect.py,
stage_b/fetch.py, tabix, bcftools, the R backend) runs unchanged and only the addresses differ.

The region: 240 biallelic variants 5 kb apart on chromosome 1, 503 EUR samples (and 20 others, so
the sample subset is exercised) whose haplotypes are a thresholded AR(1) Gaussian sequence (allele
frequency 0.3), so LD decays with distance. Every z-score vector is a combination of columns of the
EUR panel's own correlation matrix, so summary statistics and LD reference agree.

Three units, one per instrument source, with the same pQTL signals (z 12 at the sentinel, z 7
200 kb upstream):

    decode__0_0__<gene>            a deCODE folder file (collected through the download endpoint's redirect),
                                   its SMP-normalized release (S16; that endpoint answers with a JSON link),
                                   the annotation and excluded-variant files; one variant is on the excluded list
    ukbppp__OID00000__<gene>       a UKB-PPP tar and the chromosome's rsID map (collected through pre-signed links)
    interval__SYNTH.1.2.3__<gene>  an OpenGWAS regional query, GRCh37 sentinel

each for the gene GENE (unit keys are `<source>__<assay>__<gene>`, stage_b.units.unit_key),

and these outcomes:

    FINNGEN_R12_SYN_SHARED  FinnGen (tabix), signal at the sentinel           -> PP.H4, supportive
    ieu-b-9001              OpenGWAS, GRCh37, signal at the sentinel          -> PP.H4, supportive
    GCST90000001            GWAS Catalog without an index (collected once, read by two units),
                            signal 200 kb downstream                          -> PP.H3, inconclusive
    GCST90000002            GWAS Catalog with an index (tabix), 30 variants   -> fewer_than_50_shared
    GCST90000003            no study directory (404)                          -> outcome_file_unavailable
    GCST90000004            a study directory with neither a harmonised file nor
                            a summary-statistics file                         -> outcome_file_unavailable

The splicing step finds a liver sQTL signal at the sentinel and no eQTL signal (splicing_candidate).

The deCODE assay also serves a second gene, SECOND_GENE (`second_gene_unit`): a unit of its own with
the same sentinel, window and collected file. For that gene the eQTL Catalogue holds introns and no
expression trait, and VEP reports no consequence, so its unit has the same colocalization, no VEP hit
and a missing splicing flag.

`scenario` is the dry run: collect under injected faults (broken connections, an expired folder
token, a rate limit), a re-issued folder link between collect and analyze, an expired OpenGWAS token
during analyze, the checks that none of them changed a record, a checkpoint or a result, Ensembl's
error replies (`ensembl_outcomes`: only its unknown-identifier message naming the requested rsID is
an absence), a source that serves other bytes under an unchanged name, size and ETag (the finished
checkpoint is refused), and the purge of the whole files once every unit is finished.

The synthetic Ensembl answers an rsID it does not hold as rest.ensembl.org does: HTTP 400 with
`{"error": "<rsID> not found for human"}` on the variation endpoint and
`{"error": "No variant found with ID '<rsID>'"}` on the VEP endpoint (stage_b/fetch.py cites both).
"""
import gzip
import io
import json
import math
import shutil
import subprocess
import tarfile
import urllib.parse
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from fake_remote import FakeRemote, FakeSynapse
from stage_b import collect as collect_module
from stage_b import fetch as fetch_module
from stage_b import remote as remote_errors
from stage_b.assemble import build_evidence, collect_unit_dir, regional_rows
from stage_b.checkpoint import volume_collect_digest
from stage_b.collect import collect_one, collect_tasks, purge_raw, read_record
from stage_b.fetch import (SYNAPSE_RSID_MAPS, SYNAPSE_UKBPPP_EUR, Endpoints, RemoteSources, VolumeFetcher, gwas_catalog_dir,
                           gwas_catalog_study_dir)
from stage_b.schemas import (WINDOW_PRIMARY, CollectTask, HypothesisInput, InstrumentUnit, OutcomeSpec,
                             RetryableSourceError, Sentinel, SourceAbsent, SourceFile, StageBError)
from stage_b.status import error_of, status_report, volume_state
from stage_b.units import unit_key

SEED = 20261002
N_VARIANTS, N_EUR, N_OTHER, SPACING = 240, 503, 20, 5_000
SENTINEL_INDEX, SECOND_INDEX, DISTINCT_INDEX, EXCLUDED_INDEX = 120, 80, 160, 5
CHROM, CENTER_GRCH38, GRCH37_SHIFT = "1", 50_000_000, -12_345
AR1, ALLELE_FREQUENCY = 0.9, 0.3
GENE = "ENSG00000000000"
SECOND_GENE = "ENSG00000000001"     # a second gene the deCODE assay serves: introns in the eQTL Catalogue, no eQTL, no VEP hit
PQTL_N, QTL_AN, N_CASE, N_CONTROL = 35_000.0, 1000, 20_000, 80_000
TISSUES = {"blood": "Whole_Blood", "liver": "Liver", "thyroid": "Thyroid"}
DECODE_KEY, OID, SOMAMER = "0_0_SYNTH_Synthetic.txt.gz", "OID00000", "SYNTH.1.2.3"
TAR_NAME = f"SYNTH_P00000_{OID}_v1_Synthetic.tar"
MAP_NAME = "olink_rsid_map_mac5_info03_b0_7_chr1_patched_v2.tsv.gz"
SHARED, GRCH37 = "FINNGEN_R12_SYN_SHARED", "ieu-b-9001"
DISTINCT, FEW, NO_DIR, NO_FILE = "GCST90000001", "GCST90000002", "GCST90000003", "GCST90000004"
OUTCOME_SOURCE = {SHARED: "finngen", GRCH37: "opengwas", DISTINCT: "gwas_catalog", FEW: "gwas_catalog",
                  NO_DIR: "gwas_catalog", NO_FILE: "gwas_catalog"}
UNIT_OUTCOMES = {"decode": (SHARED, GRCH37, DISTINCT, FEW, NO_DIR, NO_FILE), "ukbppp": (SHARED, DISTINCT),
                 "interval": (SHARED,)}
EXPECTED_STATE = {SHARED: "supportive", GRCH37: "supportive", DISTINCT: "inconclusive", FEW: "inconclusive",
                  NO_DIR: "inconclusive", NO_FILE: "inconclusive"}
EXPECTED_NOT_RUN = {SHARED: "", GRCH37: "", DISTINCT: "", FEW: "fewer_than_50_shared", NO_DIR: "outcome_file_unavailable",
                    NO_FILE: "outcome_file_unavailable"}
TOKENS = {"decode": "folder-token-one", "decode_smp": "smp-folder-token-one", "opengwas": "opengwas-token-one"}


def rsid(i: int) -> str:
    return f"rs{9_000_000 + i}"


def panel() -> tuple[pd.DataFrame, np.ndarray]:
    """(meta with rsid/ref/alt/pos on GRCh38, alleles: variants x haplotypes), fixed by SEED. The
    first 2 x N_EUR haplotypes are the EUR samples'; allele frequency 0.3 among them."""
    rng = np.random.default_rng(SEED)
    n_hap = 2 * (N_EUR + N_OTHER)
    latent = np.empty((N_VARIANTS, n_hap))
    latent[0] = rng.standard_normal(n_hap)
    for i in range(1, N_VARIANTS):
        latent[i] = AR1 * latent[i - 1] + math.sqrt(1 - AR1 ** 2) * rng.standard_normal(n_hap)
    cut = np.quantile(latent[:, :2 * N_EUR], 1 - ALLELE_FREQUENCY, axis=1, keepdims=True)
    meta = pd.DataFrame({"rsid": [rsid(i) for i in range(N_VARIANTS)], "ref": "A", "alt": "G",
                         "pos": [CENTER_GRCH38 + SPACING * (i - SENTINEL_INDEX) for i in range(N_VARIANTS)]})
    return meta, (latent > cut).astype(int)


def eur_dosage(alleles: np.ndarray) -> np.ndarray:
    return (alleles[:, 0:2 * N_EUR:2] + alleles[:, 1:2 * N_EUR:2]).astype(float)


def region_table(meta: pd.DataFrame, dosage: np.ndarray, z: np.ndarray, n: float, case_fraction: float | None,
                 shift: int = 0) -> pd.DataFrame:
    """One row per variant (effect on the panel's alt allele) carrying z-scores `z`. The standard
    error is the one a sample of `n` gives a standardized trait, or a case-control trait with
    `case_fraction` cases, at each variant's allele frequency."""
    eaf = dosage.mean(axis=1) / 2
    scale = 1.0 if case_fraction is None else case_fraction * (1 - case_fraction)
    se = 1 / np.sqrt(2 * n * scale * eaf * (1 - eaf))
    return pd.DataFrame({"rsid": meta["rsid"], "pos": (meta["pos"] + shift).astype(int), "ea": meta["alt"], "oa": meta["ref"],
                         "eaf": eaf, "beta": z * se, "se": se, "p": [math.erfc(abs(v) / math.sqrt(2)) for v in z], "n": n})


# ---- each source's file format ---------------------------------------------------------------------

def _lines(header: list[str], rows: list[list], sep: str = "\t") -> str:
    return "".join(sep.join(str(c) for c in r) + "\n" for r in [header, *rows])


def _gz(text: str) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as f:
        f.write(text.encode())
    return buf.getvalue()


def _name(pos: int) -> str:
    return f"chr{CHROM}:{pos}:A:G"


def decode_text(t: pd.DataFrame) -> str:
    return _lines(["Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele", "Beta", "Pval", "min_log10_pval", "SE", "N",
                   "ImpMAF"],
                  [[f"chr{CHROM}", r.pos, _name(r.pos), r.rsid, r.ea, r.oa, r.beta, r.p, -math.log10(r.p), r.se, r.n,
                    min(r.eaf, 1 - r.eaf)] for r in t.itertuples()])


def decode_annotation_text(t: pd.DataFrame, excluded: bool) -> str:
    keep = t.iloc[[EXCLUDED_INDEX]] if excluded else t.drop(index=EXCLUDED_INDEX)
    header = ["Chrom", "Pos", "Name", "rsids", "effectAllele", "otherAllele"] + ([] if excluded else ["effectAlleleFreq"])
    return _lines(header, [[f"chr{CHROM}", r.pos, _name(r.pos), r.rsid, r.ea, r.oa] + ([] if excluded else [r.eaf])
                           for r in keep.itertuples()])


def ukbppp_id(pos: int) -> str:
    return f"{CHROM}:{pos}:A:G:imp:v1"


def ukbppp_tar(t: pd.DataFrame) -> bytes:
    """A per-protein tar as on Synapse: one REGENIE file per chromosome, gzipped, in one directory."""
    header = ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO", "N", "TEST", "BETA", "SE", "CHISQ", "LOG10P",
              "EXTRA"]
    chr1 = _lines(header, [[CHROM, r.pos, ukbppp_id(r.pos), r.oa, r.ea, r.eaf, 1, int(r.n), "ADD", r.beta, r.se,
                            (r.beta / r.se) ** 2, -math.log10(r.p), "NA"] for r in t.itertuples()], sep=" ")
    chr2 = _lines(header, [["2", 1000, "2:1000:A:G:imp:v1", "A", "G", 0.3, 1, 30000, "ADD", 0.0, 0.01, 0.0, 0.0, "NA"]], sep=" ")
    stem = TAR_NAME.removesuffix(".tar")
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as tar:
        for chrom, text in (("1", chr1), ("2", chr2)):
            data = _gz(text)
            info = tarfile.TarInfo(f"{stem}/discovery_chr{chrom}_SYNTH:P00000:{OID}:v1:Synthetic.gz")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def rsid_map_text(meta: pd.DataFrame) -> str:
    return _lines(["ID", "REF", "ALT", "rsid", "POS19", "POS38"],
                  [[ukbppp_id(r.pos), r.ref, r.alt, r.rsid, r.pos + GRCH37_SHIFT, r.pos] for r in meta.itertuples()])


def gwas_catalog_text(t: pd.DataFrame) -> str:
    return _lines(["hm_variant_id", "hm_rsid", "hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta",
                   "hm_odds_ratio", "hm_effect_allele_frequency", "p_value", "standard_error"],
                  [[f"{CHROM}_{r.pos}_A_G", r.rsid, CHROM, r.pos, r.oa, r.ea, r.beta, "NA", r.eaf, r.p, r.se]
                   for r in t.itertuples()])


def finngen_text(t: pd.DataFrame) -> str:
    return _lines(["#chrom", "pos", "ref", "alt", "rsids", "nearest_genes", "pval", "mlogp", "beta", "sebeta", "af_alt",
                   "af_alt_cases", "af_alt_controls"],
                  [[CHROM, r.pos, r.oa, r.ea, r.rsid, "SYNTH", r.p, -math.log10(r.p), r.beta, r.se, r.eaf, r.eaf, r.eaf]
                   for r in t.itertuples()])


def eqtl_catalogue_text(traits: list[tuple[str, str, pd.DataFrame]]) -> str:
    """`traits`: (molecular_trait_id, gene_id, table); rows in position order, as tabix needs."""
    rows = [[trait, CHROM, r.pos, r.oa, r.ea, f"chr{CHROM}_{r.pos}_A_G", 100, min(r.eaf, 1 - r.eaf), r.p, r.beta, r.se, "SNP",
             round(r.eaf * QTL_AN), QTL_AN, 1, trait, gene, 5.0, r.rsid]
            for trait, gene, t in traits for r in t.itertuples()]
    return _lines(["molecular_trait_id", "chromosome", "position", "ref", "alt", "variant", "ma_samples", "maf", "pvalue", "beta",
                   "se", "type", "ac", "an", "r2", "molecular_trait_object_id", "gene_id", "median_tpm", "rsid"],
                  sorted(rows, key=lambda r: (r[2], r[0])))


def vcf_text(meta: pd.DataFrame, alleles: np.ndarray, samples: list[str]) -> str:
    head = ("##fileformat=VCFv4.2\n" f"##contig=<ID={CHROM}>\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    rows = [[CHROM, r.pos, r.rsid, r.ref, r.alt, ".", "PASS", ".", "GT",
             *(f"{a}|{b}" for a, b in zip(alleles[i, 0::2], alleles[i, 1::2]))] for i, r in enumerate(meta.itertuples())]
    return head + _lines(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *samples], rows)


def bgzip_tabix(text: str, work: Path, name: str, tabix_args: list[str]) -> dict[str, bytes]:
    """{name: bgzipped bytes, name.tbi: its tabix index}, made with htslib's own tools."""
    work.mkdir(parents=True, exist_ok=True)
    plain = work / name.removesuffix(".gz")
    plain.write_text(text)
    subprocess.run(["bgzip", "-f", str(plain)], check=True, capture_output=True)
    made = plain.with_name(plain.name + ".gz")
    target = work / name
    if made != target:
        made.replace(target)
    subprocess.run(["tabix", "-f", *tabix_args, str(target)], check=True, capture_output=True)
    return {name: target.read_bytes(), f"{name}.tbi": target.with_name(name + ".tbi").read_bytes()}


# ---- the synthetic sources -------------------------------------------------------------------------

@dataclass
class World:
    remote: FakeRemote
    units: list[InstrumentUnit]
    endpoints: Endpoints
    folders: dict[str, dict[str, str]]       # Synapse folder -> {file name: entity id}
    annotation: Path
    excluded: Path

    def sources(self, cache: Path, decode_token: str = TOKENS["decode"], synapse: FakeSynapse | None = None) -> RemoteSources:
        syn = synapse or FakeSynapse(self.remote, self.folders)
        return RemoteSources(cache, decode_token, lambda: syn, self.endpoints, decode_smp_token=TOKENS["decode_smp"])

    def fetcher(self, root: Path, opengwas_token: str = TOKENS["opengwas"]) -> VolumeFetcher:
        return VolumeFetcher(root, opengwas_token, self.annotation, self.excluded, self.endpoints)


def _outcomes(names: tuple[str, ...]) -> tuple[OutcomeSpec, ...]:
    return tuple(OutcomeSpec(accession=a, source=OUTCOME_SOURCE[a], n_case=N_CASE, n_control=N_CONTROL, risk_coded=True)
                 for a in sorted(names))


def build_world(remote: FakeRemote, work: Path, indexed: bool = True) -> World:
    """Fill `remote` with the synthetic sources and return the units that read them. With
    `indexed=False` the tabix-indexed files are left out (no bgzip or tabix needed): enough for the
    collect phase and for the whole-file reads of the analyze phase."""
    base = remote.base
    meta, alleles = panel()
    dosage = eur_dosage(alleles)
    r = np.corrcoef(dosage)

    def signal(*weighted: tuple[float, int]) -> np.ndarray:
        return sum(w * r[:, i] for w, i in weighted)

    pqtl_z = signal((12.0, SENTINEL_INDEX), (7.0, SECOND_INDEX))
    pqtl = region_table(meta, dosage, pqtl_z, PQTL_N, None)
    case = N_CASE / (N_CASE + N_CONTROL)
    shared = region_table(meta, dosage, signal((8.0, SENTINEL_INDEX)), float(N_CASE + N_CONTROL), case)
    distinct = region_table(meta, dosage, signal((8.0, DISTINCT_INDEX)), float(N_CASE + N_CONTROL), case)
    remote.tokens.update(TOKENS)
    remote.json_links.add("decode_smp")      # one folder answers with a redirect, the other with a JSON link

    # whole files: deCODE (and its SMP release), UKB-PPP tar and rsID map, GWAS Catalog without an index
    remote.files[f"/decode/{DECODE_KEY}"] = _gz(decode_text(pqtl))
    remote.files[f"/decode_smp/{DECODE_KEY}"] = _gz(decode_text(region_table(meta, dosage, 0.9 * pqtl_z, PQTL_N, None)))
    remote.files[f"/synapse/{TAR_NAME}"] = ukbppp_tar(pqtl)
    remote.files[f"/synapse/{MAP_NAME}"] = _gz(rsid_map_text(meta))
    folders = {SYNAPSE_UKBPPP_EUR: {TAR_NAME: "syn900001", "OTHER_P11111_OID11111_v1_Synthetic.tar": "syn900002"},
               SYNAPSE_RSID_MAPS: {MAP_NAME: "syn900003", MAP_NAME.replace("chr1_", "chr2_"): "syn900004"}}
    work.mkdir(parents=True, exist_ok=True)
    annotation, excluded = work / "assocvariants.annotated.txt.gz", work / "assocvariants.excluded.txt.gz"
    annotation.write_bytes(_gz(decode_annotation_text(pqtl, excluded=False)))
    excluded.write_bytes(_gz(decode_annotation_text(pqtl, excluded=True)))
    gwascat = urllib.parse.urlsplit(gwas_catalog_dir("/gwascat", DISTINCT)).path

    def listing(*names: str) -> bytes:
        return "".join(f'<a href="{n}">{n}</a>\n' for n in names).encode()

    whole = f"1-{DISTINCT}-EFO_1.h.tsv.gz"
    remote.files[gwascat] = listing(whole, "readme.txt")
    remote.files[gwascat + whole] = _gz(gwas_catalog_text(distinct))
    remote.files[gwas_catalog_dir("/gwascat", NO_FILE)] = listing("readme.txt")
    remote.files[gwas_catalog_study_dir("/gwascat", NO_FILE)] = listing("harmonised/", "readme.txt")

    # JSON endpoints: Ensembl, OpenGWAS, GTEx
    index = {rs: i for i, rs in enumerate(meta["rsid"])}
    for build, shift in (("GRCh38", 0), ("GRCh37", GRCH37_SHIFT)):
        def variation(rest: str, q, b, h, build=build, shift=shift) -> tuple[int, object]:
            if rest not in index:
                return 400, {"error": f"{rest} not found for human"}
            return 200, {"mappings": [{"seq_region_name": CHROM, "assembly_name": build,
                                       "start": int(meta["pos"][index[rest]]) + shift}]}
        remote.json_prefixes.append(("GET", f"/ensembl/{build}/variation/human/", variation))
        remote.json_routes[("POST", f"/ensembl/{build}/vep/human/id")] = lambda q, b, h: (200, [
            {"id": json.loads(b)["ids"][0], "transcript_consequences": [{"gene_id": GENE, "consequence_terms": ["missense_variant"]}]}])
    datasets = {GRCH37: region_table(meta, dosage, signal((8.0, SENTINEL_INDEX)), float(N_CASE + N_CONTROL), case, GRCH37_SHIFT),
                "prot-a-9001": region_table(meta, dosage, pqtl_z, PQTL_N, None, GRCH37_SHIFT)}

    def associations(q, b, h) -> tuple[int, object]:
        if h.get("Authorization") != f"Bearer {remote.tokens['opengwas']}":
            return 401, {"message": "token expired"}
        req = json.loads(b)
        chrom, _, span = req["variant"][0].partition(":")
        lo, hi = (int(x) for x in span.split("-"))
        t = datasets[req["id"][0]]
        return 200, [{"id": req["id"][0], "rsid": x.rsid, "chr": chrom, "position": x.pos, "ea": x.ea, "nea": x.oa, "eaf": x.eaf,
                      "beta": x.beta, "se": x.se, "p": x.p, "n": x.n} for x in t[(t["pos"] >= lo) & (t["pos"] <= hi)].itertuples()]
    remote.json_routes[("POST", "/opengwas/associations")] = associations
    remote.json_routes[("GET", "/gtex/reference/gene")] = lambda q, b, h: (200, {"data": [{"gencodeId": f"{GENE}.1"}]})
    remote.json_routes[("GET", "/gtex/expression/medianGeneExpression")] = lambda q, b, h: (200, {"data": [
        {"median": 9.0, "tissueSiteDetailId": "Thyroid"}, {"median": 2.0, "tissueSiteDetailId": "Liver"}]})

    if indexed:   # files queried by region: FinnGen, an indexed GWAS Catalog file, the eQTL Catalogue, 1000 Genomes
        tbx = work / "tabix"
        for name, data in bgzip_tabix(finngen_text(shared), tbx, "finngen_R12_SYN_SHARED.gz", ["-s", "1", "-b", "2", "-e", "2"]).items():
            remote.files[f"/finngen/{name}"] = data
        few_dir = urllib.parse.urlsplit(gwas_catalog_dir("/gwascat", FEW)).path
        few_name = f"2-{FEW}-EFO_2.h.tsv.gz"
        few = shared.iloc[SENTINEL_INDEX - 15:SENTINEL_INDEX + 15]
        files = bgzip_tabix(gwas_catalog_text(few), tbx, few_name, ["-s", "3", "-b", "4", "-e", "4", "-S", "1"])
        remote.files[few_dir] = listing(*files)
        remote.files.update({few_dir + n: d for n, d in files.items()})
        flat = 0.8 * np.cos(np.arange(N_VARIANTS))
        paths = []
        for tissue in TISSUES:
            lead = signal((8.0, SENTINEL_INDEX)) if tissue == "liver" else flat

            def qtl(z: np.ndarray) -> pd.DataFrame:
                return region_table(meta, dosage, z, QTL_AN / 2, None)

            traits = {"ge": [(GENE, GENE, qtl(flat))],
                      "leafcutter": [(f"{tissue}_intron_a", GENE, qtl(lead)), (f"{tissue}_intron_b", GENE, qtl(flat)),
                                     (f"{tissue}_other_gene_intron", SECOND_GENE, qtl(signal((11.0, DISTINCT_INDEX))))]}
            for method, rows in traits.items():
                name = f"GTEx_{tissue}_{method}.all.tsv.gz"
                for n, d in bgzip_tabix(eqtl_catalogue_text(rows), tbx, name, ["-s", "2", "-b", "3", "-e", "3", "-S", "1"]).items():
                    remote.files[f"/eqtlcat/{n}"] = d
                paths.append(["GTEx", tissue, method, f"{base}/eqtlcat/{name}"])
        remote.files["/eqtlcat/tabix_ftp_paths.tsv"] = _lines(["study_label", "sample_group", "quant_method", "ftp_path"], paths).encode()
        samples = [f"EUR{i:04d}" for i in range(N_EUR)] + [f"AFR{i:04d}" for i in range(N_OTHER)]
        for n, d in bgzip_tabix(vcf_text(meta, alleles, samples), tbx, f"chr{CHROM}.vcf.gz", ["-p", "vcf"]).items():
            remote.files[f"/kg/{n}"] = d
        remote.files["/kg/panel"] = _lines(["sample", "pop", "super_pop", "gender"],
                                           [[s, "SYN", s[:3], "female"] for s in samples]).encode()

    endpoints = Endpoints(opengwas_api=f"{base}/opengwas", ensembl_rest={b: f"{base}/ensembl/{b}" for b in ("GRCh38", "GRCh37")},
                          gwascat_ftp=f"{base}/gwascat", finngen=f"{base}/finngen/finngen_R12_{{endpoint}}.gz",
                          kg_vcf=f"{base}/kg/chr{{chrom}}.vcf.gz", kg_panel=f"{base}/kg/panel",
                          eqtlcat_paths=f"{base}/eqtlcat/tabix_ftp_paths.tsv", gtex_api=f"{base}/gtex",
                          decode_file=f"{base}/decode/s3/download?token={{token}}&file={{key}}",
                          decode_smp_file=f"{base}/decode_smp/s3/download?token={{token}}&file={{key}}")

    def listed(source: str) -> SourceFile:
        return SourceFile(name=DECODE_KEY, size=len(remote.files[f"/{source}/{DECODE_KEY}"]), etag=remote.etag(f"/{source}/{DECODE_KEY}"))

    def unit(source: str, assay: str, locator: str, build: str, **extra) -> InstrumentUnit:
        pos = CENTER_GRCH38 + (GRCH37_SHIFT if build == "GRCh37" else 0)
        sentinel = Sentinel(source=source, assay_id=assay, rsid=rsid(SENTINEL_INDEX), chrom=CHROM, pos=pos, build=build,
                            neg_log10_p=30.0, locator=assay)
        return InstrumentUnit(unit_key=unit_key(source, assay, GENE), source=source, assay_id=assay, gene_symbol="SYNTH",
                              gene_ensembl=GENE, platform="Olink" if source == "ukbppp" else "SomaScan", sentinel=sentinel,
                              pqtl_locator=locator, outcomes=_outcomes(UNIT_OUTCOMES[source]), **extra)

    units = [unit("decode", "0_0", DECODE_KEY, "GRCh38", pqtl_listing=listed("decode"), smp_listing=listed("decode_smp")),
             unit("interval", SOMAMER, "prot-a-9001", "GRCh37"), unit("ukbppp", OID, OID, "GRCh38")]
    return World(remote=remote, units=units, endpoints=endpoints, folders=folders, annotation=annotation, excluded=excluded)


def second_gene_unit(unit: InstrumentUnit) -> InstrumentUnit:
    """The unit of `unit`'s assay for SECOND_GENE, as units.build_units forms it for an assay that
    serves two genes: the same sentinel, locator and listings, its own gene and key, and SHARED as
    its one outcome."""
    return unit.model_copy(update={"unit_key": unit_key(unit.source, unit.assay_id, SECOND_GENE), "gene_symbol": "SYNTH2",
                                   "gene_ensembl": SECOND_GENE, "outcomes": _outcomes((SHARED,))})


# ---- known answers -----------------------------------------------------------------------------------

def hypotheses(unit: InstrumentUnit) -> list[HypothesisInput]:
    """One `decrease` hypothesis per outcome of the unit, the instrument being the unit's."""
    return [HypothesisInput(hypothesis_id=f"hyp_{unit.source}_{o.accession}", gene_symbol=unit.gene_symbol,
                            gene_ensembl=unit.gene_ensembl, direction="decrease", instrument_source=unit.source,
                            instrument_assay_id=unit.assay_id, platform=unit.platform, outcome_accession=o.accession,
                            outcome_source=o.source, outcome_n_case=o.n_case, outcome_n_control=o.n_control)
            for o in unit.outcomes]


def _ld_shape(unit_dir: Path) -> tuple[int, ...]:
    with np.load(unit_dir / "ld.npz", allow_pickle=False) as z:
        return tuple(z["dosage"].shape)


def check_unit(unit: InstrumentUnit, unit_dir: Path, result: dict, collect_sha256: str) -> dict:
    """The dry-run record of one finished unit directory: what colocalization returned for each
    outcome, the evidence rows assembly forms from it, and `checks`, each True when the known
    answer of the module docstring was obtained. `collect_sha256` is the collect digest recomputed
    from the unit's collect records, as assembly recomputes it."""
    collected, metas = collect_unit_dir(unit, unit_dir, result["fingerprint"], collect_sha256)
    hyps = hypotheses(unit)
    rows = build_evidence(hyps, {h.hypothesis_id: unit.unit_key for h in hyps},
                          {unit.unit_key: (unit.source, unit.assay_id, unit.gene_ensembl)}, {unit.unit_key: collected}, {})
    state = {h.outcome_accession: r.evidence_state for h, r in zip(hyps, rows)}
    out = result["outcomes"]
    names = [o.accession for o in unit.outcomes]
    shared = out[SHARED]
    n_pqtl = N_VARIANTS - (1 if unit.source == "decode" else 0)          # deCODE drops its excluded variant
    checks = {
        "result_on_disk_equals_returned": collected == result,
        "result_carries_the_collect_digest_of_its_records": result.get("collect_sha256") == collect_sha256,
        "pqtl_available": result["pqtl_available"] is True,
        "not_run_reasons": {a: out[a]["not_run_reason"] for a in names} == {a: EXPECTED_NOT_RUN[a] for a in names},
        "evidence_states": state == {a: EXPECTED_STATE[a] for a in names},
        "shared_pp_h4_above_0.8": bool(shared["coloc_run"]) and shared["pp"][4] > 0.8,
        "shared_lead_is_sentinel_or_proxy": shared.get("lead_variant") in result["sentinel_proxies"],
        "shared_direction_risk_increasing": shared.get("genetic_direction") == 1,
        "shared_primary_window_is_500kb": shared.get("nsnps") == 2 * WINDOW_PRIMARY // SPACING + 1,
        "shared_s15c_uses_the_1mb_window": shared.get("s15c_n_shared") == n_pqtl,
        "shared_s15a_s15b_s15c_returned": all(isinstance(shared.get(k), float) for k in ("s15a_pp_h4", "s15b_pp_h4", "s15c_pp_h4")),
        "susie_found_two_pqtl_credible_sets": (shared.get("susie") or {}).get("n_cs1", 0) >= 2 and shared.get("susie_note") == "",
        "s15g_from_coloc_susie": shared.get("s15g_method") == "coloc.susie",
        "ld_panel_is_the_503_eur_samples": _ld_shape(unit_dir) == (N_VARIANTS, N_EUR),
        "vep_hit_for_the_lead": result["vep"].get(shared.get("lead_variant", ""), {}).get("hit") is True,
        "splicing_candidate": result["splicing"].get("splicing_candidate") is True,
        "regional_manifest_rows_validate": len(regional_rows({unit.unit_key: metas})) == len(metas),
    }
    if GRCH37 in out:
        checks["grch37_outcome_pp_h4_above_0.8"] = bool(out[GRCH37]["coloc_run"]) and out[GRCH37]["pp"][4] > 0.8
    if DISTINCT in out:
        checks["distinct_pp_h3_above_0.8"] = bool(out[DISTINCT]["coloc_run"]) and out[DISTINCT]["pp"][3] > 0.8
    if FEW in out:
        checks["few_counts_30_shared"] = out[FEW].get("coverage", {}).get("n_shared") == 30
    if unit.source == "decode":
        s16 = result["s16"].get(SHARED, {})
        checks["s16_pp_h4_above_0.8"] = bool(s16.get("coloc_run")) and s16.get("pp_h4", 0.0) > 0.8
    return {"checks": checks, "passed": all(checks.values()),
            "outcomes": {a: {k: rec.get(k) for k in ("coloc_run", "not_run_reason", "detail", "pp", "nsnps", "lead_variant",
                                                    "genetic_direction", "s15a_pp_h4", "s15b_pp_h4", "s15c_pp_h4",
                                                    "s15g_pp_h4", "s15g_method", "susie", "susie_note", "coverage")}
                         for a, rec in out.items()},
            "splicing": result["splicing"], "s16": result["s16"], "coloc_session": result["coloc_session"],
            "evidence_rows": [r.model_dump() for r in rows], "regional_manifest": metas,
            "unit_files": sorted(p.name for p in unit_dir.iterdir())}


# ---- the dry run -------------------------------------------------------------------------------------

def _raised(call: Callable[[], object]) -> dict:
    """The stage B error a call raised, as the status report records it; {} when it returned."""
    try:
        call()
    except StageBError as err:
        return error_of(err)
    return {}


# ---- Ensembl's error replies ---------------------------------------------------------------------------

OTHER_RSID = "rs1"      # an rsID no synthetic request asks for
# name -> (HTTP status, body) for a request about the rsID `rs`; a `bytes` body is sent as it is.
ENSEMBL_REPLIES: dict[str, Callable[[str, str], tuple[int, object]]] = {
    "unknown_id": lambda rs, genuine: (400, {"error": genuine}),
    "backend_error": lambda rs, genuine: (400, {"error": "Unknown backend error while querying the variation database"}),
    "endpoint_not_found": lambda rs, genuine: (400, {"error": "endpoint not found"}),
    "page_not_found_404": lambda rs, genuine: (404, {"error": "page not found. Please check your uri and refer to our "
                                                              "documentation https://rest.ensembl.org/"}),
    "malformed_json": lambda rs, genuine: (400, json.dumps({"error": genuine}).encode()[:-1]),
    "another_rsid": lambda rs, genuine: (400, {"error": genuine.replace(rs, OTHER_RSID)}),
    "unknown_id_and_more": lambda rs, genuine: (400, {"error": f"{genuine}; Unknown backend error"}),
    "unknown_id_beside_another_key": lambda rs, genuine: (400, {"error": genuine, "status": "backend unavailable"}),
    "post_too_large": lambda rs, genuine: (400, {"error": "POST message too large. You have submitted 201 elements but a "
                                                          "limit of 200 is in place. Request smaller regions or lists of IDs"}),
    "words_without_the_message": lambda rs, genuine: (400, {"error": f"ID '{rs}' not found"}),
    "not_an_object": lambda rs, genuine: (400, [genuine]),
    "unknown_id_with_503": lambda rs, genuine: (503, {"error": genuine}),
}
ENSEMBL_EXPECTED = {name: "absent" if name == "unknown_id" else "raised:server" if name.endswith("_503") else "raised:protocol"
                    for name in ENSEMBL_REPLIES}


def ensembl_outcomes(world: World, fetcher: VolumeFetcher, sentinel: Sentinel) -> dict[str, dict[str, str]]:
    """What the position lookup and the VEP lookup of `fetcher` make of each reply of
    ENSEMBL_REPLIES about the sentinel's rsID: "absent" (the lookup takes it as Ensembl not knowing
    the rsID: no position on the other build, SourceAbsent from VEP), "raised:<kind>" (a
    RetryableSourceError) or "answered". The routes of `world.remote` are put back afterwards."""
    remote, rs = world.remote, sentinel.rsid
    other = "GRCh37" if sentinel.build == "GRCh38" else "GRCh38"
    vep_route = ("POST", f"/ensembl/{sentinel.build}/vep/human/id")
    saved = remote.json_routes[vep_route]

    def outcome(call: Callable[[], bool]) -> str:
        try:
            return "absent" if call() else "answered"
        except SourceAbsent:
            return "absent"
        except RetryableSourceError as err:
            return f"raised:{err.kind}"

    out: dict[str, dict[str, str]] = {"position_lookup": {}, "vep": {}}
    try:
        for name, reply in ENSEMBL_REPLIES.items():
            position = reply(rs, f"{rs} not found for human")
            remote.json_prefixes.insert(0, ("GET", f"/ensembl/{other}/variation/human/", lambda rest, q, b, h, r=position: r))
            try:
                out["position_lookup"][name] = outcome(lambda: fetcher.positions(sentinel)[other] is None)
            finally:
                remote.json_prefixes.pop(0)
            vep = reply(rs, f"No variant found with ID '{rs}'")
            remote.json_routes[vep_route] = lambda q, b, h, r=vep: r
            out["vep"][name] = outcome(lambda: not fetcher.vep([rs], sentinel.build))
    finally:
        remote.json_routes[vep_route] = saved
    return out


def same_size_other_bytes(gz: bytes) -> bytes:
    """A gzip file of the same length that differs in one byte and still decompresses to its end:
    the OS field of the gzip header (byte 9), which no reader checks."""
    return gz[:9] + bytes([gz[9] ^ 0x01]) + gz[10:]


def scenario(world: World, root: Path, analyze: Callable[[InstrumentUnit, VolumeFetcher, Path], dict],
             commit: Callable[[], None] = lambda: None) -> dict:
    """Collect and analyze the synthetic units under `root` with faults injected, and report.
    `analyze(unit, fetcher, units_root)` runs one unit under its checkpoint directory and returns
    its result (modal_stage_b.checkpointed_unit in the dry run). For the duration, the wait between
    retries and between eQTL Catalogue queries is switched off and downloads are read in 512-byte
    pieces with a commit every 4 KiB, so the kilobyte-sized synthetic files go through the resume
    and checkpoint code that gigabyte files do; nothing else of the real code path is changed."""
    remote, units = world.remote, world.units
    by_source = {u.source: u for u in units}
    tasks = {(t.source, t.key): t for t in collect_tasks(units)}
    decode_path, tar_path = f"/decode/{DECODE_KEY}", f"/synapse/{TAR_NAME}"
    whole_dir = urllib.parse.urlsplit(gwas_catalog_dir("/gwascat", DISTINCT)).path
    checks: dict[str, bool] = {}
    events: dict[str, object] = {}
    saved = (remote_errors.WAIT_S, fetch_module.EQTLCAT_MIN_INTERVAL_S, collect_module.CHUNK_BYTES)
    remote_errors.WAIT_S, fetch_module.EQTLCAT_MIN_INTERVAL_S, collect_module.CHUNK_BYTES = 0.0, 0.0, 512
    try:
        cache = root / "cache"

        def collect(task: CollectTask, **kw):
            return collect_one(task, world.sources(cache, **kw), root, commit, checkpoint_bytes=4096)

        # 1. an expired folder token: the file stays pending, nothing is recorded
        decode_task = tasks[("decode", "0_0")]
        events["expired_folder_token"] = _raised(lambda: collect(decode_task, decode_token="folder-token-expired"))
        checks["expired_folder_token_raises_auth_and_records_nothing"] = (
            events["expired_folder_token"].get("kind") == "auth" and not (root / "collect" / "decode").exists())
        # 2. the connection breaks in every attempt of one call; the next call resumes with a Range request
        size = len(remote.files[decode_path])
        remote.faults[decode_path] = [("break", size // 8)] * remote_errors.ATTEMPTS
        events["broken_stream"] = _raised(lambda: collect(decode_task))
        part = root / "raw" / "decode" / f"{DECODE_KEY}.part"
        held = part.stat().st_size if part.is_file() else 0
        checks["broken_stream_raises_and_keeps_the_partial_file"] = (
            events["broken_stream"].get("kind") == "connection" and 0 < held < size and not (root / "collect" / "decode").exists())
        mark = len(remote.log)
        record = collect(decode_task)
        ranges = [rng for _m, path, rng in remote.log[mark:] if path == decode_path]
        checks["next_call_resumes_with_a_range_request"] = ranges == [f"bytes={held}-"]
        checks["resumed_file_is_the_source_file"] = (record.status, record.bytes, record.etag) == (
            "collected", size, remote.etag(decode_path))
        checks["record_carries_no_token_and_no_signed_address"] = (
            TOKENS["decode"] not in json.dumps(record.model_dump()) and "signed" not in json.dumps(record.model_dump())
            and record.source_url == f"{remote.base}/decode/s3/download?token=<token>&file={DECODE_KEY}")
        # 3. a rate limit on a listing, then every other file; a pre-signed link withdrawn mid-download is renewed
        remote.faults[whole_dir] = [("status", 429)] * remote_errors.ATTEMPTS
        events["rate_limited_listing"] = _raised(lambda: collect(tasks[("gwas_catalog", DISTINCT)]))
        checks["rate_limit_raises_and_records_nothing"] = (
            events["rate_limited_listing"].get("kind") == "rate_limit" and not (root / "collect" / "gwas_catalog").exists())
        remote.faults[tar_path] = [("break", len(remote.files[tar_path]) // 2)]
        records = {key: collect(task) for key, task in tasks.items()}
        events["collect_records"] = {f"{s}/{k}": r.status for (s, k), r in records.items()}
        checks["collect_statuses"] = events["collect_records"] == {
            "decode/0_0": "collected", "decode_smp/0_0": "collected", f"ukbppp/{OID}": "collected",
            f"ukbppp_rsid_map/{CHROM}": "collected", f"gwas_catalog/{DISTINCT}": "collected",
            f"gwas_catalog/{FEW}": "remote_indexed", f"gwas_catalog/{NO_DIR}": "absent", f"gwas_catalog/{NO_FILE}": "absent"}
        whole = [m for m, path, _ in remote.log if path.startswith(whole_dir) and path.endswith(".h.tsv.gz") and m == "GET"]
        checks["file_shared_by_two_units_is_downloaded_once"] = len(whole) == 1
        # 4. the folder link is re-issued: collecting again touches nothing and no record changes
        remote.tokens["decode"] = "folder-token-two"
        mark = len(remote.log)
        again = {key: collect(task, decode_token="folder-token-two") for key, task in tasks.items()}
        checks["reissued_link_leaves_every_record_as_it_was"] = again == records and remote.log[mark:] == []
        # 5. analyze: an expired OpenGWAS token stops the unit; nothing is recorded as unavailable
        decode_unit = by_source["decode"]
        mark = len(remote.log)
        events["expired_opengwas_token"] = _raised(
            lambda: analyze(decode_unit, world.fetcher(root, opengwas_token="opengwas-token-expired"), root / "units"))
        unit_dir = root / "units" / decode_unit.unit_key
        metas = {p.name: json.loads(p.read_text())["status"] for p in unit_dir.glob("*.meta.json")}
        checks["expired_opengwas_token_raises_auth"] = events["expired_opengwas_token"].get("kind") == "auth"
        checks["expired_token_leaves_no_result_and_no_unavailable_step"] = (
            not (unit_dir / "result.json").exists() and not (unit_dir / f"outcome__{GRCH37}.meta.json").exists()
            and "unavailable" not in {s for n, s in metas.items() if GRCH37 in n})
        results = {u.source: analyze(u, world.fetcher(root), root / "units") for u in units}
        reference = analyze(decode_unit, world.fetcher(root), root / "reference_units")
        checks["unit_resumed_after_the_token_fault_equals_a_unit_run_without_it"] = results["decode"] == reference
        read = set(remote.paths_read(mark))
        checks["analyze_reads_no_whole_file_from_the_network"] = not any(
            p.startswith(("/decode", "/synapse")) or (p.startswith(whole_dir) and p.endswith(".gz")) for p in read)
        # 6. a link re-issued and tokens renewed after the units finished: every checkpoint is still valid
        remote.tokens.update({"decode": "folder-token-three", "opengwas": "opengwas-token-two"})
        mark = len(remote.log)
        later = {u.source: analyze(u, world.fetcher(root, opengwas_token="opengwas-token-two"), root / "units") for u in units}
        checks["finished_checkpoints_stay_valid_after_reissue_and_renewal"] = later == results and remote.log[mark:] == []
        # 7. a unit whose file has no collect record stops; it is not recorded as unavailable
        events["missing_collect_record"] = _raised(lambda: analyze(decode_unit, world.fetcher(root / "empty"), root / "empty" / "units"))
        empty_dir = root / "empty" / "units" / decode_unit.unit_key
        checks["missing_collect_record_raises_and_records_nothing"] = (
            events["missing_collect_record"].get("error_class") == "CollectError" and not (empty_dir / "pqtl.meta.json").exists())
        unit_reports = {u.unit_key: check_unit(u, root / "units" / u.unit_key, results[u.source], volume_collect_digest(root, u))
                        for u in units}
        checks.update({f"{k}: {name}": ok for k, rep in unit_reports.items() for name, ok in rep["checks"].items()})
        # 7a. the deCODE assay serves a second gene: one more unit, no more files; the same extract and
        #     colocalization under another checkpoint, and the gene-specific steps run for its own gene
        second = second_gene_unit(decode_unit)
        second_dir = root / "units" / second.unit_key
        checks["second_gene_unit_adds_no_file_to_collect"] = collect_tasks([*units, second]) == list(tasks.values())
        mark = len(remote.log)
        second_result = analyze(second, world.fetcher(root, opengwas_token="opengwas-token-two"), root / "units")
        checks["second_gene_unit_reads_the_collected_file_and_none_from_the_network"] = not any(
            p.startswith(("/decode", "/synapse")) for p in remote.paths_read(mark))
        first_shared, second_shared = results["decode"]["outcomes"][SHARED], second_result["outcomes"][SHARED]
        pqtl_sha = [json.loads((d / "pqtl.meta.json").read_text())["sha256"] for d in (unit_dir, second_dir)]
        checks["second_gene_unit_has_its_own_directory_and_fingerprint"] = (
            second_dir != unit_dir and (second_dir / "result.json").is_file()
            and second_result["fingerprint"] != results["decode"]["fingerprint"])
        checks["second_gene_unit_has_the_same_pqtl_extract_and_colocalization"] = (
            pqtl_sha[0] == pqtl_sha[1] and second_shared["pp"] == first_shared["pp"]
            and second_shared["lead_variant"] == first_shared["lead_variant"])
        lead = second_shared["lead_variant"]
        checks["vep_hit_is_for_the_units_own_gene"] = (
            results["decode"]["vep"][lead]["hit"] is True and second_result["vep"][lead]["hit"] is False)
        events["second_gene_splicing"] = second_result["splicing"]
        checks["splicing_flag_is_for_the_units_own_gene"] = (
            results["decode"]["splicing"]["splicing_candidate"] is True and second_result["splicing"]["splicing_candidate"] == ""
            and {t: v["lead_event"] for t, v in second_result["splicing"]["tissues"].items()}
            == {t: f"{t}_other_gene_intron" for t in TISSUES}
            and {v["eqtl_n_shared"] for v in second_result["splicing"]["tissues"].values()} == {0})
        ids = {u.unit_key: (u.source, u.assay_id, u.gene_ensembl) for u in (decode_unit, second)}
        both = {decode_unit.unit_key: results["decode"], second.unit_key: second_result}
        own = build_evidence(hypotheses(second), {h.hypothesis_id: second.unit_key for h in hypotheses(second)}, ids, both, {})
        events["other_genes_unit"] = _raised(lambda: build_evidence(
            hypotheses(second), {h.hypothesis_id: decode_unit.unit_key for h in hypotheses(second)}, ids, both, {}))
        checks["a_hypothesis_takes_its_flags_from_its_own_genes_unit_and_is_refused_on_the_other"] = (
            [(r.evidence_state, r.protein_altering, r.splicing_candidate) for r in own] == [("supportive", False, "")]
            and events["other_genes_unit"].get("error_class") == "InputContractError")
        # 7b. Ensembl: only its unknown-identifier message naming the requested rsID is an absence
        events["ensembl"] = ensembl_outcomes(world, world.fetcher(root, opengwas_token="opengwas-token-two"), decode_unit.sentinel)
        checks["ensembl_400_is_an_absence_only_for_the_unknown_id_message_naming_the_requested_rsid"] = (
            events["ensembl"] == {"position_lookup": ENSEMBL_EXPECTED, "vep": ENSEMBL_EXPECTED})
        # 7c. the source serves other bytes under the same name, size and ETag; the file is collected
        #     again from scratch: its record differs only in sha256 and the finished checkpoint is refused
        changed, unit_name = root / "changed_bytes", decode_unit.unit_key
        original = remote.files[decode_path]
        remote.etags[decode_path], remote.files[decode_path] = remote.etag(decode_path), same_size_other_bytes(original)
        try:
            recollected = {key: collect_one(task, world.sources(changed / "cache", decode_token="folder-token-three"), changed,
                                            commit, checkpoint_bytes=4096)
                           for key, task in tasks.items() if task in collect_tasks([decode_unit])}
        finally:
            remote.files[decode_path] = original
            del remote.etags[decode_path]
        shutil.copytree(root / "units" / unit_name, changed / "units" / unit_name)
        held = {p.name: p.read_bytes() for p in sorted((changed / "units" / unit_name).iterdir())}
        mark = len(remote.log)
        events["changed_bytes"] = _raised(lambda: analyze(decode_unit, world.fetcher(changed), changed / "units"))
        old, new = records[("decode", "0_0")], recollected[("decode", "0_0")]
        differing = {k for k in old.model_dump() if getattr(old, k) != getattr(new, k)} - {"utc"}
        checks["changed_bytes_keep_name_size_and_etag_and_differ_in_sha256_only"] = (
            differing == {"sha256"} and (new.name, new.bytes, new.etag) == (old.name, old.bytes, old.etag))
        digest_before, digest_after = volume_collect_digest(root, decode_unit), volume_collect_digest(changed, decode_unit)
        checks["changed_bytes_change_the_collect_digest"] = (
            digest_before == results["decode"]["collect_sha256"] and digest_after != digest_before)
        checks["checkpoint_of_the_old_bytes_is_refused_and_left_as_it_was"] = (
            events["changed_bytes"].get("error_class") == "StaleCheckpointError" and remote.log[mark:] == []
            and {p.name: p.read_bytes() for p in sorted((changed / "units" / unit_name).iterdir())} == held)
        # 8. every unit is finished: the whole files are deleted, the records and the extracts stay
        events["purge_before_units_finish"] = _raised(lambda: purge_raw(root / "empty", units, commit))
        extracts = {p: p.read_bytes() for p in sorted((root / "units").rglob("*.tsv.gz"))}
        purged = purge_raw(root, units, commit)
        left = sorted(p.name for p in (root / "raw").rglob("*") if p.is_file())
        checks["purge_refuses_while_a_unit_is_unfinished"] = events["purge_before_units_finish"].get("error_class") == "CollectError"
        checks["purge_deletes_the_whole_files_and_keeps_records_and_extracts"] = (
            left == [] and len(purged["deleted"]) == 5 and {k: read_record(root, *k) for k in tasks} == records
            and len(extracts) >= 9 and all(p.read_bytes() == data for p, data in extracts.items()))
        status = status_report(list(tasks.values()), units, volume_state(root), {}, {"deployed": False})
        checks["status_counts"] = (
            status["files_by_source"] == {"decode": {"collected": 1}, "decode_smp": {"collected": 1},
                                          "gwas_catalog": {"absent": 2, "collected": 1, "remote_indexed": 1},
                                          "ukbppp": {"collected": 1}, "ukbppp_rsid_map": {"collected": 1}}
            and status["units_by_source"] == {s: {"done": 1} for s in ("decode", "interval", "ukbppp")})
    finally:
        remote_errors.WAIT_S, fetch_module.EQTLCAT_MIN_INTERVAL_S, collect_module.CHUNK_BYTES = saved
    absent = read_record(root, "gwas_catalog", NO_DIR)
    return {"checks": checks, "passed": all(checks.values()), "events": events,
            "collect_records": {f"{s}/{k}": r.model_dump() for (s, k), r in records.items()},
            "absent_detail": absent.detail, "status": {k: status[k] for k in ("files_by_source", "units_by_source")},
            "requests": len(remote.log), "units": unit_reports}
