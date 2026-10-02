"""V8 stage A entry point: enumeration, eligibility, mechanism labels, outcome-GWAS selection
and instrument selection (PREREG.md; output contract stages/INTERFACES.md).

Writes to the output directory: run_info.json (stage, run token, repository commit, code digest;
written first, so the directory says which run it belongs to), hypotheses.csv, funnel.csv,
outcome_gwas_selection.csv, mechanism_crosstab.csv, outcome_trait_coding.tsv (read by stage B), then
INPUTS.tsv (every input file, repo- or volume-relative, with its sha256 and source, and a
`stage_code` row carrying the code digest and the commit) and MANIFEST.tsv (one row per output file
and INPUTS.tsv), in the shared format of stages/run_guard/v8_manifest.py.

Before sealing, the run compares the outcome-blind counts PREREG §Sample size registers (S1, H1 by
class, the neurologic-or-psychiatric split, S4, H1 on blood-secreted targets) with the same counts
taken from hypotheses.csv, writes both to registered_count_check.json, and raises
RegisteredCountMismatch without writing INPUTS.tsv or MANIFEST.tsv when any differs
(stage_a/registered.py). Stage A reads no outcome, so a mismatch is investigated and the stage is
rerun under a new run token.

The output directory (or volume path) must be empty or absent, or hold the run_info.json of this
run token and no MANIFEST.tsv: a directory written by another run, by none, or already completed
is refused (`claim_output_dir`). The run records one repository commit (`v8_run_guard.run_commit`):
the value baked into the Modal image, or HEAD of a clean working tree.

Nothing runs on real study data until the OSF registration is approved and the run is logged
below the line of PREREG.md: `v8_run_guard.require_run` (stages/run_guard) refuses unless
PREREG.md passes `prereg check`, its log holds the entry `RUN_START stage=A token=<token>` for the
`--run-token` given, and OSF registration 9tzfk is public and approved. The run also refuses
while any cis-pQTL publication date is unset (stage_a.flags.PUBLICATION_DATES).

Inputs under `inputs/` and `feasibility/v2_all_indications/inputs/` are resolved from
`--inputs-root`, a directory laid out like experiments/08_mechanism_interaction; real runs use
the Modal volume `pqtl-v8-inputs` through modal_stage_a.py, which mounts it read-only. The other
inputs (the deCODE supplement, the Olink map, the EpiGraphDB list, the pilot tables) are
repository files, resolved from `--exp-dir` and `--repo-root`.

Usage:
    cd experiments/08_mechanism_interaction/stages/A
    uv run --with modal==1.4.3 modal run modal_stage_a.py --run-token <token>     # real run, on Modal
    uv run python run_stage_a.py --paths-json my_paths.json --out-dir output --run-token <token>
"""
import argparse
import json
from collections.abc import Callable
from datetime import date, datetime, timezone
from pathlib import Path

from pydantic import BaseModel, ConfigDict
from v8_manifest import (MANIFEST_NAME, RUN_INFO_NAME, InputRecord, code_record, code_sha256, guard_files,
                         portable_path, relative_files, sha256_file, write_manifest)
from v8_run_guard import require_run, run_commit

from stage_a.build import StageAInputs, StageAOutputs, run_stage_a
from stage_a.flags import PUBLICATION_DATES, karim_launched_keys, pilot_keys, pilot_label_map, \
    require_publication_dates
from stage_a.instruments import InstrumentLists, decode_lists, decode_strict_diagnostic, interval_lists, \
    ukbppp_lists
from stage_a.loaders import load_clinical_indication, load_decode, load_disease, load_epigraphdb_genes, \
    load_finngen_study, load_gwas_catalog, load_hgnc, load_hpa_blood, load_karim_launched, load_mechanisms, \
    load_molecules, load_olink_map, load_opengwas, load_pilot_tables, load_sun2018_st4, load_ukbppp_cis, \
    load_uniprot_secreted
from stage_a.registered import REGISTERED_COUNTS, RegisteredCountMismatch, count_check
from stage_a.outcome_gwas import combine_candidates, finngen_candidates, gwas_catalog_candidates, \
    label_index, opengwas_candidates
from stage_a.schemas import CROSSTAB_COLUMNS, FUNNEL_COLUMNS, HYPOTHESIS_COLUMNS, OUTCOME_SELECTION_COLUMNS, \
    TRAIT_CODING_COLUMNS

HERE = Path(__file__).resolve().parent
EXP_DIR = HERE.parents[1]
# `.parent.parent`, not `.parents[1]`: in the Modal image this file is /root/stagea/run_stage_a.py, so
# EXP_DIR is `/`, where `.parents[1]` raises on import. The image passes every root on the command line.
REPO = EXP_DIR.parent.parent
PREREG = EXP_DIR / "PREREG.md"
V2_INPUTS_REL = Path("feasibility") / "v2_all_indications" / "inputs"
GWASCAT_RELEASE = "r2026-09-13"


class InputSpec(BaseModel):
    model_config = ConfigDict(frozen=True)

    path: Path
    source: str
    release: str
    url_or_accession: str


class StageAPaths(BaseModel):
    """Every input file. Defaults are the repository locations; `--paths-json` overrides any key."""

    model_config = ConfigDict(frozen=True)

    ot_clinical_indication: InputSpec
    ot_drug_mechanism_of_action: InputSpec
    ot_drug_molecule: InputSpec
    ot_disease: InputSpec
    ot_study: InputSpec
    hgnc: InputSpec
    hpa: InputSpec
    uniprot_kw0964: InputSpec
    uniprot_sl0243: InputSpec
    gwas_catalog_studies: InputSpec
    gwas_catalog_ancestries: InputSpec
    gwas_catalog_unpublished_studies: InputSpec
    gwas_catalog_unpublished_ancestries: InputSpec
    opengwas_gwasinfo: InputSpec
    ukbppp_cis_list: InputSpec
    olink_map: InputSpec
    decode_supplement: InputSpec
    interval_supplement: InputSpec
    epigraphdb_catalog: InputSpec
    frozen_candidates_v34: InputSpec
    classification_v5: InputSpec
    classification_v5_1: InputSpec
    karim2026_workbook: InputSpec


def default_paths(inputs_root: Path = EXP_DIR, exp_dir: Path = EXP_DIR, repo: Path = REPO) -> dict[str, dict]:
    """`inputs_root` holds the files moved to the `pqtl-v8-inputs` volume (inputs/ and
    feasibility/v2_all_indications/inputs/); `exp_dir` and `repo` the repository files."""
    v2_inputs = inputs_root / V2_INPUTS_REL
    ot = v2_inputs / "ot_26.09"
    ot_url = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output"
    gc = v2_inputs / "gwas_catalog"
    gc_url = "https://www.ebi.ac.uk/gwas/api/search/downloads"

    def gcat(kind: str) -> dict:
        return {"path": gc / f"gwas-catalog-v1.0.3.1-{kind}-{GWASCAT_RELEASE}.tsv", "source": "GWAS Catalog",
                "release": GWASCAT_RELEASE, "url_or_accession": f"{gc_url}/{kind.replace('-', '_')}/v1.0.3.1"}

    return {
        "ot_clinical_indication": {"path": ot / "clinical_indication.parquet", "source": "Open Targets Platform",
                                   "release": "26.09", "url_or_accession": f"{ot_url}/clinical_indication/"},
        "ot_drug_mechanism_of_action": {"path": ot / "drug_mechanism_of_action.parquet", "source": "Open Targets Platform",
                                        "release": "26.09", "url_or_accession": f"{ot_url}/drug_mechanism_of_action/"},
        "ot_drug_molecule": {"path": ot / "drug_molecule.parquet", "source": "Open Targets Platform",
                             "release": "26.09", "url_or_accession": f"{ot_url}/drug_molecule/"},
        "ot_disease": {"path": ot / "disease.parquet", "source": "Open Targets Platform", "release": "26.09",
                       "url_or_accession": f"{ot_url}/disease/"},
        "ot_study": {"path": ot / "study.parquet", "source": "Open Targets Platform (FinnGen R12 study index)",
                     "release": "26.09", "url_or_accession": f"{ot_url}/study/"},
        "hgnc": {"path": v2_inputs / "hgnc" / "hgnc_complete_set.txt", "source": "HGNC", "release": "complete set",
                 "url_or_accession": "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"},
        "hpa": {"path": v2_inputs / "hpa" / "proteinatlas.tsv.zip", "source": "Human Protein Atlas", "release": "25.1",
                "url_or_accession": "https://www.proteinatlas.org/download/proteinatlas.tsv.zip"},
        "uniprot_kw0964": {"path": v2_inputs / "uniprot" / "secreted_KW-0964.tsv", "source": "UniProt",
                           "release": "2026_03", "url_or_accession": "https://rest.uniprot.org/uniprotkb/search keyword:KW-0964"},
        "uniprot_sl0243": {"path": v2_inputs / "uniprot" / "secreted_SL-0243.tsv", "source": "UniProt",
                           "release": "2026_03", "url_or_accession": "https://rest.uniprot.org/uniprotkb/search cc_scl_term:SL-0243"},
        "gwas_catalog_studies": gcat("studies"),
        "gwas_catalog_ancestries": gcat("ancestries"),
        "gwas_catalog_unpublished_studies": gcat("unpublished-studies"),
        "gwas_catalog_unpublished_ancestries": gcat("unpublished-ancestries"),
        "opengwas_gwasinfo": {"path": v2_inputs / "opengwas" / "gwasinfo_all.json", "source": "IEU OpenGWAS",
                              "release": "gwasinfo listing", "url_or_accession": "https://api.opengwas.io/api/gwasinfo"},
        "ukbppp_cis_list": {"path": inputs_root / "inputs" / "ukbppp" / "ukbppp_cis_pqtl_genes.json", "source": "UKB-PPP",
                            "release": "Sun et al. 2023", "url_or_accession": "syn51365303"},
        "olink_map": {"path": repo / "planning" / "perplexity_v8_design" / "F_data_readmes" / "olink_protein_map_3k_v1.tsv",
                      "source": "UKB-PPP metadata", "release": "olink_protein_map_3k_v1", "url_or_accession": "syn51396728"},
        "decode_supplement": {"path": exp_dir / "feasibility" / "inputs" / "ferkingstad2021_MOESM4_ESM.xlsx",
                              "source": "deCODE (Ferkingstad et al. 2021) Supplementary Tables ST01, ST02",
                              "release": "Nat Genet 2021", "url_or_accession": "doi:10.1038/s41588-021-00978-w"},
        "interval_supplement": {"path": v2_inputs / "interval" / "sun2018_MOESM4_ESM.xlsx",
                                "source": "INTERVAL (Sun et al. 2018) Supplementary Table ST4",
                                "release": "Nature 2018",
                                "url_or_accession": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/41586_2018_175_MOESM4_ESM.xlsx"},
        "epigraphdb_catalog": {"path": repo / "zenodo_export" / "data" / "epigraphdb" / "v34_mr_catalog.csv",
                               "source": "EpiGraphDB V3.4 gene list (gene_symbol column only; S20)", "release": "V3.4",
                               "url_or_accession": "zenodo_export/data/epigraphdb/v34_mr_catalog.csv"},
        "frozen_candidates_v34": {"path": repo / "data" / "frozen_candidates_v34.csv", "source": "V3.4 pilot indications",
                                  "release": "V3.4", "url_or_accession": "data/frozen_candidates_v34.csv"},
        "classification_v5": {"path": repo / "results" / "v5" / "classification_v5.csv", "source": "V5 pilot keys (gene, disease only)",
                              "release": "V5", "url_or_accession": "results/v5/classification_v5.csv"},
        "classification_v5_1": {"path": repo / "results" / "v5_1" / "classification_v5_1.csv",
                                "source": "V5.1 pilot keys (gene, disease only)", "release": "V5.1",
                                "url_or_accession": "results/v5_1/classification_v5_1.csv"},
        "karim2026_workbook": {"path": inputs_root / "inputs" / "karim2026" / "v2_supp_tables.xlsx",
                               "source": "Karim et al. 2026 Supplementary Table 17 (launched pQTL-supported pairs)",
                               "release": "v2 supplementary tables", "url_or_accession": "doi:10.64898/2026.02.23.26346731"},
    }


def resolve_paths(paths_json: Path | None, inputs_root: Path, exp_dir: Path, repo: Path) -> StageAPaths:
    spec = default_paths(inputs_root, exp_dir, repo)
    if paths_json is not None:
        for key, override in json.loads(paths_json.read_text()).items():
            if key not in spec:
                raise KeyError(f"{paths_json}: unknown input key {key!r}")
            spec[key] = {**spec[key], **override} if isinstance(override, dict) else {**spec[key], "path": override}
    return StageAPaths.model_validate(spec)


def load_inputs(p: StageAPaths, publication_dates: dict[str, date]) -> StageAInputs:
    ens2sym, uni2sym = load_hgnc(p.hgnc.path)
    disease = load_disease(p.ot_disease.path)

    ukb_strict, ukb_incl = ukbppp_lists(load_ukbppp_cis(p.ukbppp_cis_list.path), load_olink_map(p.olink_map.path))
    st01, st02 = load_decode(p.decode_supplement.path)
    dec_strict, dec_incl = decode_lists(st01, st02)
    itv_strict, itv_incl = interval_lists(load_sun2018_st4(p.interval_supplement.path), uni2sym,
                                          load_epigraphdb_genes(p.epigraphdb_catalog.path))

    gcat = gwas_catalog_candidates(load_gwas_catalog({
        "studies": p.gwas_catalog_studies.path, "ancestries": p.gwas_catalog_ancestries.path,
        "unpublished_studies": p.gwas_catalog_unpublished_studies.path,
        "unpublished_ancestries": p.gwas_catalog_unpublished_ancestries.path}))
    fgn = finngen_candidates(load_finngen_study(p.ot_study.path))
    ogw = opengwas_candidates(load_opengwas(p.opengwas_gwasinfo.path), gcat, fgn, label_index(disease))

    frozen, classifications = load_pilot_tables(p.frozen_candidates_v34.path,
                                                [p.classification_v5.path, p.classification_v5_1.path])
    label_map = pilot_label_map(frozen, disease)

    return StageAInputs(
        clinical_indication=load_clinical_indication(p.ot_clinical_indication.path),
        mechanisms=load_mechanisms(p.ot_drug_mechanism_of_action.path),
        molecules=load_molecules(p.ot_drug_molecule.path),
        disease=disease,
        ens2sym=ens2sym,
        hpa_blood_ensembl=load_hpa_blood(p.hpa.path),
        uniprot_secreted_symbols=load_uniprot_secreted([p.uniprot_kw0964.path, p.uniprot_sl0243.path]),
        strict_lists=InstrumentLists(ukbppp=ukb_strict, decode=dec_strict, interval=itv_strict),
        inclusive_lists=InstrumentLists(ukbppp=ukb_incl, decode=dec_incl, interval=itv_incl),
        candidates=combine_candidates(gcat, fgn, ogw),
        pilot_keys=pilot_keys(classifications, label_map),
        pilot_indications=set(label_map.values()),
        karim_keys=karim_launched_keys(load_karim_launched(p.karim2026_workbook.path), disease),
        publication_dates=publication_dates,
        decode_diagnostic=decode_strict_diagnostic(st01, st02),
    )


def input_roots(inputs_root: Path, exp_dir: Path, repo: Path) -> list[tuple[str, Path]]:
    """Where INPUTS.tsv paths are relative to: the repository, the experiment directory baked into
    the Modal image, and the `pqtl-v8-inputs` volume (mounted at inputs_root's parent)."""
    return [("", repo), ("experiments/08_mechanism_interaction", exp_dir), ("pqtl-v8-inputs:", inputs_root.parent)]


class OutputDirNotThisRun(RuntimeError):
    """The output directory holds files of another run, of no identifiable run, or a completed run."""


def stage_code_sha256() -> str:
    """Digest of the code a stage A run executes: this file, the stage_a modules and the shared guard."""
    return code_sha256({**relative_files(HERE, [Path(__file__).resolve(), *sorted((HERE / "stage_a").glob("*.py"))]),
                        **guard_files()})


def claim_output_dir(out_dir: Path, run_token: str, repo_commit: str) -> Path:
    """Refuse `out_dir` unless it is absent, empty, or an unfinished directory of this run token
    (its run_info.json names the token and no MANIFEST.tsv is there); then write run_info.json."""
    info = out_dir / RUN_INFO_NAME
    if out_dir.exists() and any(out_dir.iterdir()):
        found = json.loads(info.read_text()).get("run_token") if info.is_file() else None
        if found != run_token:
            raise OutputDirNotThisRun(f"{out_dir} is not empty and belongs to run token {found!r}, not {run_token!r}; "
                                      "stage A writes only into an empty directory or its own")
        if (out_dir / MANIFEST_NAME).exists():
            raise OutputDirNotThisRun(f"{out_dir} already holds the finished output of run token {run_token}; "
                                      "a new run needs a new RUN_START token and an empty directory")
    out_dir.mkdir(parents=True, exist_ok=True)
    info.write_text(json.dumps({"stage": "A", "run_token": run_token, "repo_commit": repo_commit,
                                "code_sha256": stage_code_sha256()}, indent=1, sort_keys=True))
    return info


COUNT_CHECK_NAME = "registered_count_check.json"


def write_outputs(out: StageAOutputs, p: StageAPaths, out_dir: Path, roots: list[tuple[str, Path]],
                  run_info: Path, registered: dict[str, int]) -> Path:
    """Write the five tables and registered_count_check.json (`registered` against the counts of
    `out.hypotheses`); raise RegisteredCountMismatch there when any count differs. Otherwise write
    INPUTS.tsv and MANIFEST.tsv (which also lists `run_info`, the run_info.json `claim_output_dir`
    wrote, and the check). Returns the manifest path."""
    info = json.loads(run_info.read_text())
    tables = {"hypotheses.csv": (out.hypotheses, HYPOTHESIS_COLUMNS),
              "funnel.csv": (out.funnel, FUNNEL_COLUMNS),
              "outcome_gwas_selection.csv": (out.outcome_gwas_selection, OUTCOME_SELECTION_COLUMNS),
              "mechanism_crosstab.csv": (out.mechanism_crosstab, CROSSTAB_COLUMNS),
              "outcome_trait_coding.tsv": (out.outcome_trait_coding, TRAIT_CODING_COLUMNS)}
    for name, (df, cols) in tables.items():
        df[cols].to_csv(out_dir / name, index=False, sep="\t" if name.endswith(".tsv") else ",")
    check = count_check(out.hypotheses, registered)
    check_path = out_dir / COUNT_CHECK_NAME
    check_path.write_text(json.dumps(check, indent=1))
    if not check["all_match"]:
        differing = {k: (q["registered"], q["obtained"]) for k, q in check["quantities"].items() if not q["match"]}
        raise RegisteredCountMismatch(f"{len(differing)} registered count(s) differ, as (registered, obtained): {differing}; "
                                      f"see {check_path}. The output is not sealed; investigate, then rerun under a new "
                                      "run token")
    inputs = [InputRecord(name=k, path=portable_path(v.path, roots), sha256=sha256_file(v.path), source=v.source,
                          release=v.release, url_or_accession=v.url_or_accession,
                          download_date=datetime.fromtimestamp(v.path.stat().st_mtime, timezone.utc).date().isoformat())
              for k, v in ((k, getattr(p, k)) for k in StageAPaths.model_fields)]
    inputs.append(code_record("A", info["code_sha256"], info["repo_commit"]))
    return write_manifest(out_dir, [run_info, check_path, *(out_dir / n for n in tables)], info["code_sha256"], inputs)


def main(argv: list[str] | None = None, osf_fetch: Callable[[], dict] | None = None,
         repo_commit: str | None = None) -> Path:
    """`osf_fetch` replaces the OSF API query and `repo_commit` the commit lookup in tests; real
    runs leave both None."""
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs-root", type=Path, default=EXP_DIR,
                    help="layout of experiments/08_mechanism_interaction holding inputs/ and the v2 inputs/")
    ap.add_argument("--exp-dir", type=Path, default=EXP_DIR)
    ap.add_argument("--repo-root", type=Path, default=REPO)
    ap.add_argument("--paths-json", type=Path, default=None,
                    help="JSON object {input key: path or {path, source, release, url_or_accession}}")
    ap.add_argument("--out-dir", type=Path, default=HERE / "output")
    ap.add_argument("--run-token", required=True,
                    help="token of the PREREG.md log entry 'RUN_START stage=A token=<token>'")
    ap.add_argument("--prereg", type=Path, default=PREREG)
    args = ap.parse_args(argv)
    require_run(args.prereg, "A", args.run_token, osf_fetch=osf_fetch)
    dates = require_publication_dates(PUBLICATION_DATES)
    paths = resolve_paths(args.paths_json, args.inputs_root, args.exp_dir, args.repo_root)
    run_info = claim_output_dir(args.out_dir, args.run_token, repo_commit or run_commit(args.repo_root))
    manifest = write_outputs(run_stage_a(load_inputs(paths, dates)), paths, args.out_dir,
                             input_roots(args.inputs_root, args.exp_dir, args.repo_root), run_info, REGISTERED_COUNTS)
    print(f"wrote {args.out_dir} (manifest {manifest})")
    return manifest


if __name__ == "__main__":
    main()
