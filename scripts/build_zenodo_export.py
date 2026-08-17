"""Assemble the Zenodo deposit for the corrected analysis.

Collects the corrected results, the full pre-registration chain with its freeze hashes,
the analysis code, and the three EpiGraphDB catalog files that live outside this
repository. Including those last three is what makes the deposit self-contained: without
them 40 of 143 pairs cannot be regenerated from anything a reader can obtain.

No manuscript PDF is included. The manuscript still asserts a mechanism dissociation that
the corrected analysis does not support, so depositing it would archive a claim the
accompanying data contradict.

Usage:
    python scripts/build_zenodo_export.py            # build into zenodo_export/
    python scripts/build_zenodo_export.py --zip      # also write the archive
"""

import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess

REPO = pathlib.Path(__file__).resolve().parent.parent
OUT = REPO / "zenodo_export"
CATALOG = pathlib.Path.home() / "Documents/GitHub/transport-wrapper/DRUGS_EXPANDED/results/v34"

# (source, destination) — sources missing on disk are reported, not silently skipped.
FILES = [
    # corrected results
    ("results/v5_1/classification_v5_1.csv", "results/classification_v5_1.csv"),
    ("results/v5_1/evaluation_v5_1.json", "results/evaluation_v5_1.json"),
    ("results/v5_1/CHECKS_V5_2.md", "results/CHECKS_V5_2.md"),
    ("results/v5_1/MECHANISM_DIFFERENCE_TEST.md", "results/MECHANISM_DIFFERENCE_TEST.md"),
    ("results/v5_1/outcome_gwas_cache_v5_1.json", "results/outcome_gwas_cache_v5_1.json"),
    ("results/v5_1/ld_proxy_map_v5_1.json", "results/ld_proxy_map_v5_1.json"),
    ("results/OUTCOME_GWAS_TRACE.md", "results/OUTCOME_GWAS_TRACE.md"),
    ("results/outcome_gwas_audit.csv", "results/outcome_gwas_audit.csv"),
    ("results/replacement_gwas_candidates.csv", "results/replacement_gwas_candidates.csv"),
    # prior analysed set, for comparison
    ("results/v5/classification_v5.csv", "results/prior/classification_v5.csv"),
    ("results/v5/evaluation_v5.json", "results/prior/evaluation_v5.json"),
    ("results/v5/outcome_gwas_cache_v5.json", "results/prior/outcome_gwas_cache_v5.json"),
    ("results/v5/decode_cis_pqtl_v5.json", "results/prior/decode_cis_pqtl_v5.json"),
    ("results/v5/ld_proxy_map_v5.json", "results/prior/ld_proxy_map_v5.json"),
    # data
    ("data/frozen_candidates_v34.csv", "data/frozen_candidates_v34.csv"),
    ("data/adjudicated_v34.csv", "data/adjudicated_v34.csv"),
    ("data/classification_v34.csv", "data/classification_v34.csv"),
    # code
    ("code/disease_gwas.py", "code/disease_gwas.py"),
    ("code/classify_v34.py", "code/classify_v34.py"),
    ("code/expand_v5.py", "code/expand_v5.py"),
    ("code/expanded_taxonomy.py", "code/expanded_taxonomy.py"),
    ("scripts/rerun_v5_1.py", "scripts/rerun_v5_1.py"),
    ("scripts/verify_outcome_gwas.py", "scripts/verify_outcome_gwas.py"),
    ("scripts/find_replacement_gwas.py", "scripts/find_replacement_gwas.py"),
    ("scripts/stratum_difference_tests.py", "scripts/stratum_difference_tests.py"),
    # deviation log
    ("protocol/DEVIATION_LOG.md", "protocol/DEVIATION_LOG.md"),
]

# The files that make the deposit self-contained; they are not in the repository.
EXTERNAL = [
    (CATALOG / "v34_mr_catalog.csv", "data/epigraphdb/v34_mr_catalog.csv"),
    (CATALOG / "outcome_gwas_cache_v34.json", "data/epigraphdb/outcome_gwas_cache_v34.json"),
    (CATALOG / "evaluation_v34.json", "data/epigraphdb/evaluation_v34.json"),
]


def copy_protocols():
    """Every PRESPEC, amendment and correction, with its sibling hash file."""
    n = 0
    for src in sorted((REPO / "protocol").rglob("*.md")):
        if src.name == "DEVIATION_LOG.md":
            continue
        rel = src.relative_to(REPO / "protocol")
        dst = OUT / "protocol" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        n += 1
        sha = src.parent / f"{src.stem}_sha256.txt"
        if sha.exists():
            shutil.copy2(sha, dst.parent / sha.name)
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", action="store_true")
    args = ap.parse_args()

    if OUT.exists():
        shutil.rmtree(OUT)
    OUT.mkdir(parents=True)

    copied, missing = [], []
    for src, dst in FILES:
        s = REPO / src
        if not s.exists():
            missing.append(src)
            continue
        d = OUT / dst
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, d)
        copied.append(dst)

    for s, dst in EXTERNAL:
        if not s.exists():
            missing.append(str(s))
            continue
        d = OUT / dst
        d.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(s, d)
        copied.append(dst + "   [external, previously undeposited]")

    n_proto = copy_protocols()

    # deposit-level documents, authored in zenodo/ and copied in
    shutil.copy2(REPO / "zenodo" / "README.md", OUT / "README.md")
    tex = REPO / "zenodo" / "analysis_report.tex"
    shutil.copy2(tex, OUT / "analysis_report.tex")
    for _ in range(2):  # twice, so the table reference settles
        subprocess.run(["pdflatex", "-interaction=nonstopmode", "analysis_report.tex"],
                       cwd=OUT, capture_output=True)
    for junk in ("analysis_report.aux", "analysis_report.log", "analysis_report.out"):
        if (OUT / junk).exists():
            (OUT / junk).unlink()
    if not (OUT / "analysis_report.pdf").exists():
        missing.append("analysis_report.pdf (pdflatex failed)")

    # manifest with a hash per file, so the deposit can be checked against itself
    manifest = []
    for p in sorted(OUT.rglob("*")):
        if p.is_file() and p.name != "MANIFEST.sha256":
            h = hashlib.sha256(p.read_bytes()).hexdigest()
            manifest.append(f"{h}  {p.relative_to(OUT)}")
    (OUT / "MANIFEST.sha256").write_text("\n".join(manifest) + "\n")

    print(f"{len(copied)} files copied, {n_proto} protocol files, {len(manifest)} in manifest")
    if missing:
        print(f"\n{len(missing)} MISSING:")
        for m in missing:
            print(f"   {m}")

    if args.zip:
        archive = shutil.make_archive(str(REPO / "zenodo_export"), "zip", OUT)
        size = pathlib.Path(archive).stat().st_size / 1e6
        print(f"\nwrote {archive}  ({size:.1f} MB)")


if __name__ == "__main__":
    main()
