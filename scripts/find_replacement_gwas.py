"""Find the OpenGWAS accession matching the source pre-specified in the V3.4 log.

Eight outcome-GWAS accessions in `code/classify_v34.py` are wrong: five point at an
unrelated trait, three are no longer in the index. The correct source consortium for each
was already named in `protocol/DEVIATION_LOG.md` in July, so this is a bug fix rather than
a post-hoc source choice: the code is being made to match what was pre-specified.

For each affected disease this lists every candidate study in the OpenGWAS index whose
trait matches, with author, year, sample size and population, so the accession named by
the deviation log can be identified. It selects nothing on its own and touches no data.

Requires the token in `.env2`:

    python scripts/find_replacement_gwas.py

Writes results/replacement_gwas_candidates.csv.
"""

import csv
import json
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
ENV = REPO / ".env2"
OUT = REPO / "results" / "replacement_gwas_candidates.csv"
API = "https://api.opengwas.io/api/gwasinfo"

# Disease -> (trait keywords to match, source named in DEVIATION_LOG.md).
# The keywords are deliberately loose; a human reads the table and picks.
TARGETS = {
    "Autism spectrum disorder": (["autism"], "iPSYCH-PGC"),
    "Chronic kidney disease": (["chronic kidney", "kidney disease", "eGFR"], "CKDGen (overlap-flagged)"),
    "Systemic lupus erythematosus": (["lupus"], "Bentham 2015"),
    "Melanoma": (["melanoma"], "GenoMEL"),
    "Anorexia nervosa": (["anorexia"], "PGC AN"),
    "Glioma": (["glioma"], "GICC"),
    "Pancreatic cancer": (["pancrea"], "PanScan"),
}


def token():
    m = re.search(r"OPEN_GWAS_TOKEN\s*=\s*(.+)", ENV.read_text())
    if not m:
        sys.exit(f"OPEN_GWAS_TOKEN not found in {ENV}")
    return m.group(1).strip().strip("'\"")


def fetch_index(tok):
    """Full study index. One call; the endpoint returns every study without an id list."""
    proc = subprocess.run(
        ["curl", "-sL", "-m", "180", API, "-H", f"Authorization: Bearer {tok}"],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        sys.exit(f"curl failed: {proc.stderr[:200]}")
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        sys.exit(f"unexpected response: {proc.stdout[:300]}")
    if isinstance(data, dict) and "message" in data:
        sys.exit(f"API error: {data['message'][:300]}")
    return data if isinstance(data, list) else list(data.values())


def main():
    studies = fetch_index(token())
    print(f"{len(studies)} studies in the index\n")

    rows = []
    for disease, (keywords, expected) in TARGETS.items():
        hits = [
            s for s in studies
            if any(k.lower() in (s.get("trait") or "").lower() for k in keywords)
        ]
        hits.sort(key=lambda s: -(s.get("sample_size") or 0))
        print("=" * 100)
        print(f"{disease}   — deviation log names: {expected}   ({len(hits)} candidates)")
        print(f"  {'accession':<22} {'trait':<34} {'author':<16} {'year':<6} {'N':>9} {'pop'}")
        for s in hits[:12]:
            print(f"  {s.get('id',''):<22} {(s.get('trait') or '')[:34]:<34} "
                  f"{(str(s.get('author') or ''))[:16]:<16} {str(s.get('year') or ''):<6} "
                  f"{str(s.get('sample_size') or ''):>9} {s.get('population') or ''}")
            rows.append({
                "disease": disease, "expected_source": expected,
                "accession": s.get("id", ""), "trait": s.get("trait", ""),
                "author": s.get("author", ""), "year": s.get("year", ""),
                "ncase": s.get("ncase", ""), "ncontrol": s.get("ncontrol", ""),
                "sample_size": s.get("sample_size", ""),
                "population": s.get("population", ""), "consortium": s.get("consortium", ""),
            })
        print()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} candidates written to {OUT}")


if __name__ == "__main__":
    main()
