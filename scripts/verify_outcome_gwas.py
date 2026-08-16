"""Verify every outcome-GWAS accession against live OpenGWAS metadata.

Two accessions are already confirmed wrong by independent literature search:

    Chronic kidney disease       -> ieu-b-4874  is BLADDER CANCER
                                    (1,279 cases / 372,016 controls, European)
    Systemic lupus erythematosus -> ieu-a-1073  is COPPER (serum copper levels)

This resolves all of them, so the audit is complete rather than anecdotal.

Requires an OpenGWAS JWT. Free registration at https://opengwas.io/.
Set it before running:

    export OPENGWAS_JWT='your-token-here'
    python scripts/verify_outcome_gwas.py

Writes results/outcome_gwas_audit.csv and prints a match/mismatch table.
No analysis is re-run and no data file is modified.
"""

import csv
import json
import os
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
SOURCE = REPO / "code" / "classify_v34.py"
OUT = REPO / "results" / "outcome_gwas_audit.csv"
API = "https://api.opengwas.io/api/gwasinfo"

# Independently confirmed by literature search, 2026-08-16, before this ran.
KNOWN_WRONG = {
    "ieu-b-4874": "Bladder cancer (1,279 cases / 372,016 controls)",
    "ieu-a-1073": "Copper (serum copper levels)",
    "ieu-a-62": "Waist circumference",
}


def disease_map():
    """Parse DISEASE_GWAS out of the analysis source rather than duplicating it."""
    text = SOURCE.read_text()
    block = text.split("DISEASE_GWAS = {", 1)[1].split("}", 1)[0]
    out = {}
    for m in re.finditer(r'"([^"]+)"\s*:\s*(?:"([^"]+)"|None)', block):
        out[m.group(1)] = m.group(2)  # None -> None
    return out


def fetch(ids, token):
    payload = json.dumps({"id": ids})
    proc = subprocess.run(
        ["curl", "-sL", "-m", "60", "-X", "POST", API,
         "-H", "Content-Type: application/json",
         "-H", f"Authorization: Bearer {token}",
         "-d", payload],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"curl failed: {proc.stderr[:200]}")
    try:
        return json.loads(proc.stdout)
    except json.JSONDecodeError:
        raise RuntimeError(f"unexpected response: {proc.stdout[:300]}")


def looks_like(label, trait):
    """Loose containment check; a human still reads the table."""
    a = re.sub(r"[^a-z]+", " ", label.lower()).split()
    b = re.sub(r"[^a-z]+", " ", (trait or "").lower())
    hits = sum(1 for w in a if len(w) > 3 and w in b)
    return hits >= 1


def main():
    token = os.environ.get("OPENGWAS_JWT") or os.environ.get("GWAS_KEY")
    dm = disease_map()
    ids = sorted({v for v in dm.values() if v})
    print(f"{len(dm)} diseases, {len(ids)} distinct accessions\n")

    if not token:
        print("No OPENGWAS_JWT set. Showing what is already known:\n")
        for disease, acc in sorted(dm.items()):
            note = KNOWN_WRONG.get(acc, "")
            flag = "  <-- WRONG" if note else ""
            print(f"  {disease:<34} {acc or '(none)':<22} {note}{flag}")
        print("\nGet a free token at https://opengwas.io/ then re-run.")
        sys.exit(1)

    meta = fetch(ids, token)
    if isinstance(meta, dict) and "message" in meta:
        print("API error:", meta["message"][:300])
        sys.exit(1)
    by_id = {m["id"]: m for m in (meta if isinstance(meta, list) else meta.values())}

    rows, bad = [], 0
    print(f"{'disease (code label)':<34} {'accession':<20} {'actual trait':<38} verdict")
    print("-" * 110)
    for disease, acc in sorted(dm.items()):
        if not acc:
            print(f"{disease:<34} {'(unmapped)':<20} {'':<38} n/a")
            continue
        m = by_id.get(acc)
        if not m:
            verdict = "NOT IN INDEX"
            trait = ""
            bad += 1
        else:
            trait = m.get("trait", "")
            verdict = "ok" if looks_like(disease, trait) else "MISMATCH"
            if verdict == "MISMATCH":
                bad += 1
        print(f"{disease:<34} {acc:<20} {trait[:38]:<38} {verdict}")
        rows.append({
            "code_label": disease, "accession": acc, "actual_trait": trait,
            "ncase": (m or {}).get("ncase", ""), "ncontrol": (m or {}).get("ncontrol", ""),
            "sample_size": (m or {}).get("sample_size", ""),
            "population": (m or {}).get("population", ""),
            "year": (m or {}).get("year", ""), "author": (m or {}).get("author", ""),
            "verdict": verdict,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n{bad} problem accession(s). Written to {OUT}")


if __name__ == "__main__":
    main()
