"""Pin more V8 inputs: hash named files, add them to modal_inputs_manifest.json, upload them.

scripts/upload_inputs_to_modal.py hashes and uploads whole directories and is run once. This adds
single files to the same manifest and the same volume, in the same order: every file is hashed and
the manifest written before any upload; each file is then uploaded with `modal volume put` and its
`uploaded_at` stamped, so an interrupted run leaves a manifest that says what reached the volume.
Entries already in the manifest are left as they are; a path already pinned with another sha256 is
refused. The manifest as it was is first copied to superseded/. Verify afterwards with
scripts/modal_verify_inputs.py, which hashes every pinned file on Modal.

Each file must already sit under experiments/08_mechanism_interaction/ in one of the manifest's
upload units, at the path it will have on the volume.

    uv run experiments/08_mechanism_interaction/scripts/add_inputs_to_modal.py \\
        inputs/ukbppp/sun2023_MOESM3_ESM.xlsx inputs/decode/ferkingstad2021_MOESM4_ESM.xlsx \\
        inputs/decode/decode_proteomics_folder_listing_2026-09-30.json
"""

import argparse
import hashlib
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "modal_inputs_manifest.json"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", help="paths relative to experiments/08_mechanism_interaction")
    rels = ap.parse_args().paths
    manifest = json.loads(MANIFEST.read_text())
    pinned = {f["path"]: f for f in manifest["files"]}
    new = []
    for rel in rels:
        p = ROOT / rel
        if not p.is_file() or p.is_symlink():
            raise RuntimeError(f"{p} is not a regular file")
        if not any(rel.startswith(unit + "/") for unit in manifest["units"]):
            raise RuntimeError(f"{rel} is in none of the upload units {manifest['units']}")
        entry = {"path": rel, "remote_path": f"{manifest['remote_prefix']}/{rel}", "bytes": p.stat().st_size,
                 "sha256": sha256(p), "uploaded_at": None}
        if rel in pinned:
            if pinned[rel]["sha256"] != entry["sha256"]:
                raise RuntimeError(f"{rel} is already pinned with sha256 {pinned[rel]['sha256']}, the file has {entry['sha256']}")
            continue
        new.append(entry)
    if not new:
        print("nothing to add: every file is already pinned")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    old = ROOT / "superseded" / f"modal_inputs_manifest_{stamp}.json"
    old.parent.mkdir(exist_ok=True)
    shutil.copyfile(MANIFEST, old)
    manifest["files"] += new
    manifest["n_files"] = len(manifest["files"])
    manifest["total_bytes"] = sum(f["bytes"] for f in manifest["files"])

    def write() -> None:
        MANIFEST.write_text(json.dumps(manifest, indent=1) + "\n")

    write()
    print(f"manifest: {len(new)} files added, {manifest['n_files']} in all; earlier manifest at {old}")
    for entry in new:
        subprocess.run(["modal", "volume", "put", manifest["volume"], str(ROOT / entry["path"]), entry["remote_path"]], check=True)
        entry["uploaded_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        write()
        print(f"{entry['uploaded_at']} uploaded {entry['path']} ({entry['bytes']} bytes, sha256 {entry['sha256']})")


if __name__ == "__main__":
    main()
