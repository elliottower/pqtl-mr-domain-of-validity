"""Hash the pinned V8 inputs, write a manifest, and upload them to a Modal volume.

Every file under the upload units below is hashed first and the manifest is written before any
upload. Each unit is then uploaded with `modal volume put`, and the manifest is rewritten with
`uploaded_at` for the files of that unit, so an interrupted run leaves a manifest that says what
reached the volume. Pattern: proteome-mr-claim-audit experiments/02_claims_meta/scripts/
upload_inputs_to_modal.py.

Remote layout: /08_mechanism_interaction/<path relative to experiments/08_mechanism_interaction>.
The volume is created beforehand with `modal volume create pqtl-v8-inputs`.

    uv run experiments/08_mechanism_interaction/scripts/upload_inputs_to_modal.py
"""

import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VOLUME = "pqtl-v8-inputs"
REMOTE_PREFIX = "/08_mechanism_interaction"
MANIFEST = ROOT / "modal_inputs_manifest.json"
UNITS = ("inputs/decode", "inputs/karim2026", "inputs/ot_26.09", "inputs/ukbppp",
         "feasibility/v2_all_indications/inputs")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def upload_units() -> list[Path]:
    units = [ROOT / u for u in UNITS]
    missing = [u for u in units if not u.is_dir()]
    if missing:
        raise RuntimeError(f"upload units missing: {missing}")
    stray = [p for p in (ROOT / "inputs").iterdir() if p.is_file()]
    if stray:
        raise RuntimeError(f"files directly under inputs/ are not handled: {stray}")
    links = [p for u in units for p in u.rglob("*") if p.is_symlink()]
    if links:
        raise RuntimeError(f"symlinks are not uploaded: {links}")
    return units


def write(manifest: dict) -> None:
    MANIFEST.write_text(json.dumps(manifest, indent=1) + "\n")


def main() -> None:
    units = upload_units()
    files = []
    for unit in units:
        for p in sorted(unit.rglob("*")):
            if p.is_file() and not p.is_symlink():
                rel = str(p.relative_to(ROOT))
                files.append({
                    "path": rel,
                    "remote_path": f"{REMOTE_PREFIX}/{rel}",
                    "bytes": p.stat().st_size,
                    "sha256": sha256(p),
                    "uploaded_at": None,
                })
    manifest = {
        "volume": VOLUME,
        "remote_prefix": REMOTE_PREFIX,
        "local_root": "experiments/08_mechanism_interaction",
        "units": list(UNITS),
        "hashed_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "n_files": len(files),
        "total_bytes": sum(f["bytes"] for f in files),
        "files": files,
    }
    write(manifest)
    print(f"manifest: {len(files)} files, {manifest['total_bytes']} bytes")

    for unit in units:
        rel = str(unit.relative_to(ROOT))
        remote = f"{REMOTE_PREFIX}/{rel}"
        subprocess.run(["modal", "volume", "put", VOLUME, str(unit), remote], check=True)
        stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
        for f in files:
            if f["path"].startswith(rel + "/"):
                f["uploaded_at"] = stamp
        write(manifest)
        print(f"{stamp} uploaded {rel}")


if __name__ == "__main__":
    main()
