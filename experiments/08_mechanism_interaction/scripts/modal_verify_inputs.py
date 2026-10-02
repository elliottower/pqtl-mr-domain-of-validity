"""Modal wrapper: hash every uploaded V8 input on the volume and compare with the local manifest.

Run from experiments/08_mechanism_interaction:
    uv run --with modal==1.4.3 modal run scripts/modal_verify_inputs.py
Writes modal_inputs_verification.json beside the manifest. The volume is mounted read-only.
Pattern: proteome-mr-claim-audit experiments/02_claims_meta/scripts/modal_hash_volume.py.
"""

import json
from datetime import datetime, timezone
from pathlib import Path

import modal

from volume_hash import compare, hash_tree

ROOT = Path(__file__).resolve().parents[1]
VOLUME = "pqtl-v8-inputs"
MOUNT = Path("/vol")

image = modal.Image.debian_slim(python_version="3.12").add_local_python_source("volume_hash")
vol = modal.Volume.from_name(VOLUME).with_mount_options(read_only=True)
app = modal.App("pqtl-v8-verify-inputs", image=image)


@app.function(volumes={str(MOUNT): vol}, timeout=3600)
def hash_remote(prefix: str) -> dict[str, dict]:
    return hash_tree(MOUNT, prefix)


@app.local_entrypoint()
def main() -> None:
    manifest = json.loads((ROOT / "modal_inputs_manifest.json").read_text())
    remote = {}
    for unit in manifest["units"]:
        remote.update(hash_remote.remote(f"{manifest['remote_prefix']}/{unit}"))
    result = compare(manifest, remote)
    result["volume"] = VOLUME
    result["verified_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    result["manifest_hashed_at"] = manifest["hashed_at"]
    out = ROOT / "modal_inputs_verification.json"
    out.write_text(json.dumps(result, indent=1) + "\n")
    print(f"{result['n_match']}/{result['n_manifest']} match; "
          f"{result['n_mismatch_or_missing']} mismatched or missing; "
          f"{len(result['remote_not_in_manifest'])} remote files not in manifest")
