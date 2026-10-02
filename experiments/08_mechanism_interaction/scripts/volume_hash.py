"""Hash every file under a directory tree (logic for scripts/modal_verify_inputs.py).

Pattern: proteome-mr-claim-audit experiments/02_claims_meta/scripts/volume_hash.py.
"""

import hashlib
from pathlib import Path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def hash_tree(mount: Path, prefix: str) -> dict[str, dict]:
    """Map remote path (/<prefix>/...) to bytes and sha256 for every file under mount/prefix."""
    base = mount / prefix.strip("/")
    out = {}
    for p in sorted(base.rglob("*")):
        if p.is_file():
            remote = "/" + str(p.relative_to(mount))
            out[remote] = {"bytes": p.stat().st_size, "sha256": sha256(p)}
    return out


def compare(manifest: dict, remote: dict[str, dict]) -> dict:
    """Compare manifest entries with remote hashes; every manifest file must match exactly."""
    rows = []
    for f in manifest["files"]:
        r = remote.get(f["remote_path"])
        rows.append({
            "path": f["path"],
            "remote_path": f["remote_path"],
            "local_sha256": f["sha256"],
            "remote_sha256": r["sha256"] if r else None,
            "local_bytes": f["bytes"],
            "remote_bytes": r["bytes"] if r else None,
            "uploaded_at": f["uploaded_at"],
            "match": bool(r) and r["sha256"] == f["sha256"] and r["bytes"] == f["bytes"],
        })
    listed = {f["remote_path"] for f in manifest["files"]}
    return {
        "n_manifest": len(rows),
        "n_remote": len(remote),
        "n_match": sum(r["match"] for r in rows),
        "n_mismatch_or_missing": sum(not r["match"] for r in rows),
        "n_not_stamped_uploaded": sum(r["uploaded_at"] is None for r in rows),
        "remote_not_in_manifest": sorted(set(remote) - listed),
        "files": rows,
    }
