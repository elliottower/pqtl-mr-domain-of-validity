"""The run fingerprint every stage D work file is bound to.

The work directory lives on a persistent volume and its file names (plan.json, designs/<id>.npz,
fits/<fit_id>/<attempt>.json, freq/<design_id>.json and the JSONL checkpoints under
freq_checkpoints/) repeat from run to run, so a name alone does not say which sealed inputs or
which code produced a file. The fingerprint is the sha256 of

    the frozen plan hash; the sealed A, B and C MANIFEST.tsv hashes; the stage D script digest
    (sha256 over run_stage_d.py, the stage_d modules and power_v9_fast.py); the installed package
    versions; the run token; the environment digest.

The environment digest is `canonical_sha256` of the runtime environment components collected in
the container (stage_d/environment.py: the Modal image id, the OS release, the machine
architecture, the installed OpenBLAS apt packages, the BLAS PyTensor and NumPy link, and every
installed Python distribution with its version). A work directory written in one image is
therefore refused in a rebuilt image whose environment differs, although the code, the package
pins and the run token are unchanged.

`claim` (prepare) writes it to <work>/FINGERPRINT.json, with the hashed components and the
environment components behind the environment digest, and refuses a non-empty work directory
that carries another fingerprint or none. `require` (every later phase) refuses unless the
directory carries this run's. plan.json, every fit attempt, final.json and every frequentist
result store the fingerprint, and the fit and frequentist results also store the sha256 of the
design file they were computed from; loaders refuse a mismatch and assembly checks each one
against the plan and the design files on disk.
"""
import hashlib
import json
from pathlib import Path

from pydantic import BaseModel, ConfigDict

FINGERPRINT_NAME = "FINGERPRINT.json"


class StaleWork(RuntimeError):
    """A work file or work directory written under another run fingerprint, or another design."""


def canonical_sha256(obj) -> str:
    """sha256 of `obj` as JSON with sorted keys and no whitespace."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


class RunFingerprint(BaseModel):
    model_config = ConfigDict(frozen=True)
    plan_sha256: str
    seals: dict[str, str]          # stage -> sealed MANIFEST.tsv sha256 (A, B, C)
    script_sha256: str
    pins: dict[str, str]           # package -> installed version
    run_token: str
    environment: dict              # runtime environment components (stage_d.environment.collect)

    @property
    def environment_sha256(self) -> str:
        return canonical_sha256(self.environment)

    def components(self) -> dict:
        """What the digest is taken over: every field, the environment as its digest."""
        return {**self.model_dump(exclude={"environment"}), "environment_sha256": self.environment_sha256}

    @property
    def digest(self) -> str:
        return canonical_sha256(self.components())

    def record(self) -> dict:
        """The digest, the hashed components and the environment components, as FINGERPRINT.json
        and results.json store them."""
        return {"fingerprint": self.digest, "components": self.components(), "environment": self.environment}


def carried(work_dir: Path) -> str | None:
    """The fingerprint <work>/FINGERPRINT.json carries, or None when the file is absent."""
    path = work_dir / FINGERPRINT_NAME
    return json.loads(path.read_text())["fingerprint"] if path.is_file() else None


def claim(work_dir: Path, fp: RunFingerprint) -> None:
    """Bind `work_dir` to `fp`: an empty or absent directory is claimed; one already carrying
    `fp` is accepted; a non-empty directory with another fingerprint, or with none, is refused."""
    found = carried(work_dir)
    if found == fp.digest:
        return
    if work_dir.is_dir() and any(work_dir.iterdir()):
        raise StaleWork(f"{work_dir} is not empty and carries fingerprint {found}, this run is {fp.digest}; "
                        "stage D does not write into another run's work directory")
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / FINGERPRINT_NAME).write_text(json.dumps(fp.record(), indent=2, sort_keys=True))


def require(work_dir: Path, fp: RunFingerprint) -> None:
    found = carried(work_dir)
    if found != fp.digest:
        raise StaleWork(f"{work_dir} carries fingerprint {found}, this run is {fp.digest}")


def check_record(record: dict, what: str, fingerprint: str, design_sha256: str) -> None:
    """A fit attempt, final fit or frequentist result must carry this run's fingerprint and the
    sha256 of the design file it is being used with."""
    got = (record.get("fingerprint"), record.get("design_sha256"))
    if got != (fingerprint, design_sha256):
        raise StaleWork(f"{what} carries (fingerprint, design sha256) {got}, this run needs "
                        f"{(fingerprint, design_sha256)}")
