"""Environment digest tests; they run in the Modal image (modal_stage_d.py::tests) and skip where PyMC is absent."""
import json
import os
import subprocess
from importlib.metadata import version
from pathlib import Path

import pytest

pytest.importorskip("pymc")

import run_stage_d  # noqa: E402
import setuptools  # noqa: E402, F401  (importing it extends sys.path with the packages it vendors)
from v8_manifest import sha256_file  # noqa: E402

from stage_d import environment  # noqa: E402
from stage_d.fingerprint import FINGERPRINT_NAME, StaleWork, canonical_sha256  # noqa: E402
from stage_d.pipeline import load_plan, prepare  # noqa: E402

TOKEN = "synthetic-token-0001"
STAGE_D = Path(__file__).resolve().parent.parent
FRESH = "import json; from stage_d.environment import collect; print(json.dumps(collect()))"
IMAGE = {
    "os_release": {"pretty_name": "Debian GNU/Linux 12 (bookworm)", "version_id": "12", "debian_version": "12.11"},
    "machine": "x86_64",
    "apt_packages": ["libopenblas-dev:amd64\t0.3.21+ds-4\tii ", "libopenblas0-pthread:amd64\t0.3.21+ds-4\tii "],
    "pytensor_blas": {"blas__ldflags": "-lopenblas", "libraries": {"openblas": {
        "soname": "libopenblas.so.0", "path": "/lib/x86_64-linux-gnu/libopenblas.so.0",
        "resolved": "/usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.21.so", "sha256": "a" * 64}}},
    "numpy_blas": {"blas": {"name": "scipy-openblas", "version": "0.3.30"}},
    "distributions": ["numpy==2.5.3", "pymc==5.28.5", "pytensor==2.38.2"],
}
REBUILT = [
    {"image_id": "im-synthetic02"},
    {"image_id": None},
    {"os_release": {**IMAGE["os_release"], "debian_version": "12.12"}},
    {"machine": "aarch64"},
    {"apt_packages": [IMAGE["apt_packages"][0], "libopenblas0-pthread:amd64\t0.3.21+ds-4+b1\tii "]},
    {"pytensor_blas": {"blas__ldflags": "-lopenblas", "libraries": {"openblas": {
        **IMAGE["pytensor_blas"]["libraries"]["openblas"], "sha256": "b" * 64}}}},
    {"numpy_blas": {"blas": {"name": "scipy-openblas", "version": "0.3.31"}}},
    {"distributions": ["numpy==2.5.3", "pymc==5.28.5", "pytensor==2.38.2", "threadpoolctl==3.6.0"]},
    {"distributions": ["numpy==2.5.3", "pymc==5.28.5", "pytensor==2.38.3"]},
]


def enter_image(monkeypatch, image_id="im-synthetic01", **changes):
    for name, value in {**IMAGE, **changes}.items():
        monkeypatch.setattr(environment, name, lambda value=value: value)
    if image_id is None:
        monkeypatch.delenv(environment.MODAL_IMAGE_ENV, raising=False)
    else:
        monkeypatch.setenv(environment.MODAL_IMAGE_ENV, image_id)


def tree(work) -> dict[str, bytes]:
    return {str(p.relative_to(work)): p.read_bytes() for p in sorted(work.rglob("*")) if p.is_file()}


def test_collect_assembles_the_collectors_and_the_image_id_only_where_it_is_set(monkeypatch):
    enter_image(monkeypatch)
    expected = {"os_release": IMAGE["os_release"], "machine": "x86_64", "apt_packages": IMAGE["apt_packages"],
                "blas": {"pytensor": IMAGE["pytensor_blas"], "numpy": IMAGE["numpy_blas"]},
                "distributions": IMAGE["distributions"]}
    assert environment.collect() == {**expected, "modal_image_id": "im-synthetic01"}
    assert environment.collect(env={}) == expected == environment.collect(env={"MODAL_IMAGE_ID": ""})
    assert environment.collect(env={"MODAL_IMAGE_ID": "im-other"})["modal_image_id"] == "im-other"


def test_a_run_resumed_in_the_same_image_is_accepted(sealed_stages, tmp_path, monkeypatch):
    root, expected = sealed_stages
    work = tmp_path / "work"
    enter_image(monkeypatch)
    fp = run_stage_d.run_fingerprint(expected, TOKEN)
    plan = prepare(root, expected, work, fp)
    stored = json.loads((work / FINGERPRINT_NAME).read_text())
    assert stored == fp.record() and stored["environment"] == environment.collect()
    assert stored["components"]["environment_sha256"] == canonical_sha256(environment.collect()) == fp.environment_sha256
    before = tree(work)
    enter_image(monkeypatch)                      # a second container of the same image
    again = run_stage_d.run_fingerprint(expected, TOKEN)
    assert again.digest == fp.digest
    assert load_plan(work, again).model_dump_json() == plan.model_dump_json()
    assert prepare(root, expected, work, again).model_dump_json() == plan.model_dump_json() and tree(work) == before


@pytest.mark.parametrize("change", REBUILT)
def test_a_run_resumed_in_a_rebuilt_image_with_another_environment_is_refused(sealed_stages, tmp_path, monkeypatch, change):
    root, expected = sealed_stages
    work = tmp_path / "work"
    enter_image(monkeypatch)
    fp = run_stage_d.run_fingerprint(expected, TOKEN)
    prepare(root, expected, work, fp)
    before = tree(work)
    enter_image(monkeypatch, **change)
    rebuilt = run_stage_d.run_fingerprint(expected, TOKEN)
    same_but_environment = {k: v for k, v in rebuilt.components().items() if k != "environment_sha256"}
    assert same_but_environment == {k: v for k, v in fp.components().items() if k != "environment_sha256"}
    assert rebuilt.environment_sha256 != fp.environment_sha256 and rebuilt.digest != fp.digest
    with pytest.raises(StaleWork, match=f"carries fingerprint {fp.digest}, this run is {rebuilt.digest}"):
        load_plan(work, rebuilt)
    with pytest.raises(StaleWork, match="not empty"):
        prepare(root, expected, work, rebuilt)
    assert tree(work) == before


def test_the_collectors_read_this_image_and_give_the_components_a_fresh_interpreter_gives():
    first = environment.collect()
    fresh = subprocess.run(["python", "-c", FRESH], cwd=STAGE_D, capture_output=True, text=True, check=True)
    assert json.loads(fresh.stdout.strip().splitlines()[-1]) == first == environment.collect()
    assert ("modal_image_id" in first) == bool(os.environ.get("MODAL_IMAGE_ID"))
    installed = {line.split("\t")[0].split(":")[0]: line.split("\t") for line in first["apt_packages"]}
    for package in ("libopenblas-dev", "libopenblas-pthread-dev", "libopenblas0-pthread"):
        assert installed[package][1] and installed[package][2].startswith("ii")
    blas = first["blas"]["pytensor"]
    assert blas["blas__ldflags"] == "-lopenblas" and set(blas["libraries"]) == {"openblas"}
    lib = blas["libraries"]["openblas"]
    assert lib["soname"].startswith("libopenblas.so") and lib["sha256"] == sha256_file(Path(lib["resolved"]))
    assert "blas" in first["blas"]["numpy"]
    for package in ("pymc", "pytensor", "numpy", "scipy"):
        assert f"{package}=={version(package)}" in first["distributions"]
