"""The runtime environment a stage D run is bound to, collected inside the container.

`collect` returns the components; `fingerprint.canonical_sha256` of them is the environment digest,
one component of the run fingerprint (stage_d/fingerprint.py). The components:

    modal_image_id  the image the container runs, from MODAL_IMAGE_ID; present only where Modal
                    sets that variable
    os_release      PRETTY_NAME and VERSION_ID of /etc/os-release, and /etc/debian_version
    machine         platform.machine()
    apt_packages    `dpkg-query -W 'libopenblas*'`: package, architecture, version and status of
                    every OpenBLAS package the image installs
    blas.pytensor   pytensor.config.blas__ldflags and, per `-l<name>` in it, the soname the linker
                    finds, the file the loader cache maps it to, that file with symlinks resolved
                    and its sha256
    blas.numpy      the "Build Dependencies" of numpy.show_config (the BLAS and LAPACK NumPy was
                    built against)
    distributions   every distribution installed in the interpreter's site-packages directories
                    (site.getsitepackages), as `name==version`: what `pip freeze` lists

Every component is a property of the image. Values that vary between containers of one image, or
within one process, are left out, so every phase of a run in one image computes the same digest:
the kernel release; the CPU count; the libraries loaded in the process (threadpoolctl lists a
library only once it is loaded, and reports the thread count and the CPU-dependent OpenBLAS
kernel); and distributions found through sys.path, which grows as modules are imported (importing
setuptools, as PyTensor does when it compiles, adds the packages setuptools vendors).
"""
import os
import platform
import re
import site
import subprocess
from collections.abc import Mapping
from ctypes.util import find_library
from importlib.metadata import distributions as installed_distributions
from pathlib import Path

import numpy as np
import pytensor
from v8_manifest import sha256_file

MODAL_IMAGE_ENV = "MODAL_IMAGE_ID"
APT_PATTERN = "libopenblas*"
LDCONFIG = "/sbin/ldconfig"
LIB_FLAG_RE = re.compile(r"(?:^|\s)-l(\S+)")


class EnvironmentUnresolved(RuntimeError):
    """A component of the runtime environment could not be read."""


def os_release() -> dict[str, str]:
    info = platform.freedesktop_os_release()
    return {"pretty_name": info["PRETTY_NAME"], "version_id": info.get("VERSION_ID", ""),
            "debian_version": Path("/etc/debian_version").read_text().strip()}


def machine() -> str:
    return platform.machine()


def apt_packages() -> list[str]:
    """`package:architecture<TAB>version<TAB>status` of every package matching APT_PATTERN."""
    res = subprocess.run(["dpkg-query", "-W", "-f", "${Package}:${Architecture}\t${Version}\t${db:Status-Abbrev}\n",
                          APT_PATTERN], capture_output=True, text=True)
    lines = sorted(line.strip() for line in res.stdout.splitlines() if line.strip())
    if res.returncode != 0 or not lines:
        raise EnvironmentUnresolved(f"dpkg-query found no package matching {APT_PATTERN}: {res.stderr.strip()[-500:]}")
    return lines


def shared_library(name: str) -> dict[str, str]:
    """The library `-l<name>` links: its soname, the file the loader cache maps the soname to, that
    file with symlinks resolved, and the sha256 of its bytes."""
    soname = find_library(name)
    if soname is None:
        raise EnvironmentUnresolved(f"no shared library found for -l{name}")
    cache = subprocess.run([LDCONFIG, "-p"], capture_output=True, text=True, check=True).stdout
    paths = [line.rsplit(" => ", 1)[1].strip() for line in cache.splitlines()
             if line.strip().startswith(f"{soname} ") and " => " in line]
    if not paths:
        raise EnvironmentUnresolved(f"{soname} is not in the loader cache ({LDCONFIG} -p)")
    resolved = Path(paths[0]).resolve()
    return {"soname": soname, "path": paths[0], "resolved": str(resolved), "sha256": sha256_file(resolved)}


def pytensor_blas() -> dict:
    ldflags = pytensor.config.blas__ldflags
    return {"blas__ldflags": ldflags, "libraries": {name: shared_library(name) for name in LIB_FLAG_RE.findall(ldflags)}}


def numpy_blas() -> dict:
    return np.show_config(mode="dicts")["Build Dependencies"]


def distributions() -> list[str]:
    return sorted(f"{d.metadata['Name']}=={d.version}" for d in installed_distributions(path=site.getsitepackages()))


def collect(env: Mapping[str, str] | None = None) -> dict:
    """The environment components (module docstring). `env` defaults to the process environment."""
    env = os.environ if env is None else env
    out = {"os_release": os_release(), "machine": machine(), "apt_packages": apt_packages(),
           "blas": {"pytensor": pytensor_blas(), "numpy": numpy_blas()}, "distributions": distributions()}
    if env.get(MODAL_IMAGE_ENV):
        out["modal_image_id"] = env[MODAL_IMAGE_ENV]
    return out
