"""Check that PyTensor links against a BLAS in the Modal images of stage D and power_v9.

Run inside the image. A fresh interpreter imports PyTensor with every warning shown, reads
`pytensor.config.blas__ldflags`, and compiles a matrix-vector and a matrix-matrix product; with
BLAS linked, PyTensor's rewrites use the C BLAS ops (CGemv, Dot22 with dgemm). PyTensor 2.38.2
warns "PyTensor could not link to a BLAS installation" when auto-detection fails; older releases
said "BLAS functions will be unavailable". Either text in stderr fails the check.
"""
import json
import os
import subprocess

PROBE = r"""
import json
import numpy as np
import pytensor
import pytensor.tensor as pt

A, v, B = pt.dmatrix("A"), pt.dvector("v"), pt.dmatrix("B")
mv = pytensor.function([A, v], pt.dot(A, v))
mm = pytensor.function([A, B], pt.dot(A, B))
rng = np.random.default_rng(0)
a, x, b = rng.normal(size=(50, 40)), rng.normal(size=40), rng.normal(size=(40, 30))
print(json.dumps({
    "pytensor": pytensor.__version__,
    "blas__ldflags": pytensor.config.blas__ldflags,
    "cxx": pytensor.config.cxx,
    "matvec_ops": [type(n.op).__name__ for n in mv.maker.fgraph.toposort()],
    "matmul_ops": [type(n.op).__name__ for n in mm.maker.fgraph.toposort()],
    "matvec_max_abs_error": float(np.abs(mv(a, x) - a @ x).max()),
    "matmul_max_abs_error": float(np.abs(mm(a, b) - a @ b).max()),
}))
"""
WARNING_TEXTS = ("could not link to a BLAS", "BLAS functions will be unavailable")


def check() -> dict:
    res = subprocess.run(["python", "-W", "always", "-c", PROBE], capture_output=True, text=True,
                         env={**os.environ, "PYTHONWARNINGS": "always"})
    if res.returncode != 0:
        raise RuntimeError(f"probe failed: {res.stderr[-2000:]}")
    out = json.loads(res.stdout.strip().splitlines()[-1])
    warned = [t for t in WARNING_TEXTS if t in res.stderr]
    dpkg = subprocess.run(["dpkg-query", "-W", "libopenblas-dev", "libopenblas-pthread-dev", "libopenblas0-pthread"],
                          capture_output=True, text=True).stdout.split("\n")
    out.update({"dpkg": [line for line in dpkg if line], "env_PYTENSOR_FLAGS": os.environ.get("PYTENSOR_FLAGS"),
                "env_OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"),
                "blas_warning_texts_found": warned, "stderr_tail": res.stderr[-2000:],
                "passes": bool(out["blas__ldflags"]) and not warned and "CGemv" in out["matvec_ops"]
                and out["matvec_max_abs_error"] < 1e-10 and out["matmul_max_abs_error"] < 1e-10})
    return out
