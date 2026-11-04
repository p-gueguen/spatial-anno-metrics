"""Run the synthetic sweep (benchmarks/synthetic/sweep.py) on the scCube world instead of the toy world.

    PYTHONPATH=src python benchmarks/worlds/sccube/run.py --pattern domains --seeds 3 --workers 16 \
        --out <dir> [--only label_noise spillover ...]

Replaces ``sim.simulate`` / ``sim.reference`` with ``world.simulate`` / ``world.reference`` (the scenario
functions look ``simulate`` up in ``sim`` at call time, and the worker processes are forked, so they inherit the
patch), points ``sim.marker_sets`` / ``sim.T_LINEAGE`` at the world's real markers (used by the panel check),
and rewires the one edge case that names a toy gene (``all_zero_marker_gene`` zeroes ``B_m0``) to zero the first
B marker instead. Every scenario is supported (see world.py for how each keyword is honoured).
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):  # before numpy loads, as sweep.py does:
    os.environ.setdefault(_v, "1")         # otherwise BLAS starts one thread per core in every forked worker

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
for p in (ROOT / "src", ROOT / "benchmarks" / "synthetic", HERE):
    sys.path.insert(0, str(p))

import numpy as np  # noqa: E402

import sim  # noqa: E402
import sweep  # noqa: E402
import world  # noqa: E402

SUPPORTED = list(sim.SCENARIOS)
_toy_edge = sweep._edge


def _edge(name):
    if name != "all_zero_marker_gene":
        return _toy_edge(name)
    a = sim.simulate(0, n_cells=600)
    X = np.asarray(a.layers["counts"].todense())
    X[:, list(a.var_names).index(a.uns["marker_sets"]["B"][0])] = 0
    a.layers["counts"] = X
    a = sim.finish(a)
    a.obs["label"] = a.obs["label"].astype(str).to_numpy()
    a.uns["marker_sets"] = {k: list(v) for k, v in a.uns["marker_sets"].items()}
    return a


def install(pattern: str):
    world.PATTERN = pattern
    os.environ["SAM_SCCUBE_PATTERN"] = pattern
    sim.simulate = world.simulate
    sim.reference = world.reference
    sim.marker_sets = world.marker_sets
    sim.T_LINEAGE = world.t_lineage()
    sweep._edge = _edge


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", choices=world.PATTERNS, default="domains")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--only", nargs="*", default=None, choices=SUPPORTED)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if mp.get_start_method(allow_none=True) not in (None, "fork"):
        raise RuntimeError("the monkeypatch reaches the workers only with the 'fork' start method")
    mp.set_start_method("fork", force=True)
    install(a.pattern)
    args = ["--seeds", str(a.seeds), "--workers", str(a.workers), "--out", a.out]
    if a.only:
        args += ["--only", *a.only]
    return sweep.main(args)


if __name__ == "__main__":
    main()
