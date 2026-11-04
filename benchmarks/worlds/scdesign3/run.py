"""Run the synthetic sweep (benchmarks/synthetic/sweep.py) on the scDesign3 world.

    python run.py [--world bands|dirichlet] [--seeds 3] [--workers 16] [--only SCEN ...] [--out DIR]

``bands`` (default): one section, samples = horizontal bands of the real tissue (as the toy world, but the
bands differ in composition because the tissue is not striped). ``dirichlet``: six separate pieces whose
composition is drawn from Dirichlet(alpha0 x natural shares) (``world.DIRICHLET_ALPHA0`` = 20 unless
``--alpha0``), the multi-sample variant for the ISC regime question.

Nothing in benchmarks/synthetic is edited: this swaps the world into the ``sim`` module before the sweep
forks its workers (``sim.simulate`` / ``sim.reference`` / ``sim.marker_sets`` / ``sim.T_LINEAGE`` /
``sim.WEIGHTS``, the ``rare_type`` levels, and the two edge cases that name toy genes).
"""
from __future__ import annotations

import argparse
import functools
import multiprocessing
import os
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks" / "synthetic"))
sys.path.insert(0, str(HERE))

import numpy as np  # noqa: E402

import sim  # noqa: E402
import sweep  # noqa: E402
import world  # noqa: E402

_toy_edge = sweep._edge


def _edge(name):
    """sweep._edge with the world's genes where the toy one names ``B_m0`` / ``NegCtrl0-1``."""
    if name not in ("all_zero_marker_gene", "extra_marker_type"):
        return _toy_edge(name)
    a = sim.simulate(0, n_cells=600)
    ms = {k: list(v) for k, v in a.uns["marker_sets"].items()}
    if name == "all_zero_marker_gene":
        X = np.asarray(a.layers["counts"].todense())
        X[:, list(a.var_names).index(ms["B"][0])] = 0
        a.layers["counts"] = X
        a = sim.finish(a)
    else:
        ms["Ghost"] = list(a.var_names[a.var["control"].to_numpy()][:2])   # an absent type: controls only
    a.obs["label"] = a.obs["label"].astype(str).to_numpy()
    a.uns["marker_sets"] = ms
    return a


def install(sample_mode: str = "bands", alpha0: float = world.DIRICHLET_ALPHA0) -> None:
    """Swap the scDesign3 world into ``sim`` (and the two toy-gene edge cases into ``sweep``)."""
    world.params()                                   # load the cache before the workers fork
    sim.simulate = functools.partial(world.simulate, sample_mode=sample_mode, alpha0=alpha0)
    sim.reference = world.reference
    sim.marker_sets = world.marker_sets
    sim.T_LINEAGE = world.lineage()
    sim.WEIGHTS = world.natural_weights()
    fn, _, kind = sim.SCENARIOS["rare_type"]
    sim.SCENARIOS["rare_type"] = (fn, world.rare_levels(), kind)
    sweep._edge = _edge


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--world", choices=["bands", "dirichlet"], default="bands")
    p.add_argument("--alpha0", type=float, default=world.DIRICHLET_ALPHA0)
    p.add_argument("--seeds", default="3")
    p.add_argument("--workers", default="16")
    p.add_argument("--only", nargs="*", default=None)
    p.add_argument("--out", default=None)
    args = p.parse_args(argv)
    # the monkeypatches reach the workers only through fork (3.14 defaults to forkserver)
    multiprocessing.set_start_method("fork", force=True)
    install(args.world, args.alpha0)
    out = args.out or str(HERE / "results" / args.world)
    sargs = ["--seeds", args.seeds, "--workers", args.workers, "--out", out]
    if args.only:
        sargs += ["--only", *args.only]
    return sweep.main(sargs)


if __name__ == "__main__":
    main()
