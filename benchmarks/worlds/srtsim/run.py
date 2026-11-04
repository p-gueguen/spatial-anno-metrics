"""Run the synthetic sweep (``benchmarks/synthetic/sweep.py``) on the SRTsim world instead of the toy world.

    python benchmarks/worlds/srtsim/run.py [--seeds 3] [--workers 16] [--out benchmarks/worlds/srtsim/results]
                                           [--only label_noise ...]

``sim.simulate`` / ``sim.reference`` / ``sim.marker_sets`` / ``sim.T_LINEAGE`` are replaced by this world's
(``world.install``); every scenario in ``sim.SCENARIOS`` then runs unchanged on the cached SRTsim sections.
Two sweep helpers assume the toy world and are adapted here, not in ``sweep.py``:

* edge case ``all_zero_marker_gene`` zeroes the toy gene ``B_m0``: here B's first marker;
* ``deconv_check`` builds 5 x 5 pseudo-spots in toy units (density 1 cell per unit^2, 25 cells per spot): here
  the coordinates (um) are rescaled to that density for the check only.

Needs the cache (``SAM_SRTSIM_CACHE``, built by prepare.py + generate.R). Workers are forked, so the patched
module state is inherited.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):  # one thread per worker, as sweep.py;
    os.environ.setdefault(_v, "1")                                         # must precede the numpy import

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks" / "synthetic"))
sys.path.insert(0, str(HERE))

import sim  # noqa: E402
import sweep  # noqa: E402
import world  # noqa: E402

SUPPORTED = list(sim.SCENARIOS)          # all 15: every simulate() knob a scenario uses is honoured (world.py)

_toy_edge = sweep._edge
_toy_deconv = sweep.deconv_check


def _edge(name):
    if name != "all_zero_marker_gene":
        return _toy_edge(name)
    a = sim.simulate(0, n_cells=600)
    X = np.asarray(a.layers["counts"].todense())
    X[:, list(a.var_names).index(a.uns["marker_sets"]["B"][0])] = 0
    a.layers["counts"] = X
    return sim.finish(a)


def _unit_density(seed, **kw):
    a = world.simulate(seed, **kw)
    xy = a.obsm["spatial"]
    area = sum(np.prod(np.ptp(xy[(a.obs["sample"] == s).to_numpy()], 0)) for s in a.obs["sample"].unique())
    a.obsm["spatial"] = xy * np.sqrt(a.n_obs / area)
    return a


def deconv_check(seeds):
    sim.simulate = _unit_density
    try:
        return _toy_deconv(seeds)
    finally:
        sim.simulate = world.simulate


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--out", default=str(HERE / "results"))
    p.add_argument("--only", nargs="*", default=None)
    a = p.parse_args(argv)
    mp.set_start_method("fork", force=True)
    world.install(sim)
    sweep._edge = _edge
    sweep.deconv_check = deconv_check
    for s in range(a.seeds):
        world._load(s)                        # load once in the parent; forked workers inherit it
    t0 = time.time()
    df, sc = sweep.main(["--seeds", str(a.seeds), "--workers", str(a.workers), "--out", a.out,
                         "--only", *(a.only or SUPPORTED)])
    secs = time.time() - t0
    Path(a.out, "runtime.json").write_text(json.dumps(
        {"sweep_seconds": round(secs, 1), "seeds": a.seeds, "workers": a.workers,
         "scenarios": a.only or SUPPORTED}, indent=1) + "\n")
    print(f"sweep: {secs:.0f} s, {len(df)} values, {len(sc)} verdicts")
    return df, sc


if __name__ == "__main__":
    main()
