"""Run the synthetic sweep (``benchmarks/synthetic/sweep.py``) on the SRTsim world instead of the toy world.

    python benchmarks/worlds/srtsim/run.py [--seeds 3] [--workers 16] [--out benchmarks/worlds/srtsim/results]
                                           [--only label_noise ...] [--fixed-crops K]

A thin wrapper around ``sweep.py --world srtsim`` that also writes ``runtime.json``. The install
(``benchmarks/worlds/__init__.py::install``) replaces ``sim.simulate`` / ``sim.reference`` / ``sim.marker_sets`` /
``sim.T_LINEAGE`` with this world's and adapts the two sweep helpers that assume the toy world:

* edge case ``all_zero_marker_gene`` zeroes the toy gene ``B_m0``: here B's first marker;
* ``deconv_check`` builds 5 x 5 pseudo-spots in toy units (density 1 cell per unit^2, 25 cells per spot): here
  the coordinates (um) are rescaled to that density for the check only.

``--fixed-crops K`` keeps the crops of world seed K for every seed and redraws only the SRTsim counts (the
control for between-crop variance; ``world.FIXED_CROPS``).

Needs the cache (``SAM_SRTSIM_CACHE``, built by prepare.py + generate.R). Workers are forked, so the patched
module state is inherited.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))                     # benchmarks/: the ``worlds`` package
import worlds  # noqa: E402

worlds.limit_threads()                                       # before numpy is imported (by sweep)
sys.path.insert(0, str(HERE.parents[1] / "synthetic"))
import sweep  # noqa: E402


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--seeds", type=int, default=3)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--out", default=None)
    p.add_argument("--only", nargs="*", default=None)
    p.add_argument("--fixed-crops", type=int, default=None)
    a = p.parse_args(argv)
    out = a.out or str(worlds.default_out("srtsim", fixed_crops=a.fixed_crops))
    sargs = ["--world", "srtsim", "--seeds", str(a.seeds), "--workers", str(a.workers), "--out", out]
    if a.fixed_crops is not None:
        sargs += ["--fixed-crops", str(a.fixed_crops)]
    if a.only:
        sargs += ["--only", *a.only]
    t0 = time.time()
    df, sc = sweep.main(sargs)
    secs = time.time() - t0
    Path(out, "runtime.json").write_text(json.dumps(
        {"sweep_seconds": round(secs, 1), "seeds": a.seeds, "workers": a.workers,
         "scenarios": a.only or list(sweep.sim.SCENARIOS), "fixed_crops": a.fixed_crops}, indent=1) + "\n")
    print(f"sweep: {secs:.0f} s, {len(df)} values, {len(sc)} verdicts")
    return df, sc


if __name__ == "__main__":
    main()
