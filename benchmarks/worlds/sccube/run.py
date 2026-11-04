"""Run the synthetic sweep (benchmarks/synthetic/sweep.py) on the scCube world instead of the toy world.

    PYTHONPATH=src python benchmarks/worlds/sccube/run.py --pattern domains --seeds 3 --workers 16 \
        [--out <dir>] [--only label_noise spillover ...]

A thin wrapper around ``sweep.py --world sccube --pattern P``. The install
(``benchmarks/worlds/__init__.py::install``) replaces ``sim.simulate`` / ``sim.reference`` with
``world.simulate`` / ``world.reference`` (the scenario functions look ``simulate`` up in ``sim`` at call time, and
the worker processes are forked, so they inherit the patch), points ``sim.marker_sets`` / ``sim.T_LINEAGE`` at the
world's real markers (used by the panel check), and rewires the one edge case that names a toy gene
(``all_zero_marker_gene`` zeroes ``B_m0``) to zero the first B marker instead. Every scenario is supported (see
world.py for how each keyword is honoured). Without ``--out`` the results go to ``results/<pattern>``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]))                     # benchmarks/: the ``worlds`` package
import worlds  # noqa: E402

worlds.limit_threads()                                       # before numpy is imported (by sweep)
sys.path.insert(0, str(HERE.parents[1] / "synthetic"))
import sweep  # noqa: E402

PATTERNS = ["stripes", "domains", "fine", "mixed"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", choices=PATTERNS, default="domains")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--only", nargs="*", default=None, choices=list(sweep.sim.SCENARIOS))
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    args = ["--world", "sccube", "--pattern", a.pattern, "--seeds", str(a.seeds), "--workers", str(a.workers),
            "--out", a.out or str(HERE / "results" / a.pattern)]
    if a.only:
        args += ["--only", *a.only]
    return sweep.main(args)


if __name__ == "__main__":
    main()
