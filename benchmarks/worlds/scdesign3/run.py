"""Run the synthetic sweep (benchmarks/synthetic/sweep.py) on the scDesign3 world.

    python run.py [--world bands|dirichlet] [--seeds 3] [--workers 16] [--only SCEN ...] [--out DIR]

A thin wrapper around ``sweep.py --world scdesign3 --sample-mode {bands,dirichlet}``; the install itself
(``sim.simulate`` / ``sim.reference`` / ``sim.marker_sets`` / ``sim.T_LINEAGE`` / ``sim.WEIGHTS``, the
``rare_type`` levels, and the two edge cases that name toy genes) is ``benchmarks/worlds/__init__.py::install``.

``bands`` (default): one section, samples = horizontal bands of the real tissue (as the toy world, but the
bands differ in composition because the tissue is not striped). ``dirichlet``: six separate pieces whose
composition is drawn from Dirichlet(alpha0 x natural shares) (``world.DIRICHLET_ALPHA0`` = 20 unless
``--alpha0``), the multi-sample variant for the ISC regime question.
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


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--world", choices=["bands", "dirichlet"], default="bands")
    p.add_argument("--alpha0", type=float, default=None)
    p.add_argument("--seeds", default="3")
    p.add_argument("--workers", default="16")
    p.add_argument("--only", nargs="*", default=None)
    p.add_argument("--out", default=None)
    args = p.parse_args(argv)
    sargs = ["--world", "scdesign3", "--sample-mode", args.world, "--seeds", args.seeds,
             "--workers", args.workers, "--out", args.out or str(HERE / "results" / args.world)]
    if args.alpha0 is not None:
        sargs += ["--alpha0", str(args.alpha0)]
    if args.only:
        sargs += ["--only", *args.only]
    return sweep.main(sargs)


if __name__ == "__main__":
    main()
