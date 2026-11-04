"""Realistic simulation worlds for the synthetic sweep, and the one place that swaps them in.

    python benchmarks/synthetic/sweep.py --world {toy,scdesign3,srtsim,sccube} [--pattern P] [--sample-mode M]

Each world lives in ``benchmarks/worlds/<name>/world.py`` and follows the contract of
``benchmarks/synthetic/sim.py`` (``simulate`` / ``reference`` / ``marker_sets``, the six type names). ``install``
patches the ``sim`` module the sweep reads, and the two sweep helpers that hard-code the toy world, before the
sweep forks its workers. The workers must be forked (``sweep.main`` asks for the ``fork`` context when a world is
installed), so they inherit the patched module state; nothing is pickled.

This module imports nothing heavy at import time: ``limit_threads`` must run before numpy is first imported.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORLDS = ("toy", "scdesign3", "srtsim", "sccube")
# one BLAS / OpenMP thread per sweep worker; unset, a 128-core node starts 128 threads in every worker
THREAD_VARS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
               "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS")


def limit_threads() -> None:
    """Set every BLAS / OpenMP thread variable to 1 unless the caller set it. Call before numpy is imported."""
    if "numpy" in sys.modules:
        import warnings
        warnings.warn("limit_threads() after numpy was imported: the BLAS thread count is already fixed",
                      RuntimeWarning, stacklevel=2)
    for v in THREAD_VARS:
        os.environ.setdefault(v, "1")


def load(name: str):
    """Import ``benchmarks/worlds/<name>/world.py`` once, as module ``sam_world_<name>``."""
    if name not in WORLDS or name == "toy":
        raise ValueError(f"unknown world {name!r}; realistic worlds: {WORLDS[1:]}")
    mod = f"sam_world_{name}"
    if mod not in sys.modules:
        spec = importlib.util.spec_from_file_location(mod, HERE / name / "world.py")
        m = importlib.util.module_from_spec(spec)
        sys.modules[mod] = m
        spec.loader.exec_module(m)
    return sys.modules[mod]


def default_out(name: str, pattern: str | None = None, sample_mode: str | None = None,
                fixed_crops: int | None = None) -> Path:
    """Where a world's sweep writes when ``--out`` is not given (never the toy world's results)."""
    if name == "toy":
        return HERE.parent / "synthetic" / "results"
    if name == "scdesign3":
        return HERE / name / "results" / (sample_mode or "bands")
    if name == "sccube":
        return HERE / name / "results" / (pattern or load("sccube").PATTERN)
    if fixed_crops is not None:
        return HERE / name / "results" / f"fixed_crops_{fixed_crops}"
    return HERE / name / "results"


def _first_marker_edge(sim, toy_edge):
    """``sweep._edge`` for a realistic world: the two edge cases that name toy genes use the world's own.

    ``all_zero_marker_gene`` zeroes B's first marker (toy: ``B_m0``); ``extra_marker_type`` gives an absent type
    the first two negative-control probes (toy: ``NegCtrl0``, ``NegCtrl1``).
    """
    import numpy as np

    def _edge(name):
        if name not in ("all_zero_marker_gene", "extra_marker_type"):
            return toy_edge(name)
        a = sim.simulate(0, n_cells=600)
        ms = {k: list(v) for k, v in a.uns["marker_sets"].items()}
        if name == "all_zero_marker_gene":
            X = np.asarray(a.layers["counts"].todense())
            X[:, list(a.var_names).index(ms["B"][0])] = 0
            a.layers["counts"] = X
            a = sim.finish(a)
        else:
            ms["Ghost"] = list(a.var_names[a.var["control"].to_numpy()][:2])
        a.obs["label"] = a.obs["label"].astype(str).to_numpy()
        a.uns["marker_sets"] = ms
        return a

    return _edge


def _unit_density_deconv(sim, w, toy_deconv):
    """``sweep.deconv_check`` with coordinates rescaled to one cell per unit^2 (the toy's 5 x 5 spot = 25 cells)."""
    import numpy as np

    def _unit_density(seed, **kw):
        a = w.simulate(seed, **kw)
        xy = a.obsm["spatial"]
        area = sum(np.prod(np.ptp(xy[(a.obs["sample"] == s).to_numpy()], 0)) for s in a.obs["sample"].unique())
        a.obsm["spatial"] = xy * np.sqrt(a.n_obs / area)
        return a

    def deconv_check(seeds):
        sim.simulate = _unit_density
        try:
            return toy_deconv(seeds)
        finally:
            sim.simulate = w.simulate

    return deconv_check


def install(name: str, sim, sweep, *, seeds: int = 3, pattern: str | None = None, sample_mode: str = "bands",
            alpha0: float | None = None, fixed_crops: int | None = None):
    """Swap world ``name`` into the ``sim`` module and the toy-specific helpers of ``sweep``.

    Parameters
    ----------
    name
        One of ``WORLDS``; ``"toy"`` changes nothing and returns ``None``.
    sim, sweep
        The module objects the sweep runs with (``sweep`` may be ``__main__``).
    seeds
        Seeds the sweep will ask for; worlds with a per-seed cache load it here, in the parent, so the forked
        workers share it.
    pattern
        scCube only: ``stripes``, ``domains``, ``fine`` or ``mixed``.
    sample_mode, alpha0
        scDesign3 only: ``bands`` (horizontal bands of one section) or ``dirichlet`` (six pieces whose
        composition is drawn from Dirichlet(``alpha0`` x natural shares)).
    fixed_crops
        SRTsim only: keep the crops of this world seed for every sweep seed and redraw only the SRTsim counts
        (the control that separates between-crop variance from the metrics).

    Returns
    -------
    The world module, or ``None`` for the toy world.
    """
    if pattern is not None and name != "sccube":
        raise ValueError("--pattern applies to the sccube world only")
    if sample_mode != "bands" and name != "scdesign3":
        raise ValueError("--sample-mode applies to the scdesign3 world only")
    if fixed_crops is not None and name != "srtsim":
        raise ValueError("--fixed-crops applies to the srtsim world only")
    if name == "toy":
        return None
    w = load(name)
    # keep the toy helpers once, so a second install does not wrap a patched one
    toy_edge = sweep.__dict__.setdefault("_toy_edge", sweep._edge)
    toy_deconv = sweep.__dict__.setdefault("_toy_deconv_check", sweep.deconv_check)
    sweep._edge = _first_marker_edge(sim, toy_edge)
    sweep.deconv_check = toy_deconv
    if name == "scdesign3":
        import functools
        w.params()                                   # load the cache before the workers fork
        a0 = w.DIRICHLET_ALPHA0 if alpha0 is None else alpha0
        sim.simulate = functools.partial(w.simulate, sample_mode=sample_mode, alpha0=a0)
        sim.reference = w.reference
        sim.marker_sets = w.marker_sets
        sim.T_LINEAGE = w.lineage()
        sim.WEIGHTS = w.natural_weights()            # None = the crop's natural composition, not the toy's
        fn, _, kind = sim.SCENARIOS["rare_type"]
        sim.SCENARIOS["rare_type"] = (fn, w.rare_levels(), kind)
    elif name == "srtsim":
        w.FIXED_CROPS = fixed_crops
        w.install(sim)
        sweep.deconv_check = _unit_density_deconv(sim, w, toy_deconv)
        for s in range(seeds):                       # load once in the parent; forked workers inherit it
            w._section(s)
    elif name == "sccube":
        pattern = pattern or w.PATTERN
        if pattern not in w.PATTERNS:
            raise ValueError(f"--pattern must be one of {w.PATTERNS}")
        w.PATTERN = pattern
        os.environ["SAM_SCCUBE_PATTERN"] = pattern
        sim.simulate = w.simulate
        sim.reference = w.reference
        sim.marker_sets = w.marker_sets
        sim.T_LINEAGE = w.t_lineage()
    return w
