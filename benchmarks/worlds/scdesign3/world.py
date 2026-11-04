"""A realistic synthetic world: scDesign3 fitted to a real Xenium breast section, same contract as sim.py.

``simulate()`` / ``reference()`` match ``benchmarks/synthetic/sim.py`` so every scenario of ``sim.SCENARIOS``
runs unchanged on it (``run.py`` swaps them in). What differs from the toy world:

* **Cells and layout are real.** A pool of 8,714 cells from a 1.3 x 1.5 mm crop of the 10x Xenium FFPE human
  breast section (replicate 1, Janesick et al. 2023, CC BY), six author cell types (GEO GSM7780153, supervised
  annotation) mapped to the toy names: CD4+ T -> T_CD4, CD8+ T -> T_CD8, B cells -> B, Macrophages_1 ->
  Myeloid, Stromal -> Fibro, Invasive_Tumor -> Epi. Cells of the other 14 author types are left out (holes in
  the tissue). A section is a random subsample of the pool (``n_cells``), at real centroids rescaled to one
  cell per unit^2 (the toy world's density), with real cell areas.
* **Counts are scDesign3 draws** (Song et al. Nat Biotechnol 2024): per gene a negative binomial GAM,
  ``log mu = cell_type + log_area + s(x, y, bs = "gp", k = 100)``, and a Gaussian copula per cell type
  (``generate.R``). The fit is cached (``params.h5ad``); this module samples from it the way
  ``scDesign3::simu_new`` does (multivariate normal with the type's copula correlation -> pnorm -> NB quantile
  at the cell's fitted mean and the gene's dispersion), so any seed is available without R.
* **Panel**: the 304 genes of the 313 detected in >= 0.5 % of crop cells, plus 28 real negative-control probes
  (``var['control']``) drawn as flat Poisson at each probe's crop mean (not fitted by scDesign3).
* **Markers** (``uns['marker_sets']``): derived from the fitted type means (``build.derive_markers``): 4 per
  type, T_CD4 and T_CD8 sharing 2 T-lineage genes as in the toy world.

Kwargs of ``sim.simulate`` and how they are honoured here:

* ``weights``: per-type shares, by stratified subsampling of the pool (``n_cells`` shrinks if a type runs
  out). ``None`` = the crop's natural composition (``natural_weights()``), not ``sim.WEIGHTS``.
* ``n_samples``: horizontal bands, as in the toy world - but the tissue is not striped, so the bands differ in
  composition. ``sample_mode="dirichlet"`` instead builds ``n_samples`` separate pieces, each a subsample of
  the pool with its own composition ~ Dirichlet(``alpha0`` x base weights) and its own count draw.
* ``infiltrate``: that fraction of T_CD8 cells moved to uniform random positions (their counts keep the
  location effect of where they came from).
* ``extra_bg``: uninformative genes ``xg*`` (level drawn from the fitted genes' mean levels, a mild per-type
  programme lognormal(0, 0.3), dispersion drawn from the fitted genes), from their own stream: counts of the
  other genes do not change.
* ``coexpr``: Myeloid cells get ``coexpr`` x Fibro's fitted mean added to their mean on Fibro's markers
  (section and reference).
* ``depth``: ``None`` keeps the real depth (median ~130 counts per cell on the fitted genes); a number rescales
  every mean so the pool mean library size equals it.

Every kwarg is honoured; nothing is silently ignored (unknown kwargs raise).
"""
from __future__ import annotations

import functools
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "synthetic"))
import sim  # noqa: E402  (finish(), TYPES, HIERARCHY - shared with the toy world)

CACHE = Path(os.environ.get("SAM_SCDESIGN3_CACHE", "/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/scdesign3"))
TYPES = sim.TYPES
HIERARCHY = sim.HIERARCHY
REF_DEPTH_FACTOR = 400 / 150           # the toy reference is 400 counts per cell vs 150 in the section
SUPPORTED_SCENARIOS = list(sim.SCENARIOS)   # all 15; rare_type with levels from rare_levels()
DIRICHLET_ALPHA0 = 20.0
_GAP = 10.0                            # units between Dirichlet pieces (kNN distances are ~1)


@functools.lru_cache(maxsize=1)
def params() -> dict:
    """The cached scDesign3 fit (``build.py pack``), loaded once per process."""
    p = ad.read_h5ad(CACHE / "params.h5ad")
    M = np.asarray(p.X, dtype=np.float64)
    types = p.obs["cell_type"].astype(str).to_numpy()
    chol = {}
    for t in TYPES:
        w, V = np.linalg.eigh(np.asarray(p.uns["copula"][t], dtype=np.float64))   # rmvnorm(method = "eigen")
        chol[t] = V * np.sqrt(np.clip(w, 0, None))
    genes = list(p.var_names)
    tmean = pd.DataFrame(M, columns=genes).groupby(types).mean().loc[TYPES]
    xy = p.obs[["x", "y"]].to_numpy(float)
    return {"M": M, "sigma": p.var["sigma"].to_numpy(float), "genes": genes, "types": types,
            "xy": xy, "lo": xy.min(0), "extent": xy.max(0) - xy.min(0),
            "area": p.obs["cell_area"].to_numpy(float), "chol": chol, "type_mean": tmean,
            "idx": {t: np.flatnonzero(types == t) for t in TYPES},
            "ctrl_names": [str(c) for c in p.uns["control_names"]],
            "ctrl_rate": np.asarray(p.uns["control_rate"], float),
            "markers": {k: [str(g) for g in v] for k, v in p.uns["marker_sets"].items()},
            "lineage": [str(g) for g in p.uns["T_lineage"]],
            "lib_mean": M.sum(1).mean()}


def marker_sets(per_type: int = 4) -> dict[str, list[str]]:
    if per_type != 4:
        raise ValueError("the scDesign3 world derives 4 markers per type (build.derive_markers)")
    return {t: list(v) for t, v in params()["markers"].items()}


def lineage() -> list[str]:
    """The two T-lineage markers shared by T_CD4 and T_CD8 (the world's ``sim.T_LINEAGE``)."""
    return list(params()["lineage"])


def natural_weights() -> np.ndarray:
    P = params()
    n = np.array([len(P["idx"][t]) for t in TYPES], float)
    return n / n.sum()


def rare_levels() -> list[float]:
    """``rare_type`` strengths: the crop's own T_CD8 share first (it cannot be raised by subsampling)."""
    return [round(float(natural_weights()[TYPES.index("T_CD8")]), 3), 0.08, 0.04, 0.02, 0.01]


# --------------------------------------------------------------------------- #
# sampling (the Python port of scDesign3::simu_new for family "nb", Gaussian copula)
# --------------------------------------------------------------------------- #
@functools.lru_cache(maxsize=16)
def _z_pool(seed: int, piece: int) -> np.ndarray:
    """Copula normals for every pool cell: z_c ~ N(0, R_type(c)). One stream per (seed, piece)."""
    P = params()
    rng = np.random.default_rng([seed, piece, 101])
    Z = np.empty(P["M"].shape, dtype=np.float32)
    for t in TYPES:
        i = P["idx"][t]
        Z[i] = rng.standard_normal((len(i), Z.shape[1])) @ P["chol"][t].T
    return Z


def _nb_quantile(U: np.ndarray, mu: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    """gamlss.dist::qNBI(U, mu, sigma): NB with size 1 / sigma (Poisson when sigma <= 1e-4); mu = 0 -> 0."""
    U = np.clip(U, 1e-12, 1 - 1e-12)
    sig = np.broadcast_to(sigma, mu.shape)
    out = np.zeros(mu.shape, dtype=np.float64)
    pos = mu > 0
    nb = pos & (sig > 1e-4)
    po = pos & ~nb
    size = 1.0 / sig[nb]
    out[nb] = stats.nbinom.ppf(U[nb], size, size / (size + mu[nb]))
    out[po] = stats.poisson.ppf(U[po], mu[po])
    return out


def _mean(cells, scale: float, coexpr: float) -> np.ndarray:
    P = params()
    M = P["M"][cells] * scale
    if coexpr:
        ms = P["markers"]["Fibro"]
        cols = [P["genes"].index(g) for g in ms]
        rows = np.flatnonzero(P["types"][cells] == "Myeloid")
        M[np.ix_(rows, cols)] += coexpr * P["type_mean"].loc["Fibro", ms].to_numpy() * scale
    return M


def _xg_params(seed: int, n: int):
    """Uninformative extra genes: base level, mild per-type programme, dispersion (own stream, as sim.world)."""
    P = params()
    r = np.random.default_rng(20_000 + seed)
    gene_mean = P["M"].mean(0)
    pick = r.integers(0, len(gene_mean), n)
    base = gene_mean[pick]
    mult = r.lognormal(0, 0.3, (len(TYPES), n))
    return pd.DataFrame(base * mult, index=TYPES, columns=[f"xg{i}" for i in range(n)]), P["sigma"][pick]


def _xg_counts(types, seed: int, n: int, stream, scale: float) -> tuple[np.ndarray, list[str]]:
    if not n:
        return np.zeros((len(types), 0)), []
    R, sig = _xg_params(seed, n)
    mu = R.loc[types].to_numpy() * scale
    return _nb_quantile(stream.random(mu.shape), mu, sig), list(R.columns)


def _depth_scale(depth) -> float:
    return 1.0 if depth is None else float(depth) / params()["lib_mean"]


def _choose(rng, n_cells: int, weights) -> np.ndarray:
    """Pool rows: a uniform subsample (natural composition), or stratified to ``weights``."""
    P = params()
    pool = len(P["types"])
    if weights is None:
        if n_cells > pool:
            raise ValueError(f"n_cells={n_cells} > the {pool} cells of the scDesign3 pool")
        return np.sort(rng.choice(pool, n_cells, replace=False))
    w = np.asarray(weights, float)
    w = w / w.sum()
    avail = np.array([len(P["idx"][t]) for t in TYPES])
    n = min(n_cells, int(np.floor(np.min(avail[w > 0] / w[w > 0]))))
    take = np.round(w * n).astype(int).clip(0, avail)
    return np.sort(np.concatenate([rng.choice(P["idx"][t], k, replace=False) for t, k in zip(TYPES, take)]))


def _piece(cells, seed: int, piece: int, scale: float, coexpr: float, extra_bg: int):
    """Counts (fitted genes + xg + controls) for pool rows ``cells`` from the (seed, piece) streams."""
    P = params()
    mu = _mean(cells, scale, coexpr)
    X = _nb_quantile(stats.norm.cdf(_z_pool(seed, piece)[cells].astype(np.float64)), mu, P["sigma"])
    Xg, xg_names = _xg_counts(P["types"][cells], seed, extra_bg, np.random.default_rng([seed, piece, 303]), scale)
    Xc = np.random.default_rng([seed, piece, 202]).poisson(P["ctrl_rate"], (len(cells), len(P["ctrl_rate"])))
    return np.hstack([X, Xg, Xc]).astype(np.float32), xg_names


def _layout(cells, scale_n: int):
    """Real centroids, shifted to 0 and rescaled to ``scale_n`` cells per crop area (1 cell / unit^2)."""
    P = params()
    s = np.sqrt(scale_n / np.prod(P["extent"]))
    return (P["xy"][cells] - P["lo"]) * s, P["extent"] * s


def simulate(seed: int = 0, n_cells: int = 2400, weights=None, n_samples: int = 6, infiltrate: float = 0.0,
             extra_bg: int = 0, coexpr: float = 0.0, depth=None, sample_mode: str = "bands",
             alpha0: float = DIRICHLET_ALPHA0, **kw) -> ad.AnnData:
    """A clean section from the scDesign3 world (see the module docstring for each kwarg)."""
    if kw:
        raise TypeError(f"scDesign3 world: unsupported simulate() kwargs {sorted(kw)}")
    P = params()
    rng = np.random.default_rng(seed)
    scale = _depth_scale(depth)
    if sample_mode == "bands":
        cells = _choose(rng, n_cells, weights)
        xy, ext = _layout(cells, len(cells))
        X, xg = _piece(cells, seed, 0, scale, coexpr, extra_bg)
        sample = (xy[:, 1] / ext[1] * n_samples).astype(int).clip(0, n_samples - 1)
        boxes = [(np.zeros(2), ext)] * len(cells)
        comp = None
    elif sample_mode == "dirichlet":
        base = natural_weights() if weights is None else np.asarray(weights, float) / np.sum(weights)
        crng = np.random.default_rng(seed + 300)
        avail = np.array([len(P["idx"][t]) for t in TYPES])
        n_s = n_cells // n_samples
        parts, comp = [], []
        for s in range(n_samples):
            pi = crng.dirichlet(alpha0 * base)
            k = np.minimum(crng.multinomial(n_s, pi), avail)
            c = np.sort(np.concatenate([rng.choice(P["idx"][t], m, replace=False) for t, m in zip(TYPES, k)]))
            xy_s, ext = _layout(c, n_s)
            off = np.array([s * (ext[0] + _GAP), 0.0])
            Xs, xg = _piece(c, seed, s, scale, coexpr, extra_bg)
            parts.append((c, xy_s + off, Xs, np.full(len(c), s), [(off, off + ext)] * len(c)))
            comp.append(dict(zip(TYPES, (k / max(k.sum(), 1)).round(4).tolist())))
        cells = np.concatenate([p[0] for p in parts])
        xy = np.vstack([p[1] for p in parts])
        X = np.vstack([p[2] for p in parts])
        sample = np.concatenate([p[3] for p in parts])
        boxes = [b for p in parts for b in p[4]]
    else:
        raise ValueError(f"sample_mode must be 'bands' or 'dirichlet', got {sample_mode!r}")
    types = P["types"][cells]
    if infiltrate:
        cd8 = np.flatnonzero(types == "T_CD8")
        ri = np.random.default_rng(seed + 7)
        move = ri.choice(cd8, int(round(infiltrate * len(cd8))), replace=False)
        for i in move:
            lo, hi = boxes[i]
            xy[i] = ri.uniform(lo, hi)
    n = len(cells)
    obs = pd.DataFrame({"truth": types, "label": types, "sample": [f"s{i}" for i in sample],
                        "cell_area": P["area"][cells]}, index=[f"c{i}" for i in range(n)])
    names = P["genes"] + xg + P["ctrl_names"]
    var = pd.DataFrame({"control": [g in set(P["ctrl_names"]) for g in names]}, index=names)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = xy
    a.uns["marker_sets"] = marker_sets()
    a.uns["hierarchy"] = dict(HIERARCHY)
    a.uns["world"] = {"seed": seed, "coexpr": coexpr, "extra_bg": extra_bg, "name": "scdesign3",
                      "sample_mode": sample_mode, "depth_scale": scale}
    if comp is not None:
        a.uns["world"]["alpha0"] = alpha0
        a.uns["world"]["sample_composition"] = pd.DataFrame(comp, index=[f"s{i}" for i in range(n_samples)])
    return sim.finish(a, seed)


def reference(seed: int = 0, n_per_type: int = 300, depth=None, coexpr: float = 0.0, extra_bg: int = 0,
              depth_factor: float = REF_DEPTH_FACTOR) -> ad.AnnData:
    """Matched dissociated reference from the same fit: ``n_per_type`` pool cells per type drawn with
    replacement, means x ``depth_factor`` (as the toy world's 400 vs 150), fresh copula draws, no controls,
    no coordinates. ``depth`` (absolute mean counts per cell) overrides ``depth_factor``."""
    P = params()
    rng = np.random.default_rng(50_000 + seed)
    scale = _depth_scale(depth) if depth is not None else depth_factor
    cells = np.concatenate([rng.choice(P["idx"][t], n_per_type, replace=True) for t in TYPES])
    types = np.repeat(TYPES, n_per_type)
    Z = np.vstack([rng.standard_normal((n_per_type, len(P["genes"]))) @ P["chol"][t].T for t in TYPES])
    X = _nb_quantile(stats.norm.cdf(Z), _mean(cells, scale, coexpr), P["sigma"])
    Xg, xg = _xg_counts(types, seed, extra_bg, np.random.default_rng([50_000 + seed, 303]), scale)
    a = ad.AnnData(X=sp.csr_matrix(np.hstack([X, Xg]).astype(np.float32)),
                   obs=pd.DataFrame({"cell_type": types}, index=[f"r{i}" for i in range(len(types))]),
                   var=pd.DataFrame(index=P["genes"] + xg))
    a.layers["counts"] = a.X.copy()
    return a
