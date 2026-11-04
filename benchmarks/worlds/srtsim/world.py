"""An SRTsim world: realistic spatial counts fitted to a real Xenium breast section, same contract as ``sim.py``.

The expensive part runs once and is cached (``prepare.py crops`` -> ``generate.R`` -> ``prepare.py cache``, see
``results/README.md``). For each world seed the cache holds six 400 x 400 um crops of the 10x Xenium FFPE human
breast section (Janesick et al. 2023, Rep1), each crop one sample. On every crop SRTsim (Zhu, Shang & Zhou,
Genome Biol 2023) fitted a per-gene, per-cell-type marginal (Poisson / NB / ZIP / ZINB, ``auto_choose``) in its
domain-based reference mode, and drew new counts that keep the rank order of the real counts at each real cell
(``srtsim_count``). Cell types are the authors' supervised annotation (GEO GSM7780153) mapped onto the six sim
types; other cells are left out. ``simulate`` and ``reference`` only read that cache and apply the scenario knobs.

``simulate(seed, n_cells, weights, n_samples, infiltrate, extra_bg, coexpr)`` -> AnnData with ``obs['truth']``,
``obs['label']``, ``obs['sample']``, ``obs['cell_area']``, ``layers['counts']`` (sparse), ``obsm['spatial']`` (um,
crops stacked in y with a 100 um gap), ``obsm['X_pca']``, ``var['control']`` (the 28 negative-control probes,
renamed ``NegCtrl0..27``; the 10x id is in ``var['feature_id']``), ``uns['marker_sets']`` (derived from the
simulated counts, ``marker_sets.json`` in the cache), ``uns['hierarchy']`` and ``uns['world']``.

How the knobs are honoured (all 15 sweep scenarios run):

* ``n_cells``: each sample keeps its ``n_cells // n_samples`` cells nearest the crop centre (a square window).
* ``weights``: read relative to ``sim.WEIGHTS``. A type's keep fraction is ``(w / WEIGHTS) / max(w / WEIGHTS)``,
  cells dropped at random, so ``weights=sim.WEIGHTS`` is the natural composition and ``rare_type`` shrinks
  T_CD8 by the same factors as in the toy (15 -> 1 % there; natural share -> share / 15 here).
* ``infiltrate``: that fraction of T_CD8 cells moved to uniform positions inside their own crop.
* ``extra_bg``: ``xg*`` genes, each drawn i.i.d. per cell from SRTsim's tissue-level fit of a random real gene
  of the same crop (no type programme, no spatial order); counts on the other genes do not change.
* ``coexpr``: Myeloid cells get Poisson(coexpr x mean Fibro count of each Fibro marker in the crop x the cell's
  relative depth) extra counts on Fibro's markers; the reference gets the same at its depth.
* ``depth``: not honoured (``None`` only): the depth is the section's. No scenario sets it.
* ``n_samples`` <= 6 (the crops). A seed outside the cache raises.

``reference(seed, n_per_type, depth, coexpr, extra_bg)``: a matched dissociated reference from the same fits:
for each type, real cells of that type are drawn with replacement across the six crops, and every gene is drawn
from the crop's fitted marginal with the mean scaled by ``rr`` (``depth`` / the section's mean depth; default
400 / 150 as in the toy) and assigned by the rank of the drawn real cells' counts - ``srtsim_count_single``
applied to resampled cells. No coordinates, no controls. Seeds outside the cache use the fit of ``seed % K``
with the draw seeded by ``seed`` (``sweep.conformal_check`` asks for ``seed + 100``).
"""
from __future__ import annotations

import functools
import json
import os
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "synthetic"))
import sim  # noqa: E402  (TYPES, WEIGHTS, HIERARCHY, finish: the toy world's contract)

CACHE = Path(os.environ.get("SAM_SRTSIM_CACHE", "/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/srtsim"))
TYPES = sim.TYPES
HIERARCHY = sim.HIERARCHY
# authors' supervised clusters (GSM7780153_Xenium_R1_Fig1-5_supervised.csv) -> sim types; the rest are left out
TYPEMAP = {"CD4+_T_Cells": "T_CD4", "CD8+_T_Cells": "T_CD8", "B_Cells": "B",
           "Macrophages_1": "Myeloid", "Macrophages_2": "Myeloid", "IRF7+_DCs": "Myeloid", "LAMP3+_DCs": "Myeloid",
           "Stromal": "Fibro",
           "Invasive_Tumor": "Epi", "Prolif_Invasive_Tumor": "Epi", "DCIS_1": "Epi", "DCIS_2": "Epi"}
N_CROPS = 6
CROP_UM = 400.0
GAP_UM = 100.0
REF_RR = 400.0 / 150.0          # the toy's reference / section depth ratio


@functools.lru_cache(maxsize=None)
def _seeds() -> tuple[int, ...]:
    return tuple(sorted(int(p.stem.split("_")[1]) for p in CACHE.glob("seed_*.h5ad")))


@functools.lru_cache(maxsize=None)
def _load(seed: int) -> ad.AnnData:
    p = CACHE / f"seed_{seed}.h5ad"
    if not p.exists():
        raise ValueError(f"SRTsim world seed {seed} is not in the cache {CACHE} (have {list(_seeds())}); "
                         "run prepare.py crops / generate.R / prepare.py cache with more seeds")
    return ad.read_h5ad(p)


@functools.lru_cache(maxsize=None)
def marker_sets() -> dict[str, list[str]]:
    return json.loads((CACHE / "marker_sets.json").read_text())["marker_sets"]


@functools.lru_cache(maxsize=None)
def _t_lineage() -> list[str]:
    return json.loads((CACHE / "marker_sets.json").read_text())["t_lineage"]


def _dense(X) -> np.ndarray:
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def _window(w: ad.AnnData, n_samples: int, n_cells: int) -> np.ndarray:
    """Indices into the cache: per crop, the ``n_cells // n_samples`` cells nearest the crop centre."""
    per = max(1, n_cells // n_samples)
    keep = []
    for c in range(n_samples):
        ids = np.flatnonzero(w.obs["crop"].to_numpy() == c)
        loc = w.obs[["x_local", "y_local"]].to_numpy()[ids]
        d = np.abs(loc - CROP_UM / 2).max(1)                      # Chebyshev: a square window
        keep.append(ids[np.argsort(d, kind="stable")[:per]])
    return np.concatenate(keep)


def _keep_by_weights(types: np.ndarray, weights, rng) -> np.ndarray:
    f = np.asarray(weights, float) / sim.WEIGHTS
    k = dict(zip(TYPES, f / f.max()))
    keep = np.ones(len(types), bool)
    for t in TYPES:
        ids = np.flatnonzero(types == t)
        drop = rng.choice(ids, len(ids) - int(round(k[t] * len(ids))), replace=False)
        keep[drop] = False
    return keep


def _draw(rng, pi0, theta, mu, n) -> np.ndarray:
    """SRTsim's draw: rbinom(n, 1, 1 - pi0) * rnbinom(n, size = theta, mu = mu); theta = Inf is Poisson."""
    if not np.isfinite(mu) or mu <= 0:
        return np.zeros(n)
    on = rng.random(n) >= pi0
    if np.isfinite(theta):
        x = rng.negative_binomial(theta, theta / (theta + mu), n)
    else:
        x = rng.poisson(mu, n)
    return on * x


def _rank_assign(rng, sim_vals: np.ndarray, real: np.ndarray) -> np.ndarray:
    """``sort(sim)[rank(real, ties.method = "random")]``."""
    order = np.lexsort((rng.random(len(real)), real))           # ascending real, random among ties
    out = np.empty_like(sim_vals)
    out[order] = np.sort(sim_vals)
    return out


def _coexpr_add(rng, X, types, crops, w, genes, s, rr=1.0):
    """Myeloid truly expresses Fibro's markers at fraction s of Fibro's mean level in the same crop."""
    col = {g: i for i, g in enumerate(genes)}
    fib = [g for g in marker_sets()["Fibro"] if g in col]
    gi = [list(w.var_names).index(g) for g in fib]
    mean_fib = w.uns["dom_mean"][:, TYPES.index("Fibro"), :][:, gi]          # crops x markers
    lib = X.sum(1)
    for c in np.unique(crops):
        ids = np.flatnonzero((types == "Myeloid") & (crops == c))
        if not len(ids):
            continue
        size = lib[ids] / max(lib[ids].mean(), 1e-9)
        lam = s * rr * mean_fib[c][None, :] * size[:, None]
        X[np.ix_(ids, [col[g] for g in fib])] += rng.poisson(lam)
    return X


def _extra_genes(w, crops, extra_bg: int, choice_seed: int, draw_seed: int, rr: float = 1.0) -> np.ndarray:
    """``extra_bg`` uninformative genes: each copies the tissue-level SRTsim marginal of a random real gene (per
    crop) and is drawn i.i.d. per cell. Nested: the first 60 of 240 are the 60 of extra_bg = 60."""
    if not extra_bg:
        return np.zeros((len(crops), 0))
    genes_ok = np.flatnonzero(~w.var["control"].to_numpy() & np.isfinite(w.uns["tis_mu"]).all(0))
    src = np.random.default_rng(20_000 + choice_seed).choice(genes_ok, extra_bg)
    out = np.zeros((len(crops), extra_bg))
    for j, g in enumerate(src):
        rng = np.random.default_rng([20_000 + draw_seed, j])
        for c in np.unique(crops):
            ids = np.flatnonzero(crops == c)
            out[ids, j] = _draw(rng, w.uns["tis_pi0"][c, g], w.uns["tis_theta"][c, g], rr * w.uns["tis_mu"][c, g],
                                len(ids))
    return out


def simulate(seed: int = 0, n_cells: int = 2400, depth=None, weights=None, n_samples: int = 6,
             infiltrate: float = 0.0, extra_bg: int = 0, coexpr: float = 0.0, **kw) -> ad.AnnData:
    if depth is not None:
        raise NotImplementedError("the SRTsim world keeps the section's depth; thin the counts instead")
    if kw:
        raise TypeError(f"simulate() got unsupported arguments {sorted(kw)}")
    if not 1 <= n_samples <= N_CROPS:
        raise ValueError(f"n_samples must be 1..{N_CROPS} (one crop per sample)")
    w = _load(seed)
    idx = _window(w, n_samples, n_cells)
    types = w.obs["truth"].astype(str).to_numpy()[idx]
    if weights is not None:
        idx = idx[_keep_by_weights(types, weights, np.random.default_rng(seed + 5))]
        types = w.obs["truth"].astype(str).to_numpy()[idx]
    crops = w.obs["crop"].to_numpy()[idx]
    genes = list(w.var_names[~w.var["control"].to_numpy()])
    ctrl = list(w.var_names[w.var["control"].to_numpy()])
    X = _dense(w.layers["counts"][idx]).astype(float)
    Xg, Xc = X[:, ~w.var["control"].to_numpy()], X[:, w.var["control"].to_numpy()]
    if coexpr:
        Xg = _coexpr_add(np.random.default_rng(seed + 31), Xg, types, crops, w, genes, coexpr)
    Xx = _extra_genes(w, crops, extra_bg, choice_seed=seed, draw_seed=seed)
    xy = w.obsm["spatial"][idx].copy()
    if infiltrate:
        ri = np.random.default_rng(seed + 7)
        cd8 = np.flatnonzero(types == "T_CD8")
        move = ri.choice(cd8, int(round(infiltrate * len(cd8))), replace=False)
        for i in move:                                          # uniform inside the cell's own crop window
            ids = np.flatnonzero(crops == crops[i])
            lo, hi = xy[ids].min(0), xy[ids].max(0)
            xy[i] = ri.uniform(lo, hi)
    xg = [f"xg{i}" for i in range(extra_bg)]
    var = pd.DataFrame({"control": [False] * (len(genes) + extra_bg) + [True] * len(ctrl),
                        "feature_id": list(w.var.loc[genes, "feature_id"]) + xg + list(w.var.loc[ctrl, "feature_id"])},
                       index=genes + xg + ctrl)
    obs = pd.DataFrame({"truth": types, "label": types, "sample": [f"s{c}" for c in crops],
                        "cell_area": w.obs["cell_area"].to_numpy()[idx]},
                       index=[f"c{i}" for i in range(len(idx))])
    a = ad.AnnData(X=sp.csr_matrix(np.hstack([Xg, Xx, Xc]).astype("float32")), obs=obs, var=var)
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = xy
    a.uns["marker_sets"] = {k: list(v) for k, v in marker_sets().items()}
    a.uns["hierarchy"] = dict(HIERARCHY)
    a.uns["world"] = {"seed": seed, "coexpr": coexpr, "extra_bg": extra_bg, "name": "srtsim"}
    return sim.finish(a, seed)


def reference(seed: int = 0, n_per_type: int = 300, depth=None, coexpr: float = 0.0, extra_bg: int = 0,
              **kw) -> ad.AnnData:
    seeds = _seeds()
    wseed = seed if seed in seeds else seeds[seed % len(seeds)]
    w = _load(wseed)
    rr = REF_RR if depth is None else float(depth) / float(w.uns["section_mean_depth"])
    rng = np.random.default_rng(50_000 + seed)
    gmask = ~w.var["control"].to_numpy()
    genes = list(w.var_names[gmask])
    real = _dense(w.layers["real"])[:, gmask]
    truth = w.obs["truth"].astype(str).to_numpy()
    crop_all = w.obs["crop"].to_numpy()
    blocks, crops_out = [], []
    for t in TYPES:
        pool = np.flatnonzero(truth == t)
        pick = rng.choice(pool, n_per_type, replace=True)
        crops = crop_all[pick]
        Xt = np.zeros((n_per_type, len(genes)))
        ti = TYPES.index(t)
        for c in np.unique(crops):
            rows = np.flatnonzero(crops == c)
            R = real[pick[rows]]
            for j, g in enumerate(np.flatnonzero(gmask)):
                mu = w.uns["dom_mu"][c, ti, g]
                if not np.isfinite(mu):                       # not fitted (< 2 non-zero cells): SRTsim gives 0
                    continue
                v = _draw(rng, w.uns["dom_pi0"][c, ti, g], w.uns["dom_theta"][c, ti, g], rr * mu, len(rows))
                Xt[rows, j] = _rank_assign(rng, v, R[:, j])
        blocks.append(Xt)
        crops_out.append(crops)
    X = np.vstack(blocks)
    types = np.repeat(TYPES, n_per_type)
    crops = np.concatenate(crops_out)
    if coexpr:
        X = _coexpr_add(np.random.default_rng(50_000 + seed + 31), X, types, crops, w, genes, coexpr, rr)
    Xx = _extra_genes(w, crops, extra_bg, choice_seed=wseed, draw_seed=50_000 + seed, rr=rr)
    names = genes + [f"xg{i}" for i in range(extra_bg)]
    a = ad.AnnData(X=sp.csr_matrix(np.hstack([X, Xx]).astype("float32")),
                   obs=pd.DataFrame({"cell_type": types}, index=[f"r{i}" for i in range(len(types))]),
                   var=pd.DataFrame(index=names))
    a.layers["counts"] = a.X.copy()
    return a


def install(sim_module=None):
    """Point ``sim`` (the toy module the sweep reads) at this world: simulate, reference, marker sets."""
    m = sim_module or sim
    m.simulate = simulate
    m.reference = reference
    m.marker_sets = lambda per_type=4: {k: list(v) for k, v in marker_sets().items()}
    m.T_LINEAGE = list(_t_lineage())
    return m
