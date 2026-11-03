"""Synthetic spatial sections with known truth, and perturbations of known strength.

One "world" = six cell types laid out as tissue stripes (T_CD4 | T_CD8 | B | Myeloid | Fibro | Epi), each
with private markers; the two T subtypes also share two lineage markers (like CD3E), so the hierarchy
T -> {T_CD4, T_CD8} is real in the data. Background genes carry a mild type programme; negative-control
probes are flat ambient. Samples are horizontal bands, so every type is present in every sample (ISC).

Every scenario is ``fn(strength, seed) -> AnnData`` with ``obs['truth']`` (true type) and ``obs['label']``
(the annotation being judged), ``layers['counts']``, ``obsm['spatial']``, ``obsm['X_pca']``,
``var['control']`` and ``uns['marker_sets']`` / ``uns['hierarchy']``. Strength 0 is always the clean world.
"""
from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

TYPES = ["T_CD4", "T_CD8", "B", "Myeloid", "Fibro", "Epi"]
WEIGHTS = np.array([0.15, 0.15, 0.12, 0.16, 0.20, 0.22])
HIERARCHY = {"T_CD4": "T", "T_CD8": "T", "B": "B", "Myeloid": "Myeloid", "Fibro": "Fibro", "Epi": "Epi"}
T_LINEAGE = ["T_lin0", "T_lin1"]


def marker_sets(per_type: int = 4) -> dict[str, list[str]]:
    ms = {t: [f"{t}_m{i}" for i in range(per_type)] for t in TYPES}
    for t in ("T_CD4", "T_CD8"):
        ms[t] = T_LINEAGE + ms[t][: per_type - len(T_LINEAGE)]
    return ms


def world(seed: int, extra_bg: int = 0, n_bg: int = 60, n_ctrl: int = 20, marker_rate: float = 10.0, off_rate: float = 0.01,
          ctrl_rate: float = 0.02, coexpr: float = 0.0):
    """Type x gene rate matrix. ``coexpr`` = Myeloid truly expresses Fibro's markers at that fraction of
    Fibro's level (biology, present in the reference too - not spillover). ``extra_bg`` uninformative genes
    ``xg*`` come from their own random stream, so adding them leaves every other rate unchanged."""
    rng = np.random.default_rng(10_000 + seed)
    ms = marker_sets()
    markers = list(dict.fromkeys(g for t in TYPES for g in ms[t]))
    bg = [f"bg{i}" for i in range(n_bg)]
    ctrl = [f"NegCtrl{i}" for i in range(n_ctrl)]
    xg = [f"xg{i}" for i in range(extra_bg)]
    genes = markers + bg + xg + ctrl
    R = pd.DataFrame(0.0, index=TYPES, columns=genes)
    R[markers] = off_rate
    for t in TYPES:
        R.loc[t, ms[t]] = marker_rate
    base = rng.lognormal(0, 1, n_bg)
    R[bg] = base * rng.lognormal(0, 0.3, (len(TYPES), n_bg))   # mild per-type background programme
    if extra_bg:
        r2 = np.random.default_rng(20_000 + seed)
        R[xg] = r2.lognormal(0, 1, extra_bg) * r2.lognormal(0, 0.3, (len(TYPES), extra_bg))
    R[ctrl] = ctrl_rate
    if coexpr:
        R.loc["Myeloid", ms["Fibro"]] = coexpr * marker_rate
    return R, ms


def _counts(R: pd.DataFrame, types: np.ndarray, depth: float, rng) -> np.ndarray:
    """Poisson counts; ``depth`` = expected counts on the base panel (``xg*`` genes add counts on top,
    as a bigger panel would, instead of diluting the base genes)."""
    base = R[[c for c in R.columns if not c.startswith("xg")]].sum(1)
    P = R.loc[types].to_numpy() / base.loc[types].to_numpy()[:, None]
    size = depth * rng.lognormal(0, 0.3, len(types))
    lam = P * size[:, None]
    xg = np.array([c.startswith("xg") for c in R.columns])
    X = np.zeros(lam.shape, dtype="float32")
    X[:, ~xg] = rng.poisson(lam[:, ~xg])               # base genes first: same draws whatever the panel size
    X[:, xg] = rng.poisson(lam[:, xg])
    return X


def _pca(X: np.ndarray, seed: int = 0) -> np.ndarray:
    from sklearn.decomposition import PCA
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    L = np.log1p(X / lib * 1e4)
    return PCA(n_components=min(20, L.shape[1] - 1, L.shape[0] - 1), random_state=seed).fit_transform(L)


def finish(a: ad.AnnData, seed: int = 0) -> ad.AnnData:
    """(Re)derive everything computed from counts / coords: X, X_pca, cleared spatial graph."""
    X = np.asarray(a.layers["counts"].todense()) if sp.issparse(a.layers["counts"]) else a.layers["counts"]
    a.layers["counts"] = sp.csr_matrix(X)
    a.X = a.layers["counts"].copy()
    a.obsm["X_pca"] = _pca(X, seed)
    for k in list(a.obsp.keys()):
        del a.obsp[k]
    return a


def simulate(seed: int = 0, n_cells: int = 2400, depth: float = 150.0, weights=None, n_samples: int = 6,
             infiltrate: float = 0.0, extra_bg: int = 0, coexpr: float = 0.0) -> ad.AnnData:
    """A clean section. ``infiltrate`` = fraction of T_CD8 cells scattered at random positions (correct
    labels, legitimately incoherent); ``extra_bg`` = uninformative genes added to the panel."""
    R, ms = world(seed, extra_bg=extra_bg, coexpr=coexpr)
    rng = np.random.default_rng(seed)                 # layout; infiltration and counts get their own streams
    rng_c = np.random.default_rng(seed + 13)
    w = WEIGHTS if weights is None else np.asarray(weights, float)
    side = np.sqrt(n_cells)
    xy = rng.uniform(0, side, (n_cells, 2))
    edges = np.cumsum(w / w.sum())[:-1] * side                  # stripes along x, widths = proportions
    types = np.array(TYPES)[np.searchsorted(edges, xy[:, 0])]
    if infiltrate:
        cd8 = np.flatnonzero(types == "T_CD8")
        ri = np.random.default_rng(seed + 7)
        move = ri.choice(cd8, int(round(infiltrate * len(cd8))), replace=False)
        xy[move] = ri.uniform(0, side, (len(move), 2))
    X = _counts(R, types, depth, rng_c)
    sample = (xy[:, 1] / side * n_samples).astype(int).clip(0, n_samples - 1)
    obs = pd.DataFrame({"truth": types, "label": types, "sample": [f"s{i}" for i in sample],
                        "cell_area": np.random.default_rng(seed + 17).lognormal(np.log(50), 0.2, n_cells)},
                       index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame({"control": [g.startswith("NegCtrl") for g in R.columns]}, index=R.columns)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    a.layers["counts"] = a.X.copy()
    a.obsm["spatial"] = xy
    a.uns["marker_sets"] = ms
    a.uns["hierarchy"] = HIERARCHY
    a.uns["world"] = {"seed": seed, "coexpr": coexpr, "extra_bg": extra_bg}   # to build the matched reference
    return finish(a, seed)


def reference(seed: int = 0, n_per_type: int = 300, depth: float = 400.0, coexpr: float = 0.0,
              extra_bg: int = 0) -> ad.AnnData:
    """Matched dissociated reference from the same world: deeper, no spatial spillover."""
    R, _ = world(seed, extra_bg=extra_bg, coexpr=coexpr)
    rng = np.random.default_rng(50_000 + seed)
    types = np.repeat(TYPES, n_per_type)
    X = _counts(R, types, depth, rng)
    keep = ~np.array([g.startswith("NegCtrl") for g in R.columns])
    a = ad.AnnData(X=sp.csr_matrix(X[:, keep]), obs=pd.DataFrame({"cell_type": types},
                   index=[f"r{i}" for i in range(len(types))]),
                   var=pd.DataFrame(index=R.columns[keep]))
    a.layers["counts"] = a.X.copy()
    return a


# --------------------------------------------------------------------------- #
# perturbations: label-side (data untouched, labels wrong)
# --------------------------------------------------------------------------- #
def _relabel(a, idx, new):
    lab = a.obs["label"].astype(str).to_numpy().copy()
    lab[idx] = new
    a.obs["label"] = lab
    return a


def label_noise(s, seed):
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    idx = rng.choice(a.n_obs, int(round(s * a.n_obs)), replace=False)
    truth = a.obs["truth"].to_numpy()
    new = [rng.choice([t for t in TYPES if t != truth[i]]) for i in idx]
    return _relabel(a, idx, new)


def sibling_confusion(s, seed):
    """A fraction of T_CD8 called T_CD4: a near miss, right lineage."""
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    cd8 = np.flatnonzero(a.obs["truth"].to_numpy() == "T_CD8")
    return _relabel(a, rng.choice(cd8, int(round(s * len(cd8))), replace=False), "T_CD4")


def type_swap(s, seed):
    """A fraction of Fibro called Epi AND the same fraction of Epi called Fibro (s = 1: names swapped)."""
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    t = a.obs["truth"].to_numpy()
    for src, dst in (("Fibro", "Epi"), ("Epi", "Fibro")):
        ids = np.flatnonzero(t == src)
        a = _relabel(a, rng.choice(ids, int(round(s * len(ids))), replace=False), dst)
    return a


def patch_error(s, seed):
    """A spatially compact patch of Fibro (fraction s, grown from the Epi border) called Epi: an error that
    is coherent in space, like a misplaced domain boundary."""
    a = simulate(seed)
    t = a.obs["truth"].to_numpy()
    fib = np.flatnonzero(t == "Fibro")
    x = a.obsm["spatial"][fib, 0]
    n = int(round(s * len(fib)))
    return _relabel(a, fib[np.argsort(-x)[:n]], "Epi")       # rightmost Fibro = next to Epi


def over_split(s, seed):
    """``s`` types (0-3) split into two random halves: over-partitioning without structure."""
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    lab = a.obs["label"].astype(str).to_numpy().copy()
    for t in ["Epi", "Fibro", "Myeloid"][: int(s)]:
        ids = np.flatnonzero(lab == t)
        half = rng.choice(ids, len(ids) // 2, replace=False)
        lab[ids] = f"{t}_a"
        lab[half] = f"{t}_b"
    a.obs["label"] = lab
    return a


def merge(s, seed):
    """``s`` pairs (0-3) merged into one label: under-partitioning (the last pair is the T siblings)."""
    a = simulate(seed)
    lab = a.obs["label"].astype(str).to_numpy().copy()
    for x, y in [("Fibro", "Epi"), ("B", "Myeloid"), ("T_CD4", "T_CD8")][: int(s)]:
        lab[np.isin(lab, [x, y])] = f"{x}+{y}"
    a.obs["label"] = lab
    return a


# --------------------------------------------------------------------------- #
# perturbations: data-side (labels = truth, data degraded)
# --------------------------------------------------------------------------- #
def _knn_mean(a, k=6):
    from sklearn.neighbors import NearestNeighbors
    xy = a.obsm["spatial"]
    idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    X = np.asarray(a.layers["counts"].todense())
    return X, X[idx].mean(1)


def spillover(s, seed):
    """Segmentation bleed: each cell keeps a (1-s) binomial share of its own counts and gains Poisson(s x the
    mean of its 6 spatial neighbours). s = 0 is the clean section, untouched."""
    a = simulate(seed)
    X, nb = _knn_mean(a)
    rng = np.random.default_rng(seed + 1)
    a.layers["counts"] = (rng.binomial(X.astype(int), 1 - s) + rng.poisson(s * nb)).astype("float32")
    return finish(a, seed)


def depth_loss(s, seed):
    """Binomial thinning: keep a (1 - s) fraction of every transcript."""
    a = simulate(seed)
    X = np.asarray(a.layers["counts"].todense()).astype(int)
    a.layers["counts"] = np.random.default_rng(seed + 1).binomial(X, 1 - s).astype("float32")
    return finish(a, seed)


def ambient(s, seed):
    """Ambient / off-target background: Poisson(s) extra counts on every gene, controls included."""
    a = simulate(seed)
    X = np.asarray(a.layers["counts"].todense())
    a.layers["counts"] = (X + np.random.default_rng(seed + 1).poisson(s, X.shape)).astype("float32")
    return finish(a, seed)


def batch(s, seed):
    """Per-sample gene-wise multiplicative effect, log-sd = s (expression shifts between replicates)."""
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    X = np.asarray(a.layers["counts"].todense())
    samples = a.obs["sample"].to_numpy()
    f = {smp: rng.lognormal(0, s, X.shape[1]) if s else np.ones(X.shape[1]) for smp in np.unique(samples)}
    F = np.vstack([f[smp] for smp in samples])
    a.layers["counts"] = (rng.poisson(X * F) if s else X).astype("float32")
    return finish(a, seed)


def coord_shuffle(s, seed):
    """Coordinates permuted among a fraction s of cells (registration / mix-up of positions)."""
    a = simulate(seed)
    rng = np.random.default_rng(seed + 1)
    idx = rng.choice(a.n_obs, int(round(s * a.n_obs)), replace=False)
    xy = a.obsm["spatial"].copy()
    xy[idx] = xy[rng.permutation(idx)]
    a.obsm["spatial"] = xy
    return finish(a, seed)


# --------------------------------------------------------------------------- #
# negative controls: labels correct, nothing wrong - a good metric stays flat
# --------------------------------------------------------------------------- #
def coexpression(s, seed):
    """Myeloid genuinely co-expresses Fibro markers at fraction s of Fibro's level (in the reference too)."""
    return simulate(seed, coexpr=s)


def rare_type(s, seed):
    """T_CD8 shrinks from 15 % of cells to ``s`` (correctly labelled)."""
    w = WEIGHTS.copy()
    w[1] = s
    w[[0, 2, 3, 4, 5]] *= (1 - s) / w[[0, 2, 3, 4, 5]].sum()
    return simulate(seed, weights=w)


def infiltration(s, seed):
    """A fraction s of T_CD8 cells scattered through the tissue (correct labels)."""
    return simulate(seed, infiltrate=s)


def panel_growth(s, seed):
    """``s`` uninformative genes added to the panel (60 -> 60 + s background genes)."""
    return simulate(seed, extra_bg=int(s))


SCENARIOS = {
    # name: (fn, levels, kind)  kind: label = labels wrong; data = data degraded; null = nothing wrong
    "label_noise": (label_noise, [0, 0.05, 0.1, 0.2, 0.4], "label"),
    "sibling_confusion": (sibling_confusion, [0, 0.1, 0.25, 0.5, 0.75], "label"),
    "type_swap": (type_swap, [0, 0.25, 0.5, 0.75, 1.0], "label"),
    "patch_error": (patch_error, [0, 0.1, 0.2, 0.35, 0.5], "label"),
    "over_split": (over_split, [0, 1, 2, 3], "label"),
    "merge": (merge, [0, 1, 2, 3], "label"),
    "spillover": (spillover, [0, 0.05, 0.1, 0.2, 0.4], "data"),
    "depth_loss": (depth_loss, [0, 0.5, 0.75, 0.88, 0.94], "data"),
    "ambient": (ambient, [0, 0.02, 0.05, 0.1, 0.2], "data"),
    "batch": (batch, [0, 0.25, 0.5, 1.0, 1.5], "data"),
    "coord_shuffle": (coord_shuffle, [0, 0.1, 0.25, 0.5, 1.0], "data"),
    "coexpression": (coexpression, [0, 0.05, 0.1, 0.2, 0.4], "null"),
    "rare_type": (rare_type, [0.15, 0.08, 0.04, 0.02, 0.01], "null"),
    "infiltration": (infiltration, [0, 0.1, 0.25, 0.5, 1.0], "null"),
    "panel_growth": (panel_growth, [0, 60, 240, 960], "null"),
}
