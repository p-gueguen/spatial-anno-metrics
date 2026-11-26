"""Generic AnnData data-augmentation and behavioral metric-contract verification.

Provides perturbation operators (label-side, count/signal-side, spatial, and negative
controls) that operate directly on real or synthetic :class:`~anndata.AnnData` objects,
plus :func:`verify_metric_contracts` to automatically verify that annotation metrics move
in the right direction on targeted failures and remain invariant on negative controls.
"""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np
import pandas as pd

__all__ = [
    "AUGMENTATIONS",
    "METRIC_CONTRACTS",
    "perturb_anndata",
    "verify_metric_contracts",
]

AUGMENTATIONS: dict[str, dict[str, Any]] = {
    "label_noise": {
        "kind": "label",
        "requires": ("labels",),
        "default_levels": [0.0, 0.05, 0.1, 0.25, 0.5],
        "description": "Randomly flip fraction `level` of labels to another observed cell type.",
    },
    "merge": {
        "kind": "label",
        "requires": ("labels",),
        "default_levels": [0, 1, 2],
        "description": "Coarsen `level` pairs of cell types into combined `<A>+<B>` labels.",
    },
    "over_split": {
        "kind": "label",
        "requires": ("labels",),
        "default_levels": [0, 1, 2],
        "description": "Randomly split `level` abundant cell types into two synthetic sublabels.",
    },
    "type_swap": {
        "kind": "label",
        "requires": ("labels",),
        "default_levels": [0, 1, 2],
        "description": "Swap labels between disjoint pairs of distinct cell types.",
    },
    "sibling_confusion": {
        "kind": "label",
        "requires": ("labels", "counts"),
        "default_levels": [0.0, 0.2, 0.5, 0.8],
        "description": "Mislabel fraction `level` of one cell type as its closest expression neighbor.",
    },
    "patch_error": {
        "kind": "label",
        "requires": ("labels", "spatial"),
        "default_levels": [0.0, 0.2, 0.4, 0.6],
        "description": "Relabel a spatially compact border patch of one type to an adjacent type.",
    },
    "spillover": {
        "kind": "data",
        "requires": ("counts", "spatial"),
        "default_levels": [0.0, 0.1, 0.2, 0.4],
        "description": "Mix each cell's expression with its `k` spatial neighbors' mean counts.",
    },
    "depth_loss": {
        "kind": "data",
        "requires": ("counts",),
        "default_levels": [0.0, 0.25, 0.5, 0.75],
        "description": "Downsample counts (binomial thinning on integers, entry dropout on floats).",
    },
    "ambient": {
        "kind": "data",
        "requires": ("counts",),
        "default_levels": [0.0, 0.05, 0.15, 0.3],
        "description": "Add uniform non-specific ambient background counts across all genes.",
    },
    "batch_shift": {
        "kind": "data",
        "requires": ("counts", "sample"),
        "default_levels": [0.0, 0.25, 0.5, 1.0],
        "description": "Apply sample-specific log-normal gene scaling across replicates.",
    },
    "coord_shuffle": {
        "kind": "spatial",
        "requires": ("spatial",),
        "default_levels": [0.0, 0.25, 0.5, 1.0],
        "description": "Permute spatial coordinates across fraction `level` of cells (within sample).",
    },
    "infiltration": {
        "kind": "spatial",
        "requires": ("labels", "spatial"),
        "default_levels": [0.0, 0.25, 0.5, 0.8],
        "description": "Scatter fraction `level` of one type through tissue preserving expr-label pairs.",
    },
    "nuisance_shuffle": {
        "kind": "null",
        "requires": (),
        "default_levels": [0, 1, 2],
        "description": "Negative control: permute cell row order without changing content.",
    },
    "panel_growth": {
        "kind": "null",
        "requires": ("counts",),
        "default_levels": [0, 60, 240],
        "description": "Append `level` uninformative Poisson background genes to the panel.",
    },
}

METRIC_CONTRACTS: dict[str, dict[str, Any]] = {
    "ext.macro_f1": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "merge": [0, 1, 2], "type_swap": [0, 1, 2]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {},
    },
    "iv.integrated": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "over_split": [0, 1, 2]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {"infiltration": "Embedding-only; invariant to spatial coordinates."},
    },
    "fid.mean_auc": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {"over_split": "AUROC can remain flat when random halves keep markers."},
    },
    "pmp": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {"spillover": "Drops only on marker-gene panels with high spillover."},
    },
    "nmp": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "type_swap": [0, 1, 2]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "Skips merged types missing from the reference."},
    },
    "rel_pw_ct": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "type_swap": [0, 1, 2]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {
            "panel_growth": "Unweighted per-gene L1 normalization drifts when adding pure-noise genes."
        },
    },
    "resid_moran.explained": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "merge": [0, 1, 2]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {
            "infiltration": "Scattering a cluster lowers raw Moran so 1 - resid/raw can wobble slightly.",
            "coord_shuffle": "Ratio 1 - resid/raw is unstable near raw <= 0; gate on raw_median_moran > 0.",
        },
    },
    "resid_moran.residual": {
        "orientation": -1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "merge": [0, 1, 2], "spillover": [0.0, 0.1, 0.2, 0.4]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {"infiltration": [0.0, 0.25, 0.5]},
        "known_caveats": {"coord_shuffle": "Raw spatial signal vanishes on shuffle so residual Moran drops near 0."},
    },
    "spatial_org.corr": {
        "orientation": 1,
        "must_detect": {"coord_shuffle": [0.0, 0.25, 0.5, 1.0]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "Diagonal self-enrichment can stay correlated across replicates under merges."},
    },
    "spatial_org.offdiag": {
        "orientation": 1,
        "must_detect": {"coord_shuffle": [0.0, 0.25, 0.5, 1.0], "label_noise": [0.0, 0.1, 0.25, 0.5]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "May shift either way when merged pairs dominate background neighbors."},
    },
    "sc.mean_coherence": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "coord_shuffle": [0.0, 0.25, 0.5, 1.0]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {
            "merge": "Inverts (rises) when neighboring cell types are merged.",
            "infiltration": "False-alarms on correctly labeled infiltrating cells.",
        },
    },
    "sc.pas": {
        "orientation": -1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "coord_shuffle": [0.0, 0.25, 0.5, 1.0]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "Inverts on adjacent-type merge.", "infiltration": "False-alarms on infiltration."},
    },
    "sc.chaos": {
        "orientation": -1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "coord_shuffle": [0.0, 0.25, 0.5, 1.0]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "Drops/inverts on merge.", "infiltration": "False-alarms on infiltration."},
    },
    "isc.integrated": {
        "orientation": 1,
        "must_detect": {"label_noise": [0.0, 0.1, 0.25, 0.5], "batch_shift": [0.0, 0.5, 1.0, 1.5]},
        "must_be_invariant": {"nuisance_shuffle": [0, 1, 2], "coord_shuffle": [0.0, 0.5, 1.0]},
        "must_not_degrade": {},
        "known_caveats": {"merge": "Can rise when noisy sibling labels are merged."},
    },
}


def _get_matrix(adata: Any, layer: str | None) -> tuple[np.ndarray, bool]:
    import scipy.sparse as sp

    if layer and layer in adata.layers:
        raw = adata.layers[layer]
    else:
        raw = adata.X
    was_sparse = sp.issparse(raw)
    arr = np.asarray(raw.todense(), dtype=np.float64) if was_sparse else np.asarray(raw, dtype=np.float64).copy()
    return arr, was_sparse


def _set_matrix(adata: Any, arr: np.ndarray, layer: str | None, was_sparse: bool) -> None:
    import scipy.sparse as sp

    mat = sp.csr_matrix(arr) if was_sparse else arr
    if layer and layer in adata.layers:
        adata.layers[layer] = mat.copy()
        adata.X = mat.copy()
    else:
        adata.X = mat.copy()
        if layer and layer in adata.layers:
            adata.layers[layer] = mat.copy()


def _recompute_pca_from(arr: np.ndarray, n_components: int = 20) -> np.ndarray:
    from sklearn.decomposition import PCA

    tot = arr.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1.0
    log_norm = np.log1p(arr / tot * 1e4)
    k = min(n_components, arr.shape[0] - 1, arr.shape[1])
    if k < 1:
        return np.zeros((arr.shape[0], 1), dtype=np.float32)
    return PCA(n_components=k, random_state=0).fit_transform(log_norm).astype(np.float32)


def _ordered_types(labels: np.ndarray) -> list[str]:
    vc = pd.Series(labels).value_counts()
    return sorted(vc.index.tolist(), key=lambda t: (-int(vc[t]), t))


def _default_pairs(labels: np.ndarray) -> list[tuple[str, str]]:
    ordered = _ordered_types(labels)
    return [(ordered[i], ordered[i + 1]) for i in range(0, len(ordered) - 1, 2)]


def _closest_type_pair(arr: np.ndarray, labels: np.ndarray, min_cells: int = 5) -> tuple[str, str] | None:
    vc = pd.Series(labels).value_counts()
    eligible = [t for t in sorted(vc.index) if int(vc[t]) >= min_cells]
    if len(eligible) < 2:
        eligible = sorted(vc.index.tolist())
    if len(eligible) < 2:
        return None
    means = []
    for t in eligible:
        m = arr[labels == t].mean(axis=0)
        s = float(m.sum())
        means.append(m / s if s > 0 else m)
    M = np.vstack(means)
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    Mn = M / norms
    sim = Mn @ Mn.T
    np.fill_diagonal(sim, -np.inf)
    i, j = np.unravel_index(int(np.argmax(sim)), sim.shape)
    return eligible[i], eligible[j]


def _closest_spatial_pair(xy: np.ndarray, labels: np.ndarray, min_cells: int = 10) -> tuple[str, str] | None:
    from scipy.spatial import cKDTree

    vc = pd.Series(labels).value_counts()
    eligible = [t for t in sorted(vc.index) if int(vc[t]) >= min_cells]
    if len(eligible) < 2:
        eligible = sorted(vc.index.tolist())
    if len(eligible) < 2:
        return None
    best_dist = np.inf
    best_pair: tuple[str, str] | None = None
    for i, t1 in enumerate(eligible):
        xy1 = xy[labels == t1]
        for t2 in eligible[i + 1 :]:
            xy2 = xy[labels == t2]
            d, _ = cKDTree(xy2).query(xy1, k=1)
            med = float(np.median(d))
            if med < best_dist:
                best_dist = med
                best_pair = (t1, t2)
    return best_pair


def _is_integer_like(arr: np.ndarray) -> bool:
    if arr.size == 0:
        return True
    if arr.min() < 0:
        return False
    return bool(np.allclose(arr, np.rint(arr), atol=1e-6))


def _clear_spatial_graphs(adata: Any) -> None:
    for k in list(adata.obsp):
        del adata.obsp[k]
    if "spatial_neighbors" in adata.uns:
        del adata.uns["spatial_neighbors"]


def perturb_anndata(
    adata: Any,
    perturbation: str,
    level: float,
    *,
    seed: int = 0,
    label_key: str = "cell_type",
    sample_key: str | None = None,
    layer: str | None = "counts",
    recompute_pca: bool = False,
    **kwargs: Any,
) -> Any:
    """Apply a controlled data, label, spatial, or null perturbation to an ``AnnData`` copy.

    Parameters
    ----------
    adata
        Input object. Never modified in place.
    perturbation
        Name in :data:`AUGMENTATIONS`.
    level
        Perturbation strength (e.g. flipped fraction, pair count, thinning fraction, or extra gene count).
    seed
        Random seed controlling stochastic perturbations.
    label_key
        Column in ``adata.obs`` holding cell-type labels.
    sample_key
        Optional sample/section column in ``adata.obs`` (used by ``batch_shift`` and within-sample shuffles).
    layer
        Count/expression layer to read and write when present; falls back to ``adata.X``.
    recompute_pca
        If ``True`` and expression counts changed, recompute ``adata.obsm["X_pca"]``.
    **kwargs
        Optional operator parameters (e.g. ``pairs``, ``k``, ``target_type``, ``src``, ``dst``).

    Returns
    -------
    anndata.AnnData
        Perturbed copy with ``uns["sam_augmentation"]`` recording the operation metadata.
    """
    if perturbation not in AUGMENTATIONS:
        raise ValueError(f"unknown perturbation {perturbation!r}; expected one of {sorted(AUGMENTATIONS)}")
    a = adata.copy()
    rng = np.random.default_rng(seed)
    lvl = float(level)
    data_changed = False

    if label_key in a.obs:
        labels = a.obs[label_key].astype(str).to_numpy().copy()
    else:
        labels = np.array([], dtype=str)

    if perturbation == "label_noise":
        if label_key not in a.obs:
            raise ValueError(f"label_key {label_key!r} not in adata.obs")
        types = sorted(set(labels))
        n_flip = int(round(lvl * len(labels)))
        if n_flip > 0 and len(types) >= 2:
            idx = rng.choice(len(labels), min(n_flip, len(labels)), replace=False)
            for i in idx:
                labels[i] = rng.choice([t for t in types if t != labels[i]])
            a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "merge":
        if label_key not in a.obs:
            raise ValueError(f"label_key {label_key!r} not in adata.obs")
        pairs = kwargs.get("pairs") or _default_pairs(labels)
        for t1, t2 in pairs[: int(round(lvl))]:
            mask = (labels == t1) | (labels == t2)
            if mask.any():
                labels[mask] = f"{t1}+{t2}"
        a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "over_split":
        if label_key not in a.obs:
            raise ValueError(f"label_key {label_key!r} not in adata.obs")
        targets = kwargs.get("types") or _ordered_types(labels)
        for t in targets[: int(round(lvl))]:
            idx = np.flatnonzero(labels == t)
            if len(idx) >= 2:
                half = rng.choice(idx, len(idx) // 2, replace=False)
                labels[idx] = f"{t}_splitA"
                labels[half] = f"{t}_splitB"
        a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "type_swap":
        if label_key not in a.obs:
            raise ValueError(f"label_key {label_key!r} not in adata.obs")
        pairs = kwargs.get("pairs") or _default_pairs(labels)
        if lvl >= 1.0:
            for t1, t2 in pairs[: int(round(lvl))]:
                m1, m2 = labels == t1, labels == t2
                labels[m1], labels[m2] = t2, t1
        elif lvl > 0 and pairs:
            t1, t2 = pairs[0]
            i1, i2 = np.flatnonzero(labels == t1), np.flatnonzero(labels == t2)
            n_pick1 = int(round(lvl * len(i1)))
            n_pick2 = int(round(lvl * len(i2)))
            p1 = rng.choice(i1, n_pick1, replace=False) if n_pick1 > 0 else np.array([], dtype=int)
            p2 = rng.choice(i2, n_pick2, replace=False) if n_pick2 > 0 else np.array([], dtype=int)
            labels[p1] = t2
            labels[p2] = t1
        a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "sibling_confusion":
        if label_key not in a.obs:
            raise ValueError(f"label_key {label_key!r} not in adata.obs")
        if lvl > 0:
            pair = kwargs.get("pair")
            if pair is None:
                arr, _ = _get_matrix(a, layer)
                pair = _closest_type_pair(arr, labels)
            if pair is not None:
                src, dst = pair
                idx = np.flatnonzero(labels == src)
                n_flip = int(round(lvl * len(idx)))
                if n_flip > 0:
                    labels[rng.choice(idx, min(n_flip, len(idx)), replace=False)] = dst
                    a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "patch_error":
        if label_key not in a.obs or "spatial" not in a.obsm:
            raise ValueError("patch_error requires obs[label_key] and obsm['spatial']")
        if lvl > 0:
            xy = np.asarray(a.obsm["spatial"], dtype=np.float64)
            src, dst = kwargs.get("src"), kwargs.get("dst")
            if src is None or dst is None:
                pair = _closest_spatial_pair(xy, labels)
                if pair is not None:
                    src, dst = pair
            if src is not None and dst is not None:
                src_idx = np.flatnonzero(labels == src)
                dst_idx = np.flatnonzero(labels == dst)
                if len(src_idx) > 0 and len(dst_idx) > 0:
                    centroid = xy[dst_idx].mean(axis=0)
                    dist = np.linalg.norm(xy[src_idx] - centroid, axis=1)
                    n_flip = int(round(lvl * len(src_idx)))
                    if n_flip > 0:
                        labels[src_idx[np.argsort(dist)[: min(n_flip, len(src_idx))]]] = dst
                        a.obs[label_key] = pd.Categorical(labels)

    elif perturbation == "spillover":
        if "spatial" not in a.obsm:
            raise ValueError("spillover requires adata.obsm['spatial']")
        if lvl > 0:
            from scipy.spatial import cKDTree

            xy = np.asarray(a.obsm["spatial"], dtype=np.float64)
            k = min(int(kwargs.get("k", 6)), max(1, a.n_obs - 1))
            _, nn = cKDTree(xy).query(xy, k=k + 1)
            nbrs = nn[:, 1:]
            arr, was_sparse = _get_matrix(a, layer)
            nb_mean = arr[nbrs].mean(axis=1)
            if _is_integer_like(arr):
                X_int = np.rint(arr).astype(np.int64)
                mixed = rng.binomial(X_int, 1.0 - lvl) + rng.poisson(np.maximum(0.0, lvl * nb_mean))
                out = mixed.astype(arr.dtype)
            else:
                out = np.maximum(0.0, (1.0 - lvl) * arr + lvl * nb_mean)
            _set_matrix(a, out, layer, was_sparse)
            data_changed = True

    elif perturbation == "depth_loss":
        if lvl > 0:
            arr, was_sparse = _get_matrix(a, layer)
            if _is_integer_like(arr):
                out = rng.binomial(np.rint(arr).astype(np.int64), 1.0 - lvl).astype(arr.dtype)
            else:
                keep = rng.random(arr.shape) >= lvl
                out = arr * keep
            _set_matrix(a, out, layer, was_sparse)
            data_changed = True

    elif perturbation == "ambient":
        if lvl > 0:
            arr, was_sparse = _get_matrix(a, layer)
            if _is_integer_like(arr):
                out = (arr + rng.poisson(lvl, size=arr.shape)).astype(arr.dtype)
            else:
                pos = arr[arr > 0]
                scale = float(np.mean(pos)) if pos.size > 0 else 1.0
                out = arr + rng.exponential(lvl * scale, size=arr.shape)
            _set_matrix(a, out, layer, was_sparse)
            data_changed = True

    elif perturbation == "batch_shift":
        if not sample_key or sample_key not in a.obs:
            raise ValueError(f"batch_shift requires valid sample_key in adata.obs (got {sample_key!r})")
        if lvl > 0:
            arr, was_sparse = _get_matrix(a, layer)
            samples = a.obs[sample_key].astype(str).to_numpy()
            factors = {s: rng.lognormal(0.0, lvl, size=arr.shape[1]) for s in pd.unique(samples)}
            F = np.vstack([factors[s] for s in samples])
            if _is_integer_like(arr):
                out = rng.poisson(np.maximum(0.0, arr * F)).astype(arr.dtype)
            else:
                out = np.maximum(0.0, arr * F)
            _set_matrix(a, out, layer, was_sparse)
            data_changed = True

    elif perturbation == "coord_shuffle":
        if "spatial" not in a.obsm:
            raise ValueError("coord_shuffle requires adata.obsm['spatial']")
        if lvl > 0:
            xy = np.asarray(a.obsm["spatial"], dtype=np.float64).copy()
            groups = (
                [np.flatnonzero(a.obs[sample_key].astype(str).to_numpy() == s)
                 for s in pd.unique(a.obs[sample_key].astype(str))]
                if sample_key and sample_key in a.obs
                else [np.arange(a.n_obs)]
            )
            for idx in groups:
                n_shuf = int(round(lvl * len(idx)))
                if n_shuf >= 2:
                    pick = rng.choice(idx, min(n_shuf, len(idx)), replace=False)
                    xy[pick] = xy[rng.permutation(pick)]
            a.obsm["spatial"] = xy
            _clear_spatial_graphs(a)

    elif perturbation == "infiltration":
        if label_key not in a.obs or "spatial" not in a.obsm:
            raise ValueError("infiltration requires obs[label_key] and obsm['spatial']")
        if lvl > 0:
            xy = np.asarray(a.obsm["spatial"], dtype=np.float64).copy()
            target = kwargs.get("target_type")
            if target is None:
                ordered = _ordered_types(labels)
                target = ordered[min(1, len(ordered) - 1)] if ordered else None
            if target is not None:
                groups = (
                    [np.flatnonzero(a.obs[sample_key].astype(str).to_numpy() == s)
                     for s in pd.unique(a.obs[sample_key].astype(str))]
                    if sample_key and sample_key in a.obs
                    else [np.arange(a.n_obs)]
                )
                for g_idx in groups:
                    t_idx = g_idx[labels[g_idx] == target]
                    o_idx = g_idx[labels[g_idx] != target]
                    n_move = min(int(round(lvl * len(t_idx))), len(o_idx))
                    if n_move > 0:
                        pick_t = rng.choice(t_idx, n_move, replace=False)
                        pick_o = rng.choice(o_idx, n_move, replace=False)
                        xy_t = xy[pick_t].copy()
                        xy[pick_t] = xy[pick_o]
                        xy[pick_o] = xy_t
                a.obsm["spatial"] = xy
                _clear_spatial_graphs(a)

    elif perturbation == "nuisance_shuffle":
        if lvl > 0:
            order = rng.permutation(a.n_obs)
            a = a[order].copy()
            _clear_spatial_graphs(a)

    elif perturbation == "panel_growth":
        n_extra = int(round(lvl))
        if n_extra > 0:
            import anndata as ad
            import scipy.sparse as sp

            arr, was_sparse = _get_matrix(a, layer)
            lam = float(kwargs.get("lam", 0.25))
            if _is_integer_like(arr):
                bg = rng.poisson(lam, size=(a.n_obs, n_extra)).astype(arr.dtype)
            else:
                bg = rng.exponential(max(lam, 1e-2), size=(a.n_obs, n_extra)).astype(arr.dtype)
            new_var = pd.DataFrame(index=[f"NegBg_{i}" for i in range(n_extra)])
            for col in a.var.columns:
                if col == "control":
                    new_var[col] = False
                elif col == "feature_type":
                    new_var[col] = "Gene Expression"
            X_full = np.hstack([arr, bg])
            X_store = sp.csr_matrix(X_full) if was_sparse else X_full
            var_full = pd.concat([a.var, new_var], axis=0)
            a_new = ad.AnnData(X=X_store, obs=a.obs.copy(), var=var_full)
            if layer and layer in a.layers:
                a_new.layers[layer] = X_store.copy()
            for k, v in a.obsm.items():
                a_new.obsm[k] = v.copy()
            for k, v in a.uns.items():
                a_new.uns[k] = v
            a = a_new
            data_changed = True

    if data_changed:
        arr_now, _ = _get_matrix(a, layer)
        if "total_counts" in a.obs:
            a.obs["total_counts"] = arr_now.sum(axis=1)
        if recompute_pca and "X_pca" in a.obsm:
            n_comp = int(a.obsm["X_pca"].shape[1]) if a.obsm["X_pca"].ndim == 2 else 20
            a.obsm["X_pca"] = _recompute_pca_from(arr_now, n_components=n_comp)

    a.uns["sam_augmentation"] = {
        "perturbation": perturbation,
        "level": level,
        "seed": seed,
    }
    return a


def _prereqs_met(adata: Any, perturbation: str, sample_key: str | None) -> tuple[bool, str | None]:
    reqs = AUGMENTATIONS[perturbation]["requires"]
    if "spatial" in reqs and "spatial" not in adata.obsm:
        return False, "missing obsm['spatial']"
    if "sample" in reqs:
        if not sample_key or sample_key not in adata.obs or adata.obs[sample_key].nunique() < 2:
            return False, "missing >=2 samples in obs[sample_key]"
    return True, None


def _score_one(
    values: pd.Series,
    levels: list[Any],
    orientation: int,
    role: str,
    *,
    min_rho: float,
    rel_tol: float,
    abs_tol: float,
    caveat: str | None,
) -> dict[str, Any]:
    from scipy.stats import spearmanr

    m = values.groupby(level=0).agg(["mean", "std"]).reindex(levels)
    if m["mean"].isna().all():
        return {
            "verdict": "skipped" if role != "must_detect" else "skipped",
            "v0": np.nan,
            "v_end": np.nan,
            "rho": np.nan,
            "threshold": np.nan,
            "first_level": None,
            "detail": caveat or "all metric values NaN",
        }
    v0 = float(m["mean"].iloc[0])
    v_end = float(m["mean"].iloc[-1]) if pd.notna(m["mean"].iloc[-1]) else np.nan
    sd0 = float(np.nan_to_num(m["std"].iloc[0]))
    o = orientation or 1
    bad = -o * (m["mean"] - v0)
    thr = max(3.0 * sd0, rel_tol * max(abs(v0), 1e-3), 1e-6)
    degrade_thr = max(thr, float(abs_tol))
    ok = ~bad.isna()
    if ok.sum() > 2 and np.nanstd(bad[ok]) > 0:
        rho = float(spearmanr(np.arange(ok.sum()), bad[ok].to_numpy())[0])
    elif ok.sum() == 2 and pd.notna(bad.iloc[-1]):
        rho = 1.0 if float(bad.iloc[-1]) > thr else (-1.0 if float(bad.iloc[-1]) < -thr else 0.0)
    else:
        rho = 0.0
    moved = bool(np.nanmax(np.abs(bad.to_numpy())) > thr)
    end_bad = float(bad.iloc[-1]) if pd.notna(bad.iloc[-1]) else 0.0
    report_thr = degrade_thr if role == "must_not_degrade" else thr
    first_lv = next((lv for lv, b in zip(levels, bad) if pd.notna(b) and float(b) > report_thr), None)

    if role == "must_detect":
        if moved and end_bad > thr and rho >= min_rho - 1e-9:
            verdict = "pass"
        elif moved and end_bad < -thr:
            verdict = "fail_inverted"
        else:
            verdict = "fail_blind"
    elif role == "must_be_invariant":
        verdict = "pass" if not moved else "fail_false_alarm"
    elif role == "must_not_degrade":
        degraded = bool(np.nanmax(bad.to_numpy()) > degrade_thr)
        verdict = "pass" if not degraded else "fail_false_alarm"
    else:
        verdict = "caveat_observed"
    return {
        "verdict": verdict,
        "v0": v0,
        "v_end": v_end,
        "rho": rho,
        "threshold": report_thr,
        "first_level": first_lv,
        "detail": caveat or "",
    }


def verify_metric_contracts(
    adata: Any,
    evaluators: Mapping[str, Callable[[Any, Any | None], Mapping[str, Any]]],
    contracts: Mapping[str, Mapping[str, Any]] | None = None,
    *,
    reference: Any | None = None,
    seeds: Sequence[int] = (0, 1),
    label_key: str = "cell_type",
    sample_key: str | None = None,
    layer: str | None = "counts",
    recompute_pca: bool = False,
    min_rho: float = 0.8,
    rel_tol: float = 0.02,
    abs_tol: float = 0.01,
    perturbation_kwargs: Mapping[str, Mapping[str, Any]] | None = None,
) -> dict[str, Any]:
    """Verify metric behaviors against controlled perturbations on any dataset.

    Parameters
    ----------
    adata
        Dataset to perturb and score.
    evaluators
        Mapping from evaluator name to ``fn(adata_perturbed, reference) -> {metric_name: value}``.
    contracts
        Contract specifications keyed by metric name. Defaults to :data:`METRIC_CONTRACTS` filtered
        to metrics returned by at least one evaluator on the unperturbed input.
    reference
        Optional matched reference ``AnnData`` passed untouched to each evaluator.
    seeds
        Random seeds per perturbation level.
    label_key
        Column in ``adata.obs`` holding cell-type labels.
    sample_key
        Optional sample/replicate column in ``adata.obs``.
    layer
        Count/expression layer forwarded to :func:`perturb_anndata`.
    recompute_pca
        Forwarded to :func:`perturb_anndata`.
    min_rho
        Minimum Spearman monotonicity required for ``"must_detect"`` checks.
    rel_tol
        Relative baseline fraction for significance threshold ``max(3*sd0, rel_tol*max(|v0|, 1e-3))``.
    perturbation_kwargs
        Optional per-perturbation extra keyword arguments forwarded to :func:`perturb_anndata`.

    Returns
    -------
    dict[str, Any]
        Result dictionary with ``status`` (``"ok"`` when ``n_failed == 0``, else ``"failed"``),
        counts (``n_checks``, ``n_passed``, ``n_failed``, ``n_caveats``, ``n_skipped``),
        ``summary`` DataFrame, and raw ``values`` DataFrame.
    """
    if not evaluators:
        return {"status": "skipped", "reason": "no evaluators provided"}

    # Probe baseline metrics cleanly so contracts apply only to active metrics
    base_metrics: dict[str, float] = {}
    for fn in evaluators.values():
        out = fn(adata, reference)
        for k, v in out.items():
            base_metrics[k] = np.nan if v is None else float(v)

    spec_map = dict(contracts if contracts is not None else METRIC_CONTRACTS)
    active_specs = {m: c for m, c in spec_map.items() if m in base_metrics}
    if not active_specs:
        return {"status": "skipped", "reason": "none of the contracted metrics were emitted by evaluators"}

    # Collect unique (perturbation, levels) needed across all active contracts
    schedule: dict[str, list[Any]] = {}
    for spec in active_specs.values():
        for role_key in ("must_detect", "must_be_invariant", "must_not_degrade"):
            for pert, lvls in spec.get(role_key, {}).items():
                cur = schedule.setdefault(pert, [])
                for lv in lvls:
                    if lv not in cur:
                        cur.append(lv)
        for pert in spec.get("known_caveats", {}):
            if pert not in schedule:
                schedule[pert] = list(AUGMENTATIONS[pert]["default_levels"])

    rows = []
    skipped_perts: dict[str, str] = {}
    p_kw = perturbation_kwargs or {}
    for pert, levels in schedule.items():
        ok, reason = _prereqs_met(adata, pert, sample_key)
        if not ok:
            skipped_perts[pert] = reason or "prerequisites not met"
            continue
        for lv in sorted(levels):
            for seed in seeds:
                aug = perturb_anndata(
                    adata,
                    pert,
                    lv,
                    seed=int(seed),
                    label_key=label_key,
                    sample_key=sample_key,
                    layer=layer,
                    recompute_pca=recompute_pca,
                    **dict(p_kw.get(pert, {})),
                )
                for fn in evaluators.values():
                    out = fn(aug, reference)
                    for m, v in out.items():
                        if m in active_specs:
                            rows.append((pert, lv, int(seed), m, np.nan if v is None else float(v)))

    val_df = pd.DataFrame(rows, columns=["perturbation", "level", "seed", "metric", "value"])
    summary_rows = []
    for metric, spec in active_specs.items():
        orient = int(spec.get("orientation", 1))
        caveats = spec.get("known_caveats", {})
        checks: list[tuple[str, str, list[Any], str | None]] = []
        for pert, lvls in spec.get("must_detect", {}).items():
            checks.append((pert, "must_detect", sorted(lvls), caveats.get(pert)))
        for pert, lvls in spec.get("must_be_invariant", {}).items():
            checks.append((pert, "must_be_invariant", sorted(lvls), caveats.get(pert)))
        for pert, lvls in spec.get("must_not_degrade", {}).items():
            checks.append((pert, "must_not_degrade", sorted(lvls), caveats.get(pert)))
        covered = set(spec.get("must_detect", {})) | set(spec.get("must_be_invariant", {})) | set(
            spec.get("must_not_degrade", {})
        )
        for pert, note in caveats.items():
            if pert not in covered:
                checks.append((pert, "known_caveat", sorted(AUGMENTATIONS[pert]["default_levels"]), note))

        for pert, role, lvls, note in checks:
            if pert in skipped_perts:
                summary_rows.append(
                    {
                        "metric": metric,
                        "perturbation": pert,
                        "role": role,
                        "verdict": "skipped",
                        "v0": np.nan,
                        "v_end": np.nan,
                        "rho": np.nan,
                        "threshold": np.nan,
                        "first_level": None,
                        "detail": skipped_perts[pert],
                    }
                )
                continue
            sub = val_df[(val_df["metric"] == metric) & (val_df["perturbation"] == pert)]
            res = _score_one(
                sub.set_index("level")["value"],
                lvls,
                orient,
                role,
                min_rho=min_rho,
                rel_tol=rel_tol,
                abs_tol=abs_tol,
                caveat=note,
            )
            summary_rows.append({"metric": metric, "perturbation": pert, "role": role, **res})

    sum_df = pd.DataFrame(summary_rows)
    n_passed = int((sum_df["verdict"] == "pass").sum()) if len(sum_df) else 0
    n_failed = int(sum_df["verdict"].str.startswith("fail_").sum()) if len(sum_df) else 0
    n_caveats = int((sum_df["verdict"] == "caveat_observed").sum()) if len(sum_df) else 0
    n_skipped = int((sum_df["verdict"] == "skipped").sum()) if len(sum_df) else 0
    return {
        "status": "ok" if n_failed == 0 else "failed",
        "n_checks": int(len(sum_df)),
        "n_passed": n_passed,
        "n_failed": n_failed,
        "n_caveats": n_caveats,
        "n_skipped": n_skipped,
        "summary": sum_df,
        "values": val_df,
    }
