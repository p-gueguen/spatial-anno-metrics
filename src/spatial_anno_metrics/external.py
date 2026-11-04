"""External harness (needs ground truth).

AvgBIO, accuracy / F1 / ARI / kappa / ECS, spatially aware ARI, hierarchical and composition accuracy, and the
deconvolution proportion metrics (OpenProblems, Spotless, Li et al. 2022) with a rank aggregate over methods.
"""
from __future__ import annotations

import numpy as np

from ._hierarchy import ancestors, lca_depth, validate_hierarchy
from ._utils import _embedding
from ._validate import check_labels, skip
from .thresholds import value as _thr

__all__ = ["avg_bio", "external_scores", "element_centric_similarity", "spatial_ari", "hierarchical_accuracy",
           "composition_accuracy", "deconvolution_metrics", "deconvolution_rank_score"]


def avg_bio(adata, label_key: str, truth_key: str, embedding: str | None = "X_pca") -> dict:
    """AvgBIO (Zhu 2026 / scIB): mean(ARI, NMI, scaled ASW), biological structure kept vs a ground truth.

    ARI/NMI compare ``label_key`` to ``truth_key``; ASW is the silhouette of the TRUE labels in ``embedding``
    (how well real types separate), scaled to [0, 1].

    Parameters
    ----------
    adata : anndata.AnnData
        Cells with both label columns in ``obs``.
    label_key : str
        ``obs`` column with the labeling to score.
    truth_key : str
        ``obs`` column with the ground truth.
    embedding : str or None
        ``obsm`` key of the embedding for the silhouette (None or missing: PCA of X).

    Returns
    -------
    dict
        ``{status, ari, nmi, asw, avg_bio}``; fewer than 2 truth labels gives a skip (silhouette undefined).
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

    pred = check_labels(adata, label_key)
    truth = check_labels(adata, truth_key)
    n_u = len(np.unique(truth))
    if n_u < 2 or n_u >= len(truth):
        return skip(f"need 2..n_cells-1 labels in adata.obs[{truth_key!r}], found {n_u} in {len(truth)} cells")
    ari = float(adjusted_rand_score(truth, pred))
    nmi = float(normalized_mutual_info_score(truth, pred))
    E = _embedding(adata, embedding)
    asw = float(silhouette_score(E, truth))
    return {"status": "ok", "ari": ari, "nmi": nmi, "asw": asw,
            "avg_bio": float(np.mean([ari, nmi, (asw + 1) / 2]))}


# --------------------------------------------------------------------------- #
# External harness (needs ground truth)
# --------------------------------------------------------------------------- #
def external_scores(*, y_true, y_pred) -> dict:
    """Reference-based scores of ``y_pred`` vs ground-truth ``y_true`` (benchmarks only).

    Thin wrapper over ``sklearn.metrics`` plus ``element_centric_similarity``.

    Parameters
    ----------
    y_true : array-like of str
        Ground-truth labels (keyword-only, sklearn order).
    y_pred : array-like of str
        Predicted labels.

    Returns
    -------
    dict
        ``{status, accuracy, balanced_accuracy, macro_f1, weighted_f1, precision_macro, recall_macro,
        per_class_f1, ari, ami, nmi, kappa, ecs}``.
    """
    from sklearn.metrics import (accuracy_score, adjusted_mutual_info_score,
                                  adjusted_rand_score, balanced_accuracy_score, cohen_kappa_score,
                                  f1_score, normalized_mutual_info_score, precision_score,
                                  recall_score)

    pred = check_labels(y_pred)
    truth = check_labels(y_true)
    if len(truth) != len(pred):
        raise ValueError(f"y_true ({len(truth)}) and y_pred ({len(pred)}) must have the same length")
    if len(truth) == 0:
        return skip("need >= 1 cell, found 0")
    labels = sorted(set(truth) | set(pred))
    per_class = f1_score(truth, pred, average=None, labels=labels, zero_division=0)
    return {
        "status": "ok", "accuracy": float(accuracy_score(truth, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, pred)),
        "macro_f1": float(f1_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "weighted_f1": float(f1_score(truth, pred, average="weighted", labels=labels, zero_division=0)),
        "precision_macro": float(precision_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "recall_macro": float(recall_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "per_class_f1": {c: float(v) for c, v in zip(labels, per_class)},
        "ari": float(adjusted_rand_score(truth, pred)),
        "ami": float(adjusted_mutual_info_score(truth, pred)),
        "nmi": float(normalized_mutual_info_score(truth, pred)),
        "kappa": float(cohen_kappa_score(truth, pred)),
        "ecs": element_centric_similarity(truth, pred)["ecs"],
    }


# --------------------------------------------------------------------------- #
# Element-Centric Similarity (Gates & Ahn 2019) - resolution-agnostic partition agreement
# --------------------------------------------------------------------------- #
def element_centric_similarity(a, b, alpha: float = _thr("external", "ecs_alpha")) -> dict:
    """Element-centric similarity (Gates et al. 2019) between two hard labelings of the same N elements.

    Exact, vectorized closed form for hard partitions (the paper's affinity-matrix definition, with the
    self-loop terms cancelling): for element k in a-cluster A (size nA) and b-cluster B (size nB) with overlap
    nAB, the L1 affinity distance is ``nAB*|α/nA − α/nB| + (nA−nAB)*α/nA + (nB−nAB)*α/nB`` and
    ``ECS_k = 1 − dist/(2α)``. Symmetric and resolution-agnostic; matches clusim 0.4 to 1e-12
    (``tests/goldens/ecs_clusim``).

    Parameters
    ----------
    a, b : array-like
        The two labelings.
    alpha : float
        Personalized-PageRank restart parameter; it cancels for hard partitions.

    Returns
    -------
    dict
        ``{status, ecs}``: the mean over elements, in [0, 1]; 1.0 iff the partitions are identical.
    """
    import pandas as pd

    if not 0 < alpha < 1:
        raise ValueError(f"alpha must be in (0, 1), got {alpha}")
    a = pd.factorize(check_labels(a))[0]
    b = pd.factorize(check_labels(b))[0]
    if len(a) != len(b):
        raise ValueError(f"a ({len(a)}) and b ({len(b)}) must have the same length")
    if len(a) == 0:
        return skip("need >= 1 element, found 0")
    _, inva, ca = np.unique(a, return_inverse=True, return_counts=True)
    _, invb, cb = np.unique(b, return_inverse=True, return_counts=True)
    nA = ca[inva].astype(float)
    nB = cb[invb].astype(float)
    pair = inva.astype(np.int64) * (invb.max() + 1) + invb
    _, invp, cp = np.unique(pair, return_inverse=True, return_counts=True)
    nAB = cp[invp].astype(float)
    dist = nAB * np.abs(alpha / nA - alpha / nB) + (nA - nAB) * (alpha / nA) + (nB - nAB) * (alpha / nB)
    return {"status": "ok", "ecs": float(np.mean(1.0 - dist / (2.0 * alpha)))}


# --------------------------------------------------------------------------- #
# spARI (Yan, Feng & Luo, Biometrics 2025) - spatially aware adjusted Rand index
# --------------------------------------------------------------------------- #
def spatial_ari(*, y_true, y_pred, coords, alpha: float = _thr("external", "spari_alpha"),
                max_cells: int | None = _thr("external", "spari_max_cells"), seed: int = 0) -> dict:
    """Spatially aware Rand index (spRI) and its chance-adjusted form (spARI) of ``y_pred`` vs ``y_true``.

    Yan, Feng & Luo, *Biometrics* 2025, doi 10.1093/biomtc/ujaf127 (preprint doi 10.1101/2025.03.25.645156).
    Like the Rand index, every pair of cells on which both labelings agree (together in both, apart in both)
    scores 1. A disagreeing pair scores by its distance ``d`` (coordinates min-max scaled to [0, 1] per axis):
    together in ``y_pred`` but apart in ``y_true`` scores ``f(d) = alpha * exp(-d^2)`` (a wrong merge of
    nearby cells is milder), apart in ``y_pred`` but together in ``y_true`` scores
    ``h(d) = alpha * (1 - exp(-d^2))`` (a wrong split of distant cells is milder). spRI is the mean pair
    score; spARI = (spRI - E) / (1 - E) with the authors' expectation
    ``E = pq + (1-p)(1-q) + (1-p) q mean(f) + p (1-q) mean(h)`` (p, q: fractions of pairs together in
    ``y_true`` / ``y_pred``; means over all pairs).
    Implemented from that definition (the authors' R package is GPL and was not copied). As in their package,
    a pair at distance 0 gets no spatial credit, and two one-cluster labelings score (1, 1). ``alpha = 0``
    gives the Rand index and ARI exactly. Asymmetric: swapping ``y_true`` and ``y_pred`` swaps f and h.

    Parameters
    ----------
    y_true : array-like of str
        Reference labels (keyword-only, sklearn order).
    y_pred : array-like of str
        Labels to score.
    coords : array-like, shape (n_cells, 2)
        Spatial coordinates, e.g. ``adata.obsm['spatial']``.
    alpha : float
        Weight scale in [0, 1); keeps a disagreeing pair strictly below an agreeing one.
    max_cells : int or None
        The pair sums are O(n^2): above this many cells, a random subset of ``max_cells`` (after scaling the
        coordinates on all cells) is scored. None scores every pair.
    seed : int
        Seed of the subsample.

    Returns
    -------
    dict
        ``{status, spri, spari, n_cells, subsampled}``; fewer than 2 cells gives a skip.
    """
    from scipy.spatial.distance import cdist

    truth, pred = check_labels(y_true), check_labels(y_pred)
    xy = np.asarray(coords, dtype=float)
    if xy.ndim != 2 or xy.shape[1] < 1 or xy.shape[0] != len(truth) or len(pred) != len(truth):
        raise ValueError(f"coords {xy.shape}, y_true ({len(truth)}) and y_pred ({len(pred)}) must have one "
                         "row per cell")
    if not np.all(np.isfinite(xy)):
        raise ValueError("coords must be finite (no NaN or infinite values)")
    if not 0 <= alpha < 1:
        raise ValueError(f"alpha must be in [0, 1), got {alpha}")
    if max_cells is not None and max_cells < 2:
        raise ValueError(f"max_cells must be >= 2 or None, got {max_cells}")
    if len(truth) < 2:
        return skip(f"need >= 2 cells, found {len(truth)}")
    span = np.ptp(xy, axis=0)
    xy = (xy - xy.min(0)) / np.where(span > 0, span, 1.0)
    n = len(truth)
    subsampled = max_cells is not None and n > max_cells
    if subsampled:
        keep = np.sort(np.random.default_rng(seed).choice(n, max_cells, replace=False))
        truth, pred, xy, n = truth[keep], pred[keep], xy[keep], max_cells
    base = {"n_cells": int(n), "subsampled": bool(subsampled)}
    r = np.unique(truth, return_inverse=True)[1]
    c = np.unique(pred, return_inverse=True)[1]
    if r.max() == 0 and c.max() == 0:
        return {"status": "ok", "spri": 1.0, "spari": 1.0, **base}

    npairs = n * (n - 1) / 2
    nij = np.zeros((r.max() + 1, c.max() + 1))
    np.add.at(nij, (r, c), 1)
    sa, sb = (nij.sum(1) ** 2).sum(), (nij.sum(0) ** 2).sum()
    agree = npairs + (nij ** 2).sum() - 0.5 * (sa + sb)       # Rand-index agreements
    s_sg = s_gs = f_all = h_all = 0.0                          # every unordered pair is seen twice below
    step = max(1, 2_000_000 // n)
    for i in range(0, n, step):
        d2 = cdist(xy[i:i + step], xy, "sqeuclidean")
        e = np.exp(-d2)
        f = np.where(d2 > 0, alpha * e, 0.0)
        h = np.where(d2 > 0, alpha * (1 - e), 0.0)
        same_r = r[i:i + step, None] == r[None, :]
        same_c = c[i:i + step, None] == c[None, :]
        s_sg += f[same_c & ~same_r].sum()
        s_gs += h[same_r & ~same_c].sum()
        f_all += f.sum()
        h_all += h.sum()
    spri = (agree + (s_sg + s_gs) / 2) / npairs
    p = 0.5 * (sa - n) / npairs
    q = 0.5 * (sb - n) / npairs
    expected = p * q + (1 - p) * (1 - q) + (1 - p) * q * f_all / 2 / npairs + p * (1 - q) * h_all / 2 / npairs
    spari = 1.0 if expected == 1 else (spri - expected) / (1 - expected)
    return {"status": "ok", "spri": float(spri), "spari": float(spari), **base}


# --------------------------------------------------------------------------- #
# Hierarchical & composition accuracy (external)
# --------------------------------------------------------------------------- #
def hierarchical_accuracy(*, y_true, y_pred, hierarchy: dict[str, str],
                          partial: float = _thr("external", "hierarchical_partial"),
                          credit: str = "lca") -> dict:
    """Accuracy with credit by distance on a cell-type tree, plus the accuracy of coarse (non-leaf) calls.

    ``credit='lca'`` (default since 0.7.0) scores each cell by the hierarchical F-measure of Kiritchenko et
    al. 2006, from the depth of the lowest common ancestor (LCA) of truth and prediction: with every label
    standing for itself plus its ancestors, ``hP = (lca + 1) / (depth(pred) + 1)`` and
    ``hR = (lca + 1) / (depth(true) + 1)`` (a root has depth 0); the credit is their harmonic mean (0 under
    different roots). An honest coarse call (the parent, 2/3 one level up) now outscores a wrong sibling
    subtype (1/2), which the old flat ``partial`` could not tell apart. ``credit='lineage'`` is the pre-0.7.0
    score: 1 exact, ``partial`` when truth and prediction share a root, else 0.

    Parameters
    ----------
    y_true : array-like of str
        Ground-truth labels (keyword-only, sklearn order).
    y_pred : array-like of str
        Predicted labels; may be coarse (an internal node of the tree).
    hierarchy : dict
        ``{child: parent}``; a label absent from it, or mapped to itself, is a root. Any depth.
    partial : float
        Credit for a shared root under ``credit='lineage'``; unused under ``'lca'``.
    credit : {'lca', 'lineage'}
        Scoring rule (see above).

    Returns
    -------
    dict
        ``{status, credit, hierarchical_accuracy, subtype_accuracy, lineage_accuracy, hierarchical_precision,
        hierarchical_recall, coarse_call_frac, coarse_call_accuracy}``. ``hierarchical_accuracy`` is the mean
        per-cell credit; ``subtype_accuracy`` exact matches; ``lineage_accuracy`` the fraction sharing a root
        (for a two-level map, the old same-parent fraction); ``hierarchical_precision`` / ``_recall`` the
        micro-averaged hP / hR. A coarse call is a prediction that is some label's parent in ``hierarchy``;
        ``coarse_call_accuracy`` is the fraction of those that are the truth or one of its ancestors, NaN when
        there are none.
    """
    import pandas as pd

    if credit not in ("lca", "lineage"):
        raise ValueError(f"credit must be 'lca' or 'lineage', got {credit!r}")
    if not 0.0 <= partial <= 1.0:
        raise ValueError(f"partial must be in [0, 1], got {partial}")
    validate_hierarchy(hierarchy)
    pred = check_labels(y_pred)
    truth = check_labels(y_true)
    if len(truth) != len(pred):
        raise ValueError(f"y_true ({len(truth)}) and y_pred ({len(pred)}) must have the same length")
    if len(truth) == 0:
        return skip("need >= 1 cell, found 0")
    codes, uniq = pd.factorize(pd.MultiIndex.from_arrays([truth, pred]))
    internal = {v for k, v in hierarchy.items() if k != v}
    rows = []
    for t, p in uniq:                                          # one tree walk per distinct (truth, pred) pair
        lca, dt, dp = lca_depth(t, p, hierarchy), len(ancestors(t, hierarchy)), len(ancestors(p, hierarchy))
        shared = lca + 1                                       # nodes on both root paths (0: different roots)
        hp, hr = shared / (dp + 1), shared / (dt + 1)
        hf = 2 * hp * hr / (hp + hr) if shared else 0.0
        score = hf if credit == "lca" else (1.0 if t == p else (partial if shared else 0.0))
        rows.append((score, t == p, shared > 0, shared, dp + 1, dt + 1, p in internal, lca == dp))
    v = np.array(rows, dtype=float)[codes]
    coarse = v[:, 6].astype(bool)
    return {
        "status": "ok", "credit": credit, "hierarchical_accuracy": float(v[:, 0].mean()),
        "subtype_accuracy": float(v[:, 1].mean()), "lineage_accuracy": float(v[:, 2].mean()),
        "hierarchical_precision": float(v[:, 3].sum() / v[:, 4].sum()),
        "hierarchical_recall": float(v[:, 3].sum() / v[:, 5].sum()),
        "coarse_call_frac": float(coarse.mean()),
        "coarse_call_accuracy": float(v[coarse, 7].mean()) if coarse.any() else float("nan"),
    }


def composition_accuracy(*, y_true, y_pred) -> dict:
    """Predicted vs true cell-type proportions from hard labels: L1, Pearson r and Jensen-Shannon divergence.

    Parameters
    ----------
    y_true : array-like of str
        Ground-truth labels (keyword-only, sklearn order).
    y_pred : array-like of str
        Predicted labels.

    Returns
    -------
    dict
        ``{status, l1, pearson, jsd, predicted, true}``; ``pearson`` is NaN when undefined (one type, or equal
        shares on either side).
    """
    import pandas as pd
    from scipy.spatial.distance import jensenshannon

    pred = pd.Series(check_labels(y_pred))
    truth = pd.Series(check_labels(y_true))
    if len(truth) != len(pred):
        raise ValueError(f"y_true ({len(truth)}) and y_pred ({len(pred)}) must have the same length")
    if len(truth) == 0:
        return skip("need >= 1 cell, found 0")
    types = sorted(set(truth) | set(pred))
    p = pred.value_counts(normalize=True).reindex(types).fillna(0).to_numpy()
    q = truth.value_counts(normalize=True).reindex(types).fillna(0).to_numpy()
    l1 = float(np.abs(p - q).sum())
    pear = float(np.corrcoef(p, q)[0, 1]) if len(types) > 1 and p.std() > 0 and q.std() > 0 else float("nan")
    jsd = float(jensenshannon(p, q) ** 2)                         # jensenshannon returns the distance (sqrt of div)
    return {"status": "ok", "l1": l1, "pearson": pear, "jsd": jsd,
            "predicted": dict(zip(types, p.tolist())), "true": dict(zip(types, q.tolist()))}


def _ssim_li(x, y) -> float:
    """SSIM of two proportion vectors as in Li et al. 2022 (each scaled by its max; NaN for an all-zero one).

    Ported from QuKunLab/SpatialBenchmarking ``Benchmarking/GenesMetrics.py`` (``cal_ssim`` + ``scale_max``),
    BSD-2-Clause, Copyright (c) 2021 Qu Kun Lab. ``M`` there is the larger max of the scaled vectors, i.e. 1.
    """
    with np.errstate(invalid="ignore", divide="ignore"):
        x, y = x / x.max(), y / y.max()
    mx, my = x.mean(), y.mean()
    sx, sy = np.sqrt(((x - mx) ** 2).mean()), np.sqrt(((y - my) ** 2).mean())
    sxy = ((x - mx) * (y - my)).mean()
    L = max(x.max(), y.max())
    c1, c2 = (0.01 * L) ** 2, (0.03 * L) ** 2
    c3 = c2 / 2
    return float((2 * mx * my + c1) / (mx ** 2 + my ** 2 + c1) * (2 * sx * sy + c2) / (sx ** 2 + sy ** 2 + c2)
                 * (sxy + c3) / (sx * sy + c3))


def _aupr_precrec(present, scores) -> float:
    """Area under the non-linearly interpolated PR curve (Davis & Goadrich 2006; R ``precrec::evalmod``).

    Spotless (Sang-aram et al. 2024) evaluates presence/absence AUPR via ``precrec::evalmod``, which applies
    Davis & Goadrich (2006, *ICML*, doi 10.1145/1143844.1143874) non-linear interpolation across tied score
    steps with ``dTP > 1`` and integrates with the trapezoidal rule. When all predicted scores are distinct,
    this equals ``sklearn.metrics.auc(recall, precision)`` and differs from ``average_precision_score``
    (right-Riemann step sum) in the last digits.
    """
    y_t = np.asarray(present, dtype=bool).ravel()
    y_s = np.asarray(scores, dtype=float).ravel()
    n_pos = int(y_t.sum())
    if n_pos == 0 or n_pos == y_t.size:
        return float("nan")
    order = np.argsort(-y_s, kind="mergesort")
    y_t, y_s = y_t[order], y_s[order]
    distinct = np.r_[np.flatnonzero(np.diff(y_s)), y_t.size - 1]
    tps = np.r_[0.0, np.cumsum(y_t)[distinct].astype(float)]
    fps = np.r_[0.0, (1 + distinct - tps[1:]).astype(float)]
    rec_pts = [0.0]
    prec_pts = [tps[1] / (tps[1] + fps[1]) if (tps[1] + fps[1]) > 0 else 1.0]
    for i in range(1, len(tps)):
        dtp = int(round(tps[i] - tps[i - 1]))
        dfp = fps[i] - fps[i - 1]
        if dtp > 1:
            xs = np.arange(1, dtp + 1, dtype=float)
            tp_x = tps[i - 1] + xs
            fp_x = fps[i - 1] + (dfp / dtp) * xs
            rec_pts.extend((tp_x / n_pos).tolist())
            prec_pts.extend((tp_x / (tp_x + fp_x)).tolist())
        else:
            tp_x, fp_x = tps[i], fps[i]
            rec_pts.append(tp_x / n_pos)
            prec_pts.append(tp_x / (tp_x + fp_x) if (tp_x + fp_x) > 0 else 1.0)
    trapz = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    return float(trapz(prec_pts, rec_pts))


def deconvolution_metrics(*, true_prop, pred_prop, aupr_method: str = "sklearn") -> dict:
    """Score predicted vs ground-truth cell-type proportion matrices (spots x cell types).

    The soft, per-spot sibling of ``composition_accuracy``: score a deconvolver's per-spot proportions
    (RCTD / TACCO / cell2location) against known mixtures such as pseudo-spots. ``r2`` and ``jsd`` are the
    OpenProblems ``task_spatial_decomposition`` metrics (``tests/goldens/deconv_openproblems``); ``aupr`` and
    ``aupr_precrec`` are the Spotless presence/absence scores (Sang-aram et al., *eLife* 2024,
    doi 10.7554/elife.88431); ``ssim`` is Li et al. 2022 (*Nat Methods*, doi 10.1038/s41592-022-01480-9).

    Parameters
    ----------
    true_prop : array-like or DataFrame, shape (n_spots, n_types)
        Ground-truth proportions (keyword-only).
    pred_prop : array-like or DataFrame, shape (n_spots, n_types)
        Predicted proportions, same rows and columns in the same order.
    aupr_method : {'sklearn', 'precrec'}
        How ``aupr`` is computed: ``'sklearn'`` (default) uses ``sklearn.metrics.average_precision_score``
        (rectangular step integration); ``'precrec'`` uses Davis & Goadrich (2006) non-linear PR interpolation
        with trapezoidal integration as in R ``precrec::evalmod`` (Spotless). ``aupr_precrec`` is always
        included in the returned dict.

    Returns
    -------
    dict
        ``{status, r2, jsd, jsd_per_spot, rmse, pearson, aupr, aupr_precrec, ssim, ssim_per_type}``:

        * ``r2`` - ``sklearn.r2_score(true, pred, multioutput='uniform_average')``, per cell type across
          spots, averaged (higher better).
        * ``jsd`` - mean over cell types of the JS distance between true and predicted columns (lower better).
          A type absent from every true spot is undefined there and left out (OpenProblems returns NaN).
        * ``jsd_per_spot`` - mean over spots of the row-wise JS distance.
        * ``rmse`` - RMSE over the whole matrix (lower better).
        * ``pearson`` - mean per-cell-type Pearson r across spots (higher better).
        * ``aupr`` - presence/absence PR score for "type present in spot" (true > 0) over all spot x type
          entries using ``aupr_method``; NaN when every entry, or none, is present.
        * ``aupr_precrec`` - Spotless / ``precrec::evalmod`` Davis-Goadrich interpolated PR AUC.
        * ``ssim_per_type`` - Li et al. SSIM of each type's column, both scaled by their max (keys: the
          DataFrame columns, else column indices); NaN for an all-zero column. ``ssim`` is their mean over the
          defined types (higher better).
    """
    from scipy.spatial.distance import jensenshannon
    from sklearn.metrics import average_precision_score, r2_score

    if aupr_method not in ("sklearn", "precrec"):
        raise ValueError(f"aupr_method must be 'sklearn' or 'precrec', got {aupr_method!r}")
    if hasattr(true_prop, "columns") and hasattr(pred_prop, "columns"):
        if list(true_prop.columns) != list(pred_prop.columns):
            raise ValueError("true_prop and pred_prop DataFrames have different columns or column order")
    true = np.asarray(true_prop, dtype=float)
    pred = np.asarray(pred_prop, dtype=float)
    if true.ndim != 2 or pred.ndim != 2 or true.shape != pred.shape or true.shape[1] < 1:
        raise ValueError(f"proportion matrices must be 2D of matching shape (n_spots, n_types >= 1): "
                         f"{true.shape} vs {pred.shape}")
    if not (np.all(np.isfinite(true)) and np.all(np.isfinite(pred))):
        raise ValueError("proportion matrices must be finite (no NaN or infinite values)")
    if (true < 0).any() or (pred < 0).any():
        raise ValueError("proportion matrices must be non-negative")
    names = list(getattr(true_prop, "columns", range(true.shape[1])))
    if true.shape[0] < 2:
        return skip(f"need >= 2 spots, found {true.shape[0]}")
    n_ct = true.shape[1]

    r2 = float(r2_score(true, pred, multioutput="uniform_average"))
    jsd_ct = float(np.nanmean([jensenshannon(true[:, k], pred[:, k]) for k in range(n_ct)]))
    jsd_spot = float(np.nanmean([jensenshannon(true[i], pred[i]) for i in range(true.shape[0])]))
    rmse = float(np.sqrt(np.mean((true - pred) ** 2)))
    pears = [np.corrcoef(true[:, k], pred[:, k])[0, 1]
             for k in range(n_ct) if true[:, k].std() > 0 and pred[:, k].std() > 0]
    pearson = float(np.mean(pears)) if pears else float("nan")
    present = (true > 0).ravel()
    aupr_sk = (float(average_precision_score(present, pred.ravel()))
               if 0 < present.sum() < present.size else float("nan"))
    aupr_pr = _aupr_precrec(present, pred.ravel())
    aupr = aupr_sk if aupr_method == "sklearn" else aupr_pr
    ssim = {name: _ssim_li(true[:, k], pred[:, k]) for k, name in enumerate(names)}
    defined = [v for v in ssim.values() if np.isfinite(v)]
    return {"status": "ok", "r2": r2, "jsd": jsd_ct, "jsd_per_spot": jsd_spot, "rmse": rmse,
            "pearson": pearson, "aupr": aupr, "aupr_precrec": aupr_pr,
            "ssim": float(np.mean(defined)) if defined else float("nan"), "ssim_per_type": ssim}


_HIGHER_IS_BETTER = {"r2": True, "pearson": True, "aupr": True, "aupr_precrec": True, "ssim": True,
                     "rmse": False, "jsd": False, "jsd_per_spot": False}


def deconvolution_rank_score(results: dict, metrics=("rmse", "jsd_per_spot", "aupr")) -> dict:
    """Rank-aggregated score of several deconvolution methods on one dataset.

    Each metric ranks the methods (1 = best, ties share their mean rank); a method's score is its mean rank
    over ``metrics``. The default metrics are the three Spotless ranks on (RMSE, JSD, AUPR; Sang-aram et al.,
    *eLife* 2024, doi 10.7554/elife.88431); Li et al. 2022's accuracy score averages ranks the same way over
    PCC, SSIM, RMSE and JSD (pass ``metrics=('pearson', 'ssim', 'rmse', 'jsd')``). Spotless further scales
    and averages over datasets; aggregate this function's per-dataset ranks for that.

    Parameters
    ----------
    results : dict
        ``{method: deconvolution_metrics(...) result}`` (any dict holding the metric keys).
    metrics : sequence of str
        Keys to rank on; each must be one of r2, pearson, aupr, aupr_precrec, ssim (higher better) or rmse,
        jsd, jsd_per_spot (lower better).

    Returns
    -------
    dict
        ``{status, metrics, mean_rank, ranks, order}``: ``mean_rank`` ``{method: rank}`` (lower better),
        ``ranks`` ``{metric: {method: rank}}``, ``order`` the methods best first. Fewer than 2 methods, or a
        NaN metric value, gives a skip naming it.
    """
    from scipy.stats import rankdata

    metrics = list(metrics)
    if not metrics:
        raise ValueError("metrics must not be empty")
    unknown = [m for m in metrics if m not in _HIGHER_IS_BETTER]
    if unknown:
        raise ValueError(f"no known direction for {unknown}; use {sorted(_HIGHER_IS_BETTER)}")
    methods = list(results)
    if len(methods) < 2:
        return skip(f"need >= 2 methods to rank, found {len(methods)}")
    missing = [(m, k) for m in methods for k in metrics
               if not isinstance(results[m], dict) or k not in results[m] or results[m][k] is None]
    if missing:
        return skip(f"undefined metric values (method, metric): {missing}")
    vals = np.array([[float(results[m][k]) for k in metrics] for m in methods])
    bad = [(m, k) for i, m in enumerate(methods) for j, k in enumerate(metrics)
           if not np.isfinite(vals[i, j])]
    if bad:
        return skip(f"undefined metric values (method, metric): {bad}")
    ranks = {k: dict(zip(methods, rankdata(-vals[:, j] if _HIGHER_IS_BETTER[k] else vals[:, j]).tolist()))
             for j, k in enumerate(metrics)}
    mean_rank = {m: float(np.mean([ranks[k][m] for k in metrics])) for m in methods}
    return {"status": "ok", "metrics": metrics, "mean_rank": mean_rank, "ranks": ranks,
            "order": sorted(methods, key=mean_rank.get)}
