"""Inter-sample consistency (reference-free, needs biological replicates): is each cell type reproducible
across samples? The ISCG / ISCL axes and per-type refinement regimes of scTypeEval
(Garnica, Andreatta & Carmona 2026)."""
from __future__ import annotations

import numpy as np

from ._utils import _dense
from ._validate import check_labels, get_counts, skip

__all__ = ["inter_sample_consistency", "isc_regime"]


# --------------------------------------------------------------------------- #
# Inter-sample consistency (ISC) - reference-free, needs biological replicates
# --------------------------------------------------------------------------- #
def _pseudobulk_profiles(adata, label_key, sample_key, *, min_cells, min_samples, layer=None):
    """One profile per (sample, type) with >= ``min_cells`` cells; a type is kept only when it has
    >= ``min_samples`` such profiles (the rest are returned in ``dropped``). With a count ``layer`` (read by
    :func:`get_counts`) the paper's recipe is used (sum raw counts, then log1p of counts-per-10k); with
    ``layer=None`` ``X`` is averaged as it is. Also returns the name of the matrix read."""
    y = check_labels(adata, label_key)
    s = check_labels(adata, sample_key)
    mat, matrix = get_counts(adata, layer) if layer else (adata.X, "X")
    X = _dense(mat).astype(float)
    profiles, labels, samples, dropped = [], [], [], {}
    for t in np.unique(y):
        rows = []
        for samp in np.unique(s):
            m = (y == t) & (s == samp)
            if m.sum() >= min_cells:
                rows.append((X[m].sum(0) if layer else X[m].mean(0), samp))
        if len(rows) < min_samples:
            dropped[t] = len(rows)
            continue
        for prof, samp in rows:
            profiles.append(prof)
            labels.append(t)
            samples.append(samp)
    P = np.vstack(profiles) if profiles else np.empty((0, X.shape[1]))
    if layer and P.size:
        P = np.log1p(P / np.maximum(P.sum(1, keepdims=True), 1.0) * 1e4)
    return P, np.asarray(labels), np.asarray(samples), dropped, matrix


def _two_label_silhouette(D, labels):
    """scTypeEval ``2label_silhouette``: a(i) = mean distance to own label (excl. self); b(i) = mean
    distance to ALL other-label observations pooled (not the nearest other cluster)."""
    n = len(labels)
    same = labels[:, None] == labels[None, :]
    np.fill_diagonal(same, False)
    other = labels[:, None] != labels[None, :]
    n_same, n_other = same.sum(1), other.sum(1)
    a = np.divide((D * same).sum(1), n_same, out=np.zeros(n), where=n_same > 0)
    b = np.divide((D * other).sum(1), n_other, out=np.zeros(n), where=n_other > 0)
    denom = np.maximum(a, b)
    return np.divide(b - a, denom, out=np.zeros(n), where=denom > 0)


def _rcm_dissimilarity(P, samples):
    """Reciprocal-classification match (RCM) distance between (sample, type) pseudobulks.

    Classifier = nearest profile by Spearman correlation (a stand-in for the paper's SingleR classic
    mode). D[i, j] = 0 when i's best hit inside j's sample is j AND j's best hit inside i's sample is
    i, else 1. A mutual cross-type swap also scores 0 - that is how RCM penalises confusable types
    (it lowers b(i) in the silhouette)."""
    from scipy.stats import rankdata

    R = rankdata(P, axis=1).astype(float)
    R -= R.mean(1, keepdims=True)
    R /= np.linalg.norm(R, axis=1, keepdims=True) + 1e-12
    C = R @ R.T
    n = len(samples)
    uniq = np.unique(samples)
    best = np.empty((n, len(uniq)), dtype=int)  # best hit of observation i inside each sample
    for k, smp in enumerate(uniq):
        cols = np.flatnonzero(samples == smp)
        best[:, k] = cols[np.argmax(C[:, cols], axis=1)]
    hit = best[:, np.searchsorted(uniq, samples)]  # hit[i, j] = best hit of i inside j's sample
    idx = np.arange(n)
    match = (hit == idx[None, :]) & (hit.T == idx[:, None])
    D = 1.0 - match.astype(float)
    np.fill_diagonal(D, 0.0)
    return D


# Garnica, Andreatta & Carmona 2026 (Fig. 5c): median - 1 IQR over a 9-atlas meta-study. Guides, not cut-offs.
ISC_THRESHOLDS = {"iscg": 0.74, "iscl": 0.45}


def isc_regime(iscg, iscl, thr_g: float = ISC_THRESHOLDS["iscg"], thr_l: float = ISC_THRESHOLDS["iscl"]) -> str:
    """Map a label's (ISCG, ISCL) to the Garnica 2026 refinement quadrant. A helper, not a metric: it
    returns the quadrant name as a ``str`` (not a result dict).

    ``consistent`` (keep) · ``remove`` (both low: multiplets, low-quality labels) · ``merge`` (high ISCG,
    low ISCL: compact but not replicable = over-partitioned) · ``split`` (low ISCG, high ISCL: matchable
    but no strong shared signature = too coarse / heterogeneous). ``not_computable`` if either is missing.
    """
    if iscg is None or iscl is None or not (np.isfinite(iscg) and np.isfinite(iscl)):
        return "not_computable"
    hi_g, hi_l = iscg >= thr_g, iscl >= thr_l
    if hi_g and hi_l:
        return "consistent"
    if not hi_g and not hi_l:
        return "remove"
    return "merge" if hi_g else "split"


def inter_sample_consistency(adata, label_key: str, sample_key: str, metric: str = "cosine", *,
                             min_cells: int = 10, min_samples: int = 5, layer: str | None = None,
                             thresholds: dict | None = None) -> dict:
    """Is each cell type reproducible across biological samples? (scTypeEval inter-sample consistency;
    Garnica, Andreatta & Carmona 2026, bioRxiv 10.64898/2026.06.10.731285.)

    One pseudobulk per (sample, type) with >= ``min_cells`` cells; a type needs >= ``min_samples`` of
    them or it is listed under ``dropped_types`` instead of being scored (the paper's defaults, 10 / 5;
    3 % of atlas labels fell below them). Pass ``layer`` naming the raw-count layer to follow the paper's
    recipe (sum counts, then log-normalise; a missing layer falls back to ``X`` only when ``X`` holds counts);
    with ``layer=None`` (default) ``X`` is averaged as-is, so hand it log-normalised HVGs. ``matrix`` names
    the matrix read. Fewer than 2 scorable types -> a skip that keeps ``n_profiles`` / ``dropped_types``.

    Returns the paper's two complementary axes plus the legacy headline:
      * ``iscg`` - GLOBAL: 2-label silhouette on pseudobulk cosine distance (compactness against all
        other types pooled; the most sensitive detector of signal degradation).
      * ``iscl`` - LOCAL: silhouette on the reciprocal-classification match distance (replicability;
        top-ranked overall in the benchmark, the detector of over-partitioning, least batch-sensitive).
      * ``integrated`` = iscg * iscl (the paper's product; R^2 > 0.98 with macro-F1 of predicted labels).
      * ``regime`` per type via :func:`isc_regime` (consistent / remove / merge / split).
      * ``consistency`` / ``per_type`` - silhouette on pseudobulk ``metric`` distance (the pre-0.5
        headline; ranked 3rd of 36 in the benchmark).
    """
    from scipy.spatial.distance import cdist
    from sklearn.metrics import silhouette_samples

    P, plabels, psamples, dropped, matrix = _pseudobulk_profiles(
        adata, label_key, sample_key, min_cells=min_cells, min_samples=min_samples, layer=layer)
    thr = {**ISC_THRESHOLDS, **(thresholds or {})}
    types = np.unique(plabels)
    if len(types) < 2 or len(plabels) <= len(types):
        return {**skip(f"need >= 2 types with >= {min_samples} samples of >= {min_cells} cells"),
                "consistency": None, "iscg": None, "iscl": None, "integrated": None,
                "n_profiles": int(len(plabels)), "dropped_types": dropped}

    def _per_type(v):
        return {t: float(v[plabels == t].mean()) for t in types}

    sil = silhouette_samples(P, plabels, metric=metric)
    D_cos = cdist(P, P, metric="cosine")
    np.fill_diagonal(D_cos, 0.0)
    glob = _two_label_silhouette(D_cos, plabels)
    loc = silhouette_samples(_rcm_dissimilarity(P, psamples), plabels, metric="precomputed")
    per_g, per_l = _per_type(glob), _per_type(loc)
    per_int = {t: per_g[t] * per_l[t] for t in types}
    return {
        "status": "ok", "consistency": float(sil.mean()), "per_type": _per_type(sil),
        "iscg": float(glob.mean()), "iscl": float(loc.mean()),
        "integrated": float(np.mean(list(per_int.values()))),
        "per_type_iscg": per_g, "per_type_iscl": per_l, "per_type_integrated": per_int,
        "regime": {t: isc_regime(per_g[t], per_l[t], thr["iscg"], thr["iscl"]) for t in types},
        "thresholds": thr, "n_profiles": int(len(plabels)), "dropped_types": dropped, "matrix": matrix,
    }
