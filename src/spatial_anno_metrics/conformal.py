"""Conformal prediction sets, calibration error and conformal novelty p-values for cell-type annotation.

Distribution-free, model-agnostic typing uncertainty (set size 1 = confident, >1 = ambiguous, 0 = novel), the
calibration check that comes before reading an annotator's probabilities as probabilities (ECE), and novelty
p-values in the reference's feature space.
"""
from __future__ import annotations

import numpy as np

from ._utils import _dense
from ._validate import check_labels, skip
from .thresholds import value as _thr

__all__ = ["conformal_prediction_sets", "calibration_error", "conformal_novelty"]
_MIN_CAL = _thr("conformal", "min_class_calibration")
_EMBED_Q = _thr("conformal", "clustered_quantiles")


def _qhat(scores, alpha):
    """Split-conformal quantile: the ceil((n + 1)(1 - alpha)) / n empirical quantile (capped at the max)."""
    n = len(scores)
    if n == 0:
        return float("nan")                                       # no calibration: no threshold
    level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
    return float(np.quantile(scores, level, method="higher"))


def _clustered_thresholds(s_cal, y_cal, classes, alpha, n_clusters, frac_clustering, seed):
    """Clustered conformal (Ding et al. 2023, arXiv 2306.09335, s. 3), written from the paper.

    Each class's calibration scores are split: a ``frac_clustering`` share embeds the class as its score
    quantiles at ``_EMBED_Q`` + {1 - alpha}; classes with fewer than n_alpha = 1/min(alpha, 0.1) - 1
    clustering scores go to the null cluster (-1). Weighted k-means (weight = sqrt(count)) groups the
    embedded classes; each cluster's threshold is the conformal quantile of its pooled held-out scores, the
    null cluster's the marginal one over all held-out scores. Defaults: n_min = rarest class count,
    n~ = max(n_min, n_alpha), K = classes with >= n~ cells, gamma = K / (75 + K), M = floor(gamma n~ / 2).
    """
    from sklearn.cluster import KMeans

    rng = np.random.default_rng(seed)
    n_alpha = int(np.ceil(1 / min(alpha, 0.1) - 1))
    counts = {c: int((y_cal == c).sum()) for c in classes}
    present = [c for c in classes if counts[c] > 0]
    n_tilde = max(min(counts[c] for c in present), n_alpha)
    k_big = sum(counts[c] >= n_tilde for c in present)
    gamma = k_big / (75 + k_big) if frac_clustering is None else float(frac_clustering)
    m = int(np.floor(gamma * n_tilde / 2)) if n_clusters is None else int(n_clusters)

    is_clu = np.zeros(len(y_cal), bool)
    for c in present:
        idx = rng.permutation(np.flatnonzero(y_cal == c))
        is_clu[idx[: int(np.floor(gamma * len(idx)))]] = True
    levels = [*_EMBED_Q, 1 - alpha]
    emb, w, embedded = [], [], []
    for c in present:
        sc = np.sort(s_cal[is_clu & (y_cal == c)])
        n = len(sc)
        if n < n_alpha:
            continue
        emb.append([np.quantile(sc, min(1.0, np.ceil((n + 1) * t) / n), method="higher") for t in levels])
        w.append(np.sqrt(n))
        embedded.append(c)
    cluster = {c: -1 for c in classes}
    m = min(max(m, 1), len(embedded))
    if m >= 1:
        km = KMeans(n_clusters=m, n_init=10, random_state=seed).fit(np.array(emb), sample_weight=np.array(w))
        cluster.update({c: int(k) for c, k in zip(embedded, km.labels_)})
    held = ~is_clu
    q_null = _qhat(s_cal[held], alpha)
    q_of = {k: _qhat(s_cal[held & np.isin(y_cal, [c for c in embedded if cluster[c] == k])], alpha)
            for k in set(cluster.values()) - {-1}}
    thr = {c: (q_null if cluster[c] == -1 or np.isnan(q_of[cluster[c]]) else q_of[cluster[c]]) for c in classes}
    return thr, cluster


# --------------------------------------------------------------------------- #
# Conformal prediction sets (catalog s3d) - distribution-free typing uncertainty
# --------------------------------------------------------------------------- #
def conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha: float = _thr("conformal", "alpha"),
                              class_conditional: bool | str = False,
                              hierarchy: dict | None = None, y_query=None, *, n_clusters: int | None = None,
                              frac_clustering: float | None = None, seed: int = 0) -> dict:
    """Distribution-free conformal prediction sets for cell-type annotation.

    Catalog s3d "conformal set size", the cleanest reference-free formalization of abstention.
    Model-agnostic: bring any annotator's per-cell class-probability matrix (RCTD weights,
    SingleR/CellTypist/scANVI posteriors, marker-score softmax). Split-conformal LAC: the
    calibration nonconformity is ``s = 1 - p(true class)``; the (1-alpha) quantile ``qhat`` sets
    the inclusion threshold and the query set is ``C(x) = {y : p(x)_y >= 1 - qhat}``. Set size
    **1 = confident, >1 = ambiguous, 0 = novel/OOD**, with the finite-sample guarantee
    ``P(true in C(x)) >= 1 - alpha`` **provided calibration and query are exchangeable**. On
    native spatial data calibrated against a dissociated reference, exchangeability is broken by
    platform shift, so the nominal coverage is approximate (recalibrate on a small platform-matched
    labeled subset, e.g. protein-confirmed cells) - the set SIZE stays a valid relative uncertainty.
    An empty set is a weak novelty signal: under depth shift it also empties for known types; use
    :func:`conformal_novelty` for a novelty test with FDR control. Depends on numpy (k-means: sklearn).

    Parameters
    ----------
    proba_query : array-like, shape (n_q, K)
        Query class probabilities.
    classes : array-like, shape (K,)
        Labels matching the probability columns.
    proba_cal : array-like, shape (n_c, K)
        Calibration probabilities.
    y_cal : array-like, shape (n_c,)
        True labels of the calibration cells.
    alpha : float
        Miscoverage rate (target coverage 1 - alpha).
    class_conditional : bool or "clustered"
        False: one marginal threshold. True: per-true-class (Mondrian) thresholds, recommended when every
        class has enough calibration cells (a class with fewer than ``min_class_calibration`` uses the
        marginal one); marginal conformal silently under-covers hard rare types. ``"clustered"``: clustered
        Mondrian conformal (Ding et al. 2023), which pools classes with similar score distributions - the
        choice when many classes have few calibration cells each.
    hierarchy : dict, optional
        ``{leaf: parent}`` map; also returns lineage-collapsed sets (scConform's ontology idea: an ambiguous
        SUBTYPE set becomes a confident LINEAGE call).
    y_query : array-like, optional
        True query labels; also returns empirical coverage, per-type coverage and CovGap.
    n_clusters, frac_clustering : int, float, optional
        Clustered mode only: number of clusters and share of each class's calibration cells used to embed
        it. None = Ding et al.'s heuristic (see ``_clustered_thresholds``).
    seed : int
        Clustered mode only: the clustering split and k-means seed.

    Returns
    -------
    dict
        ``{status, sets, set_size, summary, [clusters], [collapsed_sets, pct_ambiguous_one_parent],
        [coverage, per_type_coverage, covgap]}``. ``clusters`` maps each class to its cluster (-1 = null
        cluster, marginal threshold). ``covgap`` = mean over the query's true types of
        ``|coverage_y - (1 - alpha)|`` (Ding et al. 2023 CovGap, as a fraction, not x100). An empty
        calibration set is a skip (no threshold can be set).
    """
    if class_conditional not in (False, True, "clustered"):
        raise ValueError(f"class_conditional must be False, True or 'clustered', got {class_conditional!r}")
    if len(y_cal) == 0:
        return skip("empty calibration set: conformal sets need calibration probabilities and true labels")
    proba_query = np.asarray(proba_query, dtype=float)
    proba_cal = np.asarray(proba_cal, dtype=float)
    classes = np.asarray([str(c) for c in classes])
    y_cal = check_labels(y_cal)
    cls_idx = {c: j for j, c in enumerate(classes)}

    true_col = np.array([cls_idx[c] for c in y_cal])
    s_cal = 1.0 - proba_cal[np.arange(len(y_cal)), true_col]

    clusters = None
    if class_conditional == "clustered":
        thr, clusters = _clustered_thresholds(s_cal, y_cal, classes, alpha, n_clusters, frac_clustering, seed)
    elif class_conditional:
        thr = {c: (_qhat(s_cal[y_cal == c], alpha) if (y_cal == c).sum() >= _MIN_CAL else _qhat(s_cal, alpha))
               for c in classes}
    if class_conditional:
        incl_thresh = np.array([1.0 - thr[c] for c in classes])      # per-class inclusion threshold
        include = proba_query >= incl_thresh[None, :]
    else:
        include = proba_query >= (1.0 - _qhat(s_cal, alpha))

    sets = [[classes[j] for j in np.where(include[i])[0]] for i in range(include.shape[0])]
    size = include.sum(axis=1)
    out = {
        "status": "ok", "sets": sets, "set_size": size,
        "summary": {
            "mean_set_size": float(size.mean()) if size.size else float("nan"),
            "pct_singleton": float((size == 1).mean()) if size.size else float("nan"),
            "pct_ambiguous": float((size > 1).mean()) if size.size else float("nan"),
            "pct_empty_novel": float((size == 0).mean()) if size.size else float("nan"),
            "alpha": alpha, "class_conditional": class_conditional,
        },
    }
    if clusters is not None:
        out["clusters"] = clusters
    if hierarchy is not None:
        collapsed = [sorted({hierarchy.get(x, x) for x in s}) for s in sets]
        out["collapsed_sets"] = collapsed
        amb = size > 1
        out["pct_ambiguous_one_parent"] = float(
            sum(1 for i in range(len(sets)) if amb[i] and len(collapsed[i]) == 1) / max(1, amb.sum()))
    if y_query is not None:
        y_query = check_labels(y_query)
        covered = np.array([y_query[i] in sets[i] for i in range(len(sets))])
        out["coverage"] = float(covered.mean())
        out["per_type_coverage"] = {
            c: float(np.mean([y_query[i] in sets[i] for i in np.where(y_query == c)[0]]))
            for c in np.unique(y_query)}
        out["covgap"] = float(np.mean([abs(v - (1 - alpha)) for v in out["per_type_coverage"].values()]))
    return out


# --------------------------------------------------------------------------- #
# Calibration error (ROADMAP #24) - is the annotator's confidence a probability?
# --------------------------------------------------------------------------- #
def _binned_gap(conf, hit, n_bins):
    """sum_b |B_b|/n * |mean(hit in b) - mean(conf in b)| over equal-width (lo, hi] bins (0 joins bin 0)."""
    b = np.searchsorted(np.linspace(0, 1, n_bins + 1)[1:-1], conf, side="left")
    n = np.bincount(b, minlength=n_bins)
    ok = n > 0
    acc = np.bincount(b, weights=hit, minlength=n_bins)[ok] / n[ok]
    cf = np.bincount(b, weights=conf, minlength=n_bins)[ok] / n[ok]
    return float(np.sum(n[ok] / len(conf) * np.abs(acc - cf)))


def calibration_error(proba, classes, y_true, *, n_bins: int = _thr("conformal", "ece_n_bins")) -> dict:
    """Measure the expected calibration error (ECE) of an annotator's class probabilities.

    Check this before reading a native confidence (a softmax, an RCTD weight) as a probability. ECE (Guo et
    al. 2017, arXiv 1706.04599) bins the top-label confidence into ``n_bins`` equal-width bins and averages
    ``|accuracy - mean confidence|`` per bin, weighted by the bin's share of cells. Class-wise ECE (Kull et
    al. 2019) does the same one-vs-rest for every class's probability column and averages over classes, so
    it also sees miscalibration of the non-top classes and of rare types. Both are 0 for a calibrated
    annotator and grow with over- or under-confidence; both need true labels.

    Parameters
    ----------
    proba : array-like, shape (n, K)
        Class probabilities.
    classes : array-like, shape (K,)
        Labels matching the probability columns.
    y_true : array-like, shape (n,)
        True labels; every one must be in ``classes``.
    n_bins : int
        Number of equal-width confidence bins (Guo et al. use 15).

    Returns
    -------
    dict
        ``{status, ece, classwise_ece, per_class_ece, accuracy, mean_confidence, n}``; no cells is a skip.
        Bins are (lo, hi], as in sklearn's ``calibration_curve``.
    """
    if len(y_true) == 0:
        return skip("no cells: calibration error needs probabilities and true labels")
    y = check_labels(y_true)
    classes = np.asarray([str(c) for c in classes])
    missing = sorted(set(y) - set(classes))
    if missing:
        raise ValueError(f"y_true labels not in classes: {missing[:5]}")
    P = np.asarray(proba, dtype=float)
    top = P.argmax(1)
    conf = P[np.arange(len(y)), top]
    hit = (classes[top] == y).astype(float)
    per = {c: _binned_gap(P[:, k], (y == c).astype(float), n_bins) for k, c in enumerate(classes)}
    return {"status": "ok", "ece": _binned_gap(conf, hit, n_bins), "classwise_ece": float(np.mean(list(per.values()))),
            "per_class_ece": per, "accuracy": float(hit.mean()), "mean_confidence": float(conf.mean()),
            "n": int(len(y))}


# --------------------------------------------------------------------------- #
# Conformal novelty p-values (ROADMAP #20)
# --------------------------------------------------------------------------- #
def _conformal_pvalues(s_cal, s_query):
    """p_j = (1 + #{i : s_cal_i >= s_j}) / (n_cal + 1): higher score = more novel."""
    s_cal = np.sort(np.asarray(s_cal, float))
    n_ge = len(s_cal) - np.searchsorted(s_cal, np.asarray(s_query, float), side="left")
    return (1.0 + n_ge) / (len(s_cal) + 1.0)


def _bh(p):
    """Benjamini-Hochberg adjusted p-values (q-values)."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(max=1.0)
    out = np.empty_like(q)
    out[o] = q
    return out


def conformal_novelty(X_query, X_ref, *, method: str = "pca", fdr: float = _thr("conformal", "novelty_fdr"),
                      n_components: int = _thr("conformal", "novelty_n_components"),
                      frac_train: float = _thr("conformal", "novelty_frac_train"), seed: int = 0) -> dict:
    """Conformal novelty p-values for query cells in the reference's feature space, BH-adjusted.

    Outlier detection with conformal p-values (Bates et al. 2023, Ann Stat, doi 10.1214/22-AOS2244): a share
    ``frac_train`` of the reference cells fits a one-class nonconformity score, the rest calibrate it, and a
    query cell's p-value is the share of calibration cells at least as nonconforming. These marginal p-values
    are PRDS, so Benjamini-Hochberg at level ``fdr`` controls the false discovery rate among the cells flagged
    novel - **provided non-novel query cells are exchangeable with the reference cells**. A depth or platform
    shift between reference and section breaks that, and known cells are then flagged too (as they empty
    conformal sets): feed features in which the shift is small (e.g. a reference thinned to the section's
    depth, or a shift-robust embedding) and read a high novel share with that in mind.

    Parameters
    ----------
    X_query : array-like or sparse, shape (n_q, G)
        Query features (the same G features, same transform, as ``X_ref``).
    X_ref : array-like or sparse, shape (n_r, G)
        Reference features: cells of the known types only.
    method : {"pca", "iforest"}
        Nonconformity score. "pca": squared reconstruction error from a PCA fitted on the training split
        (``n_components`` components). "iforest": minus sklearn ``IsolationForest.score_samples``.
    fdr : float
        Benjamini-Hochberg level for the ``novel`` call.
    n_components : int
        PCA components ("pca"; capped by the training size and the feature count).
    frac_train : float
        Share of the reference that fits the score; the rest calibrates the p-values.
    seed : int
        Seed of the train/calibration split and of the isolation forest.

    Returns
    -------
    dict
        ``{status, pvalues, qvalues, novel, score, summary: {n_novel, frac_novel, fdr, method, n_train,
        n_calibration, min_attainable_p}}``. ``novel`` is ``qvalues <= fdr``. Fewer than 2 training or 1
        calibration reference cells is a skip.
    """
    if method not in ("pca", "iforest"):
        raise ValueError(f"method must be 'pca' or 'iforest', got {method!r}")
    Xq = np.asarray(_dense(X_query), dtype=float)
    Xr = np.asarray(_dense(X_ref), dtype=float)
    if Xq.shape[1] != Xr.shape[1]:
        raise ValueError(f"X_query has {Xq.shape[1]} features, X_ref {Xr.shape[1]}: use the same features")
    perm = np.random.default_rng(seed).permutation(len(Xr))
    n_tr = int(round(frac_train * len(Xr)))
    tr, cal = Xr[perm[:n_tr]], Xr[perm[n_tr:]]
    if len(tr) < 2 or len(cal) < 1:
        return skip(f"{len(Xr)} reference cells: need >= 2 to fit the score and >= 1 to calibrate it")
    if method == "pca":
        from sklearn.decomposition import PCA

        pca = PCA(n_components=max(1, min(n_components, len(tr) - 1, Xr.shape[1])), random_state=seed).fit(tr)

        def score(X):
            return ((X - pca.inverse_transform(pca.transform(X))) ** 2).sum(1)
    else:
        from sklearn.ensemble import IsolationForest

        forest = IsolationForest(random_state=seed).fit(tr)

        def score(X):
            return -forest.score_samples(X)
    s_q = score(Xq)
    p = _conformal_pvalues(score(cal), s_q)
    q = _bh(p) if len(p) else p
    novel = q <= fdr
    return {"status": "ok", "pvalues": p, "qvalues": q, "novel": novel, "score": s_q,
            "summary": {"n_novel": int(novel.sum()), "frac_novel": float(novel.mean()) if len(novel) else float("nan"),
                        "fdr": fdr, "method": method, "n_train": int(len(tr)), "n_calibration": int(len(cal)),
                        "min_attainable_p": 1.0 / (len(cal) + 1)}}
