"""Reference-free internal validity of a labeling on an embedding.

The scTypeEval family (silhouette, 2-label silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch,
average similarity, Davies-Bouldin, Calinski-Harabasz and a geometric-mean ``integrated`` score), plus the scib
rare-label scores (isolated-label ASW / F1), cLISI and the Leiden / k-means NMI and ARI of an unsupervised
clustering against the labels (Luecken et al. 2022, Nat Methods, doi 10.1038/s41592-021-01336-8).

All of them are circular on an embedding trained with the labels (scANVI, a fine-tuned foundation model):
they then score the training, not the annotation ("Metric Mirages", doi 10.1101/2024.04.02.587824,
preprint). Compute them on a label-agnostic embedding, such as a PCA of the counts.
"""
from __future__ import annotations

import numpy as np

from ._utils import _embedding, _igraph, _scale01
from ._validate import check_labels, skip
from .thresholds import value as _thr

__all__ = ["internal_validity", "isolated_label_scores", "clisi", "leiden_label_agreement",
           "kmeans_label_agreement"]
_F1_CUTS = _thr("internal", "isolated_f1_cuts_per_label")

_CIRCULAR = ("embedding_label_trained=True: internal metrics on an embedding trained with these labels are "
             "circular - they score the training, not the annotation (Metric Mirages, doi "
             "10.1101/2024.04.02.587824); use a label-agnostic embedding such as a PCA of the counts")


def _subsample(n: int, subsample: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return np.arange(n) if n <= subsample else rng.choice(n, subsample, replace=False)


# --------------------------------------------------------------------------- #
# Internal cluster-validity (scTypeEval family) - reference-free
# --------------------------------------------------------------------------- #
def internal_validity(adata, label_key: str = "cell_type", embedding: str | None = "X_pca",
                      n_neighbors: int = _thr("internal", "n_neighbors"),
                      subsample: int = _thr("internal", "subsample"), seed: int = 0, *,
                      embedding_label_trained: bool = False) -> dict:
    """Reference-free cluster-validity of a labeling on an embedding (scTypeEval family).

    Neighborhood purity, orbital-medoid and average similarity use every cell. The O(N^2) parts - the two
    silhouettes, Ward-PropMatch, Davies-Bouldin and Calinski-Harabasz - use a random ``subsample`` of cells,
    so the function runs at 1e5+ cells; with ``subsample >= n_obs`` every cell is used.

    Parameters
    ----------
    adata
        AnnData with labels in ``obs[label_key]``.
    label_key
        obs column holding the labels; NaN labels raise ``ValueError``.
    embedding
        ``obsm`` key of the embedding; when absent, a PCA of ``X`` (see ``_utils._embedding``).
    n_neighbors
        Embedding kNN size for neighborhood purity.
    subsample
        Cells drawn (without replacement, ``seed``) for the O(N^2) parts.
    seed
        Seed of the subsample (and of the PCA fallback).
    embedding_label_trained
        Say True when the embedding was trained with these labels (scANVI, a fine-tuned foundation model):
        every metric here is then circular (Metric Mirages, doi 10.1101/2024.04.02.587824), and the result
        carries a ``warning``.

    Returns
    -------
    dict
        ``{status, silhouette, silhouette_2label, neighborhood_purity, orbital_medoid, ward_propmatch,
        avg_similarity, davies_bouldin, calinski_harabasz, integrated, n_used, nan_reasons}`` (+ ``warning``).
        ``integrated`` is the geometric mean of the first six (each scaled to [0, 1]) that are defined - a
        single "is this labeling internally coherent?" number. The silhouettes, Ward-PropMatch,
        Davies-Bouldin and Calinski-Harabasz are NaN when the subsample has fewer than 2 labels or only
        singletons (the 2-label silhouette also when no label has 2 cells); ``nan_reasons`` names each NaN
        metric and why. Davies-Bouldin (lower = better) and Calinski-Harabasz are raw, not in
        ``integrated``. Fewer than 2 labels -> a skip.
    """
    from scipy.spatial.distance import cdist
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import (calinski_harabasz_score, davies_bouldin_score, pairwise_distances,
                                 silhouette_score)
    from sklearn.neighbors import NearestNeighbors

    y = check_labels(adata, label_key)
    types, code = np.unique(y, return_inverse=True)
    n, K = len(y), len(types)
    if K < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {K}")
    E = _embedding(adata, embedding, seed=seed)

    # neighborhood purity - fraction of k neighbors (excluding self) sharing the label; all cells.
    kk = int(min(n_neighbors, n - 1))
    full = NearestNeighbors(n_neighbors=kk + 1).fit(E).kneighbors(E, return_distance=False)
    hit = full == np.arange(n)[:, None]
    pos = np.where(hit.any(1), hit.argmax(1), kk)
    col = np.broadcast_to(np.arange(kk), (n, kk))
    idx = np.take_along_axis(full, col + (col >= pos[:, None]), axis=1)
    purity = float(np.mean((code[idx] == code[:, None]).mean(1)))

    # orbital-medoid - fraction of cells whose nearest type-medoid is their own type; (n, K) memory.
    M = np.empty((K, E.shape[1]))
    for j in range(K):
        pts = E[code == j]
        M[j] = pts[np.argmin(((pts - pts.mean(0)) ** 2).sum(1))]  # actual cell nearest the centroid
    orbital = float((cdist(E, M, "sqeuclidean").argmin(1) == code).mean())

    # average similarity - mean cosine of each cell to its own-type centroid (cohesion); a zero vector counts 0.
    C = np.vstack([E[code == j].mean(0) for j in range(K)])[code]
    norm = np.linalg.norm(E, axis=1) * np.linalg.norm(C, axis=1)
    cos = np.divide(np.einsum("ij,ij->i", E, C), norm, out=np.zeros(n), where=norm > 0)
    avg_sim = float(np.mean(cos))

    # O(N^2) parts on the subsample: silhouettes, Ward-PropMatch, Davies-Bouldin, Calinski-Harabasz.
    sel = _subsample(n, subsample, seed)
    Es, ys = E[sel], code[sel]
    labs = np.unique(ys)
    nan, reasons = float("nan"), {}
    sil = db = ch = ward_prop = nan
    if 2 <= len(labs) <= len(ys) - 1:
        sil = float(silhouette_score(Es, ys))
        db = float(davies_bouldin_score(Es, ys))
        ch = float(calinski_harabasz_score(Es, ys))
        ward = AgglomerativeClustering(n_clusters=len(labs), linkage="ward").fit_predict(Es)
        ward_prop = float(np.mean([np.bincount(ward[ys == t], minlength=len(labs)).max() / (ys == t).sum()
                                   for t in labs]))
    else:
        why = f"the {len(ys)}-cell subsample has {len(labs)} label(s); needs 2 to n - 1"
        reasons = dict.fromkeys(("silhouette", "ward_propmatch", "davies_bouldin", "calinski_harabasz"), why)
    D = pairwise_distances(Es)
    np.fill_diagonal(D, np.nan)
    s2 = []
    for t in labs:
        own = ys == t
        if own.sum() < 2 or (~own).sum() < 1:
            continue
        a_in = np.nanmean(D[np.ix_(own, own)], axis=1)             # mean within-type distance
        b_out = D[np.ix_(own, ~own)].mean(1)                       # mean distance to all others
        denom = np.maximum(a_in, b_out)
        s2.extend(np.divide(b_out - a_in, denom, out=np.zeros_like(denom), where=denom > 0))
    sil_2label = float(np.mean(s2)) if s2 else nan
    if not s2:
        reasons["silhouette_2label"] = (f"no label has 2 cells next to another label in the {len(ys)}-cell "
                                        f"subsample ({len(labs)} label(s))")

    parts = np.array([_scale01(sil), _scale01(sil_2label), purity, orbital, ward_prop, _scale01(avg_sim)])
    parts = parts[np.isfinite(parts)]
    integrated = float(np.exp(np.mean(np.log(np.clip(parts, 1e-9, 1.0)))))
    out = {
        "status": "ok", "silhouette": sil, "silhouette_2label": sil_2label, "neighborhood_purity": purity,
        "orbital_medoid": orbital, "ward_propmatch": ward_prop, "avg_similarity": avg_sim,
        "davies_bouldin": db, "calinski_harabasz": ch,
        "integrated": integrated, "n_used": int(len(sel)), "nan_reasons": reasons,
    }
    if embedding_label_trained:
        out["warning"] = _CIRCULAR
    return out


# --------------------------------------------------------------------------- #
# scib rare-label scores and cLISI (Luecken et al. 2022)
# --------------------------------------------------------------------------- #
def isolated_label_scores(adata, label_key: str = "cell_type", *, batch_key: str | None = None,
                          embedding: str | None = "X_pca", iso_threshold: int | None = None,
                          subsample: int = _thr("internal", "subsample"), seed: int = 0) -> dict:
    """Isolated-label ASW and F1 (scib): how well the labels seen in the fewest batches stand apart.

    A label is isolated when it occurs in at most ``iso_threshold`` batches (default: the minimum over labels).
    Without ``batch_key`` there is one batch, so every label is isolated and the scores are macro averages
    over labels: a rare label weighs as much as an abundant one, which the cell-weighted silhouette of
    :func:`internal_validity` does not do.

    * ASW, as scib-metrics ``isolated_labels``: the silhouette of every cell under all labels, rescaled to
      [0, 1] as ``(s + 1) / 2``, averaged over the cells of each isolated label.
    * F1, as scib's ``isolated_labels_f1``: the best F1 between "cell has the label" and "cell is in cluster
      c", over every cluster c of a clustering sweep. scib sweeps Leiden resolutions; here the sweep is every
      cut of one Ward tree into 2 to ``2 * n_labels`` clusters (``THRESHOLDS['internal']
      ['isolated_f1_cuts_per_label']``) (no Leiden dependency), so values are
      comparable in kind, not to the digit.

    Parameters
    ----------
    adata
        AnnData with labels in ``obs[label_key]``.
    label_key
        obs column holding the labels; NaN labels raise ``ValueError``.
    batch_key
        obs column holding the batch (sample, section); None = one batch.
    embedding
        ``obsm`` key of the embedding; when absent, a PCA of ``X``.
    iso_threshold
        Most batches an isolated label may occur in; None = the minimum over labels (scib's default).
    subsample
        Above this many cells, each label keeps at most ``max(2, subsample // n_labels)`` random cells: the
        per-label silhouette is a mean over clusters, so it is unbiased under per-label sampling, while F1 is
        then measured on that label-balanced sample. At or below it every cell is used, as in scib.
    seed
        Seed of the per-label sampling (and of the PCA fallback).

    Returns
    -------
    dict
        ``{status, asw, f1, isolated_labels, per_label: {label: {asw, f1, n}}, n_used}``; ``asw`` and ``f1``
        are means over the isolated labels, ``n`` the label's cells used. Fewer than 2 labels, or only
        singletons in the sample -> a skip.

    Notes
    -----
    Port of scib-metrics 0.5.6 ``isolated_labels`` (Copyright (c) 2022 Adam Gayoso, BSD-3-Clause,
    https://github.com/YosefLab/scib-metrics) and of scib's F1 rule (MIT); Luecken et al. 2022, Nat Methods,
    doi 10.1038/s41592-021-01336-8.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    from sklearn.metrics import silhouette_samples

    y = check_labels(adata, label_key)
    types, code = np.unique(y, return_inverse=True)
    n, K = len(y), len(types)
    if K < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {K}")
    if batch_key is None:
        n_batches = np.ones(K, dtype=int)
    else:
        b = check_labels(adata, batch_key)
        n_batches = np.array([len(np.unique(b[code == j])) for j in range(K)])
    iso = np.flatnonzero(n_batches <= (n_batches.min() if iso_threshold is None else iso_threshold))
    if len(iso) == 0:
        return skip(f"no label occurs in <= {iso_threshold} batch(es); min is {int(n_batches.min())}")
    E = _embedding(adata, embedding, seed=seed)

    if n <= subsample:
        sel = np.arange(n)
    else:
        rng, cap = np.random.default_rng(seed), max(2, subsample // K)
        sel = np.sort(np.concatenate([rng.permutation(np.flatnonzero(code == j))[:cap] for j in range(K)]))
    Es, ys = E[sel], code[sel]
    if not 2 <= len(np.unique(ys)) <= len(ys) - 1:
        return skip(f"the {len(ys)}-cell sample has {len(np.unique(ys))} label(s); the silhouette needs 2 to n - 1")

    s = (silhouette_samples(Es, ys) + 1) / 2
    Z = linkage(Es, "ward")
    cuts = [fcluster(Z, k, "maxclust") for k in range(2, min(_F1_CUTS * K, len(ys)) + 1)]
    per = {}
    for j in iso:
        own = ys == j
        f1 = max(float((2 * np.bincount(cl[own], minlength=cl.max() + 1)
                        / (np.bincount(cl) + own.sum())).max()) for cl in cuts)
        per[str(types[j])] = {"asw": float(s[own].mean()), "f1": f1, "n": int(own.sum())}
    return {"status": "ok", "asw": float(np.mean([v["asw"] for v in per.values()])),
            "f1": float(np.mean([v["f1"] for v in per.values()])), "isolated_labels": list(per),
            "per_label": per, "n_used": int(len(sel))}


def _simpson(D: np.ndarray, Idx: np.ndarray, code: np.ndarray, K: int, perplexity: float,
             tol: float = 1e-5, n_iter: int = 50) -> np.ndarray:
    """Per-cell Simpson index of the labels around each cell, Gaussian-weighted to a fixed perplexity.

    Port of scib-metrics 0.5.6 ``utils/_lisi.py`` (JAX -> numpy, vectorised over cells; Copyright (c) 2022
    Adam Gayoso, BSD-3-Clause, https://github.com/YosefLab/scib-metrics), itself after Harmony's LISI
    (Korsunsky et al. 2019). Self edges get weight 0; the kernel precision ``beta`` is bisected (at most
    ``n_iter`` steps) until the entropy is ``log(perplexity)`` within ``tol``. A cell whose weights all vanish
    is NaN (scib-metrics returns -1 there).
    """
    mask = Idx != np.arange(len(Idx))[:, None]
    target = np.log(perplexity)

    def hbeta(d, m, beta):
        P = np.exp(-d * beta[:, None]) * m
        sumP = P.sum(1)
        ok = sumP > 0
        H = np.zeros(len(d))
        H[ok] = np.log(sumP[ok]) + beta[ok] * (d[ok] * P[ok]).sum(1) / sumP[ok]
        P[ok] /= sumP[ok, None]
        return H, P

    n = len(D)
    beta, bmin, bmax = np.ones(n), np.full(n, -np.inf), np.full(n, np.inf)
    H, P = hbeta(D, mask, beta)
    Hdiff = H - target
    for _ in range(n_iter):
        act = np.flatnonzero(np.abs(Hdiff) >= tol)
        if not len(act):
            break
        up = Hdiff[act] > 0                                        # entropy too high -> sharpen the kernel
        b = beta[act]
        bmin[act[up]], bmax[act[~up]] = b[up], b[~up]
        lo, hi = bmin[act], bmax[act]
        beta[act] = np.where(up, np.where(np.isinf(hi), b * 2, (b + hi) / 2),
                             np.where(np.isinf(lo), b / 2, (b + lo) / 2))
        H[act], P[act] = hbeta(D[act], mask[act], beta[act])
        Hdiff[act] = H[act] - target
    S = np.zeros((n, K))
    np.add.at(S, (np.repeat(np.arange(n), Idx.shape[1]), code[Idx].ravel()), P.ravel())
    simpson = (S ** 2).sum(1)
    simpson[H == 0] = np.nan
    return simpson


def clisi(adata, label_key: str = "cell_type", *, embedding: str | None = "X_pca",
          n_neighbors: int = _thr("internal", "clisi_n_neighbors"), perplexity: float | None = None) -> dict:
    """Cell-type LISI (cLISI), scaled as in scib: 1 = every neighbourhood holds one label.

    Per cell, the LISI is the inverse Simpson index of the labels among its ``n_neighbors`` nearest
    neighbours, Gaussian-weighted to an effective neighbourhood of ``perplexity`` cells (Korsunsky et al.
    2019); ``clisi = (n_labels - median LISI) / (n_labels - 1)`` (Luecken et al. 2022, Nat Methods,
    doi 10.1038/s41592-021-01336-8). The median makes it blind to a small label: pair it with
    :func:`isolated_label_scores`.

    Parameters
    ----------
    adata
        AnnData with labels in ``obs[label_key]``.
    label_key
        obs column holding the labels; NaN labels raise ``ValueError``.
    embedding
        ``obsm`` key of the embedding; when absent, a PCA of ``X``.
    n_neighbors
        kNN size, the cell itself included (it gets weight 0); scib-metrics' benchmark uses 90.
    perplexity
        Effective neighbourhood size; None = ``floor(n_neighbors / 3)``, scib-metrics' rule.

    Returns
    -------
    dict
        ``{status, clisi, median_lisi, n_neighbors, perplexity}``. Fewer than 2 labels -> a skip.

    Notes
    -----
    Port of scib-metrics 0.5.6 ``clisi_knn`` / ``lisi_knn`` (Copyright (c) 2022 Adam Gayoso, BSD-3-Clause).
    scib-metrics runs in float32 (JAX); this runs in float64. On the test fixtures it agrees to 1e-6.
    """
    from sklearn.neighbors import NearestNeighbors

    if perplexity is not None and perplexity <= 1.0:
        raise ValueError(f"perplexity must be > 1, got {perplexity}")
    y = check_labels(adata, label_key)
    types, code = np.unique(y, return_inverse=True)
    K = len(types)
    if K < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {K}")
    k = int(min(n_neighbors, len(y)))
    perp = float(np.floor(k / 3)) if perplexity is None else float(perplexity)
    if perp <= 1.0 or k <= perp:
        return skip(f"need n_neighbors > perplexity > 1 (got n_neighbors={k}, perplexity={perp})")
    E = _embedding(adata, embedding)
    D, Idx = NearestNeighbors(n_neighbors=k).fit(E).kneighbors(E)
    med = float(np.nanmedian(1.0 / _simpson(D, Idx, code, K, perp)))
    return {"status": "ok", "clisi": (K - med) / (K - 1), "median_lisi": med, "n_neighbors": k,
            "perplexity": perp}


# --------------------------------------------------------------------------- #
# scib NMI / ARI of an unsupervised clustering against the labels
# --------------------------------------------------------------------------- #
def _umap_connectivities(D: np.ndarray, Idx: np.ndarray, n_iter: int = 64, tol: float = 1e-5,
                         min_scale: float = 1e-3):
    """Fuzzy-union kNN connectivities of UMAP (``local_connectivity=1``, ``set_op_mix_ratio=1``), as CSR.

    Port of umap-learn 0.5 ``fuzzy_simplicial_set`` on a precomputed kNN (``smooth_knn_dist`` +
    ``compute_membership_strengths``; numba -> numpy, vectorised over cells; Copyright (c) 2017 Leland McInnes,
    BSD-3-Clause, https://github.com/lmcinnes/umap), the graph scib-metrics 0.5.6
    ``NeighborsResults.knn_graph_connectivities`` and scanpy cluster on. Per cell, ``rho`` is the first non-zero
    distance and ``sigma`` is bisected (at most ``n_iter`` steps) until ``sum_j exp(-(d_j - rho) / sigma)`` over
    the neighbours after the first is ``log2(k)`` within ``tol``; the edge weight is ``exp(-(d - rho) / sigma)``
    (1 at or below ``rho``, 0 to self), and the graph is the fuzzy union ``A + A.T - A * A.T``. float32, as umap.
    """
    import scipy.sparse as sp

    D = np.asarray(D, dtype=np.float32)
    n, k = D.shape
    target = np.log2(float(k))
    nz = D > 0
    rho = np.where(nz.any(1), D[np.arange(n), nz.argmax(1)], 0).astype(np.float32)
    d = D[:, 1:] - rho[:, None]                                    # umap skips the first column (self)
    lo, hi, mid = np.zeros(n), np.full(n, np.inf), np.ones(n)
    act = np.arange(n)
    for _ in range(n_iter):
        if not len(act):
            break
        da = d[act]
        terms = np.where(da > 0, np.exp(-(da / mid[act, None])), 1.0)
        psum = np.cumsum(terms, axis=1)[:, -1]                    # summed left to right, as umap's loop
        go = np.abs(psum - target) >= tol
        act, psum = act[go], psum[go]
        over, under = act[psum > target], act[psum <= target]
        hi[over] = mid[over]
        mid[over] = (lo[over] + hi[over]) / 2.0
        lo[under] = mid[under]
        mid[under] = np.where(np.isinf(hi[under]), mid[under] * 2, (lo[under] + hi[under]) / 2.0)
    sigma = mid.astype(np.float32)
    floor = (min_scale * np.where(rho > 0, D.mean(1), D.mean())).astype(np.float32)
    sigma = np.where(sigma < floor, floor, sigma)

    dr = D - rho[:, None]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        val = np.where((dr <= 0) | (sigma[:, None] == 0), np.float32(1.0), np.exp(-(dr / sigma[:, None])))
    val[Idx == np.arange(n)[:, None]] = 0.0
    A = sp.coo_matrix((val.astype(np.float32).ravel(), (np.repeat(np.arange(n), k), np.asarray(Idx).ravel())),
                      shape=(n, n))
    A.eliminate_zeros()
    T = A.transpose()
    P = A.multiply(T)
    C = 1.0 * (A + T - P) + 0.0 * P                                # umap's expression at set_op_mix_ratio = 1
    C.eliminate_zeros()
    return sp.csr_matrix(C)


def _leiden_graph(ig, C):
    """Build the igraph graph of connectivities ``C`` as scib-metrics 0.5.6 ``_compute_clustering_leiden`` does.

    Directed from ``C``'s nonzeros in row-major order (igraph's ``Weighted_Adjacency``), then made undirected
    keeping each direction as its own edge (scanpy's construction). Built once and reused for every resolution
    (scib-metrics rebuilds the same graph per resolution).
    """
    C = C.tocoo()
    g = ig.Graph(C.shape[0], edges=list(zip(C.row.tolist(), C.col.tolist())), directed=True,
                 edge_attrs={"weight": C.data.tolist()})
    g.to_undirected(mode="each")
    return g


def _leiden(ig, g, resolution: float, seed: int) -> np.ndarray:
    """Leiden (modularity) on graph ``g`` at ``resolution``, as scib-metrics 0.5.6 ``_compute_clustering_leiden``.

    igraph's random generator is ``random.Random(seed)`` for the call and is restored to the ``random`` module
    (igraph's default) after it.
    """
    import random

    ig.set_random_number_generator(random.Random(seed))
    try:
        part = g.community_leiden(objective_function="modularity", weights="weight", resolution=resolution)
    finally:
        ig.set_random_number_generator(random)
    return np.asarray(part.membership)


def leiden_label_agreement(adata, label_key: str = "cell_type", embedding: str | None = "X_pca", *,
                           n_neighbors: int = _thr("internal", "leiden_n_neighbors"),
                           resolutions=_thr("internal", "leiden_resolutions"), seed: int = 0,
                           embedding_label_trained: bool = False) -> dict:
    """NMI and ARI between a Leiden clustering of the embedding and the labels (scib's Leiden NMI / ARI).

    Do the labels follow the structure an unsupervised clustering finds? The embedding's ``n_neighbors``-NN
    graph (self included) gets UMAP connectivities, Leiden (modularity) runs at every resolution in
    ``resolutions``, and the clustering with the best NMI to the labels is kept; ``ari`` is that clustering's
    ARI, not the best ARI of the sweep. Both are 1 when the clustering reproduces the labels and about 0 for
    labels unrelated to the embedding. Choosing the resolution by the labels flatters a labeling a little; the
    whole ``sweep`` is returned.

    Needs igraph, the ``leiden`` extra (``pip install "spatial-anno-metrics[leiden]"``); without it the result
    is a skip. :func:`kmeans_label_agreement` is the dependency-free sibling.

    Parameters
    ----------
    adata
        AnnData with labels in ``obs[label_key]``.
    label_key
        obs column holding the labels; NaN labels raise ``ValueError``.
    embedding
        ``obsm`` key of the embedding; when absent, a PCA of ``X`` (see ``_utils._embedding``).
    n_neighbors
        kNN size, the cell itself included (scib-metrics' benchmark uses 15); capped at ``n_obs``.
    resolutions
        Leiden resolutions to sweep; one value fixes the resolution (scib-metrics'
        ``optimize_resolution=False``).
    seed
        Seed of igraph's random generator for each Leiden run (and of the PCA fallback); scib-metrics uses 42.
    embedding_label_trained
        Say True when the embedding was trained with these labels (scANVI, a fine-tuned foundation model):
        the agreement is then circular (Metric Mirages, doi 10.1101/2024.04.02.587824), and the result
        carries a ``warning``.

    Returns
    -------
    dict
        ``{status, nmi, ari, resolution, n_clusters, n_labels, n_neighbors, clusters, sweep, n_used}``
        (+ ``warning``): ``clusters`` is the kept clustering (an int per cell, obs order), ``sweep`` a list of
        ``{resolution, nmi, ari, n_clusters}``, one per resolution. Every cell is used; no value is NaN.
        Fewer than 2 labels, or igraph missing -> a skip.

    Notes
    -----
    Mirrors scib-metrics 0.5.6 ``nmi_ari_cluster_labels_leiden`` (Copyright (c) 2022 Adam Gayoso, BSD-3-Clause,
    https://github.com/YosefLab/scib-metrics): the resolution sweep ``2 * x / 10`` for x in 1..10, UMAP
    connectivities of the 15-neighbour graph, igraph's ``community_leiden(objective_function="modularity",
    weights="weight", resolution=r)`` with its defaults (2 iterations, beta 0.01) and a fresh
    ``random.Random(seed)`` per resolution, arithmetic-mean NMI, and the first resolution with the highest
    NMI. The kNN here is exact (scikit-learn); scib-metrics' benchmark uses approximate pynndescent. Given the
    same kNN and seed, the clustering, NMI and ARI equal scib-metrics' (``tests/test_internal.py``).
    """
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    from sklearn.neighbors import NearestNeighbors

    y = check_labels(adata, label_key)
    types, code = np.unique(y, return_inverse=True)
    K = len(types)
    if K < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {K}")
    res = [float(r) for r in np.atleast_1d(resolutions)]
    if not res or min(res) <= 0:
        raise ValueError(f"resolutions must be one or more values > 0, got {resolutions!r}")
    if n_neighbors < 2:
        raise ValueError(f"n_neighbors must be >= 2 (the cell itself is one of them), got {n_neighbors}")
    try:
        ig = _igraph()
    except ImportError as e:
        return skip(f"leiden_label_agreement needs igraph, the 'leiden' extra of spatial-anno-metrics ({e})")
    E = _embedding(adata, embedding, seed=seed)
    k = int(min(n_neighbors, len(y)))
    D, Idx = NearestNeighbors(n_neighbors=k).fit(E).kneighbors(E)
    g = _leiden_graph(ig, _umap_connectivities(D, Idx))

    sweep, clusterings = [], []
    for r in res:
        cl = _leiden(ig, g, r, seed)
        clusterings.append(cl)
        sweep.append({"resolution": r,
                      "nmi": float(normalized_mutual_info_score(code, cl, average_method="arithmetic")),
                      "ari": float(adjusted_rand_score(code, cl)), "n_clusters": int(cl.max() + 1)})
    best = int(np.argmax([s["nmi"] for s in sweep]))
    out = {"status": "ok", "nmi": sweep[best]["nmi"], "ari": sweep[best]["ari"],
           "resolution": sweep[best]["resolution"], "n_clusters": sweep[best]["n_clusters"], "n_labels": K,
           "n_neighbors": k, "clusters": clusterings[best], "sweep": sweep, "n_used": int(len(y))}
    if embedding_label_trained:
        out["warning"] = _CIRCULAR
    return out


def kmeans_label_agreement(adata, label_key: str = "cell_type", embedding: str | None = "X_pca", *,
                           n_init: int = _thr("internal", "kmeans_n_init"), seed: int = 0,
                           embedding_label_trained: bool = False) -> dict:
    """NMI and ARI between a k-means clustering of the embedding (k = number of labels) and the labels.

    scib-metrics' cheaper stand-in for its Leiden NMI / ARI (:func:`leiden_label_agreement`): no graph, no
    resolution sweep, no optional dependency. With k fixed to the number of labels, k-means splits a large
    label and merges small ones when label sizes are uneven, so it is harsher than Leiden there.

    Parameters
    ----------
    adata
        AnnData with labels in ``obs[label_key]``.
    label_key
        obs column holding the labels; NaN labels raise ``ValueError``.
    embedding
        ``obsm`` key of the embedding; when absent, a PCA of ``X``.
    n_init
        k-means++ restarts (scikit-learn keeps the lowest inertia); scib-metrics uses 1.
    seed
        Seed of k-means (and of the PCA fallback).
    embedding_label_trained
        Say True when the embedding was trained with these labels: the result then carries a ``warning``
        (circular, as for :func:`internal_validity`).

    Returns
    -------
    dict
        ``{status, nmi, ari, n_clusters, clusters, n_used}`` (+ ``warning``); ``clusters`` is the k-means
        assignment (an int per cell, obs order). Every cell is used. Fewer than 2 labels -> a skip.

    Notes
    -----
    Mirrors scib-metrics 0.5.6 ``nmi_ari_cluster_labels_kmeans`` (BSD-3-Clause): k-means++ with ``n_clusters``
    = the number of labels, ``n_init=1``, ``max_iter=300``, ``tol=1e-4`` (relative to the mean feature
    variance), arithmetic-mean NMI. scib-metrics runs its own JAX k-means (float32, its own random stream);
    this runs scikit-learn's, so one run can land in another local optimum: comparable in kind, not to the
    digit.
    """
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

    y = check_labels(adata, label_key)
    types, code = np.unique(y, return_inverse=True)
    K = len(types)
    if K < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {K}")
    E = _embedding(adata, embedding, seed=seed)
    cl = KMeans(n_clusters=K, init="k-means++", n_init=n_init, max_iter=300, tol=1e-4,
                random_state=seed).fit_predict(E)
    out = {"status": "ok", "nmi": float(normalized_mutual_info_score(code, cl, average_method="arithmetic")),
           "ari": float(adjusted_rand_score(code, cl)), "n_clusters": int(len(np.unique(cl))),
           "clusters": np.asarray(cl), "n_used": int(len(y))}
    if embedding_label_trained:
        out["warning"] = _CIRCULAR
    return out
