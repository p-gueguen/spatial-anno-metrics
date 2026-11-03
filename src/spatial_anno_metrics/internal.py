"""Reference-free internal validity of a labeling on an embedding (the scTypeEval family): silhouette,
2-label silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch, average similarity,
Davies-Bouldin, Calinski-Harabasz and a geometric-mean ``integrated`` score."""
from __future__ import annotations

import numpy as np

from ._utils import _embedding, _scale01

__all__ = ["internal_validity"]


# --------------------------------------------------------------------------- #
# Internal cluster-validity (scTypeEval family) - reference-free
# --------------------------------------------------------------------------- #
def internal_validity(adata, label_key: str = "cell_type", embedding: str | None = "X_pca",
                      k: int = 30, subsample: int = 2000, seed: int = 0) -> dict:
    """Reference-free cluster-validity of a labeling on an embedding (scTypeEval family).

    Returns ``{silhouette, silhouette_2label, neighborhood_purity, orbital_medoid,
    ward_propmatch, avg_similarity, integrated, n_used}``. ``integrated`` is the geometric mean
    of all six (each scaled to [0,1]) - a single "is this labeling internally coherent?" number.
    O(N^2) parts (silhouettes) run on a ``subsample`` of cells for scalability; the rest use all.
    """
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import pairwise_distances, silhouette_score
    from sklearn.neighbors import NearestNeighbors

    E = _embedding(adata, embedding, seed=seed)
    y = adata.obs[label_key].astype(str).to_numpy()
    types = np.unique(y)
    n = len(y)
    if len(types) < 2:
        return {"status": "skipped: <2 cell types"}

    # neighborhood purity - fraction of k neighbors (excluding self) sharing the label; all cells.
    kk = int(min(k, n - 1))
    nn = NearestNeighbors(n_neighbors=kk + 1).fit(E)
    idx = nn.kneighbors(E, return_distance=False)[:, 1:]           # drop self
    purity = float(np.mean([(y[nbrs] == y[i]).mean() for i, nbrs in enumerate(idx)]))

    # orbital-medoid - fraction of cells whose nearest type-medoid is their own type.
    medoids = {}
    for t in types:
        pts = E[y == t]
        c = pts.mean(0)
        medoids[t] = pts[np.argmin(((pts - c) ** 2).sum(1))]      # actual cell nearest the centroid
    Mnames = list(types)
    M = np.vstack([medoids[t] for t in Mnames])
    d2 = ((E[:, None, :] - M[None, :, :]) ** 2).sum(2) if n * len(Mnames) < 5_000_000 else None
    if d2 is None:
        nearest = np.array([Mnames[np.argmin(((E[i] - M) ** 2).sum(1))] for i in range(n)])
    else:
        nearest = np.array(Mnames)[d2.argmin(1)]
    orbital = float((nearest == y).mean())

    # Ward-PropMatch - agreement between labels and an unsupervised Ward clustering.
    ward = AgglomerativeClustering(n_clusters=len(types), linkage="ward").fit_predict(E)
    prop = []
    for t in types:
        w = ward[y == t]
        prop.append(np.bincount(w, minlength=len(types)).max() / max(1, len(w)))
    ward_prop = float(np.mean(prop))

    # average similarity - mean cosine of each cell to its own-type centroid (cohesion).
    cents = {t: E[y == t].mean(0) for t in types}
    def _cos(u, v):
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else 0.0
    avg_sim = float(np.mean([_cos(E[i], cents[y[i]]) for i in range(n)]))

    # silhouettes - subsample the O(N^2) distance work.
    rng = np.random.default_rng(seed)
    sel = np.arange(n) if n <= subsample else rng.choice(n, subsample, replace=False)
    Es, ys = E[sel], y[sel]
    sil = float(silhouette_score(Es, ys)) if len(np.unique(ys)) > 1 else 0.0
    D = pairwise_distances(Es)
    np.fill_diagonal(D, np.nan)
    s2 = []
    for t in np.unique(ys):
        own = ys == t
        if own.sum() < 2 or (~own).sum() < 1:
            continue
        a_in = np.nanmean(D[np.ix_(own, own)], axis=1)             # mean within-type distance
        b_out = D[np.ix_(own, ~own)].mean(1)                       # mean distance to all others
        s2.extend((b_out - a_in) / np.maximum(a_in, b_out))
    sil_2label = float(np.mean(s2)) if s2 else 0.0

    # Davies-Bouldin (lower = better separated) + Calinski-Harabasz (higher = better): standard sklearn
    # internal indices, reported RAW (not folded into `integrated`) - the catalog battery-C over-cluster check.
    from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
    yenc = np.unique(ys, return_inverse=True)[1]
    db = float(davies_bouldin_score(Es, yenc)) if len(np.unique(ys)) > 1 else float("nan")
    ch = float(calinski_harabasz_score(Es, yenc)) if len(np.unique(ys)) > 1 else float("nan")

    parts = [_scale01(sil), _scale01(sil_2label), purity, orbital, ward_prop, _scale01(avg_sim)]
    integrated = float(np.exp(np.mean(np.log(np.clip(parts, 1e-9, 1.0)))))
    return {
        "silhouette": sil, "silhouette_2label": sil_2label, "neighborhood_purity": purity,
        "orbital_medoid": orbital, "ward_propmatch": ward_prop, "avg_similarity": avg_sim,
        "davies_bouldin": db, "calinski_harabasz": ch,
        "integrated": integrated, "n_used": int(len(sel)),
    }
