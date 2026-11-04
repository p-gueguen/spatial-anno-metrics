"""Inter-sample consistency: is each cell type reproducible across biological samples.

Reference-free, needs biological replicates. The ISCG / ISCL axes and per-type refinement regimes of scTypeEval
(Garnica, Andreatta & Carmona 2026).
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from ._validate import _all_finite, check_labels, get_counts, skip
from .thresholds import value as _thr

__all__ = ["inter_sample_consistency", "isc_regime"]


# --------------------------------------------------------------------------- #
# Inter-sample consistency (ISC) - reference-free, needs biological replicates
# --------------------------------------------------------------------------- #
def _pseudobulk_profiles(adata, label_key, sample_key, *, min_cells, min_samples, layer=None):
    """Build one profile per (sample, type) with >= ``min_cells`` cells.

    A type is kept only when it has >= ``min_samples`` such profiles (the rest are returned in ``dropped``). With a count ``layer`` (read by
    :func:`get_counts`) the paper's recipe is used (sum raw counts, then log1p of counts-per-10k); with
    ``layer=None`` ``X`` is averaged as it is. Also returns the name of the matrix read. Profiles are ordered by
    type, then sample; the sums are ``indicator.T @ X``, so a sparse matrix is never densified.
    """
    y = check_labels(adata, label_key)
    s = check_labels(adata, sample_key)
    if layer:
        mat, matrix = get_counts(adata, layer)
    else:
        if adata.X is None or not _all_finite(adata.X):
            raise ValueError("adata.X contains NaN or inf")
        mat, matrix = adata.X, "X"
    types, ti = np.unique(y, return_inverse=True)
    samps, si = np.unique(s, return_inverse=True)
    group = ti * len(samps) + si                                     # type-major, as the profiles are ordered
    n_groups = len(types) * len(samps)
    n_cells = np.bincount(group, minlength=n_groups)
    ok = (n_cells >= max(1, min_cells)).reshape(len(types), len(samps))
    keep_type = ok.sum(1) >= min_samples
    dropped = {t: int(k) for t, k, kt in zip(types, ok.sum(1), keep_type) if not kt}
    sel = np.flatnonzero((ok & keep_type[:, None]).ravel())
    ind = sp.csr_matrix((np.ones(len(group)), (np.arange(len(group)), group)), shape=(len(group), n_groups))
    S = ind[:, sel].T @ (mat if sp.issparse(mat) else np.asarray(mat))
    P = (S.toarray() if sp.issparse(S) else np.asarray(S)).astype(float)
    if not layer:
        P /= n_cells[sel][:, None]
    elif P.size:
        P = np.log1p(P / np.maximum(P.sum(1, keepdims=True), 1.0) * _SCALE)
    return P, types[sel // len(samps)], samps[sel % len(samps)], dropped, matrix


def _two_label_silhouette(D, labels):
    """Compute scTypeEval's ``2label_silhouette``.

    a(i) = mean distance to own label (excl. self); b(i) = mean distance to ALL other-label observations pooled
    (not the nearest other cluster).
    """
    n = len(labels)
    D_clean = np.nan_to_num(D, nan=1.0)
    same = labels[:, None] == labels[None, :]
    np.fill_diagonal(same, False)
    other = labels[:, None] != labels[None, :]
    n_same, n_other = same.sum(1), other.sum(1)
    a = np.divide(np.where(same, D_clean, 0.0).sum(1), n_same, out=np.zeros(n), where=n_same > 0)
    b = np.divide(np.where(other, D_clean, 0.0).sum(1), n_other, out=np.zeros(n), where=n_other > 0)
    denom = np.maximum(a, b)
    out = np.divide(b - a, denom, out=np.zeros(n), where=denom > 0)
    out[n_same == 0] = 0.0
    return out


def _classic_marker_n(n_labels):
    """SingleR classic ``de.n`` default: round(500 * (2/3)^log2(N)) genes per label pair."""
    return int(round(_MARKER_N * (2 / 3) ** np.log2(n_labels)))


def _classic_markers(Pk, n=None):
    """Return the union over ordered label pairs (a, b) of the top ``n`` genes by Pk[a] - Pk[b] (positive only).

    Re-implements the marker rule SingleR documents for ``de.method = "classic"`` (Aran et al. 2019; one profile
    per label here, so the median difference is the difference). No SingleR code is copied (GPL-3).
    """
    n = _classic_marker_n(len(Pk)) if n is None else n
    genes = set()
    for a in range(len(Pk)):
        for b in range(len(Pk)):
            if a != b:
                d = Pk[a] - Pk[b]
                top = np.argsort(-d, kind="stable")[:n]
                genes.update(top[d[top] > 0].tolist())
    return np.array(sorted(genes), dtype=int)


def _spearman(A, B):
    """Spearman correlation of every row of A with every row of B (average ranks for ties)."""
    from scipy.stats import rankdata

    def z(M):
        R = rankdata(M, axis=1).astype(float)
        R -= R.mean(1, keepdims=True)
        return R / (np.linalg.norm(R, axis=1, keepdims=True) + 1e-12)

    return z(A) @ z(B).T


def _rcm_dissimilarity(P, samples, labels, genes=None):
    """Reciprocal-classification match (RCM) distance between (sample, type) pseudobulks.

    Classifier = nearest profile by Spearman correlation (a stand-in for the paper's SingleR classic mode, without
    its fine-tuning and pruning), over ``genes``: None = every column of ``P``, an index array = those columns,
    ``"markers"`` = SingleR's classic marker union of the training sample (recomputed per sample). As scTypeEval's
    ``recip_classif`` + ``match_pairs``: D[i, j] = 0 when i's best hit inside j's sample is j AND j's best hit
    inside i's sample is i and both labels occur in both samples; 0.5 for such a mutual hit when either label is
    absent from the other sample; 1 otherwise, and 1 for every pair touching a sample with < 2 profiles. A mutual
    cross-type swap also scores 0 - that is how RCM penalises confusable types (it lowers b(i) in the silhouette).
    """
    n = len(samples)
    uniq, sidx = np.unique(samples, return_inverse=True)
    lab, lidx = np.unique(labels, return_inverse=True)
    C = None if isinstance(genes, str) else _spearman(*(2 * [P if genes is None else P[:, genes]]))
    best = np.empty((n, len(uniq)), dtype=int)  # best hit of observation i inside each sample
    for k in range(len(uniq)):
        cols = np.flatnonzero(sidx == k)
        if len(cols) < 2:
            best[:, k] = cols[0] if len(cols) else 0
            continue
        if C is None:
            m = _classic_markers(P[cols])
            Ck = _spearman(P[:, m], P[cols][:, m]) if len(m) >= 2 else _spearman(P, P[cols])
        else:
            Ck = C[:, cols]
        best[:, k] = cols[np.argmax(Ck, axis=1)]
    hit = best[:, sidx]  # hit[i, j] = best hit of i inside j's sample
    idx = np.arange(n)
    match = (hit == idx[None, :]) & (hit.T == idx[:, None])
    present = np.zeros((len(lab), len(uniq)), dtype=bool)
    present[lidx, sidx] = True
    shared = present[lidx][:, sidx]                         # shared[i, j]: i's label occurs in j's sample
    one = np.bincount(sidx)[sidx] < 2
    D = np.where(match, np.where(shared & shared.T, 0.0, 0.5), 1.0)
    D[one[:, None] | one[None, :]] = 1.0
    np.fill_diagonal(D, 0.0)
    return D


# Garnica, Andreatta & Carmona 2026 (Fig. 5c): median - 1 IQR over a 9-atlas meta-study. Guides, not cut-offs.
ISC_THRESHOLDS = {"iscg": _thr("isc", "iscg"), "iscl": _thr("isc", "iscl")}
_SCALE = _thr("isc", "pseudobulk_scale")
_MARKER_N = _thr("isc", "rcm_marker_n")


def isc_regime(iscg, iscl, thr_g: float = ISC_THRESHOLDS["iscg"], thr_l: float = ISC_THRESHOLDS["iscl"]) -> str:
    """Map a label's (ISCG, ISCL) to the Garnica 2026 refinement quadrant.

    A helper, not a metric: it returns the quadrant name as a ``str`` (not a result dict). Experimental: on the
    synthetic sweep only the over-split case reads right (see :func:`inter_sample_consistency`, Notes).

    Parameters
    ----------
    iscg, iscl : float or None
        The label's global and local inter-sample consistency.
    thr_g, thr_l : float
        Guide thresholds (Garnica 2026 Fig. 5c: 0.74 and 0.45).

    Returns
    -------
    str
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
                             min_cells: int = _thr("isc", "min_cells"), min_samples: int = _thr("isc", "min_samples"),
                             layer: str | None = None, thresholds: dict | None = None,
                             reduction: str | None = "pca", n_pcs: int = _thr("isc", "n_pcs"),
                             rcm_genes=None) -> dict:
    """Score how reproducible each cell type is across biological samples.

    scTypeEval inter-sample consistency (Garnica, Andreatta & Carmona 2026, bioRxiv 10.64898/2026.06.10.731285),
    checked against scTypeEval 1.0.0 on a Yerly 2022 subsample (``tests/goldens/isc_yerly2022``). One pseudobulk
    per (sample, type) with >= ``min_cells`` cells; a type needs >= ``min_samples`` of them or it is listed under
    ``dropped_types`` instead of being scored (the paper's defaults, 10 / 5; 3 % of atlas labels fell below them).

    Parameters
    ----------
    adata : AnnData
        Cells x genes.
    label_key, sample_key : str
        ``adata.obs`` columns with the labels being judged and the biological sample (replicate).
    metric : str
        Distance for the legacy ``consistency`` / ``per_type`` silhouette (any ``sklearn`` metric).
    min_cells, min_samples : int
        Cells per (sample, type) pseudobulk; pseudobulks a type needs to be scored.
    layer : str or None
        Raw-count layer: sum counts, then log1p of counts-per-10k (the paper's recipe; a missing layer falls back
        to ``X`` only when ``X`` holds counts). ``None`` (default) averages ``X`` as it is, so hand it
        log-normalised HVGs. Aggregation is ``indicator.T @ X``: a sparse matrix is never densified.
    thresholds : dict, optional
        Overrides ``{"iscg": 0.74, "iscl": 0.45}`` for the regime calls.
    reduction : {"pca", None}
        Space of the ISCG cosine. ``"pca"`` (default, scTypeEval's ``reduction = TRUE``): the first ``n_pcs``
        principal components of the centred, unscaled pseudobulk matrix (exact SVD; scTypeEval uses irlba, which
        differs by <= 1.4e-6 on the golden). ``None``: the pseudobulks themselves (scTypeEval
        ``reduction = FALSE``; the pre-0.7.0 space, which reads 0.012-0.073 lower on the golden, mean 0.03).
    n_pcs : int
        Components kept, capped at ``min(n_profiles, n_genes) - 1`` as scTypeEval does.
    rcm_genes : None, "markers" or sequence of str
        Genes of the reciprocal classifier. ``None`` (default): every gene. ``"markers"``: SingleR's classic
        marker union of each training sample (the paper's classifier, without fine-tuning). A list of var names
        (e.g. HVGs): those genes. ISCG always uses every gene.

    Returns
    -------
    dict
        ``status``; ``iscg`` - GLOBAL: 2-label silhouette on pseudobulk cosine distance, mean over labels
        (compactness against all other types pooled; the most sensitive detector of signal degradation);
        ``iscl`` - LOCAL: silhouette on the reciprocal-classification match distance, mean over labels
        (replicability; the detector of over-partitioning, least batch-sensitive); ``integrated`` - mean over
        labels of ISCG x ISCL (the paper's dataset score; R^2 > 0.98 with macro-F1); ``per_type_iscg`` /
        ``per_type_iscl`` / ``per_type_integrated``; ``iscg_profile_mean`` / ``iscl_profile_mean`` - the same
        silhouettes averaged over pseudobulks (a label weighs by its number of samples; the pre-0.7.0
        headline); ``regime`` per type via :func:`isc_regime` and ``regime_status`` = ``"experimental"``;
        ``consistency`` / ``per_type`` - silhouette on pseudobulk ``metric`` distance in the raw space (the
        pre-0.5 headline); ``thresholds``, ``reduction``, ``n_pcs`` (None without reduction), ``rcm_genes``
        (``"all"``, ``"markers"`` or the number of genes), ``n_profiles``, ``dropped_types``, ``matrix``.
        Fewer than 2 scorable types -> a skip that keeps ``n_profiles`` / ``dropped_types``.

    Notes
    -----
    Regime calls are experimental (``regime_status``). On the synthetic sweep (``benchmarks/synthetic``,
    ``check_isc_regimes.csv``, 5 seeds, 0.7.0 defaults) only over-splitting reads right (Epi halves: ``merge``
    10 / 10). Merged labels read ``consistent`` (Fibro+Epi: ISCG 0.985, ISCL 1.000; three merges: 15 / 15
    ``consistent``), 20 % random relabelling reads ``split`` for 21 / 30 labels (ISCG 0.688, ISCL 0.980) and a
    per-sample batch effect (log-sd 1.0) ``split`` for 30 / 30 (ISCG 0.617, ISCL 1.000): ISCG falls under 0.74
    while ISCL stays high (0.84 at 40 % noise). The 0.74 / 0.45 guides rest on one preprint's meta-study.

    ISCL on a toy fixture of flat, low-noise background genes is not informative: with Spearman ranks a single
    foreign cell in a pseudobulk lifts its type's genes above the flat background, so two pseudobulks that hold
    one cell of each other's type become mutual best hits (8 types x 64 genes: ISCL 0.92 -> 0.14 at 10 %
    relabelling). With a graded per-gene baseline, as real genes have, or with Pearson, ISCL stays 1.00; on real
    data (golden) our classifier makes the same mutual hits as SingleR on 1,128 / 1,128 pseudobulk pairs.
    """
    from scipy.spatial.distance import cdist
    from sklearn.metrics import silhouette_samples

    if reduction not in ("pca", None):
        raise ValueError(f"reduction must be 'pca' or None, got {reduction!r}")
    P, plabels, psamples, dropped, matrix = _pseudobulk_profiles(
        adata, label_key, sample_key, min_cells=min_cells, min_samples=min_samples, layer=layer)
    genes = rcm_genes
    if rcm_genes is not None and not isinstance(rcm_genes, str):
        rcm_list = list(rcm_genes)
        if len(rcm_list) < 2:
            raise ValueError(f"rcm_genes must contain >= 2 genes, got {len(rcm_list)}")
        genes = adata.var_names.get_indexer(rcm_list)
        if (genes < 0).any():
            raise KeyError(f"rcm_genes not in adata.var_names: {[g for g, i in zip(rcm_list, genes) if i < 0][:5]}")
    elif rcm_genes not in (None, "markers"):
        raise ValueError(f"rcm_genes must be None, 'markers' or a list of genes, got {rcm_genes!r}")
    thr = {**ISC_THRESHOLDS, **(thresholds or {})}
    types = np.unique(plabels)
    if len(types) < 2 or len(plabels) <= len(types):
        return {**skip(f"need >= 2 types with >= {min_samples} samples of >= {min_cells} cells"),
                "consistency": None, "iscg": None, "iscl": None, "integrated": None,
                "n_profiles": int(len(plabels)), "dropped_types": dropped}

    def _per_type(v):
        return {t: float(v[plabels == t].mean()) for t in types}

    sil = silhouette_samples(P, plabels, metric=metric)
    E, k = P, None
    if reduction == "pca":
        k = int(min(n_pcs, max(1, min(P.shape) - 1)))
        U, S, _ = np.linalg.svd(P - P.mean(0), full_matrices=False)
        E = U[:, :k] * S[:k]
    D_cos = np.nan_to_num(cdist(E, E, metric="cosine"), nan=1.0)
    np.fill_diagonal(D_cos, 0.0)
    glob = _two_label_silhouette(D_cos, plabels)
    D_rcm = _rcm_dissimilarity(P, psamples, plabels, genes)
    D_rcm[D_rcm == 0.5] = 1.0                                      # scTypeEval get_consistency reads 0.5 as 1
    loc = silhouette_samples(D_rcm, plabels, metric="precomputed")
    per_g, per_l = _per_type(glob), _per_type(loc)
    per_int = {t: per_g[t] * per_l[t] for t in types}
    return {
        "status": "ok", "consistency": float(sil.mean()), "per_type": _per_type(sil),
        "iscg": float(np.mean(list(per_g.values()))), "iscl": float(np.mean(list(per_l.values()))),
        "integrated": float(np.mean(list(per_int.values()))),
        "per_type_iscg": per_g, "per_type_iscl": per_l, "per_type_integrated": per_int,
        "iscg_profile_mean": float(glob.mean()), "iscl_profile_mean": float(loc.mean()),
        "regime": {t: isc_regime(per_g[t], per_l[t], thr["iscg"], thr["iscl"]) for t in types},
        "regime_status": "experimental",
        "thresholds": thr, "reduction": reduction, "n_pcs": k,
        "rcm_genes": "all" if rcm_genes is None else rcm_genes if isinstance(rcm_genes, str) else int(len(genes)),
        "n_profiles": int(len(plabels)), "dropped_types": dropped, "matrix": matrix,
    }
