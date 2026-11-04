"""Marker-program fidelity (reference-free, needs marker sets).

Do the labels track the curated marker programs? Per-type AUC-ROC / Cohen's d of the marker load, and
top-N DE overlap with the markers. Both report the labels that have no on-panel marker set
(``labels_without_markers``, ``no_dict_frac``): they are not scored, so an over-split or merged annotation
would otherwise look clean. ``marker_sets`` keys are compared as str, like the
labels. Sparse count matrices are read without densifying.
"""
from __future__ import annotations

import numpy as np

from ._validate import check_labels, dictionary_coverage, get_counts, skip
from .purity import _on_panel, _str_keys
from .thresholds import value as _thr

__all__ = ["marker_program_fidelity", "marker_gene_overlap"]


def _mean_var(M):
    """Column means and (ddof=0) variances of ``M``, dense or sparse, in float64."""
    from scipy import sparse

    if sparse.issparse(M):
        mu = np.asarray(M.mean(axis=0)).ravel()
        return mu, np.maximum(np.asarray(M.multiply(M).mean(axis=0)).ravel() - mu ** 2, 0.0)
    return M.mean(0), M.var(0)


# --------------------------------------------------------------------------- #
# Marker-program fidelity (reference-free, needs marker sets) - catalog §3g
# --------------------------------------------------------------------------- #
def marker_program_fidelity(adata, label_key: str, marker_sets: dict[str, list[str]], *,
                            layer: str | None = "counts") -> dict:
    """Do the labels track the marker programs? Per type, AUC-ROC + Cohen's d of its marker load.

    The marker load (sum of the type's on-panel marker counts per cell) should rank the type's own labeled
    cells above the rest. Reference-free given the marker sets.

    Parameters
    ----------
    adata
        Section with raw counts and labels.
    label_key
        ``obs`` column with the labels.
    marker_sets
        ``{label: [genes]}``; keys compared as str.
    layer
        Count matrix (see ``_validate.get_counts``).

    Returns
    -------
    dict
        ``{status, per_type: {auc, cohens_d, n}, mean_auc, mean_cohens_d, no_dict_frac,
        labels_without_markers, matrix}``. ``cohens_d`` is None when the marker load is constant;
        ``labels_without_markers`` = labels with no on-panel marker set (not scored), ``no_dict_frac`` their
        share of cells. Fewer than 2 labels,
        or no type scorable, -> a skip.
    """
    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    from sklearn.metrics import roc_auc_score

    ms = _on_panel(adata, _str_keys(marker_sets))
    cov = dictionary_coverage(labels, ms)
    var_index = {g: i for i, g in enumerate(map(str, adata.var_names))}
    mat, matrix = get_counts(adata, layer)
    per, aucs, ds = {}, [], []
    for t, genes in ms.items():
        idx = [var_index[g] for g in genes]
        if not idx:
            continue
        score = np.asarray(mat[:, idx].astype(np.float64).sum(axis=1)).ravel()   # marker load per cell
        y = (labels == t).astype(int)
        if y.sum() == 0 or y.sum() == len(y):
            per[t] = {"auc": None, "cohens_d": None, "n": int(y.sum())}
            continue
        try:
            auc = float(roc_auc_score(y, score))
        except Exception:
            auc = None
        pos, neg = score[y == 1], score[y == 0]
        dof = len(pos) + len(neg) - 2
        num = len(pos) * pos.var(ddof=0) + len(neg) * neg.var(ddof=0)
        sd = np.sqrt(num / dof) if dof > 0 else 0.0
        d = float((pos.mean() - neg.mean()) / sd) if sd > 0 else None     # constant load: undefined
        per[t] = {"auc": auc, "cohens_d": d, "n": int(y.sum())}
        if auc is not None:
            aucs.append(auc)
        if d is not None:
            ds.append(d)
    if not aucs and not ds:
        return {**skip("no labelled type has a marker on the panel and cells on both sides"), **cov}
    return {
        "status": "ok", "per_type": per,
        "mean_auc": float(np.mean(aucs)) if aucs else None,
        "mean_cohens_d": float(np.mean(ds)) if ds else None, **cov, "matrix": matrix,
    }


def marker_gene_overlap(adata, label_key: str, marker_sets: dict[str, list[str]],
                        n_top: int = _thr("fidelity", "overlap_n_top"), *,
                        layer: str | None = "counts") -> dict:
    """Do the labels reproduce known markers? Per type, the share of its markers among its top DE genes.

    The top-``n_top`` data-derived DE genes per label are ranked by the standardized in-vs-out mean
    difference.
    Reference-free given marker sets; the DE-overlap complement to :func:`marker_program_fidelity`'s ranking
    view.

    Parameters
    ----------
    adata
        Section with raw counts and labels.
    label_key
        ``obs`` column with the labels.
    marker_sets
        ``{label: [genes]}``; keys compared as str.
    n_top
        DE genes per label.
    layer
        Count matrix (see ``_validate.get_counts``).

    Returns
    -------
    dict
        ``{status, per_type: {overlap, n_curated}, mean_overlap, no_dict_frac, labels_without_markers,
        matrix}``; ``labels_without_markers`` = labels with no on-panel marker set, ``no_dict_frac`` their
        share of cells. Fewer than 2 labels, or no type scorable, -> a skip.
    """
    from scipy import sparse

    genes = np.array(list(map(str, adata.var_names)))
    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    ms = _on_panel(adata, _str_keys(marker_sets))
    cov = dictionary_coverage(labels, ms)
    mat, matrix = get_counts(adata, layer)
    X = mat.tocsr().astype(np.float64) if sparse.issparse(mat) else np.asarray(mat).astype(float)
    per, fracs = {}, []
    for t, cur in ms.items():
        y = labels == t
        if not cur or y.sum() < 2 or (~y).sum() < 2:
            per[t] = {"overlap": None, "n_curated": len(cur)}
            continue
        (mu_in, _), (mu_out, var_out) = _mean_var(X[y]), _mean_var(X[~y])
        score = (mu_in - mu_out) / np.sqrt(var_out + 1e-9)
        top = set(genes[np.argsort(-score)[:n_top]])
        ov = len(top & set(cur)) / len(cur)
        per[t] = {"overlap": float(ov), "n_curated": len(cur)}
        fracs.append(ov)
    if not fracs:
        return {**skip("no labelled type has a marker on the panel and >= 2 cells on both sides"), **cov}
    return {"status": "ok", "per_type": per, "mean_overlap": float(np.mean(fracs)) if fracs else None, **cov,
            "matrix": matrix}
