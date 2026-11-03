"""Marker-program fidelity (reference-free, needs marker sets): do the labels track the curated marker
programs? Per-type AUC-ROC / Cohen's d of the marker load, and top-N DE overlap with the markers."""
from __future__ import annotations

import numpy as np

from ._utils import _dense
from ._validate import check_labels, get_counts, skip

__all__ = ["marker_program_fidelity", "marker_gene_overlap"]


# --------------------------------------------------------------------------- #
# Marker-program fidelity (reference-free, needs marker sets) - catalog §3g
# --------------------------------------------------------------------------- #
def marker_program_fidelity(adata, label_key: str, marker_sets: dict[str, list[str]], *,
                            layer: str | None = "counts") -> dict:
    """Do the labels track the marker programs? Per type, AUC-ROC + Cohen's d of that type's
    marker-load score at ranking its own labeled cells above the rest (reference-free given the
    marker sets). Reads the count matrix chosen by ``layer`` (see ``get_counts``). Returns
    ``{status, per_type: {auc, cohens_d, n}, mean_auc, mean_cohens_d, matrix}``; fewer than 2 labels -> a skip."""
    from sklearn.metrics import roc_auc_score

    var_index = {g: i for i, g in enumerate(map(str, adata.var_names))}
    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    mat, matrix = get_counts(adata, layer)
    X = _dense(mat).astype(float)
    per, aucs, ds = {}, [], []
    for t, genes in marker_sets.items():
        idx = [var_index[g] for g in genes if g in var_index]
        if not idx:
            continue
        score = X[:, idx].sum(1)                                   # marker load per cell
        y = (labels == t).astype(int)
        if y.sum() == 0 or y.sum() == len(y):
            per[t] = {"auc": None, "cohens_d": None, "n": int(y.sum())}
            continue
        try:
            auc = float(roc_auc_score(y, score))
        except Exception:
            auc = None
        pos, neg = score[y == 1], score[y == 0]
        num = (len(pos) - 1) * pos.var(ddof=1) + (len(neg) - 1) * neg.var(ddof=1)
        sd = np.sqrt(num / max(1, len(pos) + len(neg) - 2))
        d = float((pos.mean() - neg.mean()) / sd) if sd > 0 else 0.0
        per[t] = {"auc": auc, "cohens_d": d, "n": int(y.sum())}
        if auc is not None:
            aucs.append(auc)
        ds.append(d)
    return {
        "status": "ok", "per_type": per,
        "mean_auc": float(np.mean(aucs)) if aucs else None,
        "mean_cohens_d": float(np.mean(ds)) if ds else None, "matrix": matrix,
    }


def marker_gene_overlap(adata, label_key: str, marker_sets: dict[str, list[str]], n_top: int = 25, *,
                        layer: str | None = "counts") -> dict:
    """Do the labels reproduce known markers? Per type, the fraction of its curated markers that land in
    the top-``n_top`` data-derived DE genes for that label (standardized in-vs-out mean difference).
    Reference-free given marker sets; the DE-overlap complement to :func:`marker_program_fidelity`'s
    ranking view. Reads the count matrix chosen by ``layer``. Returns ``{status, per_type: {overlap,
    n_curated}, mean_overlap, matrix}``; fewer than 2 labels -> a skip."""
    genes = np.array(list(map(str, adata.var_names)))
    gene_set = set(genes)
    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    mat, matrix = get_counts(adata, layer)
    X = _dense(mat).astype(float)
    per, fracs = {}, []
    for t, curated in marker_sets.items():
        cur = [g for g in map(str, curated) if g in gene_set]
        y = labels == t
        if not cur or y.sum() < 2 or (~y).sum() < 2:
            per[t] = {"overlap": None, "n_curated": len(cur)}
            continue
        score = (X[y].mean(0) - X[~y].mean(0)) / np.sqrt(X[~y].var(0) + 1e-9)
        top = set(genes[np.argsort(-score)[:n_top]])
        ov = len(top & set(cur)) / len(cur)
        per[t] = {"overlap": float(ov), "n_curated": len(cur)}
        fracs.append(ov)
    return {"status": "ok", "per_type": per, "mean_overlap": float(np.mean(fracs)) if fracs else None,
            "matrix": matrix}
