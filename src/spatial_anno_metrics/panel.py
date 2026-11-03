"""Panel adequacy: can a target gene panel resolve a reference's cell types at the panel's own depth?"""
from __future__ import annotations

import numpy as np

from ._utils import _dense
from ._validate import check_labels, skip
from .external import external_scores
from .thresholds import value as _thr

__all__ = ["panel_resolvability"]
_MIN_SHARED, _MIN_CELLS = _thr("panel", "min_shared_genes"), _thr("panel", "min_cells_per_type")
_F1_OK, _F1_WEAK = _thr("panel", "f1_resolvable"), _thr("panel", "f1_weak")
_CONF_WEAK, _CONF_NOT = _thr("panel", "confused_weak"), _thr("panel", "confused_not_resolvable")
_C, _MAX_ITER, _SCALE = _thr("panel", "logreg_C"), _thr("panel", "logreg_max_iter"), _thr("panel", "library_size_scale")


# --------------------------------------------------------------------------- #
# Panel adequacy - can a target panel resolve a reference's cell types? (§4 battery A #1)
# --------------------------------------------------------------------------- #
def panel_resolvability(reference, label_key: str, panel_genes, target_depth: float = _thr("panel", "target_depth"),
                        *, max_cells_per_type: int = _thr("panel", "max_cells_per_type"),
                        cv: int = _thr("panel", "cv_folds"), seed: int = 0) -> dict:
    """Can a target panel resolve this reference's cell types AT THE PANEL'S OWN DEPTH?

    Restrict the reference to the genes it shares with ``panel_genes``, **binomially thin each cell
    to ``target_depth`` counts** (the section's real per-cell depth), cross-validate a classifier, and
    score per-type with F1 (via :func:`external_scores`). Low per-type F1 + a dominant confuser means
    the panel cannot separate that type from a look-alike -> the reference's labels are too granular
    for the probe set (merge them).

    The thinning is load-bearing: at full scRNA depth every type separates (one-vs-rest AUC > 0.98
    even for 25 kidney types on a breast panel) - so this **supersedes**
    ``panel_check.identifiability_auc`` (one-vs-rest AUC at full depth, which inflates to ~1.0). Uses
    multiclass per-class F1, not one-vs-rest AUC, because AUC hides sibling confusion.

    IMPORTANT: this measures *separability given the panel*, not tissue relevance - it is only
    meaningful for a **tissue-matched** reference x panel pair (a wrong-tissue reference can score
    high here and still annotate nonsense).

    Returns ``{status, n_shared_genes, n_types, macro_f1, balanced_accuracy, ecs, frac_resolvable,
    frac_not_resolvable, per_type: {ct: {f1, recall, precision, confused_with, confused_frac, tier}}}``.
    Fewer than 25 shared genes or 2 types with >= 10 cells -> a skip (with ``n_shared_genes``, ``per_type``).
    """
    import scipy.sparse as sp
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import precision_recall_fscore_support
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    panel = {str(g) for g in panel_genes}
    shared = [g for g in reference.var_names if str(g) in panel]
    if len(shared) < _MIN_SHARED:
        return {**skip(f"{len(shared)} genes shared with the panel, need >= {_MIN_SHARED}"),
                "n_shared_genes": len(shared), "per_type": {}}
    sub = reference[:, shared].copy()
    xs = sub.X[:200].toarray() if sp.issparse(sub.X) else np.asarray(sub.X[:200])
    raw = bool(np.allclose(xs, np.round(xs)) and xs.max() > 1)
    if not raw and "counts" in sub.layers:
        sub.X = sub.layers["counts"]
        raw = True
    rng = np.random.default_rng(seed)
    y = check_labels(sub, label_key)
    idx: list = []
    for ct in np.unique(y):
        pos = np.where(y == ct)[0]
        idx.extend(pos if len(pos) <= max_cells_per_type
                   else rng.choice(pos, max_cells_per_type, replace=False))
    sub = sub[np.array(sorted(idx))].copy()
    y = check_labels(sub, label_key)
    keep = {ct for ct, n in zip(*np.unique(y, return_counts=True)) if n >= _MIN_CELLS}
    m = np.isin(y, list(keep))
    sub, y = sub[m].copy(), y[m]
    if len(np.unique(y)) < 2:
        return {**skip(f"fewer than 2 reference types with >= {_MIN_CELLS} cells"),
                "n_shared_genes": len(shared), "per_type": {}}
    if raw:                                                        # thin to the panel's real depth
        A = (sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)).astype(np.int64)
        tot = A.sum(1, keepdims=True)
        p = np.clip(target_depth / np.maximum(tot, 1), 0.0, 1.0)
        sub.X = sp.csr_matrix(rng.binomial(A, np.broadcast_to(p, A.shape)).astype("float32"))
    Xn = _dense(sub.X).astype(float)                             # library-size normalize + log1p (no scanpy dep)
    lib = Xn.sum(1, keepdims=True)
    lib[lib == 0] = 1.0
    X = np.log1p(Xn / lib * _SCALE)
    classes = np.unique(y)
    proba = cross_val_predict(LogisticRegression(max_iter=_MAX_ITER, C=_C), X, y,
                              cv=StratifiedKFold(cv, shuffle=True, random_state=seed),
                              method="predict_proba")
    pred = classes[np.argmax(proba, 1)]
    ext = external_scores(y_true=y, y_pred=pred)
    prec, rec, f1, _ = precision_recall_fscore_support(y, pred, labels=classes, zero_division=0)
    per = {}
    for i, ct in enumerate(classes):
        mask = y == ct
        wrong = pred[mask][pred[mask] != ct]
        cw, cf = (None, 0.0)
        if len(wrong):
            v, c = np.unique(wrong, return_counts=True)
            j = int(c.argmax())
            cw, cf = str(v[j]), float(c[j] / mask.sum())
        tier = ("not_resolvable" if f1[i] < _F1_WEAK or cf >= _CONF_NOT
                else "weak" if f1[i] < _F1_OK or cf >= _CONF_WEAK else "resolvable")
        per[str(ct)] = {"f1": float(f1[i]), "recall": float(rec[i]), "precision": float(prec[i]),
                        "confused_with": cw, "confused_frac": round(cf, 2), "tier": tier}
    tiers = [d["tier"] for d in per.values()]
    return {"status": "ok", "n_shared_genes": len(shared), "n_types": int(len(classes)),
            "macro_f1": ext["macro_f1"], "balanced_accuracy": ext["balanced_accuracy"], "ecs": ext["ecs"],
            "frac_resolvable": float(np.mean([t == "resolvable" for t in tiers])),
            "frac_not_resolvable": float(np.mean([t == "not_resolvable" for t in tiers])),
            "per_type": per}
