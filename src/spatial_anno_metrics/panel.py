"""Panel adequacy: can a target gene panel resolve a reference's cell types at the panel's own depth."""
from __future__ import annotations

import numpy as np

from ._utils import _dense
from ._validate import check_labels, get_counts, skip
from .external import external_scores
from .thresholds import value as _thr

__all__ = ["panel_resolvability"]
_MIN_SHARED, _MIN_CELLS = _thr("panel", "min_shared_genes"), _thr("panel", "min_cells_per_type")
_F1_OK, _F1_WEAK = _thr("panel", "f1_adjusted_resolvable"), _thr("panel", "f1_adjusted_weak")
_CONF_WEAK, _CONF_NOT = _thr("panel", "confused_weak"), _thr("panel", "confused_not_resolvable")
_C, _MAX_ITER, _SCALE = _thr("panel", "logreg_C"), _thr("panel", "logreg_max_iter"), _thr("panel", "library_size_scale")


# --------------------------------------------------------------------------- #
# Panel adequacy - can a target panel resolve a reference's cell types? (§4 battery A #1)
# --------------------------------------------------------------------------- #
def panel_resolvability(reference, label_key: str, panel_genes, target_depth: float = _thr("panel", "target_depth"),
                        *, max_cells_per_type: int = _thr("panel", "max_cells_per_type"),
                        cv: int = _thr("panel", "cv_folds"), seed: int = 0, layer: str | None = "counts",
                        f1_resolvable: float = _F1_OK, f1_weak: float = _F1_WEAK,
                        null_adjusted: bool = True) -> dict:
    """Score whether a target panel resolves this reference's cell types AT THE PANEL'S OWN DEPTH.

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

    The counts to thin are the matrix chosen by ``layer`` (see ``get_counts``): a log-normalised ``X``
    with no counts layer raises instead of being scored at full depth; a non-count matrix is rounded to
    counts, with one warning.

    Tiers are set on the **chance-corrected** F1, ``f1_adjusted = (f1 - f1_null) / (1 - f1_null)``, where
    ``f1_null`` is the type's share of the scored cells: the expected F1 of a classifier that ignores the
    expression and calls the type at its base rate (a label-permutation null). Raw F1 0.5 is chance for two
    balanced types and 3x chance for six, so a fixed raw cut cannot mean the same thing for both. A type is
    ``resolvable`` when ``f1_adjusted >= f1_resolvable`` and no single other type takes ``>= 25 %`` of its
    cells, ``not_resolvable`` when ``f1_adjusted < f1_weak`` or one other type takes ``>= 30 %``, else ``weak``.

    Parameters
    ----------
    reference : AnnData
        Labelled reference (scRNA-seq) with counts in ``layer`` (or count-like ``X``).
    label_key : str
        ``reference.obs`` column with the cell types.
    panel_genes : iterable of str
        The target panel's genes.
    target_depth : float
        Counts per cell to thin the reference to (the section's median depth).
    max_cells_per_type : int
        Reference cells sampled per type.
    cv : int
        Cross-validation folds.
    seed : int
        Seed for sampling, thinning and the folds.
    layer : str or None
        Count matrix, see ``get_counts``.
    f1_resolvable, f1_weak : float
        Tier cut-offs on the chance-corrected F1 (0.8, 0.5). Before 0.7 they were 0.5 and 0.3 on the raw F1,
        which called a type 'resolvable' when a quarter of its cells were misassigned (docs/ROADMAP.md 8h);
        ``f1_resolvable=0.5, f1_weak=0.3, null_adjusted=False`` gives those tiers back.
    null_adjusted : bool
        Compare the cut-offs with ``f1_adjusted`` (True) or with the raw ``f1`` (False).

    Returns
    -------
    dict
        ``{status, n_shared_genes, n_types, macro_f1, balanced_accuracy, ecs, frac_resolvable,
        frac_not_resolvable, per_type: {ct: {f1, f1_null, f1_adjusted, recall, precision, confused_with,
        confused_frac, tier}}, matrix}``. Fewer than 25 shared genes or 2 types with >= 10 cells -> a skip (with
        ``n_shared_genes``, ``per_type``).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import precision_recall_fscore_support
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    if cv < 2:
        raise ValueError(f"cv must be >= 2, got {cv}")
    if max_cells_per_type < cv:
        raise ValueError(f"max_cells_per_type ({max_cells_per_type}) must be >= cv ({cv})")
    if not np.isfinite(target_depth) or target_depth < 0:
        raise ValueError(f"target_depth must be finite and >= 0, got {target_depth}")
    panel = {str(g) for g in panel_genes}
    shared = list(dict.fromkeys(str(g) for g in reference.var_names if str(g) in panel))
    if len(shared) < _MIN_SHARED:
        return {**skip(f"{len(shared)} genes shared with the panel, need >= {_MIN_SHARED}"),
                "n_shared_genes": len(shared), "per_type": {}}
    y = check_labels(reference, label_key)
    min_req = max(_MIN_CELLS, cv)
    keep = {ct for ct, n in zip(*np.unique(y, return_counts=True)) if n >= min_req}
    var_idx = {str(g): i for i, g in enumerate(reference.var_names)}
    mat, matrix = get_counts(reference[:, [var_idx[g] for g in shared]], layer)
    rng = np.random.default_rng(seed)
    idx: list = []
    for ct in np.unique(y):
        pos = np.where(y == ct)[0]
        idx.extend(pos if len(pos) <= max_cells_per_type
                   else rng.choice(pos, max_cells_per_type, replace=False))
    idx = np.array(sorted(idx))
    y = y[idx]
    m = np.isin(y, list(keep))
    idx, y = idx[m], y[m]
    if len(np.unique(y)) < 2:
        return {**skip(f"fewer than 2 reference types with >= {min_req} cells"),
                "n_shared_genes": len(shared), "per_type": {}}
    A = np.clip(np.rint(_dense(mat[idx])), 0, None).astype(np.int64)  # thin to the panel's real depth
    tot = A.sum(1, keepdims=True)
    p = np.clip(target_depth / np.maximum(tot, 1), 0.0, 1.0)
    Xn = rng.binomial(A, np.broadcast_to(p, A.shape)).astype("float32").astype(float)
    # library-size normalize + log1p (no scanpy dep)
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
        null = float(mask.mean())                                 # base-rate (label-permutation) F1
        adj = (float(f1[i]) - null) / (1.0 - null)
        s = adj if null_adjusted else f1[i]
        tier = ("not_resolvable" if s < f1_weak or cf >= _CONF_NOT
                else "weak" if s < f1_resolvable or cf >= _CONF_WEAK else "resolvable")
        per[str(ct)] = {"f1": float(f1[i]), "f1_null": null, "f1_adjusted": adj,
                        "recall": float(rec[i]), "precision": float(prec[i]),
                        "confused_with": cw, "confused_frac": round(cf, 2), "tier": tier}
    tiers = [d["tier"] for d in per.values()]
    return {"status": "ok", "n_shared_genes": len(shared), "n_types": int(len(classes)),
            "macro_f1": ext["macro_f1"], "balanced_accuracy": ext["balanced_accuracy"], "ecs": ext["ecs"],
            "frac_resolvable": float(np.mean([t == "resolvable" for t in tiers])),
            "frac_not_resolvable": float(np.mean([t == "not_resolvable" for t in tiers])),
            "per_type": per, "matrix": matrix}
