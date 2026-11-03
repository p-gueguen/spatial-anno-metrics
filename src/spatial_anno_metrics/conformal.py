"""Conformal prediction sets for cell-type annotation: distribution-free, model-agnostic typing
uncertainty (set size 1 = confident, >1 = ambiguous, 0 = novel)."""
from __future__ import annotations

import numpy as np

from ._validate import check_labels, skip

__all__ = ["conformal_prediction_sets"]


# --------------------------------------------------------------------------- #
# Conformal prediction sets (catalog s3d) - distribution-free typing uncertainty
# --------------------------------------------------------------------------- #
def conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha: float = 0.1,
                              class_conditional: bool = False,
                              hierarchy: dict | None = None, y_query=None) -> dict:
    """Distribution-free conformal prediction SETS for cell-type annotation (catalog s3d
    "conformal set size" - the cleanest reference-free formalization of abstention).

    Model-agnostic: bring any annotator's per-cell class-probability matrix (RCTD weights,
    SingleR/CellTypist/scANVI posteriors, marker-score softmax). Split-conformal LAC: the
    calibration nonconformity is ``s = 1 - p(true class)``; the (1-alpha) quantile ``qhat`` sets
    the inclusion threshold and the query set is ``C(x) = {y : p(x)_y >= 1 - qhat}``. Set size
    **1 = confident, >1 = ambiguous, 0 = novel/OOD**, with the finite-sample guarantee
    ``P(true in C(x)) >= 1 - alpha`` **provided calibration and query are exchangeable**. On
    native spatial data calibrated against a dissociated reference, exchangeability is broken by
    platform shift, so the nominal coverage is approximate (recalibrate on a small platform-matched
    labeled subset, e.g. protein-confirmed cells) - the set SIZE stays a valid relative uncertainty.

    Parameters
    ----------
    proba_query : (n_q, K) query class probabilities. classes : (K,) labels matching the columns.
    proba_cal, y_cal : (n_c, K) calibration probabilities + their TRUE labels.
    class_conditional : per-true-class (Mondrian) quantile - **recommended**; marginal conformal
        silently under-covers rare cell types while class-conditional restores per-type coverage.
    hierarchy : optional ``{leaf: parent}`` map - also returns lineage-collapsed sets (scConform's
        ontology idea: an ambiguous SUBTYPE set becomes a confident LINEAGE call).
    y_query : optional TRUE query labels -> also returns empirical coverage + per-type coverage.

    Returns ``{status, sets, set_size, summary, [collapsed_sets, pct_ambiguous_one_parent],
    [coverage, per_type_coverage]}``; an empty calibration set -> a skip (no threshold can be set).
    Depends on numpy only.
    """
    if len(y_cal) == 0:
        return skip("empty calibration set: conformal sets need calibration probabilities and true labels")
    proba_query = np.asarray(proba_query, dtype=float)
    proba_cal = np.asarray(proba_cal, dtype=float)
    classes = np.asarray([str(c) for c in classes])
    y_cal = check_labels(y_cal)
    cls_idx = {c: j for j, c in enumerate(classes)}

    true_col = np.array([cls_idx[c] for c in y_cal])
    s_cal = 1.0 - proba_cal[np.arange(len(y_cal)), true_col]

    def _qhat(scores):
        n = len(scores)
        if n == 0:
            return float("nan")                                       # no calibration: no threshold
        level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
        return float(np.quantile(scores, level, method="higher"))

    if class_conditional:
        thr = {c: (_qhat(s_cal[y_cal == c]) if (y_cal == c).sum() >= 10 else _qhat(s_cal))
               for c in classes}
        incl_thresh = np.array([1.0 - thr[c] for c in classes])      # per-class inclusion threshold
        include = proba_query >= incl_thresh[None, :]
    else:
        q = _qhat(s_cal)
        include = proba_query >= (1.0 - q)

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
    return out
