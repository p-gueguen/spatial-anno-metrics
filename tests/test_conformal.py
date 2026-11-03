"""Tests for conformal.py - conformal prediction sets: coverage, Mondrian (class-conditional), hierarchy
collapse."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _peaked_proba(y, classes, peak=0.9, seed=0):
    """Confident-classifier probabilities: peak mass on the true class + Dirichlet-ish noise."""
    rng = np.random.default_rng(seed)
    K = len(classes)
    ci = {c: j for j, c in enumerate(classes)}
    P = np.full((len(y), K), (1 - peak) / (K - 1))
    for i, c in enumerate(y):
        P[i, ci[c]] = peak
    P = np.clip(P + rng.normal(0, 0.02, P.shape), 1e-3, None)
    return P / P.sum(1, keepdims=True)


def test_conformal_coverage_guarantee_and_singletons():
    """Split-conformal must (a) hit ~1-alpha coverage on an exchangeable cal/test split and
    (b) mostly return singletons for a confident classifier."""
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(0)
    classes = np.array(["A", "B", "C", "D"])
    y = rng.choice(classes, 1200)
    P = _peaked_proba(y, classes, peak=0.9)
    cal, te = slice(0, 800), slice(800, 1200)
    r = em.conformal_prediction_sets(P[te], classes, P[cal], y[cal], alpha=0.1, y_query=y[te])
    assert r["status"] == "ok"
    assert r["coverage"] >= 0.85                      # ~>= 1-alpha (finite-sample slack)
    assert r["summary"]["pct_singleton"] > 0.8        # confident classifier -> mostly singletons
    assert r["summary"]["mean_set_size"] < len(classes)   # not everything is in every set


def test_conformal_class_conditional_lifts_rare_type_coverage():
    """Marginal conformal under-covers a hard rare class; class-conditional restores it."""
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(1)
    classes = np.array(["maj", "rare"])
    # majority is easy (peak 0.95); the rare class is hard (peak 0.55 -> often under-covered marginally)
    y = np.array(["maj"] * 1000 + ["rare"] * 120)
    P = np.zeros((len(y), 2))
    for i, c in enumerate(y):
        peak = 0.95 if c == "maj" else 0.55
        P[i, 0 if c == "maj" else 1] = peak
        P[i, 1 if c == "maj" else 0] = 1 - peak
    P = np.clip(P + rng.normal(0, 0.03, P.shape), 1e-3, None)
    P /= P.sum(1, keepdims=True)
    idx = rng.permutation(len(y))
    cal, te = idx[:700], idx[700:]
    rm = em.conformal_prediction_sets(P[te], classes, P[cal], y[cal], alpha=0.1, y_query=y[te])
    rc = em.conformal_prediction_sets(P[te], classes, P[cal], y[cal], alpha=0.1,
                                      class_conditional=True, y_query=y[te])
    # class-conditional coverage of the rare type is at least as good as marginal (usually higher)
    assert rc["per_type_coverage"]["rare"] >= rm["per_type_coverage"]["rare"]


def test_conformal_hierarchy_collapse():
    """An ambiguous set over two subtypes of the same lineage collapses to a single parent."""
    from spatial_anno_metrics import conformal as em

    classes = np.array(["CD4", "CD8", "Bcell"])
    hierarchy = {"CD4": "T", "CD8": "T", "Bcell": "B"}
    # cal moderately confident (peak 0.5) -> inclusion threshold ~0.5, so a 0.5/0.5 query cell is a
    # 2-set while a 0.9 cell stays a singleton.
    y_cal = np.array(["CD4", "CD8", "Bcell"] * 20)
    P_cal = _peaked_proba(y_cal, classes, peak=0.5)
    P_q = np.array([[0.5, 0.5, 0.0],         # ambiguous within the T lineage
                    [0.9, 0.05, 0.05]])       # confident CD4
    r = em.conformal_prediction_sets(P_q, classes, P_cal, y_cal, alpha=0.2, hierarchy=hierarchy)
    assert len(r["sets"][0]) >= 2                              # the first cell is ambiguous
    assert r["collapsed_sets"][0] == ["T"]                     # ...but resolves to ONE lineage
    assert r["pct_ambiguous_one_parent"] == 1.0


def test_conformal_empty_calibration_is_skipped_not_all_inclusive():
    """No calibration cell -> no threshold: a skip, not qhat = 1 (which put every class in every set)."""
    from spatial_anno_metrics import conformal as em

    classes = np.array(["A", "B"])
    P = np.array([[0.9, 0.1], [0.2, 0.8]])
    out = em.conformal_prediction_sets(P, classes, np.empty((0, 2)), np.array([], dtype=str))
    assert out["status"] == "skipped" and "calibration" in out["reason"]
