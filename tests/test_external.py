"""Tests for external.py - the ground-truth harness: F1 / ARI / kappa, ECS, hierarchical and composition
accuracy, deconvolution proportions, AvgBIO."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def test_external_scores_perfect():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    out = em.external_scores(truth.copy(), truth)
    for m in ("accuracy", "balanced_accuracy", "macro_f1", "ari", "ami", "kappa", "ecs"):
        assert out[m] == pytest.approx(1.0), (m, out[m])


def test_external_scores_macro_f1_penalizes_rare_miss():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    pred = np.array(["A"] * 100)                       # miss every rare-B cell
    out = em.external_scores(pred, truth)
    assert out["accuracy"] == pytest.approx(0.9)       # inflated by the majority class
    assert out["macro_f1"] < out["accuracy"]           # macro-F1 exposes the dropped class
    assert out["balanced_accuracy"] < 0.9
    assert out["per_class_f1"]["B"] == pytest.approx(0.0)


def test_ecs_identical_is_one_and_resolution_aware():
    from spatial_anno_metrics import external as em

    a = np.array([0, 0, 1, 1, 2, 2])
    assert em.element_centric_similarity(a, a) == pytest.approx(1.0)
    allsame = np.zeros(6, dtype=int)
    refine = np.array([0, 3, 1, 1, 2, 2])              # split cluster 0 into two singletons
    assert em.element_centric_similarity(a, allsame) < 1.0
    # a near-refinement is more similar to `a` than collapsing everything to one cluster
    assert em.element_centric_similarity(a, refine) > em.element_centric_similarity(a, allsame)


def test_hierarchical_accuracy_partial_credit():
    from spatial_anno_metrics import external as em

    hier = {"CD4 T": "T", "CD8 T": "T", "B": "B"}
    truth = np.array(["CD4 T", "CD8 T", "B"])
    assert em.hierarchical_accuracy(truth.copy(), truth, hier)["hierarchical_accuracy"] == pytest.approx(1.0)

    partial = em.hierarchical_accuracy(np.array(["CD8 T", "CD8 T", "B"]), truth, hier, partial=0.5)
    assert partial["hierarchical_accuracy"] == pytest.approx((0.5 + 1 + 1) / 3)   # right lineage, wrong subtype
    assert partial["subtype_accuracy"] == pytest.approx(2 / 3)

    wrong = em.hierarchical_accuracy(np.array(["B", "CD8 T", "B"]), truth, hier, partial=0.5)
    assert wrong["hierarchical_accuracy"] == pytest.approx((0 + 1 + 1) / 3)       # wrong lineage
    assert partial["hierarchical_accuracy"] > wrong["hierarchical_accuracy"]


def test_composition_accuracy():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 60 + ["B"] * 40)
    same = em.composition_accuracy(truth.copy(), truth)
    assert same["l1"] == pytest.approx(0.0, abs=1e-9)
    assert same["jsd"] == pytest.approx(0.0, abs=1e-9)
    diff = em.composition_accuracy(np.array(["A"] * 90 + ["B"] * 10), truth)
    assert diff["l1"] > 0


def test_deconvolution_metrics_perfect_and_degraded():
    """OpenProblems spatial-decomposition metrics (R2 uniform-average + JSD axis=0) + RMSE on
    predicted vs true cell-type PROPORTION matrices (spots x types)."""
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(4), size=60)                 # 60 spots x 4 types, rows sum to 1
    perfect = em.deconvolution_metrics(true, true.copy())
    assert perfect["r2"] == pytest.approx(1.0)
    assert perfect["jsd"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["jsd_per_spot"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["rmse"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["pearson"] == pytest.approx(1.0)

    noisy = true + rng.uniform(0, 0.3, size=true.shape)
    noisy /= noisy.sum(1, keepdims=True)                      # still valid proportions
    deg = em.deconvolution_metrics(true, noisy)
    assert deg["r2"] < 1.0
    assert deg["jsd"] > 0 and deg["rmse"] > 0


def test_avg_bio_perfect_for_matching_labels(blobs, shuffle_labels):
    from spatial_anno_metrics import external as em

    a = blobs(n_per=60)
    a.obs["truth"] = a.obs["cell_type"]
    perfect = em.avg_bio(a, "cell_type", "truth", embedding="X_pca")
    assert perfect["ari"] == pytest.approx(1.0) and perfect["nmi"] == pytest.approx(1.0)
    assert perfect["avg_bio"] > em.avg_bio(shuffle_labels(a), "cell_type", "truth", embedding="X_pca")["avg_bio"]
