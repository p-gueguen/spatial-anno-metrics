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


# --------------------------------------------------------------------------- #
# calibration error (ECE / class-wise ECE), CovGap, clustered Mondrian (ROADMAP #24)
# --------------------------------------------------------------------------- #
def test_ece_hand_computed():
    """Top-label ECE on 4 cells, two sharing the (0.9, 1.0] bin: 0.5*|0.5-0.935| + 0.25*0.35 + 0.25*0.25."""
    from spatial_anno_metrics import conformal as em

    P = np.array([[0.95, 0.05], [0.92, 0.08], [0.65, 0.35], [0.25, 0.75]])
    y = np.array(["A", "B", "A", "B"])
    r = em.calibration_error(P, ["A", "B"], y, n_bins=10)
    assert r["status"] == "ok"
    assert r["ece"] == pytest.approx(0.3675)
    assert r["accuracy"] == pytest.approx(0.75)


def test_calibration_error_bin_edges_are_left_open_right_closed():
    """Bins are (lo, hi] (searchsorted side='left'): a 0.5 confidence sits in bin 0 of n_bins=2, not bin 1."""
    from spatial_anno_metrics import conformal as em

    # Cell 0: conf=0.5 (on the interior edge), hit=1; Cell 1: conf=0.9, hit=0.
    # Under (0, 0.5] and (0.5, 1.0]: bin 0 gap |1 - 0.5| = 0.5, bin 1 gap |0 - 0.9| = 0.9 -> ECE = 0.70.
    # Under [0, 0.5) and [0.5, 1.0] (side='right'): both land in bin 1 -> |0.5 - 0.7| = 0.20.
    P = np.array([[0.5, 0.5], [0.9, 0.1]])
    r = em.calibration_error(P, ["A", "B"], np.array(["A", "B"]), n_bins=2)
    assert r["ece"] == pytest.approx(0.70)


def test_classwise_ece_matches_sklearn_calibration_curve():
    """Class-wise ECE = mean over classes of the count-weighted |freq - mean p| per bin; sklearn's
    calibration_curve gives the per-bin values independently (positive control for the binning)."""
    from sklearn.calibration import calibration_curve

    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(3)
    classes = np.array(["A", "B", "C"])
    P = rng.dirichlet([1, 1, 1], 3000)
    y = np.array([rng.choice(classes, p=p ** 2 / (p ** 2).sum()) for p in P])     # miscalibrated on purpose
    r = em.calibration_error(P, classes, y, n_bins=15)
    want = []
    for k, c in enumerate(classes):
        frac, mean_p = calibration_curve(y == c, P[:, k], n_bins=15, strategy="uniform")
        n_b = np.histogram(P[:, k], bins=np.linspace(0, 1, 16))[0]
        n_b = n_b[n_b > 0]
        want.append(float(np.sum(n_b / len(y) * np.abs(frac - mean_p))))
    assert r["classwise_ece"] == pytest.approx(np.mean(want), abs=1e-12)
    assert r["per_class_ece"] == pytest.approx(dict(zip(classes, want)), abs=1e-12)


def test_ece_separates_calibrated_from_overconfident():
    """Labels drawn from p: ECE ~ 0. Same labels, probabilities sharpened (p^4): ECE large."""
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(0)
    classes = np.array(["A", "B", "C", "D"])
    P = rng.dirichlet([1, 1, 1, 1], 20000)
    y = classes[(rng.random(len(P))[:, None] > P.cumsum(1)).sum(1)]
    sharp = P ** 4 / (P ** 4).sum(1, keepdims=True)
    good, bad = em.calibration_error(P, classes, y), em.calibration_error(sharp, classes, y)
    assert good["ece"] < 0.02 and good["classwise_ece"] < 0.02
    assert bad["ece"] > 0.1 and bad["classwise_ece"] > good["classwise_ece"] + 0.02


def test_calibration_error_refuses_bad_input():
    from spatial_anno_metrics import conformal as em

    with pytest.raises(ValueError, match="not in classes"):
        em.calibration_error(np.array([[0.6, 0.4]]), ["A", "B"], ["C"])
    with pytest.raises(ValueError, match="missing"):
        em.calibration_error(np.array([[0.6, 0.4]]), ["A", "B"], np.array([None], dtype=object))
    assert em.calibration_error(np.empty((0, 2)), ["A", "B"], np.array([], dtype=str))["status"] == "skipped"


def test_covgap_is_mean_absolute_per_type_coverage_gap():
    """Ding et al. 2023 CovGap = mean_y |coverage_y - (1 - alpha)|."""
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(0)
    classes = np.array(["A", "B", "C", "D"])
    y = rng.choice(classes, 1200)
    P = _peaked_proba(y, classes, peak=0.7)
    r = em.conformal_prediction_sets(P[800:], classes, P[:800], y[:800], alpha=0.1, y_query=y[800:])
    want = np.mean([abs(v - 0.9) for v in r["per_type_coverage"].values()])
    assert r["covgap"] == pytest.approx(want)


def _two_difficulty_world(n_per_class, seed, K=40):
    """K classes; the first half easy (true-class p ~ Beta(8, 2)), the second half hard (Beta(2, 3))."""
    rng = np.random.default_rng(seed)
    classes = np.array([f"c{k:02d}" for k in range(K)])
    y = np.repeat(classes, n_per_class)
    hard = np.isin(y, classes[K // 2:])
    pt = np.where(hard, rng.beta(2, 3, len(y)), rng.beta(8, 2, len(y)))
    P = np.repeat(((1 - pt) / (K - 1))[:, None], K, axis=1)
    P[np.arange(len(y)), np.searchsorted(classes, y)] = pt
    return classes, y, P


def test_clustered_mondrian_lowers_covgap_and_finds_the_difficulty_groups():
    """Few calibration cells per class: marginal over-covers easy and under-covers hard classes; clustered
    Mondrian (Ding et al. 2023) pools similar classes and lowers CovGap."""
    from spatial_anno_metrics import conformal as em

    classes, y_cal, P_cal = _two_difficulty_world(30, seed=0)
    _, y_q, P_q = _two_difficulty_world(200, seed=1)
    kw = dict(alpha=0.1, y_query=y_q)
    rm = em.conformal_prediction_sets(P_q, classes, P_cal, y_cal, **kw)
    rk = em.conformal_prediction_sets(P_q, classes, P_cal, y_cal, class_conditional="clustered", **kw)
    assert rk["summary"]["class_conditional"] == "clustered"
    assert rk["covgap"] < rm["covgap"] - 0.03
    hard = [c for c in classes[20:] if rk["clusters"][c] >= 0]
    easy = [c for c in classes[:20] if rk["clusters"][c] >= 0]
    assert len(hard) >= 10 and len(easy) >= 10
    assert not ({rk["clusters"][c] for c in hard} & {rk["clusters"][c] for c in easy})   # never pooled together
    assert np.mean([rk["per_type_coverage"][c] for c in classes[20:]]) >= 0.85


def test_clustered_mondrian_rare_class_goes_to_the_null_cluster():
    """A class with too few calibration cells to embed falls back to the marginal threshold (cluster -1)."""
    from spatial_anno_metrics import conformal as em

    classes, y_cal, P_cal = _two_difficulty_world(30, seed=0, K=10)
    keep = (y_cal != "c09") | (np.cumsum(y_cal == "c09") <= 3)                # c09: 3 calibration cells
    r = em.conformal_prediction_sets(P_cal, classes, P_cal[keep], y_cal[keep], class_conditional="clustered")
    assert r["clusters"]["c09"] == -1
    with pytest.raises(ValueError, match="class_conditional"):
        em.conformal_prediction_sets(P_cal, classes, P_cal, y_cal, class_conditional="kmeans")


# --------------------------------------------------------------------------- #
# conformal novelty p-values (ROADMAP #20)
# --------------------------------------------------------------------------- #
def test_conformal_pvalues_and_bh_hand_and_scipy():
    """p_j = (1 + #{cal >= s_j}) / (n_cal + 1); BH matches scipy's false_discovery_control."""
    from scipy.stats import false_discovery_control

    from spatial_anno_metrics.conformal import _bh, _conformal_pvalues

    p = _conformal_pvalues(np.array([1.0, 2.0, 3.0, 4.0]), np.array([0.5, 2.0, 5.0]))
    assert p == pytest.approx([5 / 5, 4 / 5, 1 / 5])
    rng = np.random.default_rng(0)
    pv = np.concatenate([rng.uniform(size=200), rng.uniform(0, 1e-3, 20)])
    assert _bh(pv) == pytest.approx(false_discovery_control(pv, method="bh"))


@pytest.mark.parametrize("method", ["pca", "iforest"])
def test_conformal_novelty_flags_novel_and_controls_fdr(method):
    """Inliers from the reference distribution vs a shifted novel population: novel cells get small p-values,
    and the BH-flagged set keeps the inlier share near the FDR level; an all-inlier query flags (almost) none."""
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(0)
    d = 30
    W = rng.normal(size=(5, d))                              # inliers live near a 5-dim subspace
    ref = rng.normal(size=(1000, 5)) @ W + 0.3 * rng.normal(size=(1000, d))
    inl = rng.normal(size=(500, 5)) @ W + 0.3 * rng.normal(size=(500, d))
    nov = rng.normal(size=(100, 5)) @ W + 0.3 * rng.normal(size=(100, d)) + 6.0 * rng.normal(size=(1, d))
    q = np.vstack([inl, nov])
    is_nov = np.r_[np.zeros(500, bool), np.ones(100, bool)]
    r = em.conformal_novelty(q, ref, method=method, fdr=0.1, seed=0)
    assert r["status"] == "ok" and r["summary"]["method"] == method
    assert np.median(r["pvalues"][is_nov]) < 0.01 < np.median(r["pvalues"][~is_nov])
    flagged = r["novel"]
    assert flagged[is_nov].mean() > 0.8
    assert (~is_nov[flagged]).mean() <= 0.15                      # false-discovery proportion
    null = em.conformal_novelty(inl, ref, method=method, fdr=0.1, seed=0)
    assert null["novel"].mean() <= 0.02


def test_conformal_novelty_input_checks():
    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="features"):
        em.conformal_novelty(rng.normal(size=(5, 4)), rng.normal(size=(50, 3)))
    with pytest.raises(ValueError, match="method"):
        em.conformal_novelty(rng.normal(size=(5, 3)), rng.normal(size=(50, 3)), method="lof")
    bad = rng.normal(size=(5, 3))
    bad[2, 1] = np.nan
    with pytest.raises(ValueError, match="NaN or infinite"):
        em.conformal_novelty(bad, rng.normal(size=(50, 3)))
    assert em.conformal_novelty(rng.normal(size=(5, 3)), rng.normal(size=(2, 3)))["status"] == "skipped"


@pytest.mark.parametrize("method", ["pca", "iforest"])
def test_conformal_novelty_sparse_inputs_do_not_densify(monkeypatch, method):
    """Sparse CSR query and reference run without calling _dense, and PCA residuals match dense to 1e-10."""
    import scipy.sparse as sp

    from spatial_anno_metrics import conformal as em

    rng = np.random.default_rng(7)
    ref = rng.poisson(0.5, size=(120, 40)).astype(float)
    q = rng.poisson(0.5, size=(30, 40)).astype(float)
    q[-5:] += 6.0
    dense_res = em.conformal_novelty(q, ref, method=method, n_components=5, seed=0)

    monkeypatch.setattr(em, "_dense", lambda *a, **k: (_ for _ in ()).throw(AssertionError("densified")))
    sparse_res = em.conformal_novelty(sp.csr_matrix(q), sp.csr_matrix(ref), method=method, n_components=5, seed=0)
    assert sparse_res["status"] == "ok"
    if method == "pca":
        np.testing.assert_allclose(sparse_res["score"], dense_res["score"], rtol=1e-7, atol=1e-9)
        np.testing.assert_allclose(sparse_res["pvalues"], dense_res["pvalues"], rtol=0, atol=1e-12)
