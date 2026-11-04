"""Tests for external.py - the ground-truth harness: F1 / ARI / kappa, ECS, hierarchical and composition
accuracy, deconvolution proportions, AvgBIO."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def test_external_scores_perfect():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    out = em.external_scores(y_true=truth, y_pred=truth.copy())
    for m in ("accuracy", "balanced_accuracy", "macro_f1", "ari", "ami", "kappa", "ecs"):
        assert out[m] == pytest.approx(1.0), (m, out[m])
    assert out["status"] == "ok"


def test_external_scores_macro_f1_penalizes_rare_miss():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    pred = np.array(["A"] * 100)                       # miss every rare-B cell
    out = em.external_scores(y_true=truth, y_pred=pred)
    assert out["accuracy"] == pytest.approx(0.9)       # inflated by the majority class
    assert out["macro_f1"] < out["accuracy"]           # macro-F1 exposes the dropped class
    assert out["balanced_accuracy"] < 0.9
    assert out["per_class_f1"]["B"] == pytest.approx(0.0)


def test_ecs_identical_is_one_and_resolution_aware():
    from spatial_anno_metrics import external as em

    a = np.array([0, 0, 1, 1, 2, 2])
    same = em.element_centric_similarity(a, a)
    assert same["status"] == "ok" and same["ecs"] == pytest.approx(1.0)
    allsame = np.zeros(6, dtype=int)
    refine = np.array([0, 3, 1, 1, 2, 2])              # split cluster 0 into two singletons
    assert em.element_centric_similarity(a, allsame)["ecs"] < 1.0
    # a near-refinement is more similar to `a` than collapsing everything to one cluster
    assert em.element_centric_similarity(a, refine)["ecs"] > em.element_centric_similarity(a, allsame)["ecs"]


def test_hierarchical_accuracy_partial_credit():
    from spatial_anno_metrics import external as em

    hier = {"CD4 T": "T", "CD8 T": "T", "B": "B"}
    truth = np.array(["CD4 T", "CD8 T", "B"])
    exact = em.hierarchical_accuracy(y_true=truth, y_pred=truth.copy(), hierarchy=hier)
    assert exact["status"] == "ok" and exact["hierarchical_accuracy"] == pytest.approx(1.0)

    partial = em.hierarchical_accuracy(y_true=truth, y_pred=np.array(["CD8 T", "CD8 T", "B"]), hierarchy=hier,
                                       partial=0.5)
    assert partial["hierarchical_accuracy"] == pytest.approx((0.5 + 1 + 1) / 3)   # right lineage, wrong subtype
    assert partial["subtype_accuracy"] == pytest.approx(2 / 3)

    wrong = em.hierarchical_accuracy(y_true=truth, y_pred=np.array(["B", "CD8 T", "B"]), hierarchy=hier,
                                     partial=0.5)
    assert wrong["hierarchical_accuracy"] == pytest.approx((0 + 1 + 1) / 3)       # wrong lineage
    assert partial["hierarchical_accuracy"] > wrong["hierarchical_accuracy"]


def test_composition_accuracy():
    from spatial_anno_metrics import external as em

    truth = np.array(["A"] * 60 + ["B"] * 40)
    same = em.composition_accuracy(y_true=truth, y_pred=truth.copy())
    assert same["status"] == "ok" and same["l1"] == pytest.approx(0.0, abs=1e-9)
    assert same["jsd"] == pytest.approx(0.0, abs=1e-9)
    diff = em.composition_accuracy(y_true=truth, y_pred=np.array(["A"] * 90 + ["B"] * 10))
    assert diff["l1"] > 0
    flat = em.composition_accuracy(y_true=np.array(["A", "B"] * 5), y_pred=np.array(["A", "B"] * 5))
    assert np.isnan(flat["pearson"])                              # equal shares: r is undefined, not 1.0


def test_deconvolution_metrics_perfect_and_degraded():
    """OpenProblems spatial-decomposition metrics (R2 uniform-average + JSD axis=0) + RMSE on
    predicted vs true cell-type PROPORTION matrices (spots x types)."""
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(4), size=60)                 # 60 spots x 4 types, rows sum to 1
    perfect = em.deconvolution_metrics(true_prop=true, pred_prop=true.copy())
    assert perfect["status"] == "ok" and perfect["r2"] == pytest.approx(1.0)
    assert perfect["jsd"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["jsd_per_spot"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["rmse"] == pytest.approx(0.0, abs=1e-9)
    assert perfect["pearson"] == pytest.approx(1.0)

    noisy = true + rng.uniform(0, 0.3, size=true.shape)
    noisy /= noisy.sum(1, keepdims=True)                      # still valid proportions
    deg = em.deconvolution_metrics(true_prop=true, pred_prop=noisy)
    assert deg["r2"] < 1.0
    assert deg["jsd"] > 0 and deg["rmse"] > 0


def test_avg_bio_perfect_for_matching_labels(blobs, shuffle_labels):
    from spatial_anno_metrics import external as em

    a = blobs(n_per=60)
    a.obs["truth"] = a.obs["cell_type"]
    perfect = em.avg_bio(a, "cell_type", "truth", embedding="X_pca")
    assert perfect["status"] == "ok"
    assert perfect["ari"] == pytest.approx(1.0) and perfect["nmi"] == pytest.approx(1.0)
    assert perfect["avg_bio"] > em.avg_bio(shuffle_labels(a), "cell_type", "truth", embedding="X_pca")["avg_bio"]


def test_external_harness_takes_keywords_and_rejects_missing_labels():
    """(y_true, y_pred) as in sklearn, keyword-only, so a swapped positional call cannot happen silently."""
    from spatial_anno_metrics import external as em

    truth = np.array(["A", "B", "B"])
    for fn in (em.external_scores, em.composition_accuracy):
        with pytest.raises(TypeError):
            fn(truth, truth)
    with pytest.raises(TypeError):
        em.deconvolution_metrics(np.eye(2), np.eye(2))
    with pytest.raises(ValueError, match="has 1 missing"):
        em.external_scores(y_true=truth, y_pred=np.array(["A", None, "B"], dtype=object))


def test_avg_bio_skips_a_single_truth_label(blobs):
    from spatial_anno_metrics import external as em

    a = blobs(n_per=20)
    a.obs["truth"] = "T0"
    out = em.avg_bio(a, "cell_type", "truth", embedding="X_pca")
    assert out["status"] == "skipped" and "found 1" in out["reason"]


# --------------------------------------------------------------------------- #
# hierarchical_accuracy: credit by LCA depth (ROADMAP #9)
# --------------------------------------------------------------------------- #
_TREE = {"CD4 naive": "CD4 T", "CD4 mem": "CD4 T", "CD4 T": "T", "CD8 T": "T", "B naive": "B"}


def test_hierarchical_coarse_call_beats_a_wrong_sibling():
    """An honest coarse call (the parent) must score above a confident wrong sibling - the 0.7.0 defect."""
    from spatial_anno_metrics import external as em

    hier = {"CD4 T": "T", "CD8 T": "T"}
    truth = np.array(["CD4 T"])
    coarse = em.hierarchical_accuracy(y_true=truth, y_pred=np.array(["T"]), hierarchy=hier)
    sibling = em.hierarchical_accuracy(y_true=truth, y_pred=np.array(["CD8 T"]), hierarchy=hier)
    assert coarse["hierarchical_accuracy"] == pytest.approx(2 / 3)     # hP 1, hR 1/2
    assert sibling["hierarchical_accuracy"] == pytest.approx(1 / 2)    # hP 1/2, hR 1/2
    assert coarse["hierarchical_precision"] == pytest.approx(1.0)
    assert coarse["hierarchical_recall"] == pytest.approx(0.5)


def test_hierarchical_credit_is_kiritchenko_hf_from_lca_depth():
    """Per cell: hP = (lca+1)/(depth(pred)+1), hR = (lca+1)/(depth(true)+1), credit = their harmonic mean."""
    from spatial_anno_metrics import external as em

    truth = np.array(["CD4 naive", "CD4 naive", "CD4 naive", "CD4 naive", "CD4 naive", "B naive"])
    pred = np.array(["CD4 naive", "CD4 mem", "CD4 T", "CD8 T", "T", "CD4 naive"])
    # depths: CD4 naive 2, CD4 mem 2, CD4 T 1, CD8 T 1, T 0, B naive 1
    hp = np.array([3 / 3, 2 / 3, 2 / 2, 1 / 2, 1 / 1, 0])
    hr = np.array([3 / 3, 2 / 3, 2 / 3, 1 / 3, 1 / 3, 0])
    hf = np.where(hp + hr > 0, 2 * hp * hr / np.where(hp + hr > 0, hp + hr, 1), 0)
    out = em.hierarchical_accuracy(y_true=truth, y_pred=pred, hierarchy=_TREE)
    assert out["status"] == "ok" and out["credit"] == "lca"
    assert out["hierarchical_accuracy"] == pytest.approx(hf.mean())
    assert out["hierarchical_precision"] == pytest.approx((3 + 2 + 2 + 1 + 1 + 0) / (3 + 3 + 2 + 2 + 1 + 3))
    assert out["hierarchical_recall"] == pytest.approx((3 + 2 + 2 + 1 + 1 + 0) / (3 * 5 + 2))
    assert out["subtype_accuracy"] == pytest.approx(1 / 6)
    assert out["lineage_accuracy"] == pytest.approx(5 / 6)             # share the root T
    # non-leaf (coarse) calls: CD4 T and T are parents in the tree, CD8 T / CD4 mem are leaves
    assert out["coarse_call_frac"] == pytest.approx(2 / 6)
    assert out["coarse_call_accuracy"] == pytest.approx(1.0)           # both are ancestors of the truth


def test_hierarchical_coarse_call_accuracy_counts_wrong_coarse_calls_and_is_nan_without_any():
    from spatial_anno_metrics import external as em

    out = em.hierarchical_accuracy(y_true=np.array(["CD8 T", "CD4 naive"]),
                                   y_pred=np.array(["CD4 T", "CD4 T"]), hierarchy=_TREE)
    assert out["coarse_call_frac"] == pytest.approx(1.0)
    assert out["coarse_call_accuracy"] == pytest.approx(0.5)           # CD4 T is not above CD8 T
    leaf = em.hierarchical_accuracy(y_true=np.array(["CD8 T"]), y_pred=np.array(["CD8 T"]), hierarchy=_TREE)
    assert leaf["coarse_call_frac"] == 0.0 and np.isnan(leaf["coarse_call_accuracy"])


def test_hierarchical_lineage_credit_reproduces_the_old_default():
    """credit='lineage' with partial=0.5 gives the 0.7.0-and-earlier values exactly on a two-level map."""
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(1)
    hier = {"T_CD4": "T", "T_CD8": "T", "B": "B", "Myeloid": "Myeloid"}
    labels = np.array(["T_CD4", "T_CD8", "B", "Myeloid", "T"])
    truth, pred = rng.choice(labels[:4], 500), rng.choice(labels, 500)

    pp = np.array([hier.get(x, x) for x in pred])                      # the old implementation, verbatim
    pt = np.array([hier.get(x, x) for x in truth])
    old = np.where(pred == truth, 1.0, np.where(pp == pt, 0.5, 0.0))
    out = em.hierarchical_accuracy(y_true=truth, y_pred=pred, hierarchy=hier, credit="lineage", partial=0.5)
    assert out["hierarchical_accuracy"] == old.mean()
    assert out["subtype_accuracy"] == (pred == truth).mean() and out["lineage_accuracy"] == (pp == pt).mean()
    new = em.hierarchical_accuracy(y_true=truth, y_pred=pred, hierarchy=hier)
    assert new["hierarchical_accuracy"] > out["hierarchical_accuracy"]  # the coarse 'T' calls now earn 2/3
    with pytest.raises(ValueError, match="credit"):
        em.hierarchical_accuracy(y_true=truth, y_pred=pred, hierarchy=hier, credit="depth")


# --------------------------------------------------------------------------- #
# deconvolution: presence/absence AUPR, per-type SSIM, rank aggregation (ROADMAP #26)
# --------------------------------------------------------------------------- #
def _li_ssim(x, y):
    """Li et al. 2022 SSIM from the paper: columns scaled by their max, C1=(0.01)^2, C2=(0.03)^2."""
    x, y = x / x.max(), y / y.max()
    mx, my = x.mean(), y.mean()
    sx, sy, sxy = x.std(), y.std(), ((x - mx) * (y - my)).mean()
    c1, c2 = 0.01 ** 2, 0.03 ** 2
    return (2 * mx * my + c1) * (2 * sxy + c2) / ((mx ** 2 + my ** 2 + c1) * (sx ** 2 + sy ** 2 + c2))


def test_deconvolution_ssim_and_aupr():
    import pandas as pd
    from sklearn.metrics import average_precision_score

    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(4), size=60)
    true[true < 0.1] = 0                                     # absent types, so presence/absence is defined
    true /= true.sum(1, keepdims=True)
    perfect = em.deconvolution_metrics(true_prop=true, pred_prop=true.copy())
    assert perfect["aupr"] == pytest.approx(1.0) and perfect["ssim"] == pytest.approx(1.0)

    pred = 0.7 * true + 0.3 * rng.dirichlet(np.ones(4), size=60)
    out = em.deconvolution_metrics(true_prop=pd.DataFrame(true, columns=list("abcd")), pred_prop=pred)
    per = [_li_ssim(true[:, k], pred[:, k]) for k in range(4)]
    assert out["ssim_per_type"] == pytest.approx(dict(zip("abcd", per)))
    assert out["ssim"] == pytest.approx(np.mean(per))
    assert out["aupr"] == pytest.approx(average_precision_score((true > 0).ravel(), pred.ravel()))
    assert out["aupr"] < 1.0 and out["ssim"] < 1.0


def test_deconvolution_undefined_aupr_and_ssim_are_nan():
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(3), size=20)                # every type present in every spot
    out = em.deconvolution_metrics(true_prop=true, pred_prop=true[::-1].copy())
    assert np.isnan(out["aupr"])                             # no absent entry: precision-recall undefined
    true[:, 2] = 0
    true /= true.sum(1, keepdims=True)
    out = em.deconvolution_metrics(true_prop=true, pred_prop=rng.dirichlet(np.ones(3), size=20))
    assert np.isnan(out["ssim_per_type"][2]) and np.isfinite(out["ssim"])   # all-zero column: 0/0


def test_deconvolution_rank_score_orders_methods():
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(4) * 0.5, size=80)
    true[true < 0.05] = 0                                    # absent types, so AUPR is defined
    true /= true.sum(1, keepdims=True)
    noise = rng.dirichlet(np.ones(4), size=80)
    res = {f"m{s}": em.deconvolution_metrics(true_prop=true, pred_prop=(1 - s) * true + s * noise)
           for s in (0.6, 0.0, 0.3)}
    out = em.deconvolution_rank_score(res)
    assert out["status"] == "ok" and out["metrics"] == ["rmse", "jsd_per_spot", "aupr"]
    assert out["order"] == ["m0.0", "m0.3", "m0.6"]
    assert out["mean_rank"] == {"m0.0": 1.0, "m0.3": 2.0, "m0.6": 3.0}
    tie = em.deconvolution_rank_score({"a": {"rmse": 0.1}, "b": {"rmse": 0.1}, "c": {"rmse": 0.3}},
                                      metrics=("rmse",))
    assert tie["mean_rank"] == {"a": 1.5, "b": 1.5, "c": 3.0}
    assert em.deconvolution_rank_score({"a": res["m0.0"]})["status"] == "skipped"
    nan = em.deconvolution_rank_score({"a": {"aupr": float("nan")}, "b": {"aupr": 0.5}}, metrics=("aupr",))
    assert nan["status"] == "skipped" and "a" in nan["reason"]
    with pytest.raises(ValueError, match="direction"):
        em.deconvolution_rank_score({"a": {"x": 1}, "b": {"x": 2}}, metrics=("x",))


# --------------------------------------------------------------------------- #
# spARI (Yan, Feng & Luo, Biometrics 2025) - ROADMAP #27
# --------------------------------------------------------------------------- #
def _spari_bruteforce(r, c, xy, alpha):
    """The definition, one pair at a time (coordinates min-max scaled per axis; pairs at distance 0 get no
    spatial credit, as in the authors' package)."""
    xy = (xy - xy.min(0)) / (xy.max(0) - xy.min(0))
    n = len(r)
    npairs = n * (n - 1) / 2
    num = sum_f = sum_h = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(xy[i] - xy[j])
            f = alpha * np.exp(-d ** 2) if d > 0 else 0.0
            h = alpha * (1 - np.exp(-d ** 2)) if d > 0 else 0.0
            sr, sc = r[i] == r[j], c[i] == c[j]
            num += 1.0 if sr == sc else (f if sc else h)
            sum_f, sum_h = sum_f + f, sum_h + h
    p = sum((r == k).sum() * ((r == k).sum() - 1) / 2 for k in set(r)) / npairs
    q = sum((c == k).sum() * ((c == k).sum() - 1) / 2 for k in set(c)) / npairs
    e = p * q + (1 - p) * (1 - q) + (1 - p) * q * sum_f / npairs + p * (1 - q) * sum_h / npairs
    return num / npairs, (num / npairs - e) / (1 - e)


def test_spatial_ari_matches_the_pairwise_definition():
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(3)
    xy = rng.uniform(size=(70, 2)) * [5, 2]
    r = np.where(xy[:, 0] < 2.5, "L", "R")
    c = r.copy()
    c[rng.choice(70, 15, replace=False)] = rng.choice(["L", "R", "X"], 15)
    xy[1] = xy[0]                                            # one coincident pair
    spri, spari = _spari_bruteforce(r, c, xy, 0.8)
    out = em.spatial_ari(y_true=r, y_pred=c, coords=xy)
    assert out["status"] == "ok" and out["n_cells"] == 70 and out["subsampled"] is False
    assert out["spri"] == pytest.approx(spri, rel=1e-12) and out["spari"] == pytest.approx(spari, rel=1e-12)


def test_spatial_ari_is_ari_at_alpha_zero_and_rewards_near_mistakes():
    from sklearn.metrics import adjusted_rand_score

    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    xy = np.c_[np.r_[rng.uniform(0, 1, 60), rng.uniform(1.2, 2.2, 60)], rng.uniform(0, 1, 120)]
    r = np.repeat(["A", "B"], 60)
    c = r.copy()
    near = np.argsort(np.abs(xy[:60, 0] - 1.2))[:10]         # A cells closest to the B block
    far = np.argsort(xy[:60, 0])[:10]                        # A cells farthest from it
    c_near, c_far = c.copy(), c.copy()
    c_near[near], c_far[far] = "B", "B"
    assert em.spatial_ari(y_true=r, y_pred=c_near, coords=xy, alpha=0.0)["spari"] == pytest.approx(
        adjusted_rand_score(r, c_near), rel=1e-12)
    s_near = em.spatial_ari(y_true=r, y_pred=c_near, coords=xy)["spari"]
    s_far = em.spatial_ari(y_true=r, y_pred=c_far, coords=xy)["spari"]
    assert s_near > s_far                                    # same ARI, spatially milder mistake
    assert adjusted_rand_score(r, c_near) == pytest.approx(adjusted_rand_score(r, c_far))
    same = em.spatial_ari(y_true=r, y_pred=r.copy(), coords=xy)
    assert same["spri"] == pytest.approx(1.0) and same["spari"] == pytest.approx(1.0)


def test_spatial_ari_subsamples_large_inputs_and_validates():
    from spatial_anno_metrics import external as em

    rng = np.random.default_rng(0)
    xy = rng.uniform(size=(400, 2))
    r = np.where(xy[:, 0] < 0.5, "A", "B")
    out = em.spatial_ari(y_true=r, y_pred=r.copy(), coords=xy, max_cells=100)
    assert out["subsampled"] is True and out["n_cells"] == 100 and out["spari"] == pytest.approx(1.0)
    with pytest.raises(ValueError, match="alpha"):
        em.spatial_ari(y_true=r, y_pred=r, coords=xy, alpha=1.0)
    with pytest.raises(ValueError, match="coords"):
        em.spatial_ari(y_true=r, y_pred=r, coords=xy[:10])
    with pytest.raises(TypeError):
        em.spatial_ari(r, r, xy)
    assert em.spatial_ari(y_true=r[:1], y_pred=r[:1], coords=xy[:1])["status"] == "skipped"


# --------------------------------------------------------------------------- #
# Goldens: ECS vs clusim, deconvolution vs OpenProblems (tests/goldens/*/SOURCE.md)
# --------------------------------------------------------------------------- #
def _golden(name, file="expected.json"):
    import json
    from pathlib import Path

    def no_nan(token):
        raise ValueError(token)

    return json.loads((Path(__file__).parent / "goldens" / name / file).read_text(), parse_constant=no_nan)


def test_ecs_matches_clusim_golden():
    from spatial_anno_metrics import external as em

    pairs, exp = _golden("ecs_clusim", "pairs.json"), _golden("ecs_clusim")["results"]
    for key, alpha in (("alpha=0.9", 0.9), ("alpha=0.5", 0.5)):
        for name, p in pairs.items():
            ours = em.element_centric_similarity(p["a"], p["b"], alpha=alpha)["ecs"]
            assert ours == pytest.approx(exp[key][name]["ecs"], rel=1e-12, abs=1e-12), (key, name)


def test_deconvolution_matches_openproblems_golden():
    from scipy.spatial.distance import jensenshannon

    from spatial_anno_metrics import external as em

    props, exp = _golden("deconv_openproblems", "proportions.json"), _golden("deconv_openproblems")["results"]
    for name, p in props.items():
        ours = em.deconvolution_metrics(true_prop=p["true"], pred_prop=p["pred"])
        ours_pr = em.deconvolution_metrics(true_prop=p["true"], pred_prop=p["pred"], aupr_method="precrec")
        assert ours["r2"] == pytest.approx(exp[name]["r2"], rel=1e-12, abs=1e-12), name
        if exp[name]["jsd"] is not None:
            assert ours["jsd"] == pytest.approx(exp[name]["jsd"], rel=1e-12, abs=1e-12), name
        if exp[name]["aupr"] is None:
            assert np.isnan(ours["aupr"]) and np.isnan(ours["aupr_precrec"]) and np.isnan(ours_pr["aupr"])
        else:
            assert ours["aupr"] == pytest.approx(exp[name]["aupr"], rel=1e-12, abs=1e-12), name
            assert ours["aupr_precrec"] == pytest.approx(exp[name]["aupr_precrec"], rel=1e-12, abs=1e-12), name
            assert ours_pr["aupr"] == pytest.approx(exp[name]["aupr_precrec"], rel=1e-12, abs=1e-12), name
        assert ours["ssim"] == pytest.approx(exp[name]["ssim"], rel=1e-12, abs=1e-12), name
        for k, want_ssim in enumerate(exp[name]["ssim_per_type"]):
            if want_ssim is None:
                assert np.isnan(ours["ssim_per_type"][k])
            else:
                assert ours["ssim_per_type"][k] == pytest.approx(want_ssim, rel=1e-12, abs=1e-12), (name, k)
    # absent_type: OpenProblems' mean is NaN (0/0 on the all-zero true column); ours is the mean of the
    # defined cell types - a documented difference, not a match.
    t, q = np.asarray(props["absent_type"]["true"]), np.asarray(props["absent_type"]["pred"])
    assert exp["absent_type"]["jsd"] is None and not t[:, 4].any()
    defined = np.mean([jensenshannon(t[:, k], q[:, k]) for k in range(4)])
    assert em.deconvolution_metrics(true_prop=t, pred_prop=q)["jsd"] == pytest.approx(defined, rel=1e-12)
    # Davis & Goadrich (2006) non-linear PR interpolation on tied scores (dTP = 2, dFP = 2):
    # intermediate point (TP=1, FP=1) has rec=0.5, prec=0.5; end (TP=2, FP=2) has rec=1.0, prec=0.5 -> AUC = 0.5
    true_tied = np.array([[1.0, 1.0], [0.0, 0.0]])
    pred_tied = np.array([[0.5, 0.5], [0.5, 0.5]])
    tied = em.deconvolution_metrics(true_prop=true_tied, pred_prop=pred_tied, aupr_method="precrec")
    assert tied["aupr"] == pytest.approx(0.5, rel=1e-12) and tied["aupr_precrec"] == pytest.approx(0.5, rel=1e-12)
    with pytest.raises(ValueError, match="aupr_method"):
        em.deconvolution_metrics(true_prop=true_tied, pred_prop=pred_tied, aupr_method="bad")

