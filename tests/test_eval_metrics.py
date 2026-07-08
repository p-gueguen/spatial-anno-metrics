"""Tests for analysis/eval_metrics.py - the annotation-quality metric battery.

Covers the three regimes from docs/research/cell-annotation-quality-metrics.md:
  * internal (reference-free): scTypeEval cluster-validity family + inter-sample consistency;
  * marker-program fidelity (reference-free, needs marker sets): AUC-ROC / Cohen's d;
  * external (needs ground truth): F1/ARI/kappa, ECS, hierarchical & composition accuracy.

Fixtures are deterministic well-separated blobs (clean) vs shuffled labels (noise); every metric
must score the clean labeling above the shuffled one, which is the whole point of the module.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _blobs(n_per=80, k=4, dim=10, sep=8.0, seed=0):
    """k well-separated Gaussian blobs in a 10-d embedding, correctly labeled T0..T{k-1}."""
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, sep, size=(k, dim))
    X, lab = [], []
    for t in range(k):
        X.append(rng.normal(centers[t], 1.0, size=(n_per, dim)))
        lab += [f"T{t}"] * n_per
    X = np.vstack(X).astype("float32")
    a = anndata.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obsm["X_pca"] = X
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def _shuffle_labels(a, seed=1):
    a2 = a.copy()
    rng = np.random.default_rng(seed)
    a2.obs["cell_type"] = pd.Categorical(rng.permutation(a2.obs["cell_type"].to_numpy()))
    return a2


# --------------------------------------------------------------------------- #
# Internal cluster-validity (scTypeEval family)
# --------------------------------------------------------------------------- #
def test_internal_validity_high_for_clean_low_for_shuffled():
    from spatial_anno_metrics import eval_metrics as em

    a = _blobs()
    clean = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    shuf = em.internal_validity(_shuffle_labels(a), label_key="cell_type", embedding="X_pca", k=15)

    for m in ("silhouette", "silhouette_2label", "neighborhood_purity",
              "orbital_medoid", "ward_propmatch", "integrated"):
        assert clean[m] > shuf[m], (m, clean[m], shuf[m])
    assert clean["neighborhood_purity"] > 0.9      # clean blobs: neighbors share label
    assert 0.0 <= clean["integrated"] <= 1.0


def test_internal_validity_computes_pca_when_no_embedding():
    from spatial_anno_metrics import eval_metrics as em

    a = _blobs()
    del a.obsm["X_pca"]                              # force the PCA fallback
    out = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    assert out["neighborhood_purity"] > 0.9


# --------------------------------------------------------------------------- #
# Inter-sample consistency (ISC)
# --------------------------------------------------------------------------- #
def _isc_adata(consistent: bool, n_per=40, n_types=3, n_samples=4, n_genes=30, seed=0):
    """Each (sample, type) is a pseudobulk. Consistent: a type's signature block is the SAME
    across samples. Inconsistent: each sample uses a DIFFERENT block for a given type."""
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    X, types, samples = [], [], []
    for s in range(n_samples):
        for t in range(n_types):
            sig = t if consistent else rng.integers(0, n_types)  # which block is high
            base = rng.normal(0.2, 0.1, size=(n_per, n_genes))
            base[:, sig * block:(sig + 1) * block] += 5.0
            X.append(base)
            types += [f"T{t}"] * n_per
            samples += [f"S{s}"] * n_per
    a = anndata.AnnData(X=np.vstack(X).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = pd.Categorical(types)
    a.obs["sample"] = pd.Categorical(samples)
    return a


def test_inter_sample_consistency_high_for_reproducible_types():
    from spatial_anno_metrics import eval_metrics as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    assert good["consistency"] > bad["consistency"]


# --------------------------------------------------------------------------- #
# Marker-program fidelity (reference-free, needs marker sets)
# --------------------------------------------------------------------------- #
def _marker_adata(correct=True, seed=0):
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(10)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 120
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:60, gi[g]] += 10
    for g in ("MLANA", "SOX10", "DCT"):
        X[60:, gi[g]] += 10
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    labels = ["T cell"] * 60 + ["Mel"] * 60
    if not correct:
        labels = list(rng.permutation(labels))
    a.obs["cell_type"] = pd.Categorical(labels)
    return a


def test_marker_program_fidelity_high_auc_for_correct_labels():
    from spatial_anno_metrics import eval_metrics as em

    markers = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_program_fidelity(_marker_adata(correct=True), "cell_type", markers)
    bad = em.marker_program_fidelity(_marker_adata(correct=False), "cell_type", markers)
    assert good["mean_auc"] > 0.9
    assert good["mean_auc"] > bad["mean_auc"]
    assert good["per_type"]["T cell"]["auc"] > 0.9
    assert good["per_type"]["T cell"]["cohens_d"] > 0.5           # large, positive separation


# --------------------------------------------------------------------------- #
# External harness (needs ground truth)
# --------------------------------------------------------------------------- #
def test_external_scores_perfect():
    from spatial_anno_metrics import eval_metrics as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    out = em.external_scores(truth.copy(), truth)
    for m in ("accuracy", "balanced_accuracy", "macro_f1", "ari", "ami", "kappa", "ecs"):
        assert out[m] == pytest.approx(1.0), (m, out[m])


def test_external_scores_macro_f1_penalizes_rare_miss():
    from spatial_anno_metrics import eval_metrics as em

    truth = np.array(["A"] * 90 + ["B"] * 10)
    pred = np.array(["A"] * 100)                       # miss every rare-B cell
    out = em.external_scores(pred, truth)
    assert out["accuracy"] == pytest.approx(0.9)       # inflated by the majority class
    assert out["macro_f1"] < out["accuracy"]           # macro-F1 exposes the dropped class
    assert out["balanced_accuracy"] < 0.9
    assert out["per_class_f1"]["B"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# Element-Centric Similarity (ECS)
# --------------------------------------------------------------------------- #
def test_ecs_identical_is_one_and_resolution_aware():
    from spatial_anno_metrics import eval_metrics as em

    a = np.array([0, 0, 1, 1, 2, 2])
    assert em.element_centric_similarity(a, a) == pytest.approx(1.0)
    allsame = np.zeros(6, dtype=int)
    refine = np.array([0, 3, 1, 1, 2, 2])              # split cluster 0 into two singletons
    assert em.element_centric_similarity(a, allsame) < 1.0
    # a near-refinement is more similar to `a` than collapsing everything to one cluster
    assert em.element_centric_similarity(a, refine) > em.element_centric_similarity(a, allsame)


# --------------------------------------------------------------------------- #
# Hierarchical & composition accuracy
# --------------------------------------------------------------------------- #
def test_hierarchical_accuracy_partial_credit():
    from spatial_anno_metrics import eval_metrics as em

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
    from spatial_anno_metrics import eval_metrics as em

    truth = np.array(["A"] * 60 + ["B"] * 40)
    same = em.composition_accuracy(truth.copy(), truth)
    assert same["l1"] == pytest.approx(0.0, abs=1e-9)
    assert same["jsd"] == pytest.approx(0.0, abs=1e-9)
    diff = em.composition_accuracy(np.array(["A"] * 90 + ["B"] * 10), truth)
    assert diff["l1"] > 0


def test_deconvolution_metrics_perfect_and_degraded():
    """OpenProblems spatial-decomposition metrics (R2 uniform-average + JSD axis=0) + RMSE on
    predicted vs true cell-type PROPORTION matrices (spots x types)."""
    from spatial_anno_metrics import eval_metrics as em

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


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def test_annotation_quality_orchestrator():
    from spatial_anno_metrics import eval_metrics as em

    a = _blobs()
    a.var["control"] = False
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca")
    assert "internal_validity" in out
    assert out["internal_validity"]["neighborhood_purity"] > 0.9


def _panel_ref(n=120, n_genes=40, seed=0):
    """Reference where A/B have distinct on-panel markers but C/D differ ONLY off-panel (last gene)."""
    rng = np.random.default_rng(seed)
    off = n_genes - 1
    blocks, labels = [], []
    for t, hi in [("A", 0), ("B", 1), ("C", 2), ("D", 2)]:      # C,D share on-panel marker g2
        M = rng.poisson(1.0, size=(n, n_genes)).astype("float32")
        M[:, hi] += rng.poisson(20.0, size=n)
        if t == "D":
            M[:, off] += rng.poisson(20.0, size=n)              # D's only distinct gene is OFF-panel
        blocks.append(M)
        labels += [t] * n
    a = anndata.AnnData(X=np.vstack(blocks))
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = labels
    return a


def test_panel_resolvability_flags_off_panel_confusable_pair():
    """A panel that omits a type-pair's only distinguishing gene cannot resolve that pair."""
    from spatial_anno_metrics import eval_metrics as em

    a = _panel_ref()
    panel = {f"g{i}" for i in range(39)}                         # 39 shared genes; excludes g39 (D's marker)
    res = em.panel_resolvability(a, "cell_type", panel, target_depth=40)
    assert res["status"] == "ok" and res["n_types"] == 4
    assert res["per_type"]["A"]["tier"] == "resolvable"         # private on-panel markers -> resolvable
    assert res["per_type"]["B"]["tier"] == "resolvable"
    cd = {res["per_type"]["C"]["tier"], res["per_type"]["D"]["tier"]}
    assert "not_resolvable" in cd                                # identical on the panel -> confusable
    assert (res["per_type"]["C"]["confused_with"] == "D"
            or res["per_type"]["D"]["confused_with"] == "C")
    assert res["frac_resolvable"] < 1.0


def test_panel_resolvability_insufficient_overlap():
    """A reference sharing < 25 genes with the panel returns an honest skip, not a bogus score."""
    from spatial_anno_metrics import eval_metrics as em

    a = _panel_ref()
    res = em.panel_resolvability(a, "cell_type", {"zzz1", "zzz2"}, target_depth=40)
    assert res["status"] == "insufficient_overlap"


# --- conformal prediction sets (catalog s3d) ---------------------------------------------
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
    from spatial_anno_metrics import eval_metrics as em

    rng = np.random.default_rng(0)
    classes = np.array(["A", "B", "C", "D"])
    y = rng.choice(classes, 1200)
    P = _peaked_proba(y, classes, peak=0.9)
    cal, te = slice(0, 800), slice(800, 1200)
    r = em.conformal_prediction_sets(P[te], classes, P[cal], y[cal], alpha=0.1, y_query=y[te])
    assert r["coverage"] >= 0.85                      # ~>= 1-alpha (finite-sample slack)
    assert r["summary"]["pct_singleton"] > 0.8        # confident classifier -> mostly singletons
    assert r["summary"]["mean_set_size"] < len(classes)   # not everything is in every set


def test_conformal_class_conditional_lifts_rare_type_coverage():
    """Marginal conformal under-covers a hard rare class; class-conditional restores it."""
    from spatial_anno_metrics import eval_metrics as em

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
    from spatial_anno_metrics import eval_metrics as em

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
