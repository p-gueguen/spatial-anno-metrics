"""Tests for internal.py - reference-free internal validity (scTypeEval family): clean blobs must score above
shuffled labels."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")


def test_internal_validity_high_for_clean_low_for_shuffled(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    a = blobs()
    clean = em.internal_validity(a, label_key="cell_type", embedding="X_pca", n_neighbors=15)
    shuf = em.internal_validity(shuffle_labels(a), label_key="cell_type", embedding="X_pca", n_neighbors=15)

    for m in ("silhouette", "silhouette_2label", "neighborhood_purity",
              "orbital_medoid", "ward_propmatch", "integrated"):
        assert clean[m] > shuf[m], (m, clean[m], shuf[m])
    assert clean["neighborhood_purity"] > 0.9      # clean blobs: neighbors share label
    assert 0.0 <= clean["integrated"] <= 1.0
    assert clean["status"] == "ok"


def test_internal_validity_computes_pca_when_no_embedding(blobs):
    from spatial_anno_metrics import internal as em

    a = blobs()
    del a.obsm["X_pca"]                              # force the PCA fallback
    out = em.internal_validity(a, label_key="cell_type", embedding="X_pca", n_neighbors=15)
    assert out["neighborhood_purity"] > 0.9


def test_internal_validity_reports_db_ch_separating_clean_from_shuffled(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    clean = em.internal_validity(blobs(n_per=60), embedding="X_pca", n_neighbors=15)
    shuf = em.internal_validity(shuffle_labels(blobs(n_per=60)), embedding="X_pca", n_neighbors=15)
    assert "davies_bouldin" in clean and "calinski_harabasz" in clean
    assert clean["davies_bouldin"] < shuf["davies_bouldin"]        # lower = better separated
    assert clean["calinski_harabasz"] > shuf["calinski_harabasz"]  # higher = better separated


def test_internal_validity_rejects_missing_labels_and_skips_one_label(blobs):
    from spatial_anno_metrics import internal as em

    a = blobs(n_per=20)
    lab = a.obs["cell_type"].astype(str).to_numpy().astype(object)
    lab[:3] = None
    a.obs["cell_type"] = lab
    with pytest.raises(ValueError, match=r"obs\['cell_type'\] has 3 missing"):
        em.internal_validity(a, embedding="X_pca")
    a.obs["cell_type"] = "T0"
    out = em.internal_validity(a, embedding="X_pca")
    assert out["status"] == "skipped" and "found 1" in out["reason"]


def test_internal_validity_reports_avg_similarity_clean_above_shuffled(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    a = blobs()
    clean = em.internal_validity(a, embedding="X_pca", n_neighbors=15)
    shuf = em.internal_validity(shuffle_labels(a), embedding="X_pca", n_neighbors=15)
    assert -1.0 <= shuf["avg_similarity"] < clean["avg_similarity"] <= 1.0
    assert clean["avg_similarity"] > 0.9            # tight blobs far from the origin: cells point at their centroid


def test_internal_validity_runs_ward_on_the_subsample_only(blobs, monkeypatch):
    """Ward is O(N^2) memory: it must see the subsample, never every cell (40k cells took 13 GB)."""
    import sklearn.cluster

    from spatial_anno_metrics import internal as em

    seen = []
    real = sklearn.cluster.AgglomerativeClustering

    class Spy(real):
        def fit_predict(self, X, y=None):
            seen.append(len(X))
            return super().fit_predict(X)

    monkeypatch.setattr(sklearn.cluster, "AgglomerativeClustering", Spy)
    out = em.internal_validity(blobs(n_per=100), embedding="X_pca", subsample=150)
    assert seen == [150] and out["n_used"] == 150
    assert 0.0 <= out["ward_propmatch"] <= 1.0


def test_internal_validity_silhouettes_are_nan_with_a_reason_when_undefined(blobs):
    """One label in the subsample, or only singletons: the silhouettes are undefined - NaN plus a reason, never a
    fabricated 0.0 (and no sklearn ValueError for n_labels == n_samples)."""
    import numpy as np

    from spatial_anno_metrics import internal as em

    a = blobs(n_per=10, k=2)
    a.obs["cell_type"] = [f"s{i}" for i in range(a.n_obs)]         # every cell its own label
    out = em.internal_validity(a, embedding="X_pca", n_neighbors=5)
    assert out["status"] == "ok"
    for m in ("silhouette", "silhouette_2label", "ward_propmatch", "davies_bouldin", "calinski_harabasz"):
        assert np.isnan(out[m]) and m in out["nan_reasons"], m
    assert np.isfinite(out["integrated"])                          # geometric mean of the defined parts
    assert np.isfinite(out["neighborhood_purity"]) and np.isfinite(out["avg_similarity"])

    b = blobs(n_per=500, k=2)
    lab = np.array(["T0"] * b.n_obs, dtype=object)
    lab[-1] = "T1"                                                 # one T1 cell, last
    b.obs["cell_type"] = lab
    assert b.n_obs - 1 not in np.random.default_rng(0).choice(b.n_obs, 50, replace=False)   # precondition
    out = em.internal_validity(b, embedding="X_pca", n_neighbors=5, subsample=50, seed=0)
    assert np.isnan(out["silhouette"]) and "1 label" in out["nan_reasons"]["silhouette"]
    assert np.isnan(out["silhouette_2label"])


def test_internal_validity_flags_a_label_trained_embedding(blobs):
    from spatial_anno_metrics import internal as em

    a = blobs(n_per=30)
    assert "warning" not in em.internal_validity(a, embedding="X_pca")
    out = em.internal_validity(a, embedding="X_pca", embedding_label_trained=True)
    assert "circular" in out["warning"] and out["status"] == "ok"


def test_embedding_pca_of_sparse_x_matches_dense_without_densifying(blobs, monkeypatch):
    """Sparse X: PCA by implicit centering (arpack), never a dense copy; same geometry as the dense PCA."""
    import numpy as np
    import scipy.sparse as sp
    from scipy.spatial.distance import pdist

    from spatial_anno_metrics import _utils

    a = blobs(n_per=40, dim=60)
    del a.obsm["X_pca"]
    a.X = np.abs(a.X)                                              # any non-negative matrix
    dense = _utils._embedding(a, None)
    a.X = sp.csr_matrix(a.X)

    def boom(X):
        raise AssertionError("densified")

    monkeypatch.setattr(_utils, "_dense", boom)
    E = _utils._embedding(a, None)
    assert E.shape == dense.shape == (a.n_obs, 50)
    np.testing.assert_allclose(pdist(E), pdist(dense), rtol=1e-6, atol=1e-8)


# --- isolated labels + cLISI (scib; Luecken et al. 2022) -------------------------------------------------------
def _rare_section(blobs, far: bool):
    """3 blobs of 60 cells + 8 cells of a rare label 'R' seen in one batch only: far from every blob, or planted
    inside blob T0. The other labels span both batches, so 'R' is the one isolated label."""
    import anndata
    import numpy as np
    import pandas as pd

    a = blobs(n_per=60, k=3)
    X = np.asarray(a.obsm["X_pca"], dtype=float)
    rng = np.random.default_rng(5)
    centre = X.max(0) + 20.0 if far else X[:60].mean(0)
    Xr = rng.normal(centre, 1.0, size=(8, X.shape[1]))
    E = np.vstack([X, Xr])
    b = anndata.AnnData(X=E.astype("float32"))
    b.obs_names = [f"c{i}" for i in range(len(E))]
    b.obsm["X_pca"] = E
    b.obs["cell_type"] = pd.Categorical(list(a.obs["cell_type"].astype(str)) + ["R"] * 8)
    b.obs["batch"] = ["b0", "b1"] * 90 + ["b1"] * 8
    return b


def test_isolated_label_scores_see_a_buried_rare_label_that_clisi_misses(blobs):
    from spatial_anno_metrics import internal as em

    good = em.isolated_label_scores(_rare_section(blobs, True), batch_key="batch")
    bad = em.isolated_label_scores(_rare_section(blobs, False), batch_key="batch")
    assert good["status"] == "ok" and good["isolated_labels"] == ["R"]
    assert good["asw"] > 0.8 and good["f1"] > 0.95
    assert bad["asw"] < 0.55 and bad["f1"] < 0.5
    assert set(bad["per_label"]) == {"R"} and {"asw", "f1", "n"} <= set(bad["per_label"]["R"])
    # cLISI's median over cells cannot see 8 buried cells out of 188; the overall silhouette (rescaled to [0, 1])
    # drops, but stays well above the rare label's own score
    buried = _rare_section(blobs, False)
    assert em.clisi(buried)["clisi"] > 0.99
    assert bad["asw"] < (em.internal_validity(buried, embedding="X_pca")["silhouette"] + 1) / 2 - 0.2


def test_isolated_label_scores_without_batch_score_every_label(blobs):
    from spatial_anno_metrics import internal as em

    out = em.isolated_label_scores(_rare_section(blobs, True))
    assert out["isolated_labels"] == ["R", "T0", "T1", "T2"]


def test_isolated_label_scores_match_scib_metrics(blobs):
    """scib-metrics 0.5.6 ``isolated_labels`` (BSD-3, run offline, float32 JAX) on the two rare-label sections."""
    from spatial_anno_metrics import internal as em

    for far, tag in ((True, "far"), (False, "buried")):
        out = em.isolated_label_scores(_rare_section(blobs, far), batch_key="batch")
        assert out["asw"] == pytest.approx(SCIB_ISOLATED_ASW[tag], rel=1e-5), tag
    b = blobs(sep=1.5)
    b.obs["batch"] = ["b0", "b1"] * (b.n_obs // 2)                 # every label in both batches: all isolated
    out = em.isolated_label_scores(b, batch_key="batch")
    assert out["asw"] == pytest.approx(SCIB_ISOLATED_ASW["overlap"], rel=1e-5)


def test_clisi_separates_clean_from_shuffled_and_matches_scib_metrics(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    clean = em.clisi(blobs())
    shuf = em.clisi(shuffle_labels(blobs()))
    assert clean["status"] == "ok" and clean["clisi"] > 0.95 and shuf["clisi"] < 0.5
    out = em.clisi(blobs(sep=1.5))                                 # overlapping blobs: a non-trivial LISI
    assert out["clisi"] == pytest.approx(SCIB_CLISI_OVERLAP, rel=1e-5)
    assert out["median_lisi"] == pytest.approx(SCIB_MEDIAN_LISI_OVERLAP, rel=1e-5)
    assert out["n_neighbors"] == 90 and out["perplexity"] == 30.0


def test_isolated_labels_and_clisi_skip_one_label(blobs):
    from spatial_anno_metrics import internal as em

    a = blobs(n_per=20)
    a.obs["cell_type"] = "T0"
    for fn in (em.isolated_label_scores, em.clisi):
        out = fn(a)
        assert out["status"] == "skipped" and "found 1" in out["reason"]


# scib-metrics 0.5.6 (BSD-3; float32 JAX), run offline on these fixtures with the same sklearn 90-NN, self included
SCIB_ISOLATED_ASW = {"far": 0.9777132272720337, "buried": 0.5268670320510864, "overlap": 0.6455773115158081}
SCIB_CLISI_OVERLAP, SCIB_MEDIAN_LISI_OVERLAP = 0.9729409217834473, 1.0811772346496582


def test_internal_edge_cases_constant_embedding_iso_subsample_and_small_clisi(blobs):
    from spatial_anno_metrics import internal as em

    # Constant embedding: neighborhood_purity excludes self and silhouette_2label is 0.0 (not NaN)
    a = blobs(n_per=10, k=2)
    a.obsm["X_pca"][:] = 0.0
    out = em.internal_validity(a, embedding="X_pca", n_neighbors=5)
    assert out["status"] == "ok" and out["silhouette_2label"] == 0.0
    # Each cell has 9 same-label and 10 other-label peers among the 19 non-self cells
    assert 0.0 <= out["neighborhood_purity"] <= 1.0

    # isolated_label_scores per-label subsampling (n > subsample) and iso_threshold skip
    rare = _rare_section(blobs, True)
    sub_out = em.isolated_label_scores(rare, batch_key="batch", subsample=40)
    assert sub_out["status"] == "ok" and sub_out["n_used"] <= 40
    skip_iso = em.isolated_label_scores(rare, batch_key="batch", iso_threshold=0)
    assert skip_iso["status"] == "skipped" and "iso_threshold" not in skip_iso["reason"] or "batch" in skip_iso["reason"]

    # clisi skips when n_obs < 6 with default perplexity, and raises on perplexity <= 1
    small = blobs(n_per=2, k=2)
    assert em.clisi(small)["status"] == "skipped"
    with pytest.raises(ValueError, match="perplexity"):
        em.clisi(blobs(n_per=10, k=2), perplexity=1.0)
