"""Tests for internal.py - reference-free internal validity (scTypeEval family): clean blobs must score above
shuffled labels."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")


def test_internal_validity_high_for_clean_low_for_shuffled(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    a = blobs()
    clean = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    shuf = em.internal_validity(shuffle_labels(a), label_key="cell_type", embedding="X_pca", k=15)

    for m in ("silhouette", "silhouette_2label", "neighborhood_purity",
              "orbital_medoid", "ward_propmatch", "integrated"):
        assert clean[m] > shuf[m], (m, clean[m], shuf[m])
    assert clean["neighborhood_purity"] > 0.9      # clean blobs: neighbors share label
    assert 0.0 <= clean["integrated"] <= 1.0


def test_internal_validity_computes_pca_when_no_embedding(blobs):
    from spatial_anno_metrics import internal as em

    a = blobs()
    del a.obsm["X_pca"]                              # force the PCA fallback
    out = em.internal_validity(a, label_key="cell_type", embedding="X_pca", k=15)
    assert out["neighborhood_purity"] > 0.9


def test_internal_validity_reports_db_ch_separating_clean_from_shuffled(blobs, shuffle_labels):
    from spatial_anno_metrics import internal as em

    clean = em.internal_validity(blobs(n_per=60), embedding="X_pca", k=15)
    shuf = em.internal_validity(shuffle_labels(blobs(n_per=60)), embedding="X_pca", k=15)
    assert "davies_bouldin" in clean and "calinski_harabasz" in clean
    assert clean["davies_bouldin"] < shuf["davies_bouldin"]        # lower = better separated
    assert clean["calinski_harabasz"] > shuf["calinski_harabasz"]  # higher = better separated
