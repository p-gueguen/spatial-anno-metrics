"""Tests for the v0.3 eval_metrics additions: Davies-Bouldin + Calinski-Harabasz in internal_validity,
marker_gene_overlap (reference-free fidelity), and avg_bio (external structure conservation)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _blobs(n_per=60, k=4, dim=10, sep=8.0, seed=0):
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, sep, size=(k, dim))
    X, lab = [], []
    for t in range(k):
        X.append(rng.normal(centers[t], 1.0, size=(n_per, dim)))
        lab += [f"T{t}"] * n_per
    a = anndata.AnnData(X=np.vstack(X).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obsm["X_pca"] = a.X
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def _shuffle(a, seed=1):
    b = a.copy()
    b.obs["cell_type"] = pd.Categorical(np.random.default_rng(seed).permutation(b.obs["cell_type"].to_numpy()))
    return b


def test_internal_validity_reports_db_ch_separating_clean_from_shuffled():
    from spatial_anno_metrics import eval_metrics as em

    clean = em.internal_validity(_blobs(), embedding="X_pca", k=15)
    shuf = em.internal_validity(_shuffle(_blobs()), embedding="X_pca", k=15)
    assert "davies_bouldin" in clean and "calinski_harabasz" in clean
    assert clean["davies_bouldin"] < shuf["davies_bouldin"]        # lower = better separated
    assert clean["calinski_harabasz"] > shuf["calinski_harabasz"]  # higher = better separated


def test_marker_gene_overlap_high_when_labels_match_markers():
    from spatial_anno_metrics import eval_metrics as em

    rng = np.random.default_rng(3)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"BG{i}" for i in range(30)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 160
    X = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:80, gi[g]] += 15
    for g in ("MLANA", "SOX10", "DCT"):
        X[80:, gi[g]] += 15
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T"] * 80 + ["Mal"] * 80)
    ms = {"T": ["CD3D", "CD3E", "TRAC"], "Mal": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_gene_overlap(a, "cell_type", ms, n_top=10)
    assert good["mean_overlap"] == pytest.approx(1.0)              # each type's markers ARE its top DE genes
    b = a.copy()
    b.obs["cell_type"] = pd.Categorical(rng.permutation(b.obs["cell_type"].to_numpy()))
    assert em.marker_gene_overlap(b, "cell_type", ms, n_top=10)["mean_overlap"] < good["mean_overlap"]


def test_avg_bio_perfect_for_matching_labels():
    from spatial_anno_metrics import eval_metrics as em

    a = _blobs()
    a.obs["truth"] = a.obs["cell_type"]
    perfect = em.avg_bio(a, "cell_type", "truth", embedding="X_pca")
    assert perfect["ari"] == pytest.approx(1.0) and perfect["nmi"] == pytest.approx(1.0)
    assert perfect["avg_bio"] > em.avg_bio(_shuffle(a), "cell_type", "truth", embedding="X_pca")["avg_bio"]
