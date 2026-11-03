"""Tests for fidelity.py - marker-program fidelity (AUC / Cohen's d) and marker-gene overlap."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


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
    from spatial_anno_metrics import fidelity as em

    markers = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_program_fidelity(_marker_adata(correct=True), "cell_type", markers)
    bad = em.marker_program_fidelity(_marker_adata(correct=False), "cell_type", markers)
    assert good["mean_auc"] > 0.9
    assert good["mean_auc"] > bad["mean_auc"]
    assert good["per_type"]["T cell"]["auc"] > 0.9
    assert good["per_type"]["T cell"]["cohens_d"] > 0.5           # large, positive separation


def test_marker_gene_overlap_high_when_labels_match_markers():
    from spatial_anno_metrics import fidelity as em

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
