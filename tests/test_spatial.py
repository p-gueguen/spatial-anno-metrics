"""Tests for spatial.py - spatial label-quality metrics. A COHERENT labeling (types occupy distinct
spatial blocks) must score higher coherence / lower PAS than a SHUFFLED one; neighborhood_sanity must
not flag a clean block layout."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("squidpy")


def _blocks(n_side=20, seed=0):
    """A 2*n_side x n_side grid: left half labeled A, right half B (spatially coherent)."""
    rng = np.random.default_rng(seed)
    xs, ys, lab = [], [], []
    for i in range(2 * n_side):
        for j in range(n_side):
            xs.append(i + rng.normal(0, 0.01))
            ys.append(j + rng.normal(0, 0.01))
            lab.append("A" if i < n_side else "B")
    a = anndata.AnnData(X=np.zeros((len(xs), 3), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(len(xs))]
    a.obsm["spatial"] = np.column_stack([xs, ys]).astype(float)
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def test_spatial_coherence_coherent_vs_shuffled():
    from spatial_anno_metrics import spatial as sp

    a = _blocks()
    coh = sp.spatial_coherence(a, label_key="cell_type", k=8)
    assert "spatial_coherence" in a.obs
    assert coh["mean_coherence"] > 0.8 and coh["pas"] < 0.15      # clean blocks: neighbors share label

    b = a.copy()
    b.obs["cell_type"] = pd.Categorical(np.random.default_rng(1).permutation(b.obs["cell_type"].to_numpy()))
    del b.obsp["spatial_connectivities"]                          # rebuild graph for the shuffled copy
    shuf = sp.spatial_coherence(b, label_key="cell_type", k=8)
    assert shuf["mean_coherence"] < coh["mean_coherence"]
    assert shuf["pas"] > coh["pas"]


def test_neighborhood_sanity_clean_layout_not_flagged():
    from spatial_anno_metrics import spatial as sp

    out = sp.neighborhood_sanity(_blocks(), label_key="cell_type", n_neighs=8)
    assert set(out["categories"]) == {"A", "B"}
    # a clean 2-block layout: each type's strongest partner is itself -> nothing suspicious
    assert out["suspicious"] == []
    assert len(out["self_enrichment"]) == 2
