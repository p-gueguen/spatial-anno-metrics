"""Tests for confidence.py - per-cell marker margin/entropy + subsampling label stability. Clean,
well-separated cells => high margin, low entropy, low flip rate; ambiguous cells the opposite."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

MARKERS = {"T": ["CD3D", "CD3E", "TRAC"], "Mal": ["MLANA", "SOX10", "DCT"]}
GENES = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"BG{i}" for i in range(20)]
_GI = {g: i for i, g in enumerate(GENES)}


def _section(n=200, ambiguous=False, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.3, size=(n, len(GENES))).astype("float32")
    half = n // 2
    for g in MARKERS["T"]:
        X[:half, _GI[g]] += 12
    for g in MARKERS["Mal"]:
        X[half:, _GI[g]] += 12
    if ambiguous:                                   # give every cell BOTH lineages -> low margin
        for g in MARKERS["T"] + MARKERS["Mal"]:
            X[:, _GI[g]] += 12
    a = anndata.AnnData(X=X)
    a.var_names = GENES
    a.obs_names = [f"c{i}" for i in range(n)]
    a.layers["counts"] = X.copy()
    a.obs["cell_type"] = pd.Categorical(["T"] * half + ["Mal"] * (n - half))
    return a


def test_per_cell_confidence_margin_and_entropy():
    from spatial_anno_metrics import confidence as cf

    clean = cf.per_cell_confidence(_section(), MARKERS)
    amb = cf.per_cell_confidence(_section(ambiguous=True, seed=1), MARKERS)
    assert clean["mean_margin"] > amb["mean_margin"]        # clean cells: one lineage dominates
    assert clean["mean_entropy"] < amb["mean_entropy"]      # ambiguous cells: flatter softmax


def _marker_argmax(adata):
    """A tiny marker-argmax labeller for the stability callback (mean expr per set)."""
    import numpy as np
    panel = set(adata.var_names)
    names, cols = [], []
    for name, genes in MARKERS.items():
        gi = [list(adata.var_names).index(g) for g in genes if g in panel]
        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        names.append(name)
        cols.append(X[:, gi].mean(1))
    S = np.column_stack(cols)
    return np.array(names)[S.argmax(1)]


def test_label_stability_low_flip_for_separated_cells():
    from spatial_anno_metrics import confidence as cf

    a = _section(seed=2)
    out = cf.label_stability(a, _marker_argmax, drop_frac=0.2, reps=4, seed=3)
    assert "label_stability" in a.obs
    assert 0.0 <= out["mean_flip_rate"] <= 1.0
    assert out["mean_flip_rate"] < 0.1                     # well-separated markers -> stable under 20% dropout
