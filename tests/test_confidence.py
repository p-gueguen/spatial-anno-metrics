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
    assert clean["status"] == "ok" and clean["matrix"] == "layers['counts']"
    assert clean["mean_margin"] > amb["mean_margin"]        # clean cells: one lineage dominates
    assert clean["mean_entropy"] < amb["mean_entropy"]      # ambiguous cells: flatter softmax


def test_per_cell_confidence_writes_two_named_columns():
    from spatial_anno_metrics import confidence as cf

    a = _section()
    cf.per_cell_confidence(a, MARKERS)
    assert {"sam_marker_margin", "sam_marker_entropy"} <= set(a.obs)
    cf.per_cell_confidence(a, MARKERS, key_added=("m", "e"))
    assert {"m", "e"} <= set(a.obs) and a.uns["sam"]["per_cell_confidence"]["key_added"] == ["m", "e"]
    with pytest.raises(TypeError, match="two column names"):
        cf.per_cell_confidence(a, MARKERS, key_added="me")


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
    assert "sam_label_stability" in a.obs and a.uns["sam"]["label_stability"]["reps"] == 4
    assert 0.0 <= out["mean_flip_rate"] <= 1.0
    assert out["mean_flip_rate"] < 0.1                     # well-separated markers -> stable under 20% dropout
    assert out["status"] == "ok" and out["matrix"] == "layers['counts']"


def test_per_cell_confidence_skips_with_fewer_than_two_marker_sets():
    """One scorable set has no runner-up: no margin to report (was a fabricated margin 0 / entropy 1)."""
    from spatial_anno_metrics import confidence as cf

    a = _section()
    out = cf.per_cell_confidence(a, {"T": MARKERS["T"], "Ghost": ["not_on_panel"]})
    assert out["status"] == "skipped" and "found 1" in out["reason"]
    assert "sam_marker_margin" not in a.obs and "per_cell_confidence" not in a.uns.get("sam", {})


def test_label_stability_does_not_truncate_float_values():
    """Thinning rounds a non-integer matrix (with a warning) instead of flooring it: 0.6 -> 1, not 0."""
    import warnings

    from spatial_anno_metrics import confidence as cf

    a = _section(n=20)
    del a.layers["counts"]
    a.X = np.full(a.X.shape, 0.6, dtype="float32")
    seen = []

    def total(x):
        seen.append(float(np.asarray(x.X).sum()))
        return np.array(["T"] * x.n_obs)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        cf.label_stability(a, total, drop_frac=1e-9, reps=1, layer=None)
    assert seen[1] == a.n_obs * a.n_vars                    # every 0.6 kept as 1 count (astype(int) gave 0)
    assert sum(issubclass(x.category, UserWarning) for x in w) == 1
