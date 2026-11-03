"""Tests for decomposition.py - decomposition_purity from deconvolution (e.g. RCTD full-mode) weights."""
from __future__ import annotations

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")


def _adata(W):
    a = anndata.AnnData(X=np.zeros((len(W), 1), dtype="float32"))
    a.obsm["rctd_weights"] = np.asarray(W, dtype=float)
    return a


def test_decomposition_purity_values():
    import pandas as pd

    from spatial_anno_metrics.decomposition import decomposition_purity

    W = [[1, 0, 0, 0], [0.25, 0.25, 0.25, 0.25], [2, 1, 1, 0],      # pure, even, unnormalised (max 2 / 4)
         [np.nan, 1, 0, 0], [0, 0, 0, 0], [1, -0.1, 0, 0]]          # undefined: NaN, zero total, negative
    a = _adata(W)
    out = decomposition_purity(a, "rctd_weights")
    v = a.obs["sam_decomposition_purity"].to_numpy()
    assert np.allclose(v[:3], [0.0, 0.75, 0.5]) and np.isnan(v[3:]).all()
    assert (out["n_cells"], out["n_undefined"], out["n_types"]) == (6, 3, 4)
    assert out["median_decomposition_purity"] == pytest.approx(0.5)
    assert out["mean_decomposition_purity"] == pytest.approx(1.25 / 3)
    assert a.uns["sam"]["decomposition_purity"]["key_added"] == "sam_decomposition_purity"
    b = anndata.AnnData(X=np.zeros((2, 1), dtype="float32"))
    b.obsm["w"] = pd.DataFrame([[3.0, 1.0], [1.0, 1.0]], index=b.obs_names, columns=["T", "B"])
    decomposition_purity(b, "w", key_added="mix")
    assert np.allclose(b.obs["mix"], [0.25, 0.5])


def test_decomposition_purity_skips():
    from spatial_anno_metrics.decomposition import decomposition_purity

    a = _adata([[0, 0], [np.nan, 1]])
    out = decomposition_purity(a, "rctd_weights")
    assert out["status"] == "skipped" and "sam_decomposition_purity" not in a.obs   # all undefined
    assert decomposition_purity(a, "missing")["status"] == "skipped"
