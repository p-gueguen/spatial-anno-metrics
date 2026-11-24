"""Core-dependency tests for squidpy-free spatial metrics (v0.8.0)."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _structured_expression_blocks(n_side: int = 12, seed: int = 0):
    """3-band tissue (A, B, C) with distinct spatially localized gene programs + 1 NegCtrl probe."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(n_side), np.arange(n_side))
    xy = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    x = xy[:, 0]
    labels = np.where(x < n_side // 3, "A", np.where(x < 2 * n_side // 3, "B", "C"))
    n = len(labels)
    X = np.zeros((n, 7), dtype=np.float32)
    X[:, 0:2] = rng.poisson(np.where(labels[:, None] == "A", 20.0, 0.2))
    X[:, 2:4] = rng.poisson(np.where(labels[:, None] == "B", 20.0, 0.2))
    X[:, 4:6] = rng.poisson(np.where(labels[:, None] == "C", 20.0, 0.2))
    X[:, 6] = rng.poisson(0.1, size=n)
    a = anndata.AnnData(
        X=X,
        obs=pd.DataFrame({"cell_type": labels}, index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(
            {"feature_type": ["Gene Expression"] * 6 + ["Negative Control Probe"]},
            index=[f"g{i}" for i in range(6)] + ["NegCtrl0"],
        ),
        obsm={"spatial": xy},
    )
    a.layers["counts"] = X.copy()
    return a


def test_residual_spatial_autocorrelation_catches_merge_and_ignores_coords_only_error():
    from spatial_anno_metrics import spatial as sp

    a = _structured_expression_blocks(n_side=12, seed=0)
    clean = sp.residual_spatial_autocorrelation(a)
    assert clean["status"] == "ok"
    assert clean["n_genes"] == 6
    assert "NegCtrl0" not in clean["per_gene"]
    assert clean["raw_median_moran"] > 0.5
    assert clean["residual_median_moran"] < 0.1
    assert clean["explained_fraction"] > 0.85

    # Merging B and C inflates spatial coherence, yet leaves high within-label residual Moran's I
    merged = a.copy()
    merged.obs["cell_type"] = np.where(merged.obs["cell_type"].isin(["B", "C"]), "B+C", "A")
    m_res = sp.residual_spatial_autocorrelation(merged)
    assert m_res["status"] == "ok"
    assert m_res["residual_median_moran"] > clean["residual_median_moran"] + 0.3
    assert m_res["explained_fraction"] < clean["explained_fraction"] - 0.3
    assert m_res["per_type_residual_moran"]["B+C"] > m_res["per_type_residual_moran"]["A"] + 0.3

    # Shuffling spatial coordinates collapses raw spatial autocorrelation
    shuf = a.copy()
    shuf.obsm["spatial"] = shuf.obsm["spatial"][np.random.default_rng(1).permutation(shuf.n_obs)]
    s_res = sp.residual_spatial_autocorrelation(shuf)
    assert abs(s_res["raw_median_moran"]) < 0.1

    # Missing spatial coords skips cleanly
    no_xy = a.copy()
    del no_xy.obsm["spatial"]
    assert sp.residual_spatial_autocorrelation(no_xy)["status"] == "skipped"


def test_spatial_organization_consistency_across_replicates():
    from spatial_anno_metrics import spatial as sp

    rng = np.random.default_rng(7)
    parts = []
    # 3 biological samples with reproducible niche architecture:
    # A always touches B (adjacent strips), C is separate
    for s_idx in range(3):
        b = _structured_expression_blocks(n_side=12, seed=s_idx)
        b.obs["sample"] = f"S{s_idx}"
        b.obs_names = [f"S{s_idx}_{c}" for c in b.obs_names]
        parts.append(b)
    cohort = anndata.concat(parts, join="outer")

    clean = sp.spatial_organization_consistency(cohort, label_key="cell_type", sample_key="sample", n_neighs=6)
    assert clean["status"] == "ok"
    assert clean["n_samples"] == 3 and clean["n_types"] == 3
    assert clean["mean_corr"] > 0.9
    assert clean["offdiag_mean_corr"] > 0.8

    # Shuffling coordinates in one sample destroys cross-sample spatial organization consistency
    bad = cohort.copy()
    m_s0 = (bad.obs["sample"] == "S0").to_numpy()
    xy = bad.obsm["spatial"].copy()
    xy[m_s0] = xy[m_s0][rng.permutation(m_s0.sum())]
    bad.obsm["spatial"] = xy
    shuf = sp.spatial_organization_consistency(bad, label_key="cell_type", sample_key="sample", n_neighs=6)
    assert shuf["mean_corr"] < clean["mean_corr"] - 0.3

    # Single sample skips cleanly
    one = cohort[cohort.obs["sample"] == "S0"].copy()
    assert sp.spatial_organization_consistency(one, label_key="cell_type", sample_key="sample")["status"] == "skipped"
