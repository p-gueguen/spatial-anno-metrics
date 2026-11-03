"""Tests for analysis/signal_qc.py - annotation-independent signal/QC metrics (SpatialQM port).

Fixtures build a section where the panel genes are spatially structured + well above the
negative-control probes, so each metric must separate real signal from background noise.
"""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")


def _spatial_adata(n=300, n_gene=6, n_ctrl=4, seed=0):
    """Panel genes spatially structured (each high in a spatial band) + low random controls."""
    rng = np.random.default_rng(seed)
    coords = rng.uniform(0, 100, size=(n, 2))
    genes = rng.poisson(0.5, size=(n, n_gene)).astype("float32")
    for g in range(n_gene):
        genes[coords[:, 0] < (g + 1) * 100 / (n_gene + 1), g] += 8.0   # spatial gradient
    ctrl = rng.poisson(0.1, size=(n, n_ctrl)).astype("float32")        # random, low
    X = np.hstack([genes, ctrl])
    a = anndata.AnnData(X=sp.csr_matrix(X))
    a.var_names = [f"g{i}" for i in range(n_gene)] + [f"NegControl{i}" for i in range(n_ctrl)]
    a.var["control"] = [False] * n_gene + [True] * n_ctrl
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["spatial"] = coords
    return a


def test_moran_signal_separates_genes_from_controls():
    pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    out = sq.moran_signal(_spatial_adata())
    assert out["status"] == "ok" and out["matrix"] == "X"
    assert out["n_genes"] == 6 and out["n_controls"] == 4
    # real genes are spatially autocorrelated; controls are spatially random.
    assert out["median_gene_moran"] > out["median_control_moran"]
    assert out["separation"] > 0.1
    assert out["frac_genes_above_control_p95"] > 0.8


def test_moran_signal_survives_without_controls():
    pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.var["control"] = False                                   # controls dropped from the panel
    out = sq.moran_signal(a)
    assert out["n_controls"] == 0
    assert out["median_control_moran"] is None
    assert not np.isnan(out["median_gene_moran"])


def test_signal_to_noise_positive_when_signal_above_control():
    from spatial_anno_metrics import signal_qc as sq

    out = sq.signal_to_noise(_spatial_adata())
    assert out["status"] == "ok" and out["matrix"] == "X"
    assert out["median_mean_signal_ratio"] > 0                 # genes > controls (log2 ratio)
    assert out["median_max_signal_ratio"] > out["median_mean_signal_ratio"]
    a = _spatial_adata()
    sq.signal_to_noise(a)
    assert {"sam_mean_signal_ratio", "sam_max_signal_ratio"} <= set(a.obs)
    sq.signal_to_noise(a, key_added=("mean_r", "max_r"))
    assert {"mean_r", "max_r"} <= set(a.obs)


def test_signal_to_noise_skips_without_controls():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.var["control"] = False
    out = sq.signal_to_noise(a)
    assert out["status"] == "skipped" and "control" in out["reason"]


def test_sparsity():
    from spatial_anno_metrics import signal_qc as sq

    X = np.zeros((100, 10), dtype="float32")
    X[:, :5] = 1.0                                             # exactly half the entries nonzero
    a = anndata.AnnData(X=sp.csr_matrix(X))
    out = sq.sparsity(a)
    assert out["zero_fraction"] == pytest.approx(0.5)
    assert "sam_sparsity" in a.obs and a.uns["sam"]["sparsity"]["zero_fraction"] == pytest.approx(0.5)


def test_detection_entropy_low_for_dominance_high_for_uniform():
    from spatial_anno_metrics import signal_qc as sq

    X = np.zeros((2, 10), dtype="float32")
    X[0, 0] = 10.0                                             # one gene dominates -> entropy 0
    X[1, :5] = 1.0                                             # uniform over 5 genes -> entropy 1
    a = anndata.AnnData(X=sp.csr_matrix(X))
    sq.detection_entropy(a, key_added="ent")
    assert a.obs["ent"].iloc[0] == pytest.approx(0.0, abs=1e-9)
    assert a.obs["ent"].iloc[1] == pytest.approx(1.0, abs=1e-6)


def test_tx_per_area():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.obs["cell_area"] = 10.0
    out = sq.tx_per_area(a)
    assert out["median_tx_per_area"] > 0
    assert "sam_tx_per_area" in a.obs


def test_run_signal_qc_headline():
    pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.obs["cell_area"] = 10.0
    head = sq.run_signal_qc(a)
    assert head["status"] == "ok"
    for key in ("sparsity", "detection_entropy", "signal_to_noise", "tx_per_area", "moran"):
        assert head[key]["status"] == "ok", (key, head[key])      # a skipped moran fails here


def test_run_signal_qc_reports_skipped_parts():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    del a.obsm["spatial"]
    head = sq.run_signal_qc(a)                                     # no cell_area, no coordinates
    assert head["tx_per_area"]["status"] == "skipped" and head["moran"]["status"] == "skipped"
    assert head["sparsity"]["status"] == "ok"


def test_signal_qc_refuses_log_X_without_a_counts_layer():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.X = a.X.log1p()
    with pytest.raises(ValueError, match="layer="):
        sq.sparsity(a)
    a.layers["counts"] = a.X.expm1()
    a.layers["counts"].data = np.round(a.layers["counts"].data)
    assert sq.sparsity(a)["matrix"] == "layers['counts']"
