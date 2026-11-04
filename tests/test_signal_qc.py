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
    for method in ("squidpy", "spatialqm"):
        assert sq.moran_signal(a, method=method)["status"] == "skipped"


def test_signal_qc_refuses_log_X_without_a_counts_layer():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.X = a.X.log1p()
    with pytest.raises(ValueError, match="layer="):
        sq.sparsity(a)
    a.layers["counts"] = a.X.expm1()
    a.layers["counts"].data = np.round(a.layers["counts"].data)
    assert sq.sparsity(a)["matrix"] == "layers['counts']"


# --------------------------------------------------------------------------- #
# 0.7: control rule, NaN-safe Moran, SpatialQM recipes, eps-free SNR, #21
# --------------------------------------------------------------------------- #
def test_controls_follow_feature_type_not_the_name():
    """var['feature_type'] flags controls when var['control'] is absent; a name never does."""
    from spatial_anno_metrics.signal_qc import _control_mask

    a = _spatial_adata()
    del a.var["control"]
    a.var_names = ["g0", "NegControlLookalike", "g2", "g3", "g4", "g5", "antisense_X", "BLANK_1", "c2", "c3"]
    a.var["feature_type"] = (["Gene Expression"] * 6 + ["Negative Control Probe", "Blank Codeword",
                                                        "Negative Control Probe", "Negative Control Codeword"])
    assert _control_mask(a).tolist() == [False] * 6 + [True] * 4
    assert _control_mask(a, probes_only=True).tolist() == [False] * 6 + [True, False, True, False]
    a.var["control"] = [False] * 9 + [True]                    # an explicit flag wins over the type
    assert _control_mask(a).sum() == 1
    b = _spatial_adata()
    del b.var["control"]
    assert _control_mask(b).sum() == 0                         # 'NegControl*' names alone flag nothing


def test_moran_signal_is_nan_safe_and_counts_excluded_features():
    pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    X = a.X.toarray()
    X[:, [0, 7]] = 0                                           # one gene and one control with no count
    a.X = sp.csr_matrix(X)
    out = sq.moran_signal(a)
    assert np.isfinite([out["median_gene_moran"], out["median_control_moran"], out["separation"]]).all()
    assert (out["n_genes"], out["n_genes_excluded"], out["n_controls"], out["n_controls_excluded"]) == (6, 1, 4, 1)
    assert out["per_feature"].isna().sum() == 2


@pytest.mark.parametrize("method", ["squidpy", "spatialqm"])
def test_moran_shuffled_coordinates_null(method):
    """Shuffling the coordinates keeps every expression value and removes the spatial structure: the
    gene-vs-control separation must collapse."""
    if method == "squidpy":
        pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata(n=600)
    real = sq.moran_signal(a, method=method)["separation"]
    b = _spatial_adata(n=600)
    b.obsm["spatial"] = np.random.default_rng(1).permutation(b.obsm["spatial"])
    null = sq.moran_signal(b, method=method)["separation"]
    assert real > 0.5 and abs(null) < 0.05, (real, null)


def test_moran_spatialqm_needs_no_squidpy_and_writes_no_graph():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    out = sq.moran_signal(a, method="spatialqm")
    assert out["status"] == "ok" and out["method"] == "spatialqm" and out["separation"] > 0.3
    assert "spatial_connectivities" not in a.obsp and "moranI" not in a.uns
    with pytest.raises(ValueError, match="method"):
        sq.moran_signal(a, method="voyager")


def test_moran_spatialqm_skips_duplicated_coordinates():
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    a.obsm["spatial"][1] = a.obsm["spatial"][0]
    out = sq.moran_signal(a, method="spatialqm")
    assert out["status"] == "skipped" and "coordinates" in out["reason"]


def _two_cell_snr():
    """Cell 0: genes 4 and 2, control 1 -> mean 3 vs 1. Cell 1: genes 2, 0, no control count."""
    X = np.array([[4.0, 2.0, 1.0], [2.0, 0.0, 0.0]], dtype="float32")
    a = anndata.AnnData(X=sp.csr_matrix(X))
    a.var["control"] = [False, False, True]
    return a


def test_signal_to_noise_has_no_eps_floor():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()
    out = sq.signal_to_noise(a)
    assert out["n_cells_with_controls"] == 1 and out["frac_cells_with_controls"] == 0.5
    assert out["median_mean_signal_ratio"] == pytest.approx(np.log2(3.0))
    assert out["median_max_signal_ratio"] == pytest.approx(np.log2(4.0))
    assert out["pooled_mean_signal_ratio"] == pytest.approx(np.log2((8 / 4) / (1 / 2)))   # all cells pooled
    assert out["frac_cells_signal_below_control"] == 0.0
    assert np.isnan(a.obs["sam_mean_signal_ratio"].iloc[1])     # undefined, not log2(1/eps)


def test_signal_to_noise_skips_when_no_cell_has_a_control_count():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()
    a.X = sp.csr_matrix(np.array([[4.0, 2.0, 0.0], [2.0, 0.0, 0.0]], dtype="float32"))
    out = sq.signal_to_noise(a)
    assert out["status"] == "skipped" and out["n_cells_with_controls"] == 0
    assert "sam_mean_signal_ratio" not in a.obs


def test_signal_to_noise_old_eps_reproduces_0_6():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()
    out = sq.signal_to_noise(a, eps=1e-9)
    old_mean = np.log2((np.array([3.0, 1.0]) + 1e-9) / (np.array([1.0, 0.0]) + 1e-9))
    old_max = np.log2((np.array([4.0, 2.0]) + 1e-9) / (np.array([1.0, 0.0]) + 1e-9))
    assert out["median_mean_signal_ratio"] == pytest.approx(float(np.median(old_mean)), rel=1e-15)
    assert out["median_max_signal_ratio"] == pytest.approx(float(np.median(old_max)), rel=1e-15)
    assert np.allclose(a.obs["sam_mean_signal_ratio"], old_mean, rtol=0, atol=0)


def test_signal_to_noise_spatialqm_by_hand():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()                                         # gene means 3, 1; control mean 0.5
    out = sq.signal_to_noise(a, method="spatialqm")
    want = np.mean(np.log10([3.1, 1.1])) - np.log10(0.6)
    assert out["mean_signal_ratio_log10"] == pytest.approx(want)
    assert out["max_signal_ratio_log10"] == pytest.approx(np.log10(3.0) - np.log10(0.5))
    assert "sam_mean_signal_ratio" not in a.obs


def test_detection_entropy_spatialqm_by_hand():
    from spatial_anno_metrics import signal_qc as sq

    a = anndata.AnnData(X=sp.csr_matrix(np.array([[0.0, 1.0], [1.0, 2.0]], dtype="float32")))
    out = sq.detection_entropy(a, method="spatialqm")             # values 0, 1, 1, 2 -> p = 1/4, 1/2, 1/4
    assert out["entropy_bits"] == pytest.approx(1.5) and out["n_distinct_values"] == 3


def test_tx_per_area_reports_mean_and_median():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()
    a.obs["cell_area"] = [1.0, 2.0]                             # totals 7 and 2 -> 7 and 1
    out = sq.tx_per_area(a)
    assert out["mean_tx_per_area"] == pytest.approx(4.0) and out["median_tx_per_area"] == pytest.approx(4.0)


def test_tx_per_cell_counts_genes_only():
    from spatial_anno_metrics import signal_qc as sq

    a = _two_cell_snr()
    out = sq.tx_per_cell(a)
    assert out["mean_tx_per_cell"] == pytest.approx(4.0) and out["median_tx_per_cell"] == pytest.approx(4.0)
    assert out["n_genes"] == 2 and a.obs["sam_tx_per_cell"].tolist() == [6.0, 2.0]


def test_transcript_qc_by_hand():
    import pandas as pd

    from spatial_anno_metrics import signal_qc as sq

    tx = pd.DataFrame({"feature_name": ["A", "A", "B", "NEG1", "NEG2", "A", "C", "A"],
                       "cell_id": [1, 1, 2, -1, 3, 3, 4, -1]})
    out = sq.transcript_qc(tx, ["NEG1", "NEG2", "NEG3"])
    assert out["frac_tx_in_cells"] == pytest.approx(6 / 8)
    assert out["fdr"] == pytest.approx(2 / 8 * 3 / 2)            # 3 genes and 2 negatives detected
    assert (out["n_genes_detected"], out["n_negatives_detected"]) == (3, 2)
    none = sq.transcript_qc(tx[~tx.feature_name.str.startswith("NEG")], ["NEG1"])
    assert none["fdr"] is None and none["n_negative_tx"] == 0   # undefined, not a made-up 0
    assert sq.transcript_qc(tx.iloc[:0], ["NEG1"])["status"] == "skipped"


@pytest.mark.parametrize("n_ctrl_cells", [0, 1, 2, 20])
def test_moran_spatialqm_small_control_group_is_nan_not_a_crash(n_ctrl_cells):
    """spdep's knearneigh(k = 20) needs more than 20 points, so a group with <= 20 cells holding a count has no
    defined I: its features are NaN and counted as excluded, and the gene result survives (0 / 1 cell crashed or
    skipped everything; 2 cells gave the degenerate I = -1)."""
    from spatial_anno_metrics import signal_qc as sq

    a = _spatial_adata()
    X = a.X.toarray()
    X[n_ctrl_cells:, 6:] = 0
    X[:n_ctrl_cells, 6:] = 1 + (np.arange(n_ctrl_cells)[:, None] + np.arange(4)) % 3   # varying counts
    a.X = sp.csr_matrix(X)
    a.var["feature_type"] = ["Gene Expression"] * 6 + ["Negative Control Probe"] * 4
    out = sq.moran_signal(a, method="spatialqm")
    assert out["status"] == "ok" and np.isfinite(out["median_gene_moran"])
    assert out["n_controls"] == out["n_controls_excluded"] == 4 and out["median_control_moran"] is None
    assert out["per_feature"].iloc[6:].isna().all() and out["n_cells_controls"] == n_ctrl_cells


def test_moran_spatialqm_feature_subset_matches_full_run():
    """``genes=`` picks what is scored, not the cells or library sizes the group is normalised on."""
    from spatial_anno_metrics import signal_qc as sq

    full = sq.moran_signal(_spatial_adata(), method="spatialqm")
    sub = ["g0", "g3", "NegControl1"]
    part = sq.moran_signal(_spatial_adata(), method="spatialqm", genes=sub)
    assert list(part["per_feature"].index) == sub
    assert np.allclose(part["per_feature"], full["per_feature"][sub], rtol=0, atol=1e-12)
    assert (part["n_cells_genes"], part["n_cells_controls"]) == (full["n_cells_genes"], full["n_cells_controls"])


def test_transcript_qc_reads_bytes_and_string_ids():
    """Xenium transcripts.parquet stores feature_name as parquet binary: pandas gives bytes."""
    import pandas as pd

    from spatial_anno_metrics import signal_qc as sq

    names = ["A", "A", "B", "NEG1", "NEG2", "A", "C", "A"]
    ids = ["1", "1", "2", "-1", "3", "3", "4", "UNASSIGNED"]
    want = sq.transcript_qc(pd.DataFrame({"feature_name": names, "cell_id": [1, 1, 2, -1, 3, 3, 4, -1]}),
                            ["NEG1", "NEG2"])
    for tx in (pd.DataFrame({"feature_name": names, "cell_id": ids}),                       # str '-1' too
               pd.DataFrame({"feature_name": [n.encode() for n in names], "cell_id": [i.encode() for i in ids]})):
        out = sq.transcript_qc(tx, ["NEG1", "NEG2"])
        assert out["n_unassigned"] == 2 and out["n_negative_tx"] == 2
        assert out["fdr"] == pytest.approx(want["fdr"]) and out["fdr"] == pytest.approx(2 / 8 * 3 / 2)


def test_signal_qc_stored_zeros_and_edge_skips():
    import pandas as pd

    from spatial_anno_metrics import signal_qc as sq

    # Non-canonical CSR with explicit stored zeros (0.0 in X.data) must not produce NaN entropy or wrong zero_fraction
    dense = np.array([[0.0, 2.0, 2.0], [1.0, 0.0, 3.0]], dtype="float32")
    csr_canon = sp.csr_matrix(dense)
    csr_zeros = sp.csr_matrix((np.array([0.0, 2.0, 2.0, 1.0, 0.0, 3.0], dtype="float32"),
                               np.array([0, 1, 2, 0, 1, 2]), np.array([0, 3, 6])), shape=(2, 3))
    assert csr_zeros.nnz == 6 and csr_canon.nnz == 4
    a1, a2 = anndata.AnnData(X=csr_canon), anndata.AnnData(X=csr_zeros)
    assert sq.sparsity(a1)["zero_fraction"] == pytest.approx(sq.sparsity(a2)["zero_fraction"]) == pytest.approx(2 / 6)
    e1, e2 = sq.detection_entropy(a1), sq.detection_entropy(a2)
    assert np.isfinite(e2["median_detection_entropy"])
    assert e1["median_detection_entropy"] == pytest.approx(e2["median_detection_entropy"])

    # All-control panel skips signal_to_noise and tx_per_cell
    a_ctrl = anndata.AnnData(X=csr_canon.copy())
    a_ctrl.var["control"] = True
    assert sq.signal_to_noise(a_ctrl)["status"] == "skipped"
    assert sq.tx_per_cell(a_ctrl)["status"] == "skipped"

    # Non-positive cell_area skips tx_per_area
    a_area = anndata.AnnData(X=csr_canon.copy())
    a_area.obs["cell_area"] = [0.0, -1.0]
    assert sq.tx_per_area(a_area)["status"] == "skipped"

    # transcript_qc with leading None and bytes in subsequent rows
    tx_na = pd.DataFrame({"feature_name": [None, b"A", b"NEG1"], "cell_id": [None, b"c1", b"UNASSIGNED"]})
    out_na = sq.transcript_qc(tx_na, ["NEG1"])
    assert out_na["n_unassigned"] == 2 and out_na["n_negative_tx"] == 1
