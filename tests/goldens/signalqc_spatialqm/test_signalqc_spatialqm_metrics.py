"""signal_qc against SpatialQM (expected.json, see SOURCE.md) on the Xenium breast crop.

`method="spatialqm"` reproduces getMorans, getMaxRatio and getEntropy to float precision; the mean of
`tx_per_area`, `sparsity`, `tx_per_cell` and `transcript_qc` reproduce getTxPerArea, getSparsity, getTxPerCell,
getCellTxFraction and getGlobalFDR. SpatialQM picks controls by name (its regexes miss the 8 `antisense_*`
negative-control probes); our rule is the feature type, so each test also pins what that rule gives.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

anndata = pytest.importorskip("anndata")

HERE = Path(__file__).parent
E = json.loads((HERE / "expected.json").read_text())
REL = 1e-12                                   # expected.json holds 15 significant digits
SQM_NOISE = r"Neg*|SystemControl*|Blank*|BLANK*|Unassigned*"            # getMeanSignalRatio / getMaxRatio
SQM_TPC = r"Unassigned*|NegControl*|BLANK*|Blank*|SystemControl*|NegPrb*"  # getTxPerCell
SQM_FDR = r"(?i)^(Neg|Blank|BLANK)"                                      # getGlobalFDR


def _crop(control: str | None = None):
    """AnnData of the crop: X = counts of all 541 features, var['feature_type'], obs cell_area/total_counts,
    obsm['spatial']. ``control`` = a SpatialQM name regex for var['control']; None leaves the feature-type rule."""
    c = np.load(HERE / "counts.npz", allow_pickle=False)
    X = sp.csr_matrix((c["data"], c["indices"], c["indptr"]), shape=tuple(c["shape"]))
    cells = pd.read_csv(HERE / "cells.csv")
    names = c["feature_name"].astype(str)
    a = anndata.AnnData(X=X.astype(np.float32),
                        var=pd.DataFrame({"feature_type": c["feature_type"].astype(str)}, index=names))
    if control is not None:
        a.var["control"] = [bool(re.search(control, n)) for n in names]
    a.obs_names = cells["cell_id"].astype(str).to_numpy()
    a.obs["cell_area"] = cells["cell_area"].to_numpy()
    a.obs["total_counts"] = cells["total_counts"].to_numpy()
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()
    return a


def _transcripts():
    c = np.load(HERE / "counts.npz", allow_pickle=False)
    t = np.load(HERE / "transcripts.npz", allow_pickle=False)
    names = c["feature_name"].astype(str)
    return pd.DataFrame({"feature_name": names[t["feature_idx"]], "cell_id": t["cell_id"]}), names


def test_sparsity_and_tx_per_area_mean():
    from spatial_anno_metrics import signal_qc as sq

    a = _crop()
    assert sq.sparsity(a)["zero_fraction"] == pytest.approx(E["getSparsity"], rel=REL)
    out = sq.tx_per_area(a)
    assert out["mean_tx_per_area"] == pytest.approx(E["getTxPerArea"], rel=REL)
    assert out["median_tx_per_area"] == pytest.approx(1.056031914269972, rel=REL)   # unchanged since 0.6


def test_detection_entropy_spatialqm():
    from spatial_anno_metrics import signal_qc as sq

    a = _crop()
    out = sq.detection_entropy(a, method="spatialqm")
    assert out["entropy_bits"] == pytest.approx(E["getEntropy"], rel=REL)
    assert "sam_detection_entropy" not in a.obs                  # a dataset scalar: nothing per cell


def test_signal_to_noise_spatialqm():
    from spatial_anno_metrics import signal_qc as sq

    out = sq.signal_to_noise(_crop(SQM_NOISE), method="spatialqm")
    assert out["max_signal_ratio_log10"] == pytest.approx(E["getMaxRatio"], rel=REL)
    # SpatialQM recycles the 220 noise means against the 321 gene means; that pairing moves its value within
    # 0.49696-0.49725 over noise-row orders (SOURCE.md finding 7). Ours pairs every gene with every control.
    assert out["mean_signal_ratio_log10"] == pytest.approx(E["getMeanSignalRatio"], abs=5e-4)
    assert (out["n_genes"], out["n_controls"]) == (321, 220)
    ft = sq.signal_to_noise(_crop(), method="spatialqm")          # feature-type rule: antisense_* are controls
    assert (ft["n_genes"], ft["n_controls"]) == (313, 228)
    assert ft["max_signal_ratio_log10"] == pytest.approx(4.001, abs=5e-4)


def test_moran_spatialqm_reproduces_getmorans():
    from spatial_anno_metrics import signal_qc as sq

    a = _crop()
    out = sq.moran_signal(a, method="spatialqm")
    ref = pd.Series(E["getMorans"]["I"], index=E["getMorans"]["feature"], dtype=float)
    ours = out["per_feature"].reindex(ref.index)
    assert ours.notna().sum() == ref.notna().sum() == 334
    assert (ours.isna() == ref.isna()).all()                     # NaN at the same 7 features
    assert np.nanmax(np.abs(ours - ref)) < 1e-12                  # 5e-15 measured
    types = pd.Series(E["getMorans"]["type"], index=ref.index)
    assert out["median_gene_moran"] == pytest.approx(ref[types == "Gene"].median(), abs=1e-12)
    assert out["median_control_moran"] == pytest.approx(ref[types == "Control"].median(), abs=1e-12)
    assert (out["n_genes"], out["n_controls"]) == (313, 28)      # negative-control probes only, antisense included
    assert (out["n_genes_excluded"], out["n_controls_excluded"]) == (1, 6)
    assert (out["n_cells_genes"], out["n_cells_controls"]) == (2000, 150)


def test_moran_default_summaries_are_nan_free():
    """0.6 returned NaN medians here: 139 features have no count in the crop and a NaN Moran's I."""
    pytest.importorskip("squidpy")
    from spatial_anno_metrics import signal_qc as sq

    out = sq.moran_signal(_crop())
    assert out["median_gene_moran"] == pytest.approx(0.10637, abs=5e-6)
    assert out["median_control_moran"] == pytest.approx(-0.00058, abs=5e-6)
    assert out["separation"] == pytest.approx(0.10695, abs=5e-6)
    assert out["frac_genes_above_control_p95"] == pytest.approx(0.929, abs=5e-4)
    assert (out["n_genes"], out["n_controls"]) == (313, 228)
    assert (out["n_genes_excluded"], out["n_controls_excluded"]) == (1, 138)


def test_signal_to_noise_cell_does_not_sit_on_the_eps_floor():
    """88.8 % of crop cells have no control count; 0.6 put their ratio at log2(1/1e-9) and reported 27.956."""
    from spatial_anno_metrics import signal_qc as sq

    a = _crop()
    out = sq.signal_to_noise(a)
    assert out["n_cells_with_controls"] == 224
    assert out["frac_cells_with_controls"] == pytest.approx(0.112)
    assert out["median_mean_signal_ratio"] == pytest.approx(6.491, abs=5e-4)
    assert out["pooled_mean_signal_ratio"] == pytest.approx(9.511, abs=5e-4)
    assert a.obs["sam_mean_signal_ratio"].isna().sum() == 2000 - 224
    old = sq.signal_to_noise(_crop(), eps=1e-9)                   # the 0.6 recipe, on request
    assert old["median_mean_signal_ratio"] == pytest.approx(27.95603501634586, rel=1e-12)
    assert old["median_max_signal_ratio"] == pytest.approx(32.70470777624996, rel=1e-12)


def test_tx_per_cell():
    from spatial_anno_metrics import signal_qc as sq

    assert sq.tx_per_cell(_crop(SQM_TPC))["mean_tx_per_cell"] == pytest.approx(E["getTxPerCell"], rel=REL)
    assert sq.tx_per_cell(_crop())["mean_tx_per_cell"] == pytest.approx(128.1905, rel=REL)  # antisense_* out


def test_transcript_qc():
    from spatial_anno_metrics import signal_qc as sq

    tx, names = _transcripts()
    neg = [n for n in names if re.search(SQM_FDR, n)]
    out = sq.transcript_qc(tx, neg)
    assert out["frac_tx_in_cells"] == pytest.approx(E["getCellTxFraction"], rel=REL)
    assert out["fdr"] / 100 == pytest.approx(E["getGlobalFDR"], rel=REL)     # SpatialQM divides by 100
    assert (out["n_transcripts"], out["n_unassigned"], out["n_negative_tx"]) == (321637, 65, 1223)
    ftype = np.load(HERE / "counts.npz", allow_pickle=False)["feature_type"].astype(str)
    ft = sq.transcript_qc(tx, names[ftype != "Gene Expression"])             # antisense_* as negatives
    assert ft["fdr"] / 100 == pytest.approx(6.64377e-05, rel=1e-5)
    xoa13 = tx.assign(cell_id=tx["cell_id"].astype(str).replace("-1", "UNASSIGNED"))   # Xenium >= 1.3 ids
    assert sq.transcript_qc(xoa13, neg)["frac_tx_in_cells"] == out["frac_tx_in_cells"]
