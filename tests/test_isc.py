"""Tests for isc.py - inter-sample consistency: ISCG / ISCL, regimes, over-partitioning, dropped types."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _isc_adata(consistent: bool, n_per=40, n_types=3, n_samples=5, n_genes=30, seed=0):
    """Each (sample, type) is a pseudobulk. Consistent: a type's signature block is the SAME
    across samples. Inconsistent: each sample uses a DIFFERENT block for a given type."""
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    X, types, samples = [], [], []
    for s in range(n_samples):
        for t in range(n_types):
            sig = t if consistent else rng.integers(0, n_types)  # which block is high
            base = rng.normal(0.2, 0.1, size=(n_per, n_genes))
            base[:, sig * block:(sig + 1) * block] += 5.0
            X.append(base)
            types += [f"T{t}"] * n_per
            samples += [f"S{s}"] * n_per
    a = anndata.AnnData(X=np.vstack(X).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = pd.Categorical(types)
    a.obs["sample"] = pd.Categorical(samples)
    return a


def test_inter_sample_consistency_high_for_reproducible_types():
    from spatial_anno_metrics import isc as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    assert good["consistency"] > bad["consistency"]


def test_inter_sample_consistency_two_axes_and_regime():
    """Garnica 2026: ISCG (2-label silhouette | pseudobulk cosine) and ISCL (silhouette | reciprocal
    classification match) both separate reproducible from non-reproducible types; product + regime."""
    from spatial_anno_metrics import isc as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    for k in ("iscg", "iscl", "integrated"):
        assert good[k] > bad[k], k
    assert good["iscl"] > 0.9 and set(good["regime"].values()) == {"consistent"}
    assert good["dropped_types"] == {} and good["n_profiles"] == 15
    assert good["thresholds"] == {"iscg": 0.74, "iscl": 0.45}
    assert good["status"] == "ok" and good["matrix"] == "X"          # layer=None: X averaged as it is


def test_inter_sample_consistency_over_partitioning_drops_local_not_global():
    """Benchmark task 2 (the paper's positive control for ISCL): split one type into two identical
    halves in every sample. The halves cannot be reciprocally matched -> ISCL collapses; they stay
    compact -> ISCG stays high; regime = 'merge'. Untouched types stay 'consistent'."""
    from spatial_anno_metrics import isc as em

    base = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    a = _isc_adata(consistent=True)
    lab = a.obs["cell_type"].astype(str).to_numpy()
    lab[np.flatnonzero(lab == "T0")[::2]] = "T0.b"  # every 2nd T0 cell, i.e. 20/20 per sample
    a.obs["cell_type"] = pd.Categorical(lab)
    split = em.inter_sample_consistency(a, "cell_type", "sample")
    for t in ("T0", "T0.b"):
        assert split["per_type_iscl"][t] < base["per_type_iscl"]["T0"] - 0.5, t
        assert split["per_type_iscg"][t] > 0.74, t
        assert split["regime"][t] == "merge", t
    assert split["regime"]["T1"] == "consistent" and split["regime"]["T2"] == "consistent"


def test_inter_sample_consistency_drops_underpowered_types_and_discloses():
    """A type present in fewer than min_samples samples is reported under dropped_types, not scored
    (the paper's 3 % 'not computable'); lowering min_samples brings it back."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True)
    lab, smp = a.obs["cell_type"].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    a = a[~((lab == "T2") & np.isin(smp, ["S0", "S1", "S2"]))].copy()  # T2 survives in 2 samples only
    out = em.inter_sample_consistency(a, "cell_type", "sample")
    assert out["dropped_types"] == {"T2": 2} and "T2" not in out["per_type_iscl"]
    relaxed = em.inter_sample_consistency(a, "cell_type", "sample", min_samples=2)
    assert relaxed["dropped_types"] == {} and "T2" in relaxed["per_type_iscl"]
    none = em.inter_sample_consistency(a, "cell_type", "sample", min_samples=6)    # only 5 samples
    assert none["status"] == "skipped" and ">= 6 samples" in none["reason"]
    assert none["iscg"] is None and none["consistency"] is None and len(none["dropped_types"]) == 3


def test_inter_sample_consistency_counts_layer_recipe():
    """layer= names a raw-count layer (was counts_layer=): pseudobulks are summed counts, log1p CP10k."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True)
    a.layers["counts"] = np.round(np.clip(a.X, 0, None) * 10)
    out = em.inter_sample_consistency(a, "cell_type", "sample", layer="counts")
    assert out["status"] == "ok" and out["matrix"] == "layers['counts']" and out["iscl"] > 0.9
    with pytest.raises(TypeError):
        em.inter_sample_consistency(a, "cell_type", "sample", counts_layer="counts")


def test_isc_regime_quadrants():
    from spatial_anno_metrics import isc as em

    assert em.isc_regime(0.9, 0.9) == "consistent"
    assert em.isc_regime(0.9, 0.1) == "merge"
    assert em.isc_regime(0.1, 0.9) == "split"
    assert em.isc_regime(0.1, 0.1) == "remove"
    assert em.isc_regime(None, 0.5) == "not_computable"
    assert isinstance(em.isc_regime(0.9, 0.9), str)                  # a helper: returns a str, not a dict
