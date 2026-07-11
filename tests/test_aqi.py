"""Tests for the AQI headline (annotation_quality_index): the single 0-1 quality score.

AQI = w_coh * soft-min_{p=-4}({A adequacy, C contamination, G agreement, M marker-fidelity}), an
INDEX (not P(correct)), with any undefined component DROPPED (never 0) and coherence a bounded
multiplier. These pin the math (bounds, soft-min <= arithmetic mean, drop-undefined, regime flip).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _counts_fixture(n_per=70, n_bg=16, seed=0):
    """A count section + a matched raw-count reference, each with per-type marker blocks so purity,
    marker-fidelity, panel-adequacy and cross-method agreement all resolve. 3 well-separated types."""
    from sklearn.decomposition import PCA

    rng = np.random.default_rng(seed)
    types = ["A", "B", "C"]
    blocks = {"A": list(range(0, 8)), "B": list(range(8, 16)), "C": list(range(16, 24))}
    genes = [f"g{i}" for i in range(24 + n_bg)]

    def make(n):
        X, lab = [], []
        for t in types:
            M = rng.poisson(0.5, size=(n, len(genes))).astype(float)
            M[:, blocks[t]] += rng.poisson(15, size=(n, len(blocks[t])))
            X.append(M)
            lab += [t] * n
        return np.vstack(X), lab

    Xs, labs = make(n_per)
    a = anndata.AnnData(X=Xs.copy())
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = pd.Categorical(labs)
    a.layers["counts"] = Xs.copy()
    logn = np.log1p(Xs / Xs.sum(1, keepdims=True) * np.median(Xs.sum(1)))
    a.obsm["X_pca"] = PCA(n_components=8, random_state=0).fit_transform(logn).astype("float32")

    Xr, labr = make(120)
    ref = anndata.AnnData(X=Xr.copy())
    ref.var_names = genes
    ref.obs["cell_type"] = pd.Categorical(labr)

    for m in ("m1", "m2", "m3"):                                   # 3 voters, ~12% noise each
        v = np.array(labs, dtype=object)
        flip = rng.random(len(v)) < 0.12
        v[flip] = rng.choice(types, size=int(flip.sum()))
        a.obs[m] = pd.Categorical(v)

    marker_sets = {t: [genes[i] for i in blocks[t]] for t in types}
    return a, ref, marker_sets, genes


def test_aqi_full_regime_bounds_and_softmin():
    from spatial_anno_metrics import eval_metrics as em

    a, ref, msets, genes = _counts_fixture()
    depth = float(np.median(np.asarray(a.layers["counts"]).sum(1)))
    out = em.annotation_quality(
        a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
        reference=ref, panel_genes=genes, median_depth=depth,
        method_label_cols=["m1", "m2", "m3"], platform="xenium", n_panel_genes=len(genes),
        normalization="log1p-median", marker_set_name="curated")
    aqi = out["aqi"]

    assert set(aqi["active_set"]) == {"A", "C", "G", "M"}          # all four defined on clean data
    comps = [aqi["components"][k] for k in ("A", "C", "G", "M")]
    for v in comps:
        assert v is not None and 0.0 <= v <= 1.0
    assert 0.0 <= aqi["aqi"] <= 1.0 and 0.0 <= aqi["aqi_core"] <= 1.0
    assert aqi["regime"] == "agreement_anchored"                  # >=3 voters
    assert aqi["argmin"] in ("panel/depth adequacy", "contamination",
                             "cross-method agreement", "marker fidelity")
    # soft-min: the power mean (p=-4) sits between the weakest link and the arithmetic mean.
    assert min(comps) - 1e-6 <= aqi["aqi_core"] <= float(np.mean(comps)) + 1e-6
    assert aqi["aqi"] <= aqi["aqi_core"] + 1e-9                    # coherence multiplier <= 1
    assert aqi["n_ge50_types"] == 3
    assert aqi["provenance"]["marker_set"] == "curated" and aqi["provenance"]["n_voters"] == 3


def test_aqi_drops_undefined_components_and_flips_regime():
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, msets, _genes = _counts_fixture()
    # <3 voters -> agreement UNDEFINED; no reference -> adequacy UNDEFINED. Both must DROP, not become 0.
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
                                method_label_cols=["m1"])
    aqi = out["aqi"]
    assert aqi["regime"] == "coherence_index"
    assert aqi["components"]["G"] is None and "G" not in aqi["active_set"]
    assert aqi["components"]["A"] is None and aqi["flags"]["reference_unanchored"] is True
    assert set(aqi["active_set"]) == {"C", "M"}                    # only the reference-free terms survive
    assert 0.0 <= aqi["aqi"] <= 1.0


def test_aqi_empty_active_set_is_none_not_zero():
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, _m, _g = _counts_fixture()
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca")  # no markers/ref/methods
    aqi = out["aqi"]
    assert aqi["active_set"] == [] and aqi["aqi"] is None and aqi["regime"] == "coherence_index"
