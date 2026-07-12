"""Tests for the AQI headline (annotation_quality_index): the single 0-1 quality score.

AQI = w_coh * min(A ceiling, soft-min_{p=-4}({C contamination, M marker-fidelity})), an INDEX (not
P(correct)), with any undefined term DROPPED (never 0) and coherence a bounded multiplier. Cross-method
agreement G is a WITHIN-section abstention signal (reported, not a cross-section term - it does not
transfer across sections per the 2026-07 validation). These pin the math (bounds, soft-min <= arithmetic
mean of {C,M}, A caps, drop-undefined, regime = abstention availability).
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

    # C, M drive the soft-min; A joins as the ceiling. G is reported under abstention, NOT active_set.
    assert set(aqi["active_set"]) == {"C", "M", "A"}
    for k in ("A", "C", "M"):
        v = aqi["components"][k]
        assert v is not None and 0.0 <= v <= 1.0
    assert 0.0 <= aqi["aqi"] <= 1.0 and 0.0 <= aqi["aqi_core"] <= 1.0
    assert aqi["regime"] == "with_agreement"                       # >=3 voters -> abstention available
    assert aqi["abstention"]["available"] is True and aqi["abstention"]["n_voters"] == 3
    assert aqi["abstention"]["signal"] is not None
    assert aqi["argmin"] in ("contamination", "marker fidelity", "panel/depth ceiling")
    # soft-min core: the power mean (p=-4) of {C, M} sits between the weaker link and their mean.
    cm = [aqi["components"]["C"], aqi["components"]["M"]]
    assert min(cm) - 1e-6 <= aqi["aqi_core"] <= float(np.mean(cm)) + 1e-6
    assert aqi["aqi"] <= aqi["aqi_core"] + 1e-9                     # A-cap and coherence multiplier <= core
    assert aqi["n_ge50_types"] == 3
    assert aqi["provenance"]["marker_set"] == "curated" and aqi["provenance"]["n_voters"] == 3


def test_aqi_a_is_a_ceiling_not_a_bottleneck():
    """A LOW adequacy caps the index (min), but a HIGH adequacy never lifts it above the {C,M} soft-min:
    that is the whole point of the 2026-07 correction (A anti-orders accuracy, so it may only cap)."""
    from spatial_anno_metrics import eval_metrics as em

    a, ref, msets, genes = _counts_fixture()
    depth = float(np.median(np.asarray(a.layers["counts"]).sum(1)))
    out = em.annotation_quality(
        a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
        reference=ref, panel_genes=genes, median_depth=depth,
        method_label_cols=["m1", "m2", "m3"], platform="xenium", n_panel_genes=len(genes))
    aqi = out["aqi"]
    A, core = aqi["components"]["A"], aqi["aqi_core"]
    # aqi (pre-w_coh) = min(A, core): capped by A only when A < core, else equals core. Never exceeds core.
    assert aqi["aqi"] <= min(A, core) * aqi["w_coh"] + 1e-9
    if A >= core:                                                  # rich panel must not inflate the score
        assert abs(aqi["aqi"] - aqi["w_coh"] * core) < 1e-9


def test_aqi_drops_undefined_components_and_flips_regime():
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, msets, _genes = _counts_fixture()
    # <3 voters -> agreement UNDEFINED (abstention unavailable); no reference -> adequacy UNDEFINED. Both DROP.
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
                                method_label_cols=["m1"])
    aqi = out["aqi"]
    assert aqi["regime"] == "index_only"
    assert aqi["components"]["G"] is None and aqi["abstention"]["available"] is False
    assert aqi["components"]["A"] is None and aqi["flags"]["adequacy_unknown"] is True
    assert set(aqi["active_set"]) == {"C", "M"}                    # only the reference-free terms survive
    assert 0.0 <= aqi["aqi"] <= 1.0


def test_aqi_empty_active_set_is_none_not_zero():
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, _m, _g = _counts_fixture()
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca")  # no markers/ref/methods
    aqi = out["aqi"]
    assert aqi["active_set"] == [] and aqi["aqi"] is None and aqi["regime"] == "index_only"


def test_aqi_coverage_flags_a_collapsed_lineage():
    """The advisory coverage flag must catch a lineage the markers strongly support but the annotation
    dropped ENTIRELY - the AQI index blind spot confirmed on colon. Collapse type C into A: C's markers
    stay strong in those cells, no C label survives, so coverage.missing must name C."""
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, msets, _genes = _counts_fixture()
    lab = a.obs["cell_type"].astype(str).to_numpy()
    lab[lab == "C"] = "A"                                        # collapse C -> A
    a.obs["cell_type"] = pd.Categorical(lab)
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca")
    aqi = out["aqi"]
    cov = aqi["coverage"]
    assert cov is not None
    assert "C" in cov["missing"], cov                           # markers say C is here, no C label
    assert "A" in cov["annotated_on_axis"]
    assert aqi["flags"]["missing_lineages"] is True
    assert cov["missing_frac"] > 0.1                            # a real fraction of cells look like C


def test_aqi_coverage_empty_when_all_lineages_present():
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, msets, _genes = _counts_fixture()                  # A, B, C all annotated
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca")
    cov = out["aqi"]["coverage"]
    assert cov is not None and cov["missing"] == []
    assert out["aqi"]["flags"]["missing_lineages"] is False


def test_aqi_cov_ref_flags_dropped_reference_lineage():
    """Reference-anchored coverage catches a lineage ABUNDANT in the reference but absent from the labels -
    the fix for the marker-based blind spot on off-panel-marker panels (UM). Collapse C -> A."""
    from spatial_anno_metrics import eval_metrics as em

    a, ref, msets, genes = _counts_fixture()
    lab = a.obs["cell_type"].astype(str).to_numpy()
    lab[lab == "C"] = "A"
    a.obs["cell_type"] = pd.Categorical(lab)
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
                                reference=ref, ref_key="cell_type", panel_genes=genes)
    cov = out["aqi"]["coverage"]
    assert "C" in cov["missing_vs_reference"], cov
    assert out["aqi"]["flags"]["missing_vs_reference"] is True


def test_aqi_cov_ref_clears_when_reference_lineages_all_present():
    from spatial_anno_metrics import eval_metrics as em

    a, ref, msets, genes = _counts_fixture()                    # A, B, C all annotated == ref types
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
                                reference=ref, ref_key="cell_type", panel_genes=genes)
    assert out["aqi"]["coverage"]["missing_vs_reference"] == []
    assert out["aqi"]["flags"]["missing_vs_reference"] is False


def test_aqi_cohort_signal_surfaces_with_sample_key():
    """When a cohort exists (sample_key), inter-sample consistency is reported under `cohort` (not scored)."""
    from spatial_anno_metrics import eval_metrics as em

    a, _ref, msets, _genes = _counts_fixture()
    a.obs["sample"] = pd.Categorical(np.tile(["s1", "s2"], a.n_obs // 2 + 1)[:a.n_obs])
    out = em.annotation_quality(a, label_key="cell_type", marker_sets=msets, embedding="X_pca",
                                sample_key="sample")
    coh = out["aqi"]["cohort"]
    assert coh is not None and "consistency" in coh and coh["n_profiles"] is not None
