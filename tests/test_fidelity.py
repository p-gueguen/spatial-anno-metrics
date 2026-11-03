"""Tests for fidelity.py - marker-program fidelity (AUC / Cohen's d) and marker-gene overlap."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _marker_adata(correct=True, seed=0):
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"N{i}" for i in range(10)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 120
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.2, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:60, gi[g]] += 10
    for g in ("MLANA", "SOX10", "DCT"):
        X[60:, gi[g]] += 10
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    labels = ["T cell"] * 60 + ["Mel"] * 60
    if not correct:
        labels = list(rng.permutation(labels))
    a.obs["cell_type"] = pd.Categorical(labels)
    return a


def test_marker_program_fidelity_high_auc_for_correct_labels():
    from spatial_anno_metrics import fidelity as em

    markers = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_program_fidelity(_marker_adata(correct=True), "cell_type", markers)
    bad = em.marker_program_fidelity(_marker_adata(correct=False), "cell_type", markers)
    assert good["mean_auc"] > 0.9
    assert good["mean_auc"] > bad["mean_auc"]
    assert good["per_type"]["T cell"]["auc"] > 0.9
    assert good["per_type"]["T cell"]["cohens_d"] > 0.5           # large, positive separation
    assert good["status"] == "ok" and good["matrix"] == "X"          # no counts layer, X holds counts


def test_marker_gene_overlap_high_when_labels_match_markers():
    from spatial_anno_metrics import fidelity as em

    rng = np.random.default_rng(3)
    genes = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"BG{i}" for i in range(30)]
    gi = {g: i for i, g in enumerate(genes)}
    n = 160
    X = rng.poisson(0.3, size=(n, len(genes))).astype("float32")
    for g in ("CD3D", "CD3E", "TRAC"):
        X[:80, gi[g]] += 15
    for g in ("MLANA", "SOX10", "DCT"):
        X[80:, gi[g]] += 15
    a = anndata.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T"] * 80 + ["Mal"] * 80)
    ms = {"T": ["CD3D", "CD3E", "TRAC"], "Mal": ["MLANA", "SOX10", "DCT"]}
    good = em.marker_gene_overlap(a, "cell_type", ms, n_top=10)
    assert good["mean_overlap"] == pytest.approx(1.0)              # each type's markers ARE its top DE genes
    b = a.copy()
    b.obs["cell_type"] = pd.Categorical(rng.permutation(b.obs["cell_type"].to_numpy()))
    assert em.marker_gene_overlap(b, "cell_type", ms, n_top=10)["mean_overlap"] < good["mean_overlap"]


def test_fidelity_reads_counts_and_skips_a_single_label():
    from spatial_anno_metrics import fidelity as em

    markers = {"T cell": ["CD3D", "CD3E", "TRAC"], "Mel": ["MLANA", "SOX10", "DCT"]}
    a = _marker_adata()
    a.layers["counts"] = a.X.copy()
    a.X = np.log1p(a.X)                                             # the common scanpy layout
    assert em.marker_program_fidelity(a, "cell_type", markers)["matrix"] == "layers['counts']"
    assert em.marker_gene_overlap(a, "cell_type", markers)["matrix"] == "layers['counts']"
    del a.layers["counts"]
    with pytest.raises(ValueError, match="layer="):
        em.marker_program_fidelity(a, "cell_type", markers)        # log X is not silently read as counts
    a.obs["cell_type"] = "Mel"
    for fn in (em.marker_program_fidelity, em.marker_gene_overlap):
        out = fn(a, "cell_type", markers, layer=None)
        assert out["status"] == "skipped" and "found 1" in out["reason"]


def test_fidelity_skips_when_no_marker_is_on_the_panel():
    from spatial_anno_metrics import fidelity as fd

    a = _marker_adata()
    off = {"T cell": ["NOPE1"], "Mel": ["NOPE2"]}
    for out in (fd.marker_program_fidelity(a, "cell_type", off), fd.marker_gene_overlap(a, "cell_type", off)):
        assert out["status"] == "skipped" and "on the panel" in out["reason"]


def test_cohens_d_is_undefined_not_zero_for_a_constant_marker_load():
    from spatial_anno_metrics import fidelity as fd

    a = _marker_adata()
    a.X[:, a.var_names.get_loc("N0")] = 0                        # a marker set that is all zeros
    out = fd.marker_program_fidelity(a, "cell_type", {"T cell": ["CD3D"], "Mel": ["N0"]})
    assert out["per_type"]["Mel"]["cohens_d"] is None
    assert out["mean_cohens_d"] == pytest.approx(out["per_type"]["T cell"]["cohens_d"])
