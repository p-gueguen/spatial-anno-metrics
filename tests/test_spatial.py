"""Tests for spatial.py - spatial label-quality metrics. A COHERENT labeling (types occupy distinct
spatial blocks) must score higher coherence / lower PAS than a SHUFFLED one; neighborhood_sanity must
not flag a clean block layout."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")
pytest.importorskip("squidpy")


def _blocks(n_side=20, seed=0):
    """A 2*n_side x n_side grid: left half labeled A, right half B (spatially coherent)."""
    rng = np.random.default_rng(seed)
    xs, ys, lab = [], [], []
    for i in range(2 * n_side):
        for j in range(n_side):
            xs.append(i + rng.normal(0, 0.01))
            ys.append(j + rng.normal(0, 0.01))
            lab.append("A" if i < n_side else "B")
    a = anndata.AnnData(X=np.zeros((len(xs), 3), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(len(xs))]
    a.obsm["spatial"] = np.column_stack([xs, ys]).astype(float)
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def test_spatial_coherence_coherent_vs_shuffled():
    from spatial_anno_metrics import spatial as sp

    a = _blocks()
    coh = sp.spatial_coherence(a, label_key="cell_type", n_neighs=8)
    assert "sam_spatial_coherence" in a.obs and a.uns["sam"]["spatial_coherence"]["n_neighs"] == 8
    assert coh["mean_coherence"] > 0.8 and coh["pas"] < 0.15      # clean blocks: neighbors share label
    assert coh["status"] == "ok" and coh["n_neighs"] == 8

    b = a.copy()
    b.obs["cell_type"] = pd.Categorical(np.random.default_rng(1).permutation(b.obs["cell_type"].to_numpy()))
    del b.obsp["spatial_connectivities"]                          # rebuild graph for the shuffled copy
    shuf = sp.spatial_coherence(b, label_key="cell_type", n_neighs=8)
    assert shuf["mean_coherence"] < coh["mean_coherence"]
    assert shuf["pas"] > coh["pas"]


def test_neighborhood_sanity_clean_layout_not_flagged():
    from spatial_anno_metrics import spatial as sp

    out = sp.neighborhood_sanity(_blocks(), label_key="cell_type", n_neighs=8)
    assert out["status"] == "ok" and set(out["categories"]) == {"A", "B"}


def test_neighborhood_sanity_leaves_the_label_column_alone():
    """It used to cast the caller's label column to category and write squidpy's result into its uns; an
    unused category level changed the result."""
    from spatial_anno_metrics import spatial as sp

    a = _blocks()
    a.obs["cell_type"] = a.obs["cell_type"].astype(str).to_numpy().astype(object)
    before = a.obs["cell_type"].copy()
    out = sp.neighborhood_sanity(a, label_key="cell_type", n_neighs=8)
    pd.testing.assert_series_equal(a.obs["cell_type"], before)
    assert not any("nhood_enrichment" in k for k in a.uns)
    b = _blocks()
    b.obs["cell_type"] = pd.Categorical(b.obs["cell_type"].astype(str), categories=["A", "B", "Ghost"])
    ghost = sp.neighborhood_sanity(b, label_key="cell_type", n_neighs=8)
    assert ghost["categories"] == ["A", "B"] and ghost["self_enrichment"] == out["self_enrichment"]
    # a clean 2-block layout: each type's strongest partner is itself -> nothing suspicious
    assert out["suspicious"] == []
    assert len(out["self_enrichment"]) == 2


SPILL_MARKERS = {"T": ["CD3E"], "Fib": ["COL1A1", "ACTA2"], "SMC": ["ACTA2", "MYH11"]}


def _spill_tissue():
    """10 cells on an explicit contact graph: c0 T touches c1 Fib | c2-c3 T pair, c2 carries COL1A1 with no
    Fib neighbor | c4 T touches c5, a Fib cell detecting no Fib marker | c6 T touches c7 Fib + c8 SMC,
    whose marker sets share ACTA2 | c9 an empty (0-count) T cell touching c1."""
    from scipy.sparse import csr_matrix

    X = np.array([                                       # CD3E COL1A1 ACTA2 MYH11 BG
        [9, 1, 0, 0, 0],                                 # c0 T
        [0, 10, 0, 0, 0],                                # c1 Fib
        [9, 1, 0, 0, 0],                                 # c2 T
        [10, 0, 0, 0, 0],                                # c3 T
        [9, 1, 0, 0, 0],                                 # c4 T
        [0, 0, 0, 0, 5],                                 # c5 Fib
        [8, 0, 2, 0, 0],                                 # c6 T
        [0, 5, 5, 0, 0],                                 # c7 Fib
        [0, 0, 5, 5, 0],                                 # c8 SMC
        [0, 0, 0, 0, 0],                                 # c9 T
    ], dtype="float32")
    i, j = zip(*[(0, 1), (2, 3), (4, 5), (6, 7), (6, 8), (1, 9)])
    a = anndata.AnnData(X=X)
    a.var_names = ["CD3E", "COL1A1", "ACTA2", "MYH11", "BG"]
    a.obs_names = [f"c{k}" for k in range(10)]
    a.obs["cell_type"] = ["T", "Fib", "T", "T", "T", "Fib", "T", "Fib", "SMC", "T"]
    a.obsp["spatial_connectivities"] = csr_matrix((np.ones(12), (i + j, j + i)), shape=(10, 10))
    return a


def test_neighbor_contamination_per_cell():
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    out = sp.neighbor_contamination(a, SPILL_MARKERS)
    # c0 1/10 next to a COL1A1+ Fib; c2/c3 no other type adjacent and c4's Fib neighbor detects no Fib marker
    # -> not evaluable; c6 ACTA2 2/10 is a marker of BOTH adjacent sources -> counted once; c9 has no counts
    want = [0.1, 0.0, np.nan, np.nan, np.nan, 0.0, 0.2, 0.0, 0.0, np.nan]
    np.testing.assert_allclose(a.obs["sam_neighbor_contamination"].to_numpy(), want)
    assert out["mean_contamination"] == pytest.approx(0.05) and out["n_evaluable"] == 6
    assert out["status"] == "ok" and out["matrix"] == "X"


def test_neighbor_contamination_without_neighbor_expression_check():
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    sp.neighbor_contamination(a, SPILL_MARKERS, require_neighbor_expression=False)
    assert a.obs["sam_neighbor_contamination"].iloc[4] == pytest.approx(0.1)   # any adjacent Fib cell suffices


def test_neighbor_contamination_ignores_own_markers():
    """Fib and SMC sets share ACTA2: a Fib cell's own ACTA2 is not spillover from the SMC it touches."""
    from scipy.sparse import csr_matrix

    from spatial_anno_metrics import spatial as sp

    a = anndata.AnnData(X=np.array([[5, 5, 0], [0, 5, 5]], dtype="float32"))     # COL1A1 ACTA2 MYH11
    a.var_names = ["COL1A1", "ACTA2", "MYH11"]
    a.obs_names = ["fib", "smc"]
    a.obs["cell_type"] = ["Fib", "SMC"]
    a.obsp["spatial_connectivities"] = csr_matrix(np.array([[0, 1], [1, 0]], dtype=float))
    sp.neighbor_contamination(a, SPILL_MARKERS)
    np.testing.assert_allclose(a.obs["sam_neighbor_contamination"].to_numpy(), [0.0, 0.0])


def test_neighbor_contamination_source_target_matrix():
    from spatial_anno_metrics import spatial as sp

    out = sp.neighbor_contamination(_spill_tissue(), SPILL_MARKERS)
    m, n = out["source_target"], out["n_source_target"]  # rows = source type, columns = target type
    assert m.loc["Fib", "T"] == pytest.approx(0.15) and n.loc["Fib", "T"] == 2    # c0 0.1, c6 0.2
    assert m.loc["SMC", "T"] == pytest.approx(0.2) and n.loc["SMC", "T"] == 1
    assert m.loc["T", "Fib"] == pytest.approx(0.0) and n.loc["T", "Fib"] == 3     # c1, c5, c7
    assert np.isnan(m.loc["Fib", "SMC"]) and n.loc["Fib", "SMC"] == 0             # never adjacent


def test_spatial_label_metrics_skip_a_single_label():
    from spatial_anno_metrics import spatial as sp

    a = _blocks(n_side=6)
    a.obs["cell_type"] = "A"
    for out in (sp.spatial_coherence(a, label_key="cell_type"), sp.neighborhood_sanity(a, label_key="cell_type")):
        assert out["status"] == "skipped" and "found 1" in out["reason"]


def test_spatial_key_added():
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    sp.neighbor_contamination(a, SPILL_MARKERS, key_added="nc")
    assert "nc" in a.obs and "sam_neighbor_contamination" not in a.obs
    b = _blocks(n_side=6)
    sp.spatial_coherence(b, label_key="cell_type", key_added="coh")
    assert "coh" in b.obs


def test_neighbor_contamination_skips_without_two_typed_labels_or_contact():
    """One label, every marker off the panel, or no cell touching another typed cell: a skip, no obs column."""
    from scipy.sparse import csr_matrix

    from spatial_anno_metrics import spatial as sp

    one = _spill_tissue()
    one.obs["cell_type"] = "T"
    off = _spill_tissue()
    apart = _spill_tissue()
    apart.obsp["spatial_connectivities"] = csr_matrix((10, 10))
    for a, markers, why in ((one, SPILL_MARKERS, "found 1"), (off, {"T": ["X1"], "Fib": ["X2"]}, "found 0"),
                            (apart, SPILL_MARKERS, "no cell")):
        out = sp.neighbor_contamination(a, markers)
        assert out["status"] == "skipped" and why in out["reason"]
        assert "sam_neighbor_contamination" not in a.obs and "sam" not in a.uns
