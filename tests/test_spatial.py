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
    assert np.isnan(m.loc["Fib", "SMC"]) and np.isnan(n.loc["Fib", "SMC"])        # never adjacent: NaN, as SegTraQ


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


# --------------------------------------------------------------------------- #
# 0.7.0: PAS to the cited definition, CHAOS, the graph cache
# --------------------------------------------------------------------------- #
def _scatter(n=600, seed=0):
    """Random points, labels from 6 noisy patches of 3 types (a mix of clean and mixed neighbourhoods)."""
    from scipy.spatial import cKDTree

    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, 30, (n, 2))
    lab = cKDTree(rng.uniform(0, 30, (6, 2))).query(xy)[1] % 3
    flip = rng.random(n) < 0.3
    lab[flip] = rng.integers(0, 3, flip.sum())
    a = anndata.AnnData(X=np.zeros((n, 2), dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obsm["spatial"] = xy
    a.obs["cell_type"] = pd.Categorical(np.array(["A", "B", "C"])[lab])
    return a


def test_pas_defaults_are_shang_zhou_2022():
    """Abnormal = the label differs from >= 6 of the 10 nearest neighbours (fx_PAS in SpatialPCA)."""
    import inspect

    from scipy.spatial import cKDTree

    from spatial_anno_metrics import spatial as sp

    sig = inspect.signature(sp.spatial_coherence).parameters
    assert sig["n_neighs"].default == 10 and sig["pas_threshold"].default == 0.4
    a = _scatter()
    out = sp.spatial_coherence(a)
    xy, lab = a.obsm["spatial"], a.obs["cell_type"].to_numpy()
    nn = cKDTree(xy).query(xy, k=11)[1][:, 1:]                  # 10 nearest, self excluded
    want = np.mean((lab[nn] != lab[:, None]).sum(axis=1) >= 6)
    assert 0.05 < want < 0.95 and out["pas"] == pytest.approx(want)


def test_pas_threshold_is_inclusive():
    """A cell with exactly 4 of 10 same-label neighbours (0.4) is abnormal; 5 of 10 is not."""
    from scipy.sparse import csr_matrix

    from spatial_anno_metrics import spatial as sp

    for n_same, want in ((4, 1 / 11), (5, 0.0)):
        a = anndata.AnnData(X=np.zeros((11, 1), dtype="float32"))
        a.obs_names = [f"c{i}" for i in range(11)]
        a.obs["cell_type"] = ["A"] + ["A"] * n_same + ["B"] * (10 - n_same)
        a.obsp["spatial_connectivities"] = csr_matrix((np.ones(10), (np.zeros(10, int), np.arange(1, 11))),
                                                      shape=(11, 11))   # only c0 has neighbours
        out = sp.spatial_coherence(a)
        assert a.obs["sam_spatial_coherence"].iloc[0] == pytest.approx(n_same / 10)
        assert out["pas"] == pytest.approx(want)


def test_spatial_coherence_honours_n_neighs_over_a_cached_graph():
    """A graph another call cached (another n_neighs) used to be reused silently."""
    import squidpy as sq

    from spatial_anno_metrics import spatial as sp

    fresh = sp.spatial_coherence(_scatter(), n_neighs=10)
    a = _scatter()
    sp.neighborhood_sanity(a, n_neighs=6)                       # caches a 6-NN graph
    assert sp.spatial_coherence(a, n_neighs=10)["pas"] == pytest.approx(fresh["pas"])
    b = _scatter()
    sq.gr.spatial_neighbors(b, coord_type="generic", n_neighs=6)  # a squidpy graph built elsewhere
    assert sp.spatial_coherence(b, n_neighs=10)["mean_coherence"] == pytest.approx(fresh["mean_coherence"])
    assert b.uns["spatial_neighbors"]["params"]["n_neighbors"] == 6   # ... and left in place
    assert sp.spatial_coherence(b, n_neighs=10, graph="obsp")["mean_coherence"] != pytest.approx(
        fresh["mean_coherence"])                                 # graph="obsp" uses it as it is


def test_spatial_coherence_dense_graph():
    from spatial_anno_metrics import spatial as sp

    a = _scatter(n=200)
    ref = sp.spatial_coherence(a)
    b = _scatter(n=200)
    b.obsp["spatial_connectivities"] = a.obsp["spatial_connectivities"].toarray()
    assert sp.spatial_coherence(b)["pas"] == pytest.approx(ref["pas"])


def _chaos_by_hand(xy, lab):
    """SpatialPCA's CHAOS definition (Shang & Zhou 2022, Methods): coordinates scaled per axis, each cell's
    distance to its nearest same-label cell, summed over labels and divided by N (a singleton adds 0)."""
    z = (xy - xy.mean(0)) / xy.std(0, ddof=1)
    tot = 0.0
    for t in np.unique(lab):
        p = z[lab == t]
        if len(p) > 1:
            d = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1))
            np.fill_diagonal(d, np.inf)
            tot += d.min(1).sum()
    return tot / len(lab)


def test_chaos_coherent_vs_shuffled():
    from spatial_anno_metrics import spatial as sp

    a = _blocks(n_side=10)
    coh = sp.spatial_coherence(a)
    b = _blocks(n_side=10)
    b.obs["cell_type"] = pd.Categorical(np.random.default_rng(1).permutation(b.obs["cell_type"].to_numpy()))
    shuf = sp.spatial_coherence(b)
    assert coh["chaos"] < shuf["chaos"]
    for x in (a, b):
        want = _chaos_by_hand(x.obsm["spatial"], x.obs["cell_type"].to_numpy())
        assert sp.spatial_coherence(x)["chaos"] == pytest.approx(want)
    c = _scatter(n=100)
    c.obs["cell_type"] = c.obs["cell_type"].astype(str)
    c.obs.iloc[0, 0] = "lonely"                                  # a singleton label adds 0, N still counts it
    assert sp.spatial_coherence(c)["chaos"] == pytest.approx(
        _chaos_by_hand(c.obsm["spatial"], c.obs["cell_type"].to_numpy()))


def test_chaos_nan_without_coordinates():
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()                                          # a hand graph, no obsm['spatial']
    assert np.isnan(sp.spatial_coherence(a)["chaos"])


def test_neighborhood_sanity_flags_two_types_forced_adjacent():
    """Positive control: A and B alternate column by column (every A touches B), C is its own block."""
    from spatial_anno_metrics import spatial as sp

    a = _blocks()
    i = np.round(a.obsm["spatial"][:, 0]).astype(int)
    a.obs["cell_type"] = pd.Categorical(np.where(i >= 20, "C", np.where(i % 2 == 0, "A", "B")))
    out = sp.neighborhood_sanity(a, n_neighs=8)
    flagged = {d["type"]: d["stronger_with"] for d in out["suspicious"]}
    assert flagged == {"A": "B", "B": "A"}


# --------------------------------------------------------------------------- #
# neighbor_contamination: pinned behaviour (WF1 review) and the 0.7.0 additions
# --------------------------------------------------------------------------- #
def _tiny(X, var, labels, edges, n=None, directed=False):
    from scipy.sparse import csr_matrix

    n = n or len(labels)
    i, j = zip(*edges)
    r, c = (list(i), list(j)) if directed else (list(i + j), list(j + i))
    a = anndata.AnnData(X=np.asarray(X, dtype="float32"))
    a.var_names = var
    a.obs_names = [f"c{k}" for k in range(n)]
    a.obs["cell_type"] = labels
    a.obsp["spatial_connectivities"] = csr_matrix((np.ones(len(r)), (r, c)), shape=(n, n))
    return a


def test_neighbor_contamination_or_across_sources():
    """A T cell touching only an ACTA2+ Fib, with an SMC elsewhere (SMC also lists ACTA2): still 0.2."""
    from spatial_anno_metrics import spatial as sp

    var = ["CD3E", "COL1A1", "ACTA2", "MYH11"]
    X = [[8, 0, 2, 0], [0, 5, 5, 0], [0, 0, 5, 5], [0, 0, 5, 5]]   # T, Fib, SMC, SMC
    a = _tiny(X, var, ["T", "Fib", "SMC", "SMC"], [(0, 1), (2, 3)])
    sp.neighbor_contamination(a, SPILL_MARKERS)
    assert a.obs["sam_neighbor_contamination"].iloc[0] == pytest.approx(0.2)


def test_neighbor_contamination_graph_direction():
    """Row i of the graph lists i's neighbours (SegTraQ's convention), on an asymmetric graph."""
    from spatial_anno_metrics import spatial as sp

    var = ["CD3E", "COL1A1"]
    X = [[9, 1], [0, 10]]
    a = _tiny(X, var, ["T", "Fib"], [(0, 1)], directed=True)     # c0 lists c1; c1 lists nobody
    sp.neighbor_contamination(a, SPILL_MARKERS)
    np.testing.assert_allclose(a.obs["sam_neighbor_contamination"].to_numpy(), [0.1, np.nan])
    b = _tiny(X, var, ["T", "Fib"], [(1, 0)], directed=True)     # c1 lists c0
    sp.neighbor_contamination(b, SPILL_MARKERS)
    np.testing.assert_allclose(b.obs["sam_neighbor_contamination"].to_numpy(), [np.nan, 0.0])


def test_neighbor_contamination_matrix_is_the_mean():
    """Three T cells next to a Fib: 0.1, 0.1, 0.4 -> 0.2 (the median would be 0.1)."""
    from spatial_anno_metrics import spatial as sp

    X = [[9, 1], [9, 1], [6, 4], [0, 10]]
    a = _tiny(X, ["CD3E", "COL1A1"], ["T", "T", "T", "Fib"], [(0, 3), (1, 3), (2, 3)])
    out = sp.neighbor_contamination(a, SPILL_MARKERS)
    assert out["source_target"].loc["Fib", "T"] == pytest.approx(0.2)


def test_neighbor_contamination_skips_nested_sets_like_mecr():
    """T = {CD3E} sits inside CD4 T = {CD3E, CD4}: hierarchical, not spillover (mecr skips the pair too)."""
    from spatial_anno_metrics import spatial as sp

    X = [[8, 2, 0], [8, 2, 0], [0, 0, 10]]                      # CD3E CD4 COL1A1
    a = _tiny(X, ["CD3E", "CD4", "COL1A1"], ["T", "CD4 T", "Fib"], [(0, 1), (1, 2)])
    out = sp.neighbor_contamination(a, {"T": ["CD3E"], "CD4 T": ["CD3E", "CD4"], "Fib": ["COL1A1"]})
    assert np.isnan(a.obs["sam_neighbor_contamination"].iloc[0])   # its only neighbour is its subtype
    assert np.isnan(out["source_target"].loc["CD4 T", "T"])


def test_neighbor_contamination_int_keys():
    """Integer labels with an int-keyed marker_sets used to give all-NaN and a 0x0 matrix."""
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    want = sp.neighbor_contamination(a, SPILL_MARKERS)["mean_contamination"]
    code = {"T": 0, "Fib": 1, "SMC": 2}
    b = _spill_tissue()
    b.obs["cell_type"] = [code[x] for x in b.obs["cell_type"]]
    out = sp.neighbor_contamination(b, {code[k]: v for k, v in SPILL_MARKERS.items()})
    assert out["mean_contamination"] == pytest.approx(want) and out["n_evaluable"] == 6


def test_neighbor_contamination_dense_graph_and_sparse_counts():
    from scipy.sparse import csr_matrix

    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    sp.neighbor_contamination(a, SPILL_MARKERS)
    b = _spill_tissue()
    b.obsp["spatial_connectivities"] = b.obsp["spatial_connectivities"].toarray()
    b.X = csr_matrix(b.X)
    sp.neighbor_contamination(b, SPILL_MARKERS)
    np.testing.assert_allclose(b.obs["sam_neighbor_contamination"], a.obs["sam_neighbor_contamination"])


def test_neighbor_contamination_explicit_negatives():
    """SegTraQ's marker format: only the target's listed negatives that a neighbour type lists as positive."""
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    markers = {"T": {"positive": ["CD3E"], "negative": ["COL1A1"]},
               "Fib": {"positive": ["COL1A1", "ACTA2"]}, "SMC": {"positive": ["ACTA2", "MYH11"]}}
    sp.neighbor_contamination(a, markers)
    v = a.obs["sam_neighbor_contamination"].to_numpy()
    assert v[0] == pytest.approx(0.1)                            # COL1A1 from c1 still counts
    assert v[6] == pytest.approx(0.0)                            # its 2 ACTA2 are not a listed T negative (0.2 implicit)
    assert v[1] == pytest.approx(0.0)                            # Fib without negatives: implicit rule


def test_neighbor_contamination_reports_dictionary_coverage():
    from spatial_anno_metrics import spatial as sp

    a = _spill_tissue()
    a.obs["cell_type"] = ["T", "Fib", "T", "T", "T", "Fib", "T", "Fib", "Mystery", "T"]
    out = sp.neighbor_contamination(a, {"T": ["CD3E"], "Fib": ["COL1A1", "ACTA2"]})
    assert out["labels_without_markers"] == ["Mystery"] and out["no_dict_frac"] == pytest.approx(0.1)


def test_neighbor_contamination_delaunay_graph():
    import squidpy as sq

    from spatial_anno_metrics import spatial as sp

    a = _scatter(n=300)
    rng = np.random.default_rng(3)
    X = rng.poisson(1.0, (300, 4)).astype("float32")
    a = anndata.AnnData(X=X, obs=a.obs, obsm={"spatial": a.obsm["spatial"]})
    a.var_names = ["a1", "b1", "c1", "bg"]
    ms = {"A": ["a1"], "B": ["b1"], "C": ["c1"]}
    out = sp.neighbor_contamination(a, ms, graph="delaunay")
    assert out["graph"] == "delaunay"
    b = a.copy()
    sq.gr.spatial_neighbors(b, coord_type="generic", delaunay=True)
    ref = sp.neighbor_contamination(b, ms, graph="obsp")
    assert out["mean_contamination"] == pytest.approx(ref["mean_contamination"])
    knn = sp.neighbor_contamination(a.copy(), ms)
    assert knn["graph"] == "knn(n_neighs=6)" and knn["mean_contamination"] != pytest.approx(
        out["mean_contamination"])


def _two_block_probe(spill: bool, seed=0, n_rows=20, coexpr_border=37.0):
    """A 20-column grid (rook contacts): columns 0-9 A, 10-19 B. B carries its marker 'b'; A cells in column 9
    touch B. co-expression: every A cell expresses b (~3.4 % of counts, 3.7 % at the border); spillover: only
    the border A cells carry b (~4 %), the others ~0.1 %."""
    from scipy.sparse import csr_matrix

    rng = np.random.default_rng(seed)
    cells = [(c, r) for c in range(20) for r in range(n_rows)]
    idx = {p: k for k, p in enumerate(cells)}
    e = [(idx[(c, r)], idx[(c + dc, r + dr)]) for c, r in cells for dc, dr in ((1, 0), (0, 1))
         if (c + dc, r + dr) in idx]
    i, j = zip(*e)
    n = len(cells)
    col = np.array([c for c, _ in cells])
    is_a = col < 10
    border = col == 9
    mb = np.where(is_a, np.where(border, 40.0, 1.0) if spill else np.where(border, coexpr_border, 34.0), 500.0)
    X = np.column_stack([rng.poisson(np.where(is_a, 1000 - mb - 60, 0)),
                         rng.poisson(np.where(is_a, 0, 500)) + rng.poisson(np.where(is_a, mb, 0)),
                         rng.poisson(60, n)]).astype("float32")
    a = anndata.AnnData(X=X)
    a.var_names = ["a", "b", "bg"]
    a.obs_names = [f"c{k}" for k in range(n)]
    a.obs["cell_type"] = np.where(is_a, "A", "B")
    a.obsm["spatial"] = np.array(cells, dtype=float)
    a.obsp["spatial_connectivities"] = csr_matrix((np.ones(2 * len(i)), (i + j, j + i)), shape=(n, n))
    return a


def test_neighbor_contamination_baseline_separates_spillover_from_coexpression():
    """docs/ROADMAP.md #19: raw B->A 0.037 (co-expression) vs 0.040 (spillover) look alike; the same-type
    cells away from the border (0.034 vs 0.001) and the permutation p separate them."""
    from spatial_anno_metrics import spatial as sp

    ms = {"A": ["a"], "B": ["b"]}
    co = sp.neighbor_contamination(_two_block_probe(spill=False), ms, n_perm=999, seed=0)
    sp_ = sp.neighbor_contamination(_two_block_probe(spill=True), ms, n_perm=999, seed=0)
    for out, raw, base in ((co, 0.037, 0.034), (sp_, 0.040, 0.001)):
        assert out["source_target"].loc["B", "A"] == pytest.approx(raw, abs=0.002)
        assert out["baseline"].loc["B", "A"] == pytest.approx(base, abs=0.001)
        assert out["excess"].loc["B", "A"] == pytest.approx(
            out["source_target"].loc["B", "A"] - out["baseline"].loc["B", "A"])
        assert out["n_baseline"].loc["B", "A"] == 180
    assert sp_["excess"].loc["B", "A"] > 0.035 and co["excess"].loc["B", "A"] < 0.005
    assert sp_["p_value"].loc["B", "A"] == pytest.approx(1 / 1000)
    assert co["p_value"].loc["B", "A"] > sp_["p_value"].loc["B", "A"]
    pure = sp.neighbor_contamination(_two_block_probe(spill=False, coexpr_border=34.0), ms, n_perm=999, seed=0)
    assert pure["p_value"].loc["B", "A"] > 0.05                  # co-expression alone, no border excess
    off = sp.neighbor_contamination(_two_block_probe(spill=True), ms, n_perm=0)
    assert np.isnan(off["p_value"].loc["B", "A"]) and off["excess"].loc["B", "A"] > 0.035


def test_contamination_p_value_is_one_sided():
    """A DEPLETION next to the source (fewer source-marker counts than far away) is not spillover: p ~ 1."""
    from spatial_anno_metrics.spatial import _perm_p

    rng = np.random.default_rng(0)
    far, near = rng.normal(0.03, 0.005, 100), rng.normal(0.01, 0.005, 20)
    assert _perm_p(near, far, 199, np.random.default_rng(0)) > 0.9
    assert _perm_p(far[:20], near, 199, np.random.default_rng(0)) == pytest.approx(1 / 200)


def test_neighbor_contamination_split_and_cellspa_distance_readouts():
    """Roadmap #19 follow-ups: SPLIT (contaminant share vs source neighbour share) and CellSPA (contaminant
    share vs distance to nearest source cell) separate boundary spillover from uniform co-expression."""
    from spatial_anno_metrics import spatial as sp

    ms = {"A": ["a"], "B": ["b"]}
    sp_ = sp.neighbor_contamination(_two_block_probe(spill=True), ms, n_perm=99, seed=0)
    co = sp.neighbor_contamination(_two_block_probe(spill=False, coexpr_border=34.0), ms, n_perm=99, seed=0)

    # SPLIT readout: under boundary spillover from B into A, contaminant share rises with B's neighbour share
    assert sp_["split_slope"].loc["B", "A"] > 0.10
    assert sp_["split_corr"].loc["B", "A"] > 0.80
    assert abs(co["split_slope"].loc["B", "A"]) < 0.02
    assert abs(co["split_corr"].loc["B", "A"]) < 0.15

    # CellSPA distance readout: under spillover, B marker share in A cells decays with distance from B
    assert sp_["distance_corr"].loc["B", "A"] < -0.20
    assert abs(co["distance_corr"].loc["B", "A"]) < 0.15
    curve = sp_["distance_curve"]
    ba = curve[(curve["source"] == "B") & (curve["target"] == "A")].sort_values("bin")
    assert len(ba) == 4 and int(ba["n_cells"].sum()) == 200
    assert ba.iloc[0]["mean_contamination"] > ba.iloc[-1]["mean_contamination"] + 0.005


def test_spatial_metrics_skip_when_no_coordinates_or_graph():
    """Without obsm['spatial'] or obsp['spatial_connectivities'], all three spatial label metrics skip cleanly."""
    from spatial_anno_metrics import spatial as sp

    a = _blocks(n_side=6)
    del a.obsm["spatial"]
    a.X = np.ones((a.n_obs, 3), dtype="float32")
    a.var_names = ["g0", "g1", "g2"]
    for out in (sp.spatial_coherence(a), sp.neighborhood_sanity(a),
                sp.neighbor_contamination(a, {"A": ["g0"], "B": ["g1"]})):
        assert out["status"] == "skipped" and "spatial" in out["reason"]


def test_neighborhood_sanity_ignores_nan_from_isolated_type_and_zero_degree_skips():
    from scipy.sparse import csr_matrix

    from spatial_anno_metrics import spatial as sp

    # Add 1 isolated cell of type 'C' (degree 0) to a clean two-block section: A and B must not be flagged suspicious
    a = _blocks(n_side=6)
    A_orig = sp._spatial_graph(a, 4)[0].tocsr()
    n = a.n_obs
    A_ext = csr_matrix((A_orig.data, A_orig.indices, np.append(A_orig.indptr, A_orig.indptr[-1])), shape=(n + 1, n + 1))
    b = anndata.AnnData(X=np.zeros((n + 1, 2), dtype="float32"))
    b.obs["cell_type"] = list(a.obs["cell_type"].astype(str)) + ["C"]
    b.obsp["spatial_connectivities"] = A_ext
    out = sp.neighborhood_sanity(b, graph="obsp")
    assert out["status"] == "ok" and out["suspicious"] == []

    # All-zero graph skips spatial_coherence instead of returning mean_coherence = NaN
    b.obsp["spatial_connectivities"] = csr_matrix((n + 1, n + 1))
    assert sp.spatial_coherence(b, graph="obsp")["status"] == "skipped"

    # Small n_obs <= n_neighs clamps kNN instead of crashing
    small = _blocks(n_side=2)
    assert sp.spatial_coherence(small, n_neighs=20)["status"] == "ok"
