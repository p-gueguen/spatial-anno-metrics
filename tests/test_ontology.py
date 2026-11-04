"""Tests for ontology.py - do label centroid distances follow Cell Ontology distances?"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

anndata = pytest.importorskip("anndata")

OBO = Path(__file__).parent / "data" / "cl_tiny.obo"
CL = {"T": "CL:0000084", "B": "CL:0000236", "mono": "CL:0000576", "macrophage": "CL:0000235",
      "keratinocyte": "CL:0000312", "epithelial": "CL:0000066", "fibroblast": "CL:0000057"}
# centroids laid out like the ontology: lymphoid pair, myeloid pair, then the non-immune trio far away
CENTRES = {"T": (0, 0, 0), "B": (1, 0, 0), "mono": (4, 0, 0), "macrophage": (3.5, 1, 0),
           "keratinocyte": (0, 8, 0), "epithelial": (1, 8, 0), "fibroblast": (2.5, 8, 0)}


def _section(centres=CENTRES, n=30, seed=0):
    rng = np.random.default_rng(seed)
    labels = [t for t in centres for _ in range(n)]
    E = np.vstack([rng.normal(centres[t], 0.3, size=(n, 3)) for t in centres])
    a = anndata.AnnData(X=np.zeros((len(labels), 2), dtype="float32"))
    a.obs["cell_type"] = labels
    a.obsm["X_pca"] = E
    return a


def test_read_obo_parses_is_a_only_and_skips_obsolete_and_typedefs():
    from spatial_anno_metrics.ontology import _read_obo

    parents = _read_obo(OBO)
    assert len(parents) == 13                                     # 14 [Term] stanzas, one obsolete
    assert parents["CL:0000235"] == ["CL:0000763", "CL:0000145"]  # two is_a; trailing {..} and ! dropped
    assert parents["CL:0000000"] == []
    assert "CL:0000999" not in parents and "develops_from" not in parents


def test_ontology_distance_is_the_shortest_path_through_a_common_ancestor():
    from spatial_anno_metrics.ontology import _ontology_distances, _read_obo

    d = _ontology_distances(["CL:0000084", "CL:0000236", "CL:0000235", "CL:0000312"], _read_obo(OBO))
    assert d[0, 1] == 2                                           # T - lymphocyte - B
    assert d[0, 2] == 4                                           # 2nd parent: T-lymph-leuk-APC-mac, not 5 via myeloid
    assert d[0, 3] == 6                                           # T ... cell ... keratinocyte
    assert d[2, 3] == 5                                           # macrophage - myeloid - hemat. - cell - epi - kerat.
    assert (d == d.T).all() and (np.diag(d) == 0).all()


def test_consistent_layout_scores_high_with_a_small_permutation_p():
    import spatial_anno_metrics as sam

    r = sam.ontology_consistency(_section(), "cell_type", OBO, CL, n_perm=199)
    assert r["status"] == "ok" and r["n_labels"] == 7 and r["n_perm"] == 199
    assert r["rho"] > 0.6 and r["p_value"] <= 0.01
    assert r["per_label"]["T"]["nearest_embedding"] == "B" and r["per_label"]["T"]["nearest_ontology"] == ["B"]


def test_planted_macrophage_inside_epithelium_is_named():
    """A 'macrophage' centroid sitting in the epithelium lowers rho and is the label whose removal helps most."""
    import spatial_anno_metrics as sam

    good = sam.ontology_consistency(_section(), "cell_type", OBO, CL, n_perm=199)
    bad = sam.ontology_consistency(_section({**CENTRES, "macrophage": (0.5, 8.5, 0)}), "cell_type", OBO, CL,
                                   n_perm=199)
    assert bad["rho"] < good["rho"] - 0.2
    assert bad["most_inconsistent"] == "macrophage"
    mac = bad["per_label"]["macrophage"]
    assert mac["nearest_embedding"] in {"keratinocyte", "epithelial"} and mac["nearest_ontology"] == ["mono"]
    assert mac["loo_delta"] == max(d["loo_delta"] for d in bad["per_label"].values()) > 0


def test_unstructured_embedding_is_not_significant():
    """Negative control: centroids placed at random give a permutation p above 0.05. This layout has a small
    POSITIVE rho, so only the permutation null (not the sign of rho) keeps p high."""
    import spatial_anno_metrics as sam

    rng = np.random.default_rng(21)
    r = sam.ontology_consistency(_section({t: tuple(rng.uniform(0, 8, 3)) for t in CENTRES}), "cell_type", OBO, CL,
                                 n_perm=199)
    assert r["status"] == "ok" and 0.1 < r["rho"] < 0.4 and r["p_value"] > 0.05


def test_min_cells_excludes_a_label_and_nearest_ontology_keeps_ties():
    """A label below min_cells is excluded; with the monocyte gone, four labels tie as the macrophage's nearest
    ontology term (4 is_a steps each) and all four are listed."""
    import spatial_anno_metrics as sam

    a = _section()
    a = a[~((a.obs["cell_type"] == "mono").to_numpy() & (np.arange(a.n_obs) % 30 >= 5))].copy()   # 5 mono cells
    r = sam.ontology_consistency(a, "cell_type", OBO, CL, n_perm=19, min_cells=10)
    assert r["n_labels"] == 6 and r["excluded"] == {"mono": "5 cells < min_cells=10"}
    assert sorted(r["per_label"]["macrophage"]["nearest_ontology"]) == ["B", "T", "epithelial", "fibroblast"]
    assert sam.ontology_consistency(a, "cell_type", OBO, CL, n_perm=19, min_cells=5)["n_labels"] == 7


def test_excluded_labels_and_skips():
    import spatial_anno_metrics as sam

    a = _section()
    cl = {**CL, "B": "CL:9999999"}
    del cl["fibroblast"]
    r = sam.ontology_consistency(a, "cell_type", OBO, cl, n_perm=19)
    assert r["status"] == "ok" and r["n_labels"] == 5
    assert set(r["excluded"]) == {"B", "fibroblast"}
    assert "not in the ontology" in r["excluded"]["B"] and "no CL id" in r["excluded"]["fibroblast"]
    s = sam.ontology_consistency(a, "cell_type", OBO, {k: CL[k] for k in ("T", "B", "mono")})
    assert s["status"] == "skipped" and "need >= 4" in s["reason"] and len(s["excluded"]) == 4
    with pytest.raises(ValueError, match="missing"):
        b = a.copy()
        b.obs["cell_type"] = b.obs["cell_type"].astype(object)
        b.obs.iloc[0, 0] = None
        sam.ontology_consistency(b, "cell_type", OBO, CL)


def test_p_value_matches_the_exact_mantel_permutation_distribution():
    """With 5 labels the Mantel null has 5! = 120 relabellings: the sampled p converges on the exact one, which
    permutes rows AND columns of the ontology distance matrix together (rows only would give 0.83 here)."""
    from itertools import permutations

    from scipy.stats import spearmanr

    import spatial_anno_metrics as sam
    from spatial_anno_metrics.ontology import _ontology_distances, _read_obo

    labs = ["T", "B", "mono", "macrophage", "keratinocyte"]
    a = _section(dict(zip(labs, map(tuple, np.random.default_rng(3).uniform(0, 8, (5, 3))))))
    r = sam.ontology_consistency(a, "cell_type", OBO, {k: CL[k] for k in labs}, n_perm=5000)
    E = np.asarray(a.obsm["X_pca"])
    C = np.vstack([E[(a.obs["cell_type"] == lab).to_numpy()].mean(0) for lab in sorted(labs)])
    De = np.sqrt(((C[:, None] - C[None]) ** 2).sum(-1))
    Do = _ontology_distances([CL[k] for k in sorted(labs)], _read_obo(OBO))
    iu = np.triu_indices(5, 1)
    rho = spearmanr(De[iu], Do[iu])[0]
    null = [spearmanr(De[iu], Do[np.ix_(p, p)][iu])[0] for p in map(list, permutations(range(5)))]
    exact = np.mean(np.asarray(null) >= rho - 1e-12)
    assert r["rho"] == pytest.approx(rho)
    assert 0.4 < exact < 0.65 and r["p_value"] == pytest.approx(exact, abs=0.03)


def test_ontology_consistency_on_yerly2022_bcc_real_data():
    """Real scRNA-seq (Yerly et al. 2022 human BCC, 2,620 cells x 2,000 HVGs): lymphoid (Tcell, Bcell), myeloid
    (MoMacDC), epithelial (Normal_Kerat) and stromal (Fibroblast) PCA centroids follow the Cell Ontology, and a
    planted cross-lineage swap (Tcell <-> Normal_Kerat) flips rho negative and names the swapped type."""
    import pandas as pd
    import scipy.sparse as sp

    import spatial_anno_metrics as sam

    gdir = Path(__file__).parent / "goldens" / "isc_yerly2022"
    with np.load(gdir / "counts.npz", allow_pickle=False) as z:
        X = sp.csr_matrix((z["data"].astype(np.float64), z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.read_csv(gdir / "obs.csv")
    obs.index = obs["cell_id"].astype(str)
    tot = np.asarray(X.sum(axis=1)).ravel()
    scale = np.divide(1e4, tot, out=np.zeros_like(tot), where=tot > 0)
    X.data *= np.repeat(scale, np.diff(X.indptr))
    np.log1p(X.data, out=X.data)
    a = anndata.AnnData(X=X, obs=obs)

    cl_map = {
        "Tcell": "CL:0000084",
        "Bcell": "CL:0000236",
        "MoMacDC": "CL:0000235",
        "Normal_Kerat": "CL:0000312",
        "Fibroblast": "CL:0000057",
    }
    good = sam.ontology_consistency(a, "annotation", OBO, cl_map, n_perm=499, seed=0)
    assert good["status"] == "ok" and good["n_labels"] == 5
    assert good["rho"] > 0.45
    assert good["per_label"]["Tcell"]["nearest_embedding"] == "Bcell"
    assert good["per_label"]["Bcell"]["nearest_embedding"] == "Tcell"
    assert good["per_label"]["Normal_Kerat"]["nearest_embedding"] == "Fibroblast"

    # On the 4-type immune + keratinocyte subset, rho > 0.75
    four = sam.ontology_consistency(a, "annotation", OBO, {k: v for k, v in cl_map.items() if k != "Fibroblast"},
                                    n_perm=199, seed=0)
    assert four["rho"] > 0.75

    # Swapping Tcell and Normal_Kerat flips rho negative and flags Normal_Kerat
    swapped = a.copy()
    swapped.obs["annotation"] = swapped.obs["annotation"].replace(
        {"Tcell": "Normal_Kerat", "Normal_Kerat": "Tcell"})
    bad = sam.ontology_consistency(swapped, "annotation", OBO, cl_map, n_perm=499, seed=0)
    assert bad["rho"] < -0.35 and bad["rho"] < good["rho"] - 0.8
    assert bad["most_inconsistent"] in {"Normal_Kerat", "Tcell"}


def test_ontology_rwr_distance_resolves_multi_parent_ties_and_flags_swaps():
    """distance='rwr' (scGraph-OntoRWR) uses Random-Walk-with-Restart stationary distributions to break ties
    where shortest-path hop counts are equal (T -> mono vs T -> epithelial both 5 hops on cl_tiny.obo)."""
    import spatial_anno_metrics as sam
    from spatial_anno_metrics.ontology import _ontology_distances, _ontology_rwr_distances, _read_obo

    parents = _read_obo(OBO)
    ids = ["CL:0000084", "CL:0000236", "CL:0000576", "CL:0000235", "CL:0000312", "CL:0000066", "CL:0000057"]
    d_path = _ontology_distances(ids, parents)
    d_rwr = _ontology_rwr_distances(ids, parents, restart=0.5)
    # Under shortest-path hop count, T (0) is 5 hops from mono (2), epithelial (5), and fibroblast (6)
    assert d_path[0, 2] == d_path[0, 5] == d_path[0, 6] == 5.0
    # Under RWR, T shares hematopoietic/leukocyte paths with mono, so d_rwr(T, mono) < d_rwr(T, epithelial)
    assert d_rwr[0, 2] < d_rwr[0, 5] < d_rwr[0, 4]

    good = sam.ontology_consistency(_section(), "cell_type", OBO, CL, distance="rwr", n_perm=199)
    assert good["status"] == "ok" and good["distance"] == "rwr" and good["rho"] > 0.6 and good["p_value"] <= 0.03
    bad = sam.ontology_consistency(_section({**CENTRES, "macrophage": (0.5, 8.5, 0)}), "cell_type", OBO, CL,
                                   distance="rwr", n_perm=199)
    assert bad["rho"] < good["rho"] - 0.2 and bad["most_inconsistent"] == "macrophage"
    with pytest.raises(ValueError, match="distance"):
        sam.ontology_consistency(_section(), "cell_type", OBO, CL, distance="bad")
    with pytest.raises(ValueError, match="rwr_restart"):
        sam.ontology_consistency(_section(), "cell_type", OBO, CL, distance="rwr", rwr_restart=1.0)
