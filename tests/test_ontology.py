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
    """Negative control: centroids placed at random give a permutation p well above 0.05."""
    import spatial_anno_metrics as sam

    rng = np.random.default_rng(3)
    r = sam.ontology_consistency(_section({t: tuple(rng.uniform(0, 8, 3)) for t in CENTRES}), "cell_type", OBO, CL,
                                 n_perm=199)
    assert r["status"] == "ok" and r["p_value"] > 0.05


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
