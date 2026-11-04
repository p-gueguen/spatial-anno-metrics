"""Tests for _validate.py - label checks, the count-matrix rule, the skip shape and marker-dictionary coverage."""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _adata(X, labels=None):
    a = anndata.AnnData(X=np.asarray(X, dtype="float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.var_names = [f"g{j}" for j in range(a.n_vars)]
    if labels is not None:
        a.obs["cell_type"] = labels
    return a


def _user_warnings(fn):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out = fn()
    return out, [w for w in caught if issubclass(w.category, UserWarning)]


def test_check_labels_rejects_missing_labels_naming_column_and_count():
    from spatial_anno_metrics._validate import check_labels

    a = _adata(np.ones((5, 2)), pd.Series(["A", None, "B", np.nan, "A"], dtype=object).to_numpy())
    with pytest.raises(ValueError, match=r"obs\['cell_type'\] has 2 missing"):
        check_labels(a, "cell_type")
    with pytest.raises(ValueError, match="the label array has 1 missing"):
        check_labels(pd.Series(["A", None, "B"]))


def test_check_labels_returns_str_labels_of_the_used_levels():
    from spatial_anno_metrics._validate import check_labels

    a = _adata(np.ones((4, 2)), pd.Categorical([1, 2, 2, 1], categories=[1, 2, 3]))   # level 3 unused
    y = check_labels(a, "cell_type")
    assert list(y) == ["1", "2", "2", "1"] and all(isinstance(v, str) for v in y)
    assert list(np.unique(y)) == ["1", "2"]
    assert list(check_labels([0, 1, 1])) == ["0", "1", "1"]
    with pytest.raises(ValueError, match="at least 3"):
        check_labels(a, "cell_type", min_labels=3)                  # 2 used labels, not 3 categories


def test_get_counts_picks_the_layer_or_count_like_X_and_names_it():
    import scipy.sparse as sp

    from spatial_anno_metrics._validate import get_counts

    a = _adata(np.arange(12).reshape(4, 3))
    a.layers["counts"] = a.X.copy()
    m, which = get_counts(a)
    assert which == "layers['counts']" and m is a.layers["counts"]
    del a.layers["counts"]
    m, which = get_counts(a)                                         # no counts layer, X holds counts
    assert which == "X" and m is a.X
    assert get_counts(a, layer=None)[1] == "X"
    b = _adata(np.zeros((3, 2)))
    b.X = sp.csr_matrix(np.array([[0, 2], [1, 0], [0, 0]], dtype="float32"))
    assert get_counts(b)[1] == "X"


def test_get_counts_refuses_a_missing_layer_when_X_is_not_counts():
    from spatial_anno_metrics._validate import get_counts

    a = _adata(np.log1p(np.arange(12).reshape(4, 3)))
    with pytest.raises(ValueError, match="layer="):
        get_counts(a)                                                # default 'counts' absent, X log-normalised
    with pytest.raises(ValueError, match="'countz'"):
        get_counts(a, layer="countz")                                # a typo is not silently replaced by X
    b = _adata(-np.ones((3, 2)))
    with pytest.raises(ValueError):
        get_counts(b)                                                # negative integers are not counts either


def test_get_counts_warns_once_on_a_non_integer_matrix():
    from spatial_anno_metrics._validate import get_counts

    a = _adata(np.log1p(np.arange(12).reshape(4, 3)))
    (_, which), warned = _user_warnings(lambda: get_counts(a, layer=None))
    assert which == "X" and len(warned) == 1 and "not raw counts" in str(warned[0].message)
    a.layers["lognorm"] = a.X.copy()
    assert len(_user_warnings(lambda: get_counts(a, layer="lognorm"))[1]) == 1
    b = _adata(np.arange(12).reshape(4, 3))
    assert _user_warnings(lambda: get_counts(b, layer=None))[1] == []


def test_skip_shape():
    from spatial_anno_metrics._validate import skip

    assert skip("why") == {"status": "skipped", "reason": "why"}


def test_dictionary_coverage():
    from spatial_anno_metrics._validate import dictionary_coverage

    out = dictionary_coverage(np.array(["A", "A", "B", "C"]), {"A": ["g1"], "B": [], "D": ["g2"]})
    assert out == {"no_dict_frac": 0.5, "labels_without_markers": ["B", "C"]}   # B empty, C absent
    assert dictionary_coverage(["A"], {"A": ["g1"]}) == {"no_dict_frac": 0.0, "labels_without_markers": []}


def test_get_counts_checks_rows_beyond_the_first():
    """Counts in row 0, a non-integer further down: still not raw counts."""
    from spatial_anno_metrics._validate import get_counts

    X = np.ones((10, 3))
    X[7, 1] = 0.5
    a = _adata(X)
    with pytest.raises(ValueError, match="layer="):
        get_counts(a, "counts")
    _, w = _user_warnings(lambda: get_counts(a, None))
    assert len(w) == 1


def test_get_counts_refuses_nan_or_inf_in_dense_or_sparse_matrix():
    """NaN or +/-inf in the count matrix raises ValueError immediately instead of leaking into metrics."""
    import scipy.sparse as sp

    from spatial_anno_metrics._validate import get_counts

    X = np.ones((2000, 3), dtype="float32")
    X[1234, 1] = np.nan                                              # row 1234 is not in the 1000-row sample
    a = _adata(X)
    with pytest.raises(ValueError, match="NaN or infinite"):
        get_counts(a, layer=None)
    with pytest.raises(ValueError, match="NaN or infinite"):
        get_counts(a, layer="counts")
    a.layers["counts"] = sp.csr_matrix(np.array([[1.0, np.inf], [0.0, 2.0]], dtype="float32")[:a.n_obs]
                                       if a.n_obs == 2 else X)
    with pytest.raises(ValueError, match="NaN or infinite"):
        get_counts(a, layer="counts")


def test_layer_is_keyword_only_in_every_public_function():
    import inspect

    import spatial_anno_metrics as sam

    kinds = {}
    for name in sam.__all__:
        fn = getattr(sam, name)
        if callable(fn) and "layer" in inspect.signature(fn).parameters:
            kinds[name] = inspect.signature(fn).parameters["layer"].kind
    assert {"panel_resolvability", "per_cell_confidence", "label_stability", "inter_sample_consistency"} <= set(kinds)
    assert len(kinds) >= 18
    assert [n for n, k in kinds.items() if k is not inspect.Parameter.KEYWORD_ONLY] == []


def _nan_label_section():
    """40 cells x 30 genes of counts with coordinates, an embedding, a sample column and 2 NaN labels."""
    rng = np.random.default_rng(0)
    X = rng.poisson(2.0, size=(40, 30)).astype("float32")
    a = _adata(X, pd.Series(["A"] * 19 + [None] + ["B"] * 19 + [np.nan], dtype=object).to_numpy())
    a.layers["counts"] = a.X.copy()
    a.obsm["X_pca"] = rng.normal(size=(40, 5))
    a.obsm["spatial"] = rng.uniform(0, 10, size=(40, 2))
    a.obs["sample"] = ["s1", "s2"] * 20
    a.obs["truth"] = ["A"] * 20 + ["B"] * 20
    return a


_MARKERS = {"A": ["g0", "g1"], "B": ["g2", "g3"]}
_NAN_LABEL_CALLS = {
    "internal_validity": lambda sam, a: sam.internal_validity(a, "cell_type"),
    "inter_sample_consistency": lambda sam, a: sam.inter_sample_consistency(a, "cell_type", "sample"),
    "marker_program_fidelity": lambda sam, a: sam.marker_program_fidelity(a, "cell_type", _MARKERS),
    "marker_gene_overlap": lambda sam, a: sam.marker_gene_overlap(a, "cell_type", _MARKERS),
    "avg_bio": lambda sam, a: sam.avg_bio(a, "cell_type", "truth"),
    "panel_resolvability": lambda sam, a: sam.panel_resolvability(a, "cell_type", list(a.var_names)),
    "pmp": lambda sam, a: sam.pmp(a, _MARKERS, "cell_type"),
    "nmp": lambda sam, a: sam.nmp(a, a.copy(), "cell_type", "truth"),
    "spatial_coherence": lambda sam, a: sam.spatial_coherence(a, "cell_type"),
    "neighborhood_sanity": lambda sam, a: sam.neighborhood_sanity(a, "cell_type"),
    "neighbor_contamination": lambda sam, a: sam.neighbor_contamination(a, _MARKERS, "cell_type"),
    "isolated_label_scores": lambda sam, a: sam.isolated_label_scores(a, "cell_type"),
    "clisi": lambda sam, a: sam.clisi(a, "cell_type"),
    "coexpression_similarity": lambda sam, a: sam.coexpression_similarity(a, a.copy(), label_key="cell_type",
                                                                          ref_label_key="truth"),
    "consensus_entropy": lambda sam, a: sam.consensus_entropy(a, ["cell_type", "cell_type"]),
    "ontology_consistency": lambda sam, a: sam.ontology_consistency(
        a, "cell_type", Path(__file__).parent / "data" / "cl_tiny.obo", {"A": "CL:0000084", "B": "CL:0000236"}),
}


@pytest.mark.parametrize("name", sorted(_NAN_LABEL_CALLS))
def test_every_label_metric_refuses_nan_labels(name):
    """Each label metric routes its labels through check_labels (a bypass scores NaN as a type called 'nan')."""
    import spatial_anno_metrics as sam

    with pytest.raises(ValueError, match=r"obs\['cell_type'\] has 2 missing"):
        _NAN_LABEL_CALLS[name](sam, _nan_label_section())
