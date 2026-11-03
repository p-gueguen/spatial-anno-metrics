"""Tests for _validate.py - label checks, the count-matrix rule, the skip shape and marker-dictionary coverage."""
from __future__ import annotations

import warnings

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
    with pytest.raises(ValueError, match="has 1 missing"):
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
