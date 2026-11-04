"""Tests for _utils.py - the guarded squidpy import."""
from __future__ import annotations

import sys

import pytest

anndata = pytest.importorskip("anndata")


def test_squidpy_import_failure_becomes_an_import_error(monkeypatch):
    """squidpy can fail with a non-ImportError (a dask NotImplementedError seen in the wild): re-raise it as an
    ImportError naming squidpy, so optional-dependency handling catches it."""
    from spatial_anno_metrics._utils import _squidpy

    class Broken:
        def find_spec(self, name, path=None, target=None):
            if name == "squidpy":
                raise NotImplementedError("dask says no")
            return None

    monkeypatch.delitem(sys.modules, "squidpy", raising=False)
    monkeypatch.setattr(sys, "meta_path", [Broken(), *sys.meta_path])
    with pytest.raises(ImportError, match="squidpy failed to import: dask says no") as e:
        _squidpy()
    assert isinstance(e.value.__cause__, NotImplementedError)


def test_embedding_rejects_nan_or_bad_shape(blobs):
    import numpy as np

    from spatial_anno_metrics._utils import _embedding

    a = blobs(n_per=10, k=2)
    a.obsm["X_pca"][0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN or inf"):
        _embedding(a, "X_pca")
    del a.obsm["X_pca"]
    a.X[0, 0] = np.inf
    with pytest.raises(ValueError, match="NaN or inf"):
        _embedding(a, None)
