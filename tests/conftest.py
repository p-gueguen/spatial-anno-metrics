"""Shared fixtures: well-separated labelled blobs and a label shuffler, as factories so a test picks the size."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_SRC = str(Path(__file__).resolve().parents[1] / "src")
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)


def _blobs(n_per=80, k=4, dim=10, sep=8.0, seed=0):
    """k well-separated Gaussian blobs in a 10-d embedding, correctly labeled T0..T{k-1}."""
    anndata = pytest.importorskip("anndata")
    rng = np.random.default_rng(seed)
    centers = rng.normal(0, sep, size=(k, dim))
    X, lab = [], []
    for t in range(k):
        X.append(rng.normal(centers[t], 1.0, size=(n_per, dim)))
        lab += [f"T{t}"] * n_per
    X = np.vstack(X).astype("float32")
    a = anndata.AnnData(X=X)
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obsm["X_pca"] = X
    a.obs["cell_type"] = pd.Categorical(lab)
    return a


def _shuffle_labels(a, seed=1):
    a2 = a.copy()
    rng = np.random.default_rng(seed)
    a2.obs["cell_type"] = pd.Categorical(rng.permutation(a2.obs["cell_type"].to_numpy()))
    return a2


@pytest.fixture
def blobs():
    """``blobs(n_per=80, k=4, dim=10, sep=8.0, seed=0)`` -> AnnData labelled ``obs['cell_type']``, X_pca = X."""
    return _blobs


@pytest.fixture
def shuffle_labels():
    """``shuffle_labels(a, seed=1)`` -> a copy of ``a`` with ``obs['cell_type']`` permuted."""
    return _shuffle_labels
