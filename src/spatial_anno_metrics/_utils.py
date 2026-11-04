"""Private helpers shared by the metric modules.

Dense conversion, [-1, 1] -> [0, 1] scaling, the cell embedding (``obsm`` or a quick PCA), the guarded squidpy
import, and the ``uns['sam']`` summary.
"""
from __future__ import annotations

import numpy as np

from .thresholds import value as _thr

__all__ = []
_PCA_N = _thr("internal", "pca_n_components")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _dense(X):
    return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)


def _scale01(x: float) -> float:
    """Map a silhouette-like value in [-1, 1] to [0, 1]."""
    return float((np.clip(x, -1.0, 1.0) + 1.0) / 2.0)


def _embedding(adata, embedding: str | None, seed: int = 0) -> np.ndarray:
    """Return the cell embedding: ``obsm[embedding]`` if present, else a quick PCA of ``X``.

    A sparse ``X`` is never densified (1e6 cells x 5k genes would be 40 GB): PCA by implicit mean-centering
    (arpack, scikit-learn >= 1.4), the same components as the dense PCA up to sign and solver tolerance.
    It densifies as before when arpack cannot run (scikit-learn < 1.4, or ``n_components >= min(X.shape)``).
    """
    if embedding and embedding in adata.obsm:
        return np.asarray(adata.obsm[embedding], dtype=float)
    import scipy.sparse as sp
    from sklearn.decomposition import PCA

    X = adata.X
    n_comp = int(min(_PCA_N, X.shape[1], max(2, X.shape[0] - 1)))
    if sp.issparse(X) and n_comp < min(X.shape):
        try:
            return PCA(n_components=n_comp, svd_solver="arpack", random_state=seed).fit_transform(X.astype(float))
        except TypeError:                                    # scikit-learn < 1.4: no sparse PCA
            pass
    return PCA(n_components=n_comp, random_state=seed).fit_transform(_dense(X).astype(float))


def _squidpy():
    """Import squidpy, re-raising any failure as ImportError.

    The import can fail with errors other than ImportError (e.g. a dask NotImplementedError), which would
    escape ``except ImportError`` handlers.
    """
    import contextlib
    import warnings

    with contextlib.suppress(Exception), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import dask.dataframe  # noqa: F401 - prevent spatialdata<0.4 + dask>=2025.1 query-planning error
    try:
        import squidpy
    except Exception as e:  # noqa: BLE001 - deliberately broad, see docstring
        raise ImportError(f"squidpy failed to import: {e}") from e
    return squidpy


def _store(adata, name: str, result: dict, key_added) -> None:
    """Keep a metric's scalar results and the obs column(s) it wrote in ``adata.uns['sam'][name]``.

    Frames, per-type dicts and arrays stay in the returned dict only, so the summary is small and
    h5ad-writable.
    """
    summary = {k: v for k, v in result.items() if isinstance(v, (str, int, float, np.number, np.bool_))}
    summary["key_added"] = key_added if isinstance(key_added, str) else list(key_added)
    adata.uns.setdefault("sam", {})[name] = summary
