"""Private helpers shared by the metric modules: dense conversion, [-1, 1] -> [0, 1] scaling, the cell
embedding (``obsm`` or a quick PCA), and the guarded squidpy import."""
from __future__ import annotations

import numpy as np

__all__ = []


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _dense(X):
    return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)


def _scale01(x: float) -> float:
    """Map a silhouette-like value in [-1, 1] to [0, 1]."""
    return float((np.clip(x, -1.0, 1.0) + 1.0) / 2.0)


def _embedding(adata, embedding: str | None, seed: int = 0) -> np.ndarray:
    """Return the cell embedding: ``obsm[embedding]`` if present, else a quick PCA of ``X``."""
    if embedding and embedding in adata.obsm:
        return np.asarray(adata.obsm[embedding], dtype=float)
    from sklearn.decomposition import PCA

    X = _dense(adata.X).astype(float)
    n_comp = int(min(50, X.shape[1], max(2, X.shape[0] - 1)))
    return PCA(n_components=n_comp, random_state=seed).fit_transform(X)


def _squidpy():
    """Import squidpy. It can fail with errors other than ImportError (e.g. a dask NotImplementedError), which
    would escape ``except ImportError`` handlers: re-raise any failure as ImportError."""
    try:
        import squidpy
    except Exception as e:  # noqa: BLE001 - deliberately broad, see docstring
        raise ImportError(f"squidpy failed to import: {e}") from e
    return squidpy
