"""Input checks shared by the metrics: labels, the count-matrix rule, the skip shape, marker-dictionary coverage.

Every metric routes its labels through :func:`check_labels` and its count matrix through :func:`get_counts`, so
the rules live in one place.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

__all__ = []


def skip(reason: str) -> dict:
    """Return the one "could not score this" shape: ``{"status": "skipped", "reason": reason}``.

    A success carries ``"status": "ok"``. Callers may add context keys (counts, dropped types) next to these two.
    """
    return {"status": "skipped", "reason": reason}


def check_labels(adata_or_values, key: str | None = None, *, min_labels: int = 1) -> np.ndarray:
    """Labels as strings, refusing missing ones.

    ``check_labels(adata, "cell_type")`` reads ``adata.obs["cell_type"]``; ``check_labels(values)`` takes any
    1-D array-like. A NaN / None label raises ``ValueError`` naming the column and the count: as a string it
    would become a cell type called "nan" (pandas < 3) or crash a sort (pandas >= 3). Returns the labels of
    the cells as str, so unused categorical levels do not appear; fewer than ``min_labels`` distinct labels
    raises. A metric that needs >= 2 labels checks that itself and returns :func:`skip` instead of raising.
    """
    raw, where = ((adata_or_values, "the label array") if key is None
                  else (adata_or_values.obs[key], f"adata.obs[{key!r}]"))
    n_na = int(np.sum(pd.isna(raw)))
    if n_na:
        raise ValueError(f"{where} has {n_na} missing (NaN/None) labels of {len(raw)}: drop those cells or "
                         "give them a label (e.g. 'Unassigned') before scoring")
    # same conversions as before 0.7.0 (obs column -> str objects; array-like -> numpy str), so results keep
    # their key types
    y = np.asarray(raw).astype(str) if key is None else raw.astype(str).to_numpy()
    n = len(np.unique(y))
    if n < min_labels:
        raise ValueError(f"{where} has {n} distinct label(s); need at least {min_labels}")
    return y


def _count_like(M) -> bool:
    """Non-negative integers (NaN-free) on up to 1000 evenly spaced rows: the "raw counts" test."""
    import scipy.sparse as sp

    n = M.shape[0]
    rows = np.unique(np.linspace(0, n - 1, num=min(n, 1000)).astype(int))
    sub = M[rows]
    v = np.asarray(sub.data if sp.issparse(sub) else sub)
    return bool(np.all(np.isfinite(v)) and np.all(v >= 0) and np.all(v == np.round(v)))


def matrix_name(adata, layer: str | None) -> str:
    """Return the name :func:`get_counts` reports for ``layer``, without reading the matrix.

    ``"layers['<layer>']"`` when the layer exists, else ``"X"`` (``get_counts`` raises instead when that ``X`` is
    not raw counts).
    """
    return f"layers[{layer!r}]" if layer is not None and layer in adata.layers else "X"


def get_counts(adata, layer: str | None = "counts"):
    """``(matrix, which)``: the count matrix a count-based metric reads, and its name for the result.

    * ``layer=None`` -> ``X``.
    * ``layer`` present -> ``adata.layers[layer]`` (``which = "layers['<layer>']"``).
    * ``layer`` absent -> ``X`` if ``X`` holds raw counts (non-negative integers on a sample of rows),
      ``which = "X"``; otherwise ``ValueError`` asking for ``layer=``. A typo is never silently replaced by
      log-normalised data.

    A chosen matrix that is not raw counts (non-integer or negative values) is used, with one ``UserWarning``.
    Every count-based metric returns ``which`` as ``"matrix"`` in its result.
    """
    which = matrix_name(adata, layer)
    if layer is not None and which == "X":                         # the named layer is missing
        if not _count_like(adata.X):
            named = sorted(k for k in adata.layers if k is not None) or "none"   # anndata >= 0.12 lists X as None
            raise ValueError(f"layer {layer!r} is not in adata.layers ({named}) and X "
                             "is not raw counts (negative or non-integer values): pass layer= naming the "
                             "raw-count layer, or layer=None to read X as it is")
        return adata.X, which
    mat = adata.X if which == "X" else adata.layers[layer]
    if not _count_like(mat):
        warnings.warn(f"{which} is not raw counts (negative or non-integer values) but this metric reads it as "
                      "counts; pass layer= naming the raw-count layer", UserWarning, stacklevel=3)
    return mat, which


def dictionary_coverage(labels, marker_sets: dict) -> dict:
    """Report which labels have no marker set.

    Returns ``{"no_dict_frac": fraction of cells whose label has no entry (or an empty one) in ``marker_sets``,
    "labels_without_markers": [sorted labels]}``. A marker metric silently skips such labels, so report this
    beside it. Pass ``marker_sets`` already restricted to the panel to count a set whose genes are all off-panel
    as missing too.
    """
    y = check_labels(labels)
    missing = sorted(str(t) for t in np.unique(y) if not marker_sets.get(t))
    return {"no_dict_frac": float(np.isin(y, missing).mean()), "labels_without_markers": missing}
