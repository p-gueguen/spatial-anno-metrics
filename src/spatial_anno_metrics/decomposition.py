"""Decomposition metrics: scores read from per-cell cell-type weights of a deconvolution (e.g. RCTD full mode).

  * ``decomposition_purity`` - 1 - the largest normalised weight per cell: the share of a cell's expression
    the deconvolution assigns outside its dominant type. Spatial Touchstone (Plummer et al., Nat Biotechnol
    2025) used it as a segmentation-error readout; the computation follows SegTraQ ``sp.rctd_contamination``
    (github.com/LazDaria/SegTraQ, commit 05e90b6, MIT), which runs RCTD itself; here the weights are an input.

Depends on: numpy (pandas for DataFrame weights).
"""
from __future__ import annotations

import numpy as np

from ._utils import _store
from ._validate import skip

__all__ = ["decomposition_purity"]


def decomposition_purity(adata, weights_key: str, *, key_added: str = "sam_decomposition_purity") -> dict:
    """Per-cell share of the deconvolution weight outside the dominant type.

    ``value = 1 - max_t w_t / sum_t w_t``: 0 for a cell fully assigned to one type, ``1 - 1/T`` for even
    weights over T types. Higher means more mixed, despite the name (Spatial Touchstone's), so read it as
    contamination: a high value flags a cell whose segmentation merged neighbours, or a type the reference
    lacks. It depends on the reference and on the deconvolution, not on the labels. Ported from SegTraQ
    ``rctd_contamination`` (MIT): rows with a non-finite or negative weight, or a non-positive total, are
    undefined (NaN); set a row to NaN to drop a non-converged fit.

    Parameters
    ----------
    adata : AnnData
        Cells with the weights in ``obsm[weights_key]``.
    weights_key : str
        ``obsm`` key holding cells x types weights (array or DataFrame), e.g. RCTD full-mode weights.
    key_added : str
        Obs column for the per-cell value.

    Returns
    -------
    dict
        ``status``, ``median_decomposition_purity``, ``mean_decomposition_purity`` (over defined cells),
        ``n_cells``, ``n_undefined``, ``n_types``. A skip, writing nothing, when ``obsm[weights_key]`` is absent
        or no cell is defined.
    """
    if weights_key not in adata.obsm:
        return skip(f"no obsm[{weights_key!r}]")
    raw = adata.obsm[weights_key]
    W = np.asarray(raw.toarray() if hasattr(raw, "toarray") else raw, dtype=float)
    if W.ndim != 2 or W.shape[0] == 0 or W.shape[1] == 0:
        return {**skip(f"obsm[{weights_key!r}] must be a non-empty 2D matrix"), "n_cells": int(len(W))}
    with np.errstate(invalid="ignore", over="ignore"):
        total, top = W.sum(1), W.max(1)
        ok = np.isfinite(W).all(1) & (W >= 0).all(1) & np.isfinite(total) & (total > 0)
    val = np.full(len(W), np.nan)
    val[ok] = 1.0 - top[ok] / total[ok]
    if not ok.any():
        return {**skip(f"no cell has valid weights in obsm[{weights_key!r}]"), "n_cells": len(W)}
    adata.obs[key_added] = val
    out = {"status": "ok", "median_decomposition_purity": float(np.median(val[ok])),
           "mean_decomposition_purity": float(val[ok].mean()), "n_cells": int(len(W)),
           "n_undefined": int((~ok).sum()), "n_types": int(W.shape[1])}
    _store(adata, "decomposition_purity", out, key_added)
    return out
