"""Per-cell confidence / uncertainty (reference-free, given marker sets).

"How sure are we about THIS cell's label?" - the per-cell complement to the dataset-level validity
metrics. All reference-free:

* :func:`per_cell_confidence` - marker-score **margin** (top1 - top2; the margin, not the absolute
  score, is the confidence) + softmax **entropy** (near 1 => the classifier is guessing).
* :func:`label_stability` - annotation-robustness under transcript subsampling: drop a fraction of
  each cell's counts, re-annotate via a caller-supplied ``annotate_fn``, measure the label-flip rate
  (MERFISH: ~20% dropout flips 10-15% of labels). Annotator-agnostic (pass any labeller).

Per-cell values go to ``obs[key_added]`` (default ``sam_<name>``) and the scalar summary to
``uns['sam'][<function name>]``. Source: SpatialScribe QC Layer-5. Depends on: numpy (scipy for sparse). See
``docs/cell-annotation-quality-metrics.md`` §3d.
"""
from __future__ import annotations

import numpy as np

from ._utils import _store
from ._validate import get_counts, skip
from .thresholds import value as _thr

__all__ = ["per_cell_confidence", "label_stability"]


def per_cell_confidence(adata, marker_sets: dict[str, list[str]], *, layer: str | None = "counts",
                        key_added: tuple[str, str] = ("sam_marker_margin", "sam_marker_entropy")) -> dict:
    """Per-cell marker-score **margin** + **entropy** from panel-restricted lineage marker sets.

    Per set, the score is the cell's mean expression over that set's on-panel genes; the margin is
    ``top1 - top2`` across sets (high = one lineage clearly wins) and the entropy is the normalized
    Shannon entropy of the softmax over set scores (0 = decisive, 1 = uniform / guessing). Reference-free
    (needs marker sets).

    Parameters
    ----------
    adata
        AnnData with the counts.
    marker_sets
        ``{lineage: [genes]}``; genes off the panel are ignored.
    layer
        Count matrix to read (see ``get_counts``).
    key_added
        The two obs columns written: (margin, entropy).

    Returns
    -------
    dict
        ``{status, mean_margin, mean_entropy, n_sets, matrix}``; the per-cell values go to ``obs[key_added]``.
        Fewer than 2 marker sets on the panel -> a skip (no runner-up, so no margin).
    """
    if isinstance(key_added, str) or len(key_added) != 2:
        raise TypeError("key_added takes two column names: (margin column, entropy column)")
    if adata.n_obs == 0:
        return skip("no cells in adata")
    panel = set(adata.var_names)
    mat, matrix = get_counts(adata, layer)
    var_idx = {g: i for i, g in enumerate(adata.var_names)}

    names, cols = [], []
    for name, genes in marker_sets.items():
        gi = [var_idx[g] for g in genes if g in panel]
        if not gi:
            continue
        names.append(name)
        cols.append(np.asarray(mat[:, gi].mean(axis=1), dtype=float).ravel())
    if len(names) < 2:
        return skip(f"need >= 2 marker sets with on-panel genes, found {len(names)}")

    S = np.column_stack(cols)                                  # cells x sets
    srt = np.sort(S, axis=1)
    margin = srt[:, -1] - srt[:, -2]                           # top1 - top2
    Z = S - S.max(1, keepdims=True)
    P = np.exp(Z)
    P /= P.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(P * np.log(np.clip(P, 1e-12, 1.0))).sum(1) / np.log(P.shape[1])
    adata.obs[key_added[0]] = margin
    adata.obs[key_added[1]] = ent
    out = {"status": "ok", "mean_margin": float(np.mean(margin)), "mean_entropy": float(np.mean(ent)),
           "n_sets": len(names), "matrix": matrix}
    _store(adata, "per_cell_confidence", out, key_added)
    return out


def _subsample_counts(counts, drop_frac: float, rng):
    """Binomial transcript thinning: keep each count with prob ``1 - drop_frac`` (preserves CSR sparsity).

    Non-integer values are rounded to the nearest count first (``astype(int)`` floored them, so 0.6 became 0).
    """
    import scipy.sparse as sp

    if sp.issparse(counts):
        csr = counts.tocsr().astype("float32", copy=True)
        if drop_frac <= 0:
            return csr
        csr.data = rng.binomial(np.rint(csr.data).astype(np.int64), 1.0 - drop_frac).astype("float32")
        csr.eliminate_zeros()
        return csr
    arr = np.asarray(counts, dtype="float32")
    if drop_frac <= 0:
        return arr.copy()
    return rng.binomial(np.rint(arr).astype(np.int64), 1.0 - drop_frac).astype("float32")


def label_stability(adata, annotate_fn, drop_frac: float = _thr("confidence", "stability_drop_frac"),
                    reps: int = _thr("confidence", "stability_reps"), seed: int = 0, *,
                    layer: str | None = "counts", key_added: str = "sam_label_stability") -> dict:
    """Annotation robustness under transcript subsampling (reference-free).

    Baseline labels come from ``annotate_fn(adata)``; then for ``reps`` reps, drop ``drop_frac`` of each
    cell's counts, re-run ``annotate_fn`` on the thinned object, and record the label-flip rate.

    The thinned object is a new AnnData, not ``adata.copy()``: it holds the thinned counts in ``X`` and, if
    ``layer`` exists, in ``layers[layer]``, plus copies of ``obs`` and ``var`` and shallow copies of ``uns``,
    ``obsm`` and ``obsp`` (their arrays are shared with ``adata``: a labeller must not modify them in place).
    Other layers and ``raw`` are left out - they were the memory cost of the full copy, and a labeller
    reading them got unthinned data, so its labels could not flip.

    Parameters
    ----------
    adata
        AnnData with the counts.
    annotate_fn
        Any labeller ``f(AnnData) -> array of n_obs labels`` (marker argmax, a classifier, ...).
    drop_frac
        Share of each cell's transcripts dropped per rep (binomial thinning).
    reps
        Thinning repetitions.
    seed
        Seed of the thinning.
    layer
        Count matrix to thin (see ``get_counts``).
    key_added
        obs column for the per-cell flip rate.

    Returns
    -------
    dict
        ``{status, mean_flip_rate, reps, drop_frac, matrix}``; ``obs[key_added]`` holds each cell's flip rate
        in [0, 1] (higher = LESS stable).
    """
    import anndata

    if not np.isfinite(drop_frac) or not (0.0 <= drop_frac <= 1.0):
        raise ValueError(f"drop_frac must be in [0, 1], got {drop_frac}")
    if reps < 1:
        raise ValueError(f"reps must be >= 1, got {reps}")
    if adata.n_obs == 0:
        return skip("no cells in adata")
    counts, matrix = get_counts(adata, layer)
    base = np.asarray(annotate_fn(adata)).astype(str)
    if base.ndim != 1 or len(base) != adata.n_obs:
        raise ValueError(f"annotate_fn must return a 1D array of {adata.n_obs} labels, got shape {base.shape}")
    rng = np.random.default_rng(seed)
    flips = np.zeros(adata.n_obs, dtype=float)
    for _ in range(reps):
        thinned = _subsample_counts(counts, drop_frac, rng)
        tmp = anndata.AnnData(X=thinned, obs=adata.obs.copy(), var=adata.var.copy(), uns=dict(adata.uns),
                              obsm=dict(adata.obsm), obsp=dict(adata.obsp),
                              layers={layer: thinned} if layer and layer in adata.layers else None)
        lab = np.asarray(annotate_fn(tmp)).astype(str)
        if lab.ndim != 1 or len(lab) != adata.n_obs:
            raise ValueError(f"annotate_fn must return a 1D array of {adata.n_obs} labels, got shape {lab.shape}")
        flips += (lab != base).astype(float)
    stability = flips / max(1, reps)
    adata.obs[key_added] = stability
    out = {"status": "ok", "mean_flip_rate": float(np.mean(stability)), "reps": reps, "drop_frac": drop_frac,
           "matrix": matrix}
    _store(adata, "label_stability", out, key_added)
    return out
