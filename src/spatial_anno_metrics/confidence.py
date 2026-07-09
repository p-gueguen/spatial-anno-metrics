"""Per-cell confidence / uncertainty (reference-free, given marker sets).

"How sure are we about THIS cell's label?" - the per-cell complement to the dataset-level validity
metrics. All reference-free:

* :func:`per_cell_confidence` - marker-score **margin** (top1 - top2; the margin, not the absolute
  score, is the confidence) + softmax **entropy** (near 1 => the classifier is guessing).
* :func:`label_stability` - annotation-robustness under transcript subsampling: drop a fraction of
  each cell's counts, re-annotate via a caller-supplied ``annotate_fn``, measure the label-flip rate
  (MERFISH: ~20% dropout flips 10-15% of labels). Annotator-agnostic (pass any labeller).

Source: SpatialScribe QC Layer-5. Depends on: numpy (scipy for sparse). See
``docs/cell-annotation-quality-metrics.md`` §3d.
"""
from __future__ import annotations

import numpy as np


def _dense(x):
    return x.toarray() if hasattr(x, "toarray") else np.asarray(x)


def per_cell_confidence(adata, marker_sets: dict[str, list[str]], layer: str | None = None) -> dict:
    """Per-cell marker-score **margin** + **entropy** from panel-restricted lineage marker sets.

    Per set, the score is the cell's mean expression over that set's on-panel genes; the margin is
    ``top1 - top2`` across sets (high = one lineage clearly wins) and the entropy is the normalized
    Shannon entropy of the softmax over set scores (0 = decisive, 1 = uniform / guessing). Writes
    ``obs['marker_margin']`` and ``obs['marker_entropy']``; returns ``{mean_margin, mean_entropy,
    n_sets}``. Reference-free (needs marker sets). ``layer`` picks a counts/normalised layer (else X)."""
    panel = set(adata.var_names)
    mat = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    var_idx = {g: i for i, g in enumerate(adata.var_names)}

    names, cols = [], []
    for name, genes in marker_sets.items():
        gi = [var_idx[g] for g in genes if g in panel]
        if not gi:
            continue
        names.append(name)
        cols.append(_dense(mat[:, gi]).mean(1).ravel().astype(float))
    if len(names) < 2:
        adata.obs["marker_margin"] = 0.0
        adata.obs["marker_entropy"] = 1.0
        return {"mean_margin": 0.0, "mean_entropy": 1.0, "n_sets": len(names)}

    S = np.column_stack(cols)                                  # cells x sets
    srt = np.sort(S, axis=1)
    margin = srt[:, -1] - srt[:, -2]                           # top1 - top2
    Z = S - S.max(1, keepdims=True)
    P = np.exp(Z); P /= P.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(P * np.log(np.clip(P, 1e-12, 1.0))).sum(1) / np.log(P.shape[1])
    adata.obs["marker_margin"] = margin
    adata.obs["marker_entropy"] = ent
    return {"mean_margin": float(np.mean(margin)), "mean_entropy": float(np.mean(ent)), "n_sets": len(names)}


def _subsample_counts(counts, drop_frac: float, rng) -> "np.ndarray":
    """Binomial transcript thinning: keep each count with prob ``1 - drop_frac`` (dense float32)."""
    arr = _dense(counts).astype("float32")
    if drop_frac <= 0:
        return arr
    return rng.binomial(arr.astype(int), 1.0 - drop_frac).astype("float32")


def label_stability(adata, annotate_fn, drop_frac: float = 0.2, reps: int = 5, seed: int = 0,
                    layer: str | None = "counts") -> dict:
    """Annotation robustness under transcript subsampling (reference-free). Baseline labels come from
    ``annotate_fn(adata) -> array[str] of length n_obs``; then for ``reps`` reps, drop ``drop_frac`` of
    each cell's counts, re-run ``annotate_fn`` on the thinned copy, and record the label-flip rate.
    Writes ``obs['label_stability']`` in [0,1] (higher = LESS stable); returns ``{mean_flip_rate,
    reps, drop_frac}``. ``annotate_fn`` is any labeller (marker argmax, a classifier, ...); the thinned
    copy carries the reduced counts in both ``.X`` and (if ``layer`` present) that layer."""
    counts = adata.layers[layer] if (layer and layer in adata.layers) else adata.X
    base = np.asarray(annotate_fn(adata)).astype(str)
    rng = np.random.default_rng(seed)
    flips = np.zeros(adata.n_obs, dtype=float)
    for _ in range(reps):
        tmp = adata.copy()
        thinned = _subsample_counts(counts, drop_frac, rng)
        tmp.X = thinned
        if layer and layer in tmp.layers:
            tmp.layers[layer] = thinned
        lab = np.asarray(annotate_fn(tmp)).astype(str)
        flips += (lab != base).astype(float)
    stability = flips / max(1, reps)
    adata.obs["label_stability"] = stability
    return {"mean_flip_rate": float(np.mean(stability)), "reps": reps, "drop_frac": drop_frac}
