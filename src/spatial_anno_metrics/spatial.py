"""Spatial label-quality metrics - "does the cell-type label make SPATIAL sense?"

The reference-free spatial regime the package was named for but only covered on the *signal* side
(``signal_qc.moran_signal`` scores the data, not the labels). These score the LABELS against the
tissue geometry:

* :func:`spatial_coherence` - per-cell fraction of spatial neighbors sharing the label + the
  dataset PAS (Proportion of Abnormal Spots).
* :func:`neighborhood_sanity` - flags cell types whose strongest neighborhood-enrichment partner is
  a DIFFERENT type (a systematic-mislabeling smell), via a squidpy permutation null.
* :func:`neighbor_contamination` - per-cell share of counts on markers of the ADJACENT cell types
  (spillover from neighbors) + a directed source x target matrix (SegTraQ ``contamination_strength``).

**Caveat (load-bearing):** rare infiltrating cells are legitimately spatially incoherent. Use these
to DOWN-WEIGHT / flag, never to hard-delete. Source: SpatialScribe QC Layer-6; grounded in Salas
et al. Nat Methods 2025. Depends on: numpy, scipy; squidpy (builds/uses the spatial graph).
See ``docs/cell-annotation-quality-metrics.md`` §3e.
"""
from __future__ import annotations

import numpy as np

__all__ = ["spatial_coherence", "neighborhood_sanity", "neighbor_contamination"]


def _ensure_spatial_graph(adata, n_neighs: int) -> None:
    """Build ``obsp['spatial_connectivities']`` via squidpy kNN if absent (needs ``obsm['spatial']``)."""
    if "spatial_connectivities" in adata.obsp:
        return
    import squidpy as sq
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)


def spatial_coherence(adata, label_key: str = "cell_type", k: int = 15,
                      pas_threshold: float = 0.2) -> dict:
    """Per-cell fraction of same-label spatial neighbors (writes ``obs['spatial_coherence']``), plus the
    dataset **PAS** (Proportion of Abnormal Spots = fraction of cells with < ``pas_threshold`` same-label
    neighbors). Reuses ``obsp['spatial_connectivities']`` (built via squidpy kNN if absent). Reference-free
    (labels + coords). Returns ``{pas, mean_coherence, k, pas_threshold}``. DOWN-WEIGHT ONLY - rare
    infiltrating cells are legitimately incoherent, so this never drops cells."""
    from scipy.sparse import csr_matrix

    _ensure_spatial_graph(adata, n_neighs=k)
    A = adata.obsp["spatial_connectivities"].tocsr()

    labels = adata.obs[label_key].astype(str).to_numpy()
    cats = {c: i for i, c in enumerate(np.unique(labels))}
    lab_idx = np.array([cats[lb] for lb in labels])
    onehot = csr_matrix((np.ones(adata.n_obs), (np.arange(adata.n_obs), lab_idx)),
                        shape=(adata.n_obs, len(cats)))
    same = np.asarray((A @ onehot)[np.arange(adata.n_obs), lab_idx]).ravel()
    deg = np.asarray(A.sum(1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(deg > 0, same / deg, np.nan)
    adata.obs["spatial_coherence"] = frac

    pas = float(np.mean(np.nan_to_num(frac, nan=1.0) < pas_threshold))
    return {"pas": pas, "mean_coherence": float(np.nanmean(frac)), "k": k, "pas_threshold": pas_threshold}


def neighborhood_sanity(adata, label_key: str = "cell_type", n_neighs: int = 6, seed: int = 0) -> dict:
    """A healthy annotation self-associates spatially. Runs squidpy's neighborhood-enrichment
    (permutation z-score) and flags any type whose STRONGEST partner is a different type - a mild
    systematic-mislabeling smell (e.g. two labels that are really one, or a swap). Reference-free
    (labels + coords). Returns ``{categories, self_enrichment (diagonal), suspicious[...]}``."""
    import pandas as pd
    import squidpy as sq

    _ensure_spatial_graph(adata, n_neighs=n_neighs)
    if not isinstance(adata.obs[label_key].dtype, pd.CategoricalDtype):
        adata.obs[label_key] = adata.obs[label_key].astype("category")
    sq.gr.nhood_enrichment(adata, cluster_key=label_key, seed=seed, show_progress_bar=False)
    res = adata.uns.get(f"{label_key}_nhood_enrichment", {})
    z = np.asarray(res.get("zscore", np.array([])))
    cats = list(adata.obs[label_key].cat.categories)
    diag = np.diag(z) if z.ndim == 2 and z.shape[0] == z.shape[1] else np.array([])

    suspicious = []
    if z.ndim == 2 and z.shape[0] > 1:
        for i, c in enumerate(cats):
            j = int(np.argmax(z[i]))
            if j != i:
                suspicious.append({"type": str(c), "stronger_with": str(cats[j]),
                                   "self_z": float(z[i][i]), "cross_z": float(z[i][j])})
    return {"categories": [str(c) for c in cats],
            "self_enrichment": [float(x) for x in diag], "suspicious": suspicious}


def neighbor_contamination(adata, marker_sets: dict[str, list[str]], label_key: str = "cell_type",
                           n_neighs: int = 6, require_neighbor_expression: bool = True) -> dict:
    """Spillover from ADJACENT cells (reference-free: labels + coords + curated ``marker_sets``). For a
    cell of type t, a gene is a local contaminant when it is a marker of a type s touching the cell and not
    a marker of t; with ``require_neighbor_expression`` a touching s cell must also detect it. Writes
    ``obs['neighbor_contamination']`` = the cell's counts on its local contaminants / its total counts (NaN
    when no other marker-set type touches it, or its own type has no marker set). The spatial counterpart
    to :func:`~spatial_anno_metrics.purity.pmp` / ``nmp``, which cannot tell a neighbor's spillover from a
    mislabel. Returns ``{mean_contamination, n_evaluable, matrix, n_matrix}``: ``matrix`` is the source
    (rows) x target (columns) mean per-cell contamination over the ``n_matrix`` evaluable target cells.
    Re-implements SegTraQ's ``contamination_strength`` (github.com/LazDaria/SegTraQ), vectorized over
    cells; reuses ``obsp['spatial_connectivities']`` (squidpy kNN if absent)."""
    import pandas as pd

    _ensure_spatial_graph(adata, n_neighs=n_neighs)
    A = adata.obsp["spatial_connectivities"].tocsc()

    panel = {g: i for i, g in enumerate(adata.var_names)}
    sets = {ct: {g for g in mk if g in panel} for ct, mk in marker_sets.items()}
    genes = sorted(set().union(*sets.values()))
    col = {g: j for j, g in enumerate(genes)}
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    sub = mat[:, [panel[g] for g in genes]]
    Xm = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)        # cells x marker genes
    total = np.asarray(mat.sum(axis=1)).ravel().astype(float)
    labels = adata.obs[label_key].astype(str).to_numpy()
    types = [str(t) for t in np.unique(labels) if sets.get(t)]

    local = np.zeros(Xm.shape, dtype=bool)             # cell x gene: a local contaminant of that cell
    M = pd.DataFrame(np.nan, index=types, columns=types)
    N = pd.DataFrame(0, index=types, columns=types)
    for s in types:
        src = labels == s
        if require_neighbor_expression:                # an s neighbor detects the gene
            nb = np.asarray((A[:, src] @ (Xm[src] > 0)) > 0)
        else:                                          # any s neighbor
            nb = np.broadcast_to(np.asarray(A[:, src].sum(axis=1)).reshape(-1, 1) > 0, Xm.shape)
        for t in types:
            g = [col[x] for x in sorted(sets[s] - sets[t])]
            if t == s or not g:
                continue
            tgt = np.flatnonzero((labels == t) & (total > 0))
            hit = nb[np.ix_(tgt, g)]
            ok = hit.any(axis=1)
            if ok.any():
                M.loc[s, t] = ((Xm[np.ix_(tgt, g)] * hit).sum(axis=1)[ok] / total[tgt][ok]).mean()
                N.loc[s, t] = int(ok.sum())
            local[np.ix_(tgt, g)] |= hit               # a gene two sources share is counted once per cell

    evaluable = local.any(axis=1)
    frac = np.full(adata.n_obs, np.nan)
    frac[evaluable] = (Xm * local).sum(axis=1)[evaluable] / total[evaluable]
    adata.obs["neighbor_contamination"] = frac
    return {"mean_contamination": float(np.nanmean(frac)) if evaluable.any() else float("nan"),
            "n_evaluable": int(evaluable.sum()), "matrix": M, "n_matrix": N}
