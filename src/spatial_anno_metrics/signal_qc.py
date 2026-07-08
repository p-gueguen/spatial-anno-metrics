"""Annotation-independent signal / QC metrics, squidpy-compatible.

Ported from **SpatialQM** (Center for Spatial OMICs, github.com/Center-for-Spatial-OMICs/SpatialQM,
an R package) - the metrics that score the *signal quality* of an imaging-ST section (is the
measured expression real structured biology or background noise?) independent of any cell-type
labels. They complement the contamination metrics in ``purity.py``, the section QC in ``qc.py``,
and the label-quality battery in ``eval_metrics.py``.

  * ``moran_signal``      - Moran's I spatial autocorrelation per gene vs per negative-control probe
                            (SpatialQM ``getMorans``). Real genes are spatially structured (high I);
                            controls are spatially random (I~0). The separation is a reference-free,
                            spatial signal-quality score. Uses squidpy's neighbor graph +
                            ``sq.gr.spatial_autocorr`` (fast, parallel).
  * ``signal_to_noise``   - per-cell log2 ratio of gene expression to the negative-control mean
                            (SpatialQM ``getMeanSignalRatio`` / ``getMaxRatio``).
  * ``sparsity``          - zero fraction of the count matrix, dataset + per cell (``getSparsity``).
  * ``detection_entropy`` - per-cell normalized Shannon entropy over detected genes
                            (``getEntropy``); low = one gene dominates (rRNA / artifact).
  * ``tx_per_area``       - transcript density per unit cell area (``getTxPerArea``; needs cell_area).
  * ``run_signal_qc``     - orchestrator -> one headline dict.

squidpy-compatible: ``moran_signal`` builds / reuses ``adata.obsp['spatial_connectivities']``
(``sq.gr.spatial_neighbors``) and delegates to ``sq.gr.spatial_autocorr``. Everything else is
vectorized sparse ops that scale to 1e5-1e6 cells. Control probes are read from ``var['control']``
(set by ``io.load``); metrics degrade gracefully when controls / coords / areas are absent.

Depends on: numpy, scipy; squidpy (``moran_signal`` only). Heavy imports live inside the functions.
"""
from __future__ import annotations

import numpy as np


def _control_mask(adata) -> np.ndarray:
    if "control" in adata.var.columns:
        return np.asarray(adata.var["control"], dtype=bool)
    return np.zeros(adata.n_vars, dtype=bool)


def _csr(X):
    import scipy.sparse as sp
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


# --------------------------------------------------------------------------- #
# Moran's I gene-vs-control (SpatialQM getMorans) - squidpy-native
# --------------------------------------------------------------------------- #
def moran_signal(adata, n_neighs: int = 6, genes: list[str] | None = None,
                 n_perms: int | None = None, seed: int = 0) -> dict:
    """Moran's I spatial autocorrelation split into panel genes vs negative-control probes.

    Real genes carry spatial structure (I well above 0); control probes do not (I~0), so the
    **separation** between the two distributions is a reference-free spatial signal-quality score.
    Builds ``obsp['spatial_connectivities']`` via ``sq.gr.spatial_neighbors`` if absent, then calls
    the fast, parallel ``sq.gr.spatial_autocorr``. Returns medians, the separation, and the fraction
    of genes whose I exceeds the controls' 95th percentile (needs coords in ``obsm['spatial']``).
    """
    import squidpy as sq

    if "spatial_connectivities" not in adata.obsp:
        sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)
    use_genes = list(genes) if genes is not None else list(map(str, adata.var_names))
    sq.gr.spatial_autocorr(adata, mode="moran", genes=use_genes, n_perms=n_perms, seed=seed)
    mi = adata.uns["moranI"]["I"]

    ctrl = _control_mask(adata)
    names = np.asarray(adata.var_names, dtype=str)
    gene_set, ctrl_set = set(names[~ctrl]), set(names[ctrl])
    gI = mi[[g for g in mi.index if g in gene_set]].to_numpy(dtype=float)
    cI = mi[[g for g in mi.index if g in ctrl_set]].to_numpy(dtype=float)

    out: dict = {
        "median_gene_moran": float(np.median(gI)) if gI.size else float("nan"),
        "median_control_moran": float(np.median(cI)) if cI.size else None,
        "n_genes": int(gI.size), "n_controls": int(cI.size),
    }
    if cI.size:
        out["separation"] = out["median_gene_moran"] - float(np.median(cI))
        out["frac_genes_above_control_p95"] = float((gI > np.quantile(cI, 0.95)).mean())
    return out


# --------------------------------------------------------------------------- #
# Signal-to-noise vs controls (SpatialQM getMeanSignalRatio / getMaxRatio)
# --------------------------------------------------------------------------- #
def signal_to_noise(adata) -> dict:
    """Per-cell log2 ratio of gene expression to the negative-control mean.

    ``mean_signal_ratio`` = log2(mean gene expr / mean control expr); ``max_signal_ratio`` =
    log2(max gene expr / mean control expr). Writes both to ``obs``; returns their medians and the
    fraction of cells whose signal does not clear the control floor. Skips (no raise) without
    control probes in the matrix.
    """
    ctrl = _control_mask(adata)
    if ctrl.sum() == 0:
        return {"status": "skipped: no control probes in the matrix"}
    X = _csr(adata.X).astype(float)
    Xg, Xc = X[:, ~ctrl], X[:, ctrl]
    eps = 1e-9
    gmean = np.asarray(Xg.mean(1)).ravel()
    cmean = np.asarray(Xc.mean(1)).ravel()
    gmax = np.asarray(Xg.max(1).todense()).ravel()
    msr = np.log2((gmean + eps) / (cmean + eps))
    maxr = np.log2((gmax + eps) / (cmean + eps))
    adata.obs["mean_signal_ratio"] = msr
    adata.obs["max_signal_ratio"] = maxr
    return {
        "median_mean_signal_ratio": float(np.median(msr)),
        "median_max_signal_ratio": float(np.median(maxr)),
        "frac_cells_signal_below_control": float((msr <= 0).mean()),
    }


# --------------------------------------------------------------------------- #
# Sparsity (SpatialQM getSparsity)
# --------------------------------------------------------------------------- #
def sparsity(adata) -> dict:
    """Zero fraction of the count matrix (dataset) + per-cell sparsity (written to ``obs``)."""
    X = _csr(adata.X)
    n, g = X.shape
    zero_frac = 1.0 - X.nnz / (n * g)
    per_cell = 1.0 - np.asarray((X > 0).sum(1)).ravel() / g
    adata.obs["sparsity"] = per_cell
    return {"zero_fraction": float(zero_frac), "median_cell_sparsity": float(np.median(per_cell))}


# --------------------------------------------------------------------------- #
# Detection entropy (SpatialQM getEntropy)
# --------------------------------------------------------------------------- #
def detection_entropy(adata) -> dict:
    """Per-cell normalized Shannon entropy over the cell's detected genes (0..1).

    Vectorized: ``H = ln(total) - (1/total) * sum(x*ln x)`` over nonzero counts, normalized by
    ``ln(n_detected)``. Low (→0) means one gene dominates the cell (rRNA / probe artifact / extreme
    specialization); high (→1) means counts are spread across many genes. Writes ``obs['detection_entropy']``.
    """
    X = _csr(adata.X).astype(float)
    total = np.asarray(X.sum(1)).ravel()
    M = X.copy()
    M.data = X.data * np.log(X.data)
    sum_xlogx = np.asarray(M.sum(1)).ravel()
    safe_total = np.where(total > 0, total, 1.0)
    ent = np.log(safe_total) - sum_xlogx / safe_total          # Shannon entropy in nats
    n_det = np.asarray((X > 0).sum(1)).ravel()
    norm = np.log(np.maximum(n_det, 1))
    norm = np.where(norm > 0, norm, 1.0)
    ent_norm = np.clip(ent / norm, 0.0, 1.0)
    adata.obs["detection_entropy"] = ent_norm
    return {"median_detection_entropy": float(np.median(ent_norm)),
            "frac_cells_single_gene_dominated": float((ent_norm < 0.1).mean())}


# --------------------------------------------------------------------------- #
# Transcript density per area (SpatialQM getTxPerArea)
# --------------------------------------------------------------------------- #
def tx_per_area(adata) -> dict:
    """Per-cell transcript density = total counts / cell area (needs ``obs['cell_area']``)."""
    if "cell_area" not in adata.obs.columns:
        return {"status": "skipped: no cell_area column"}
    if "total_counts" in adata.obs.columns:
        counts = np.asarray(adata.obs["total_counts"], dtype=float)
    else:
        counts = np.asarray(_csr(adata.X).sum(1)).ravel().astype(float)
    area = np.asarray(adata.obs["cell_area"], dtype=float)
    dens = np.divide(counts, area, out=np.full_like(counts, np.nan), where=area > 0)
    adata.obs["tx_per_area"] = dens
    return {"median_tx_per_area": float(np.nanmedian(dens))}


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def run_signal_qc(adata, do_moran: bool = True) -> dict:
    """Run the annotation-independent signal-QC battery and return one headline dict. Each metric
    is guarded (returns ``{'status': 'skipped: ...'}`` rather than raising) so this is safe on any
    section. ``do_moran=False`` skips the (heavier) Moran's I pass."""
    out: dict = {}
    for name, fn in (("sparsity", sparsity), ("detection_entropy", detection_entropy),
                     ("signal_to_noise", signal_to_noise), ("tx_per_area", tx_per_area)):
        try:
            out[name] = fn(adata)
        except Exception as e:  # pragma: no cover - defensive
            out[name] = {"status": f"skipped: {e}"}
    if do_moran and "spatial" in adata.obsm:
        try:
            out["moran"] = moran_signal(adata)
        except Exception as e:  # pragma: no cover
            out["moran"] = {"status": f"skipped: {e}"}
    return out
