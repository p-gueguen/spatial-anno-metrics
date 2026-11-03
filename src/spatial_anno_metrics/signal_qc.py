"""Annotation-independent signal / QC metrics, squidpy-compatible.

Ported from **SpatialQM** (Center for Spatial OMICs, github.com/Center-for-Spatial-OMICs/SpatialQM, MIT, an R
package published with the Spatial Touchstone study, Plummer et al., Nat Biotechnol 2025) - the metrics that
score the *signal quality* of an imaging-ST section (is the measured expression real structured biology or
background noise?) independent of any cell-type labels. They complement the contamination metrics in
``purity.py`` and the label-quality battery in ``battery.py``.

  * ``moran_signal``      - Moran's I per gene vs per negative control (SpatialQM ``getMorans``).
  * ``signal_to_noise``   - gene expression over the negative-control level (``getMeanSignalRatio`` /
                            ``getMaxRatio`` with ``method="spatialqm"``).
  * ``sparsity``          - zero fraction of the count matrix, dataset + per cell (``getSparsity``).
  * ``detection_entropy`` - per-cell normalised entropy, or SpatialQM's ``getEntropy`` (``method="spatialqm"``).
  * ``tx_per_area``       - transcripts per unit cell area, mean (``getTxPerArea``) and median.
  * ``tx_per_cell``       - gene transcripts per cell, mean (``getTxPerCell``) and median.
  * ``transcript_qc``     - fraction of transcripts in cells (``getCellTxFraction``) and the false discovery
                            rate from negative codewords (``getGlobalFDR``), from a transcripts DataFrame.
  * ``run_signal_qc``     - orchestrator -> one headline dict.

Where a function has ``method="spatialqm"``, that method reproduces the SpatialQM function on the parity golden
(``tests/goldens/signalqc_spatialqm``); the other methods are this package's recipes and say how they differ.

**Which features are controls.** Never the feature name: in the 10x Xenium breast panel 8 of the 28 negative-
control probes are named ``antisense_<GENE>``, and SpatialQM's name regexes count them as genes. The rule:

1. ``var['control']`` (bool) when present - an explicit flag always wins;
2. else ``var['feature_type']`` (or ``var['feature_types']``, as ``scanpy.read_10x_h5`` names it): every feature
   whose type is not ``"Gene Expression"`` is a control (Xenium: negative-control probes, negative-control and
   blank / unassigned codewords);
3. else no controls (the control-based metrics skip).

``moran_signal(method="spatialqm")`` uses only the negative-control *probes*, as SpatialQM does:
``feature_type == "Negative Control Probe"`` when a feature-type column exists, else rule 1.

Every count-based metric reads the matrix chosen by ``layer=`` (default ``layers['counts']``, or ``X`` when it
holds raw counts; see ``_validate.get_counts``) and names it under ``"matrix"``. Per-cell values go to
``obs[key_added]`` (default ``sam_<name>``) and the scalar summary to ``uns['sam'][<function name>]``;
``moran_signal``, ``transcript_qc`` and every ``method="spatialqm"`` call only return their result.

Depends on: numpy, scipy, pandas; squidpy (``moran_signal(method="squidpy")`` only).
"""
from __future__ import annotations

import numpy as np

from ._utils import _squidpy, _store
from ._validate import get_counts, skip
from .thresholds import value as _thr

__all__ = ["moran_signal", "signal_to_noise", "sparsity", "detection_entropy", "tx_per_area", "tx_per_cell",
           "transcript_qc", "run_signal_qc"]


_ONE_GENE = _thr("signal_qc", "single_gene_entropy")
_CTRL_Q = _thr("signal_qc", "control_quantile")
_SQM_K = _thr("signal_qc", "spatialqm_moran_n_neighs")
_SQM_PSEUDO = _thr("signal_qc", "spatialqm_pseudocount")
_GENE_TYPE, _PROBE_TYPE = "Gene Expression", "Negative Control Probe"


def _feature_types(adata):
    for col in ("feature_type", "feature_types"):
        if col in adata.var.columns:
            return np.asarray(adata.var[col], dtype=str)
    return None


def _control_mask(adata, probes_only: bool = False) -> np.ndarray:
    """Flag controls by the rule in the module docstring; ``probes_only`` keeps negative-control probes only."""
    ft = _feature_types(adata)
    if probes_only and ft is not None:
        return ft == _PROBE_TYPE
    if "control" in adata.var.columns:
        return np.asarray(adata.var["control"], dtype=bool)
    if ft is not None:
        return ft != _GENE_TYPE
    return np.zeros(adata.n_vars, dtype=bool)


def _csr(X):
    import scipy.sparse as sp
    return X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)


def _check_method(method, allowed):
    if method not in allowed:
        raise ValueError(f"method must be one of {allowed}, got {method!r}")


# --------------------------------------------------------------------------- #
# Moran's I gene-vs-control (SpatialQM getMorans)
# --------------------------------------------------------------------------- #
def _knn_idw(xy, k):
    """Build spdep ``knearneigh(k)`` + ``nb2listwdist(type="idw", style="W")``: weights 1/d, rows summing to 1."""
    import scipy.sparse as sp
    from scipy.spatial import cKDTree

    n = len(xy)
    d, j = cKDTree(xy).query(xy, k + 1)
    if (d[:, 1:] == 0).any():
        raise ValueError("duplicated cell coordinates: inverse-distance weights are undefined")
    w = 1.0 / d[:, 1:]
    w /= w.sum(1, keepdims=True)
    return sp.csr_matrix((w.ravel(), j[:, 1:].ravel(), np.arange(0, n * k + 1, k)), shape=(n, n))


def _moran_columns(V, W) -> np.ndarray:
    """Moran's I of every column of sparse ``V`` on weights ``W``; NaN for a constant column."""
    n = V.shape[0]
    s0 = W.sum()
    m = np.asarray(V.mean(0)).ravel()
    vWv = np.asarray(V.multiply(W @ V).sum(0)).ravel()
    cV = np.asarray(W.sum(0)).ravel() @ V                       # 1'W v
    rV = np.asarray(W.sum(1)).ravel() @ V                       # v'W 1
    zWz = vWv - m * cV - m * rV + m * m * s0
    vv = np.asarray(V.multiply(V).sum(0)).ravel()
    zz = vv - n * m * m
    with np.errstate(invalid="ignore", divide="ignore"):
        out = (n / s0) * zWz / zz
    out[zz <= 1e-12 * vv] = np.nan                              # constant column (all-zero included)
    return out


def _moran_spatialqm(C, xy, score):
    """Score one SpatialQM group; return (I of the ``score`` columns, n cells kept).

    Keep cells with a count anywhere in the group, ``scater::logNormCounts`` over the whole group (library-size
    factors centred at 1, log2(x / sf + 1)), Moran's I on the 20-NN 1/d graph. With <= 20 kept cells spdep's
    ``knearneigh(k = 20)`` stops (SpatialQM gives no value): every column is NaN.
    """
    C = C.astype(np.float64)                                    # float32 library sizes cost 1e-9 in I
    lib = np.asarray(C.sum(1)).ravel()
    keep = lib > 0
    n = int(keep.sum())
    if n <= _SQM_K:
        return np.full(int(score.sum()), np.nan), n
    C = C[keep].tocsr()
    sf = lib[keep] / lib[keep].mean()
    C.data = np.log2(C.data / np.repeat(sf, np.diff(C.indptr)) + 1.0)
    return _moran_columns(C[:, score], _knn_idw(xy[keep], _SQM_K)), n


def moran_signal(adata, n_neighs: int = _thr("signal_qc", "moran_n_neighs"), genes: list[str] | None = None,
                 n_perms: int | None = None, seed: int = 0, *, method: str = "squidpy",
                 layer: str | None = "counts") -> dict:
    """Moran's I spatial autocorrelation split into panel genes vs negative controls.

    Real genes carry spatial structure (I well above 0); controls do not (I near 0), so the **separation**
    between the two medians is a reference-free spatial signal-quality score. Summaries ignore features whose
    I is undefined (no count, or constant, in the cells used) and report how many were left out.

    Parameters
    ----------
    adata : AnnData
        Section with coordinates in ``obsm['spatial']`` and controls flagged as in the module docstring.
    n_neighs : int
        Spatial kNN size of the graph ``method="squidpy"`` builds (binary weights). Ignored by
        ``method="spatialqm"``, which uses SpatialQM's 20 neighbours weighted 1/d.
    genes : list of str, optional
        Features to score (default: all). With ``method="spatialqm"`` the cells and library sizes still come
        from the whole group, so a feature's I does not depend on which others are requested.
    n_perms : int, optional
        Permutations for squidpy's p-values (``method="squidpy"`` only).
    seed : int
        Seed for those permutations.
    method : {"squidpy", "spatialqm"}
        ``"squidpy"`` (default, the 0.6 recipe): raw counts of every cell, ``sq.gr.spatial_neighbors(n_neighs)``
        reused from ``obsp['spatial_connectivities']`` when present, ``sq.gr.spatial_autocorr``; controls = every
        control feature. Writes squidpy's graph cache (``obsp['spatial_connectivities']``,
        ``obsp['spatial_distances']``, ``uns['spatial_neighbors']``) and ``uns['moranI']``.
        ``"spatialqm"``: SpatialQM ``getMorans`` - per group (genes; negative-control probes), the cells with a
        count in that group, ``scater::logNormCounts``, 20 nearest neighbours weighted 1/d and row-standardised.
        Needs no squidpy and writes nothing to ``adata``. Its control median sits near the null
        E[I] = -1/(n - 1) of the few cells with a probe count. A group with <= 20 such cells (spdep's
        ``knearneigh(k = 20)`` stops there) gets NaN for all its features, counted as excluded.
    layer : str or None
        Count matrix (see ``_validate.get_counts``).

    Returns
    -------
    dict
        ``status``, ``median_gene_moran``, ``median_control_moran`` (None without controls), ``n_genes``,
        ``n_controls``, ``n_genes_excluded`` / ``n_controls_excluded`` (undefined I), ``separation`` and
        ``frac_genes_above_control_p95`` (with controls), ``per_feature`` (pandas Series of I, NaN = undefined),
        ``method``, ``matrix``; ``"spatialqm"`` adds ``n_cells_genes`` and ``n_cells_controls``. A skip when
        ``method="spatialqm"`` meets duplicated coordinates.
    """
    import pandas as pd

    _check_method(method, ("squidpy", "spatialqm"))
    mat, matrix = get_counts(adata, layer)
    names = np.asarray(adata.var_names, dtype=str)
    use = np.isin(names, list(genes)) if genes is not None else np.ones(len(names), dtype=bool)
    gene = ~_control_mask(adata)
    extra: dict = {}
    if method == "squidpy":
        sq = _squidpy()
        if "spatial_connectivities" not in adata.obsp:
            sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=n_neighs)
        sq.gr.spatial_autocorr(adata, mode="moran", genes=list(names[use]), n_perms=n_perms, seed=seed,
                               layer=None if matrix == "X" else layer)
        per = adata.uns["moranI"]["I"].reindex(names).astype(float)
        ctrl = ~gene
    else:
        ctrl = _control_mask(adata, probes_only=True)
        X, xy = _csr(mat), np.asarray(adata.obsm["spatial"], dtype=float)[:, :2]
        per = pd.Series(np.nan, index=names)
        try:
            for key, grp in (("n_cells_genes", gene), ("n_cells_controls", ctrl)):
                if (grp & use).any():
                    per[grp & use], extra[key] = _moran_spatialqm(X[:, grp], xy, use[grp])
        except ValueError as e:
            return skip(str(e))
    per = per[use]
    gI, cI = per[gene[use]].to_numpy(), per[ctrl[use]].to_numpy()
    gok, cok = gI[~np.isnan(gI)], cI[~np.isnan(cI)]
    out: dict = {
        "status": "ok", "median_gene_moran": float(np.median(gok)) if gok.size else float("nan"),
        "median_control_moran": float(np.median(cok)) if cok.size else None,
        "n_genes": int(gI.size), "n_controls": int(cI.size),
        "n_genes_excluded": int(gI.size - gok.size), "n_controls_excluded": int(cI.size - cok.size),
    }
    if cok.size:
        out["separation"] = out["median_gene_moran"] - out["median_control_moran"]
        out["frac_genes_above_control_p95"] = (float((gok > np.quantile(cok, _CTRL_Q)).mean()) if gok.size
                                               else float("nan"))
    out.update(extra)
    out.update({"per_feature": per, "method": method, "matrix": matrix})
    return out


# --------------------------------------------------------------------------- #
# Signal-to-noise vs controls (SpatialQM getMeanSignalRatio / getMaxRatio)
# --------------------------------------------------------------------------- #
def signal_to_noise(adata, *, method: str = "cell", eps: float | None = None, layer: str | None = "counts",
                    key_added: tuple[str, str] = ("sam_mean_signal_ratio", "sam_max_signal_ratio")) -> dict:
    """Gene expression over the negative-control level.

    Parameters
    ----------
    adata : AnnData
        Section with controls flagged as in the module docstring.
    method : {"cell", "spatialqm"}
        ``"cell"`` (default): per cell, ``log2(mean gene count / mean control count)`` and
        ``log2(max gene count / mean control count)``, written to the two ``key_added`` columns. The ratio is
        defined only for cells with at least one control count (NaN otherwise; -inf when such a cell has no
        gene count). On imaging panels most cells have no control count (88.8 % on the Xenium breast golden),
        so the result also reports the pooled ratio over all cells and how many cells the medians rest on.
        ``"spatialqm"``: SpatialQM's dataset scalars on per-feature means over cells,
        ``mean_signal_ratio_log10 = mean_g log10(mean_g + 0.1) - mean_c log10(mean_c + 0.1)`` and
        ``max_signal_ratio_log10 = log10(max_g mean_g) - log10(mean_c mean_c)``. ``getMaxRatio`` is reproduced
        exactly; ``getMeanSignalRatio`` pairs gene i with control i mod n_controls (R recycling), which this
        port replaces by the difference of the two means (0.49702 vs ours within 5e-4 on the golden). Writes
        nothing to ``adata``.
    eps : float, optional
        ``method="cell"`` only. None (default) leaves undefined cells NaN. A number restores the 0.6 recipe,
        ``log2((gene + eps) / (control + eps))`` for every cell (0.6 used 1e-9, whose floor log2(1/1e-9) = 29.9
        set the median whenever most cells had no control count, docs/ROADMAP.md 8g).
    layer : str or None
        Count matrix (see ``_validate.get_counts``).
    key_added : tuple of str
        Obs columns for the per-cell mean and max ratios (``method="cell"``).

    Returns
    -------
    dict
        ``"cell"``: ``status``, ``median_mean_signal_ratio`` and ``median_max_signal_ratio`` (log2, over the
        cells with a control count; over all cells when ``eps`` is set), ``pooled_mean_signal_ratio`` (log2 of
        mean gene count over mean control count, all cells), ``n_cells_with_controls``,
        ``frac_cells_with_controls``, ``frac_cells_signal_below_control`` (mean gene <= mean control, among cells
        with a control count), ``matrix``. A skip (with ``n_cells_with_controls``: 0) when no cell has a control
        count and ``eps`` is None. ``"spatialqm"``: ``status``, ``mean_signal_ratio_log10``,
        ``max_signal_ratio_log10``, ``n_genes``, ``n_controls``, ``method``, ``matrix``; a skip when the
        controls hold no count. Both skip without control features.
    """
    _check_method(method, ("cell", "spatialqm"))
    if isinstance(key_added, str) or len(key_added) != 2:
        raise TypeError("key_added takes two column names: (mean-ratio column, max-ratio column)")
    ctrl = _control_mask(adata)
    if ctrl.sum() == 0:
        return skip("no control features (var['control'] or var['feature_type'])")
    mat, matrix = get_counts(adata, layer)
    X = _csr(mat).astype(float)
    Xg, Xc = X[:, ~ctrl], X[:, ctrl]
    if method == "spatialqm":
        gm, cm = np.asarray(Xg.mean(0)).ravel(), np.asarray(Xc.mean(0)).ravel()
        if cm.sum() == 0:
            return skip("the control features hold no count")
        return {"status": "ok",
                "mean_signal_ratio_log10": float(np.log10(gm + _SQM_PSEUDO).mean() - np.log10(cm + _SQM_PSEUDO).mean()),
                "max_signal_ratio_log10": float(np.log10(gm.max()) - np.log10(cm.mean())),
                "n_genes": int(gm.size), "n_controls": int(cm.size), "method": method, "matrix": matrix}
    gmean = np.asarray(Xg.mean(1)).ravel()
    cmean = np.asarray(Xc.mean(1)).ravel()
    gmax = np.asarray(Xg.max(1).todense()).ravel()
    has = cmean > 0
    n_has = int(has.sum())
    if eps is None and n_has == 0:
        return {**skip("no cell has a control count: the ratio is undefined"), "n_cells_with_controls": 0}
    with np.errstate(divide="ignore", invalid="ignore"):
        if eps is None:
            msr = np.where(has, np.log2(gmean / cmean), np.nan)
            maxr = np.where(has, np.log2(gmax / cmean), np.nan)
            pooled = float(np.log2(gmean.mean() / cmean.mean()))
        else:
            msr = np.log2((gmean + eps) / (cmean + eps))
            maxr = np.log2((gmax + eps) / (cmean + eps))
            pooled = float(np.log2(gmean.mean() / cmean.mean())) if n_has else float("nan")
    adata.obs[key_added[0]] = msr
    adata.obs[key_added[1]] = maxr
    sel = has if eps is None else np.ones_like(has)
    out = {
        "status": "ok", "median_mean_signal_ratio": float(np.median(msr[sel])),
        "median_max_signal_ratio": float(np.median(maxr[sel])), "pooled_mean_signal_ratio": pooled,
        "n_cells_with_controls": n_has, "frac_cells_with_controls": n_has / adata.n_obs,
        "frac_cells_signal_below_control": (float((msr <= 0).mean()) if eps is not None
                                            else float((gmean[has] <= cmean[has]).mean())),
        "matrix": matrix,
    }
    _store(adata, "signal_to_noise", out, key_added)
    return out


# --------------------------------------------------------------------------- #
# Sparsity (SpatialQM getSparsity)
# --------------------------------------------------------------------------- #
def sparsity(adata, *, layer: str | None = "counts", key_added: str = "sam_sparsity") -> dict:
    """Zero fraction of the count matrix (SpatialQM ``getSparsity``) and per cell.

    Parameters
    ----------
    adata : AnnData
        Any section; every feature counts, controls included.
    layer : str or None
        Count matrix (see ``_validate.get_counts``).
    key_added : str
        Obs column for the per-cell zero fraction.

    Returns
    -------
    dict
        ``status``, ``zero_fraction`` (dataset), ``median_cell_sparsity``, ``matrix``.
    """
    mat, matrix = get_counts(adata, layer)
    X = _csr(mat)
    n, g = X.shape
    zero_frac = 1.0 - X.nnz / (n * g)
    per_cell = 1.0 - np.asarray((X > 0).sum(1)).ravel() / g
    adata.obs[key_added] = per_cell
    out = {"status": "ok", "zero_fraction": float(zero_frac), "median_cell_sparsity": float(np.median(per_cell)),
           "matrix": matrix}
    _store(adata, "sparsity", out, key_added)
    return out


# --------------------------------------------------------------------------- #
# Detection entropy (SpatialQM getEntropy)
# --------------------------------------------------------------------------- #
def detection_entropy(adata, *, method: str = "cell", layer: str | None = "counts",
                      key_added: str = "sam_detection_entropy") -> dict:
    """Shannon entropy of the counts.

    Parameters
    ----------
    adata : AnnData
        Any section.
    method : {"cell", "spatialqm"}
        ``"cell"`` (default): per cell, the entropy of its counts over its detected features divided by
        ln(n detected) (0..1), written to ``obs[key_added]``; low means one feature dominates the cell (rRNA,
        probe artefact). ``"spatialqm"``: SpatialQM ``getEntropy`` (``BioQC::entropy`` of the whole matrix) -
        the entropy in bits of the histogram of count values over every entry, zeros included; one dataset
        scalar, nothing written to ``adata``.
    layer : str or None
        Count matrix (see ``_validate.get_counts``).
    key_added : str
        Obs column for the per-cell value (``method="cell"``).

    Returns
    -------
    dict
        ``"cell"``: ``status``, ``median_detection_entropy``, ``frac_cells_single_gene_dominated``, ``matrix``.
        ``"spatialqm"``: ``status``, ``entropy_bits``, ``n_distinct_values``, ``method``, ``matrix``.
    """
    _check_method(method, ("cell", "spatialqm"))
    mat, matrix = get_counts(adata, layer)
    X = _csr(mat).astype(float)
    if method == "spatialqm":
        vals, cnt = np.unique(X.data, return_counts=True)
        n_zero = X.shape[0] * X.shape[1] - X.nnz
        cnt = np.concatenate([cnt[vals != 0], [n_zero + cnt[vals == 0].sum()]])
        p = cnt[cnt > 0] / cnt.sum()
        return {"status": "ok", "entropy_bits": float(-(p * np.log2(p)).sum()), "n_distinct_values": int(p.size),
                "method": method, "matrix": matrix}
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
    adata.obs[key_added] = ent_norm
    out = {"status": "ok", "median_detection_entropy": float(np.median(ent_norm)),
           "frac_cells_single_gene_dominated": float((ent_norm < _ONE_GENE).mean()), "matrix": matrix}
    _store(adata, "detection_entropy", out, key_added)
    return out


# --------------------------------------------------------------------------- #
# Transcripts per area / per cell (SpatialQM getTxPerArea / getTxPerCell)
# --------------------------------------------------------------------------- #
def tx_per_area(adata, *, layer: str | None = "counts", key_added: str = "sam_tx_per_area") -> dict:
    """Per-cell transcript density, total counts / cell area (SpatialQM ``getTxPerArea`` is the mean).

    Parameters
    ----------
    adata : AnnData
        Section with ``obs['cell_area']``. Total counts come from ``obs['total_counts']`` when present (Xenium:
        every feature, controls included), else from the matrix.
    layer : str or None
        Count matrix, read only without ``obs['total_counts']``.
    key_added : str
        Obs column for the per-cell density (NaN where the area is not positive).

    Returns
    -------
    dict
        ``status``, ``mean_tx_per_area`` (= ``getTxPerArea``), ``median_tx_per_area``, ``matrix``; a skip without
        ``obs['cell_area']``.
    """
    if "cell_area" not in adata.obs.columns:
        return skip("no obs['cell_area']")
    if "total_counts" in adata.obs.columns:
        counts, matrix = np.asarray(adata.obs["total_counts"], dtype=float), "obs['total_counts']"
    else:
        mat, matrix = get_counts(adata, layer)
        counts = np.asarray(_csr(mat).sum(1)).ravel().astype(float)
    area = np.asarray(adata.obs["cell_area"], dtype=float)
    dens = np.divide(counts, area, out=np.full_like(counts, np.nan), where=area > 0)
    adata.obs[key_added] = dens
    out = {"status": "ok", "mean_tx_per_area": float(np.nanmean(dens)), "median_tx_per_area": float(np.nanmedian(dens)),
           "matrix": matrix}
    _store(adata, "tx_per_area", out, key_added)
    return out


def tx_per_cell(adata, *, layer: str | None = "counts", key_added: str = "sam_tx_per_cell") -> dict:
    """Gene transcripts per cell (SpatialQM ``getTxPerCell`` is the mean).

    Parameters
    ----------
    adata : AnnData
        Any section; controls (module docstring rule) are left out of the count.
    layer : str or None
        Count matrix (see ``_validate.get_counts``).
    key_added : str
        Obs column for the per-cell gene count.

    Returns
    -------
    dict
        ``status``, ``mean_tx_per_cell``, ``median_tx_per_cell``, ``n_genes``, ``matrix``.
    """
    mat, matrix = get_counts(adata, layer)
    gene = ~_control_mask(adata)
    per_cell = np.asarray(_csr(mat)[:, gene].sum(1)).ravel().astype(float)
    adata.obs[key_added] = per_cell
    out = {"status": "ok", "mean_tx_per_cell": float(per_cell.mean()), "median_tx_per_cell": float(np.median(per_cell)),
           "n_genes": int(gene.sum()), "matrix": matrix}
    _store(adata, "tx_per_cell", out, key_added)
    return out


# --------------------------------------------------------------------------- #
# Transcript-level QC (SpatialQM getCellTxFraction / getGlobalFDR)
# --------------------------------------------------------------------------- #
def _decoded(col):
    """Decode a bytes column (Xenium stores ``feature_name`` as parquet binary) to str."""
    if col.dtype == object and len(col) and isinstance(col.iloc[0], bytes):
        return col.str.decode("utf-8")
    return col


def transcript_qc(transcripts, negative_features, *, feature_key: str = "feature_name",
                  cell_key: str = "cell_id") -> dict:
    """Fraction of transcripts in cells and the negative-codeword false discovery rate, per transcript.

    Every row counts (filter by ``qv`` beforehand to match a q-filtered count matrix; SpatialQM does not).

    Parameters
    ----------
    transcripts : pandas.DataFrame
        One row per transcript (e.g. Xenium ``transcripts.parquet``), with ``feature_key`` and ``cell_key``;
        bytes values (parquet ``binary``, as Xenium writes ``feature_name``) are decoded as UTF-8.
    negative_features : iterable of str
        Names of the negative features. Pick them by feature type, not by name, e.g.
        ``adata.var_names[adata.var["feature_type"] != "Gene Expression"]`` (negative-control probes, negative-
        control and blank codewords); SpatialQM's ``^(Neg|Blank|BLANK)`` misses the ``antisense_*`` probes.
    feature_key : str
        Column with the feature name.
    cell_key : str
        Column with the cell id; ``-1`` (Xenium <= 1.2) and ``"UNASSIGNED"`` (>= 1.3) mean no cell.

    Returns
    -------
    dict
        ``status``, ``frac_tx_in_cells`` (= ``getCellTxFraction``), ``fdr`` = negative transcripts / all
        transcripts x (genes detected / negatives detected), the share of gene calls expected to be false
        (SpatialQM ``getGlobalFDR`` returns ``fdr / 100``); None when no negative feature is detected.
        Also ``n_transcripts``, ``n_unassigned``, ``n_negative_tx``, ``n_genes_detected``,
        ``n_negatives_detected``. A skip for an empty table.
    """
    import pandas as pd

    n = len(transcripts)
    if n == 0:
        return skip("no transcripts")
    cid, feat = _decoded(transcripts[cell_key]), _decoded(transcripts[feature_key])
    if pd.api.types.is_numeric_dtype(cid):
        unassigned = int((cid == -1).sum())
    else:
        unassigned = int(cid.astype(str).isin(["UNASSIGNED", "-1"]).sum())
    neg = feat.isin(set(map(str, negative_features))).to_numpy()
    n_neg = int(neg.sum())
    n_neg_det, n_gene_det = int(feat[neg].nunique()), int(feat[~neg].nunique())
    return {"status": "ok", "frac_tx_in_cells": (n - unassigned) / n,
            "fdr": n_neg / n * (n_gene_det / n_neg_det) if n_neg_det else None,
            "n_transcripts": n, "n_unassigned": unassigned, "n_negative_tx": n_neg,
            "n_genes_detected": n_gene_det, "n_negatives_detected": n_neg_det}


# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
def run_signal_qc(adata, do_moran: bool = True, *, layer: str | None = "counts") -> dict:
    """Run the annotation-independent signal-QC battery and return one headline dict.

    Each metric is guarded (a ``{'status': 'skipped', 'reason': ...}`` entry rather than an exception), so this
    is safe on any section. Every metric runs with its default method and writes its default ``sam_*`` columns.

    Parameters
    ----------
    adata : AnnData
        Section to score.
    do_moran : bool
        False skips the (heavier) Moran's I pass.
    layer : str or None
        Count matrix, passed to every metric.

    Returns
    -------
    dict
        ``status`` plus one entry per metric: ``sparsity``, ``detection_entropy``, ``signal_to_noise``,
        ``tx_per_area``, ``tx_per_cell``, ``moran``.
    """
    out: dict = {"status": "ok"}
    for name, fn in (("sparsity", sparsity), ("detection_entropy", detection_entropy),
                     ("signal_to_noise", signal_to_noise), ("tx_per_area", tx_per_area),
                     ("tx_per_cell", tx_per_cell)):
        try:
            out[name] = fn(adata, layer=layer)
        except Exception as e:  # pragma: no cover - defensive
            out[name] = skip(f"{type(e).__name__}: {e}")
    if not do_moran:
        out["moran"] = skip("do_moran=False")
    elif "spatial" not in adata.obsm:
        out["moran"] = skip("no obsm['spatial']")
    else:
        try:
            out["moran"] = moran_signal(adata, layer=layer)
        except Exception as e:  # pragma: no cover
            out["moran"] = skip(f"{type(e).__name__}: {e}")
    return out
