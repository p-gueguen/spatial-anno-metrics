"""Input adapters: get the AnnData the metrics read out of other containers."""
from __future__ import annotations

import numpy as np

from ._validate import skip

__all__ = ["from_spatialdata", "segtraq_covariates"]

SEGTRAQ_OBS_COLUMNS: tuple[str, ...] = (
    "border_admixture_score",
    "border_admixture_p_value",
    "similarity_nucleus_cell",
    "similarity_nucleus_cytoplasm",
    "similarity_top_bottom",
    "similarity_top_bottom_p_value",
    "heterotypic_overlap_area",
    "heterotypic_overlap_fraction",
    "vertical_signal_integrity",
    "rctd_contamination",
)


def from_spatialdata(sdata, table_key: str = "table"):
    """Return the AnnData table ``table_key`` of a SpatialData object, for the metrics in this package.

    Returns the table itself, not a copy: obs columns and ``uns['sam']`` summaries the metrics write land in
    ``sdata.tables[table_key]``, and any pre-computed transcript-level admixture covariates written by SegTraQ
    (``border_admixture_score``, ``similarity_nucleus_cytoplasm``, ``similarity_top_bottom``,
    ``heterotypic_overlap_fraction``, ``vertical_signal_integrity``, ``rctd_contamination``; see
    :func:`segtraq_covariates`) are preserved in ``obs``. Spatial metrics need cell coordinates in
    ``obsm['spatial']``; the spatialdata-io readers (Xenium, MERSCOPE, CosMx, ...) put the cell centroids
    there. Needs the optional ``spatialdata`` package.
    """
    import contextlib
    import warnings

    with contextlib.suppress(Exception), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import dask.dataframe  # noqa: F401 - prevent spatialdata<0.4 + dask>=2025.1 query-planning error
    try:
        import spatialdata  # noqa: F401 - only to fail early with a clear message
    except Exception as e:  # noqa: BLE001 - optional dependency, any import failure means "not usable"
        raise ImportError(f"from_spatialdata needs spatialdata, which failed to import: {e}") from e
    if table_key not in sdata.tables:
        raise KeyError(f"no table {table_key!r} in sdata.tables (available: {list(sdata.tables)})")
    return sdata.tables[table_key]


def segtraq_covariates(adata_or_sdata, table_key: str = "table") -> dict:
    """Summarize pre-computed SegTraQ transcript-level admixture columns in ``adata.obs``.

    Transcript-coordinate admixture metrics (border-vs-interior molecular neighborhood similarity,
    nucleus-vs-cytoplasm similarity, z-plane top-vs-bottom similarity, polygon heterotypic overlap) require
    subcellular molecule coordinates and cell-boundary polygons in a ``SpatialData`` object and are computed
    by ``segtraq`` (Lazarenko et al. 2025) into ``sdata.tables[table_key].obs``. This helper detects those
    columns on an ``AnnData`` or ``SpatialData`` container and returns their per-column summary statistics.

    Parameters
    ----------
    adata_or_sdata : AnnData or SpatialData
        Annotated cells or a ``SpatialData`` container.
    table_key : str
        Table name when ``adata_or_sdata`` is a ``SpatialData`` object.

    Returns
    -------
    dict
        ``{status, columns, summary: {col: {mean, median, n_defined}}}``; returns a skip when no SegTraQ
        transcript-admixture column is present in ``obs``.
    """
    adata = (from_spatialdata(adata_or_sdata, table_key=table_key)
             if hasattr(adata_or_sdata, "tables") else adata_or_sdata)
    cols = [c for c in SEGTRAQ_OBS_COLUMNS if c in adata.obs.columns]
    if not cols:
        return skip("no SegTraQ transcript-admixture columns found in adata.obs")
    summary = {}
    for c in cols:
        vals = np.asarray(adata.obs[c], dtype=float)
        fin = vals[np.isfinite(vals)]
        summary[c] = {
            "mean": float(np.mean(fin)) if fin.size else float("nan"),
            "median": float(np.median(fin)) if fin.size else float("nan"),
            "n_defined": int(fin.size),
        }
    return {"status": "ok", "columns": cols, "summary": summary}
