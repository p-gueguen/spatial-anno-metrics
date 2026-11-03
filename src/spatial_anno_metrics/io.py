"""Input adapters: get the AnnData the metrics read out of other containers."""
from __future__ import annotations

__all__ = ["from_spatialdata"]


def from_spatialdata(sdata, table_key: str = "table"):
    """Return the AnnData table ``table_key`` of a SpatialData object, for the metrics in this package.

    Returns the table itself, not a copy: obs columns and ``uns['sam']`` summaries the metrics write land in
    ``sdata.tables[table_key]``. Spatial metrics need cell coordinates in ``obsm['spatial']``; the
    spatialdata-io readers (Xenium, MERSCOPE, CosMx, ...) put the cell centroids there. Needs the optional
    ``spatialdata`` package.
    """
    try:
        import spatialdata  # noqa: F401 - only to fail early with a clear message
    except Exception as e:  # noqa: BLE001 - optional dependency, any import failure means "not usable"
        raise ImportError(f"from_spatialdata needs spatialdata, which failed to import: {e}") from e
    if table_key not in sdata.tables:
        raise KeyError(f"no table {table_key!r} in sdata.tables (available: {list(sdata.tables)})")
    return sdata.tables[table_key]
