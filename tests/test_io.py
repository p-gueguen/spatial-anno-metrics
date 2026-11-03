"""Tests for io.py - SpatialData input."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def test_from_spatialdata_returns_the_table_itself():
    spatialdata = pytest.importorskip("spatialdata")
    from spatialdata.models import TableModel

    from spatial_anno_metrics import io

    t = anndata.AnnData(X=np.ones((3, 2), dtype="float32"), obs=pd.DataFrame(index=["a", "b", "c"]))
    t.obsm["spatial"] = np.zeros((3, 2))
    sdata = spatialdata.SpatialData(tables={"cells": TableModel.parse(t)})
    a = io.from_spatialdata(sdata, table_key="cells")
    assert a is sdata.tables["cells"]                  # obs columns the metrics write land in the SpatialData
    with pytest.raises(KeyError, match="'table'.*cells"):
        io.from_spatialdata(sdata)                     # default key absent: name what is there
