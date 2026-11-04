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


def test_segtraq_covariates_passthrough_and_summary():
    from spatial_anno_metrics import io

    t = anndata.AnnData(X=np.ones((4, 2), dtype="float32"), obs=pd.DataFrame(index=["a", "b", "c", "d"]))
    assert io.segtraq_covariates(t)["status"] == "skipped"
    t.obs["border_admixture_score"] = [0.1, 0.2, 0.3, np.nan]
    t.obs["similarity_top_bottom"] = [0.9, 0.8, 0.7, 0.6]
    out = io.segtraq_covariates(t)
    assert out["status"] == "ok" and out["columns"] == ["border_admixture_score", "similarity_top_bottom"]
    assert out["summary"]["border_admixture_score"]["mean"] == pytest.approx(0.2)
    assert out["summary"]["border_admixture_score"]["n_defined"] == 3
    assert out["summary"]["similarity_top_bottom"]["median"] == pytest.approx(0.75)

    t.obs["rctd_contamination"] = [np.nan, np.nan, np.nan, np.nan]
    out_nan = io.segtraq_covariates(t)
    assert out_nan["summary"]["rctd_contamination"]["n_defined"] == 0
    assert np.isnan(out_nan["summary"]["rctd_contamination"]["mean"])
    assert np.isnan(out_nan["summary"]["rctd_contamination"]["median"])

    spatialdata = pytest.importorskip("spatialdata")
    from spatialdata.models import TableModel

    sdata = spatialdata.SpatialData(tables={"cells": TableModel.parse(t)})
    out_sd = io.segtraq_covariates(sdata, table_key="cells")
    assert out_sd["status"] == "ok" and "border_admixture_score" in out_sd["columns"]
