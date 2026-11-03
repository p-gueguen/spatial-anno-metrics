"""Tests for battery.py - the reference-free orchestrator (AQI behaviour is in test_aqi.py)."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")


def test_annotation_quality_orchestrator(blobs):
    from spatial_anno_metrics import battery as em

    a = blobs()
    a.var["control"] = False
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca")
    assert "internal_validity" in out
    assert out["internal_validity"]["neighborhood_purity"] > 0.9
