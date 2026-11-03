"""Tests for battery.py - the reference-free orchestrator (AQI behaviour is in test_aqi.py)."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")


def test_annotation_quality_orchestrator(blobs):
    from spatial_anno_metrics import battery as em

    a = blobs()
    a.var["control"] = False
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca", n_neighbors=15)
    assert out["status"] == "ok" and out["internal_validity"]["status"] == "ok"
    assert out["internal_validity"]["neighborhood_purity"] > 0.9


def test_annotation_quality_reports_a_failed_sub_metric_as_a_skip(blobs):
    from spatial_anno_metrics import battery as em

    a = blobs(n_per=20)                                # X is a float embedding, not counts
    out = em.annotation_quality(a, label_key="cell_type", embedding="X_pca", marker_sets={"T0": ["0"]},
                                compute_aqi=False)
    fid = out["marker_fidelity"]                       # count-based: refuses X instead of crashing the battery
    assert fid["status"] == "skipped" and "layer=" in fid["reason"]
    a.obs["cell_type"] = "T0"
    iv = em.annotation_quality(a, label_key="cell_type", embedding="X_pca", compute_aqi=False)["internal_validity"]
    assert iv == {"status": "skipped", "reason": iv["reason"]} and "found 1" in iv["reason"]
