"""Tests for _hierarchy.py - ancestors and lowest-common-ancestor depth on a {child: parent} map."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")

PARENT = {"CD4 T": "T", "CD8 T": "T", "T": "Lymphoid", "B": "Lymphoid", "Epi": "Epi"}   # Epi: its own parent


def test_ancestors_nearest_first():
    from spatial_anno_metrics._hierarchy import ancestors, validate_hierarchy

    assert ancestors("CD4 T", PARENT) == ["T", "Lymphoid"]
    assert ancestors("Lymphoid", PARENT) == []                   # a root
    assert ancestors("Epi", PARENT) == [] and ancestors("unknown", PARENT) == []
    assert ancestors("root_none", {"root_none": None, "child": "root_none"}) == []
    assert ancestors("child", {"root_none": None, "child": "root_none"}) == ["root_none"]
    with pytest.raises(ValueError, match="cycle"):
        ancestors("x", {"x": "y", "y": "x"})
    with pytest.raises(ValueError, match="cycle"):
        validate_hierarchy({"ok": "root", "x": "y", "y": "x"})


def test_lca_depth():
    from spatial_anno_metrics._hierarchy import lca_depth

    assert lca_depth("CD4 T", "CD8 T", PARENT) == 1              # siblings meet at T, one below the root
    assert lca_depth("CD4 T", "B", PARENT) == 0                  # only the root in common
    assert lca_depth("CD4 T", "CD4 T", PARENT) == 2              # a label is its own LCA
    assert lca_depth("CD4 T", "T", PARENT) == lca_depth("T", "CD4 T", PARENT) == 1
    assert lca_depth("CD4 T", "Epi", PARENT) == -1               # different roots
