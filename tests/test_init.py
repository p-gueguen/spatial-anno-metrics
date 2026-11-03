"""Tests for the package __init__ - version and public names."""
from __future__ import annotations

import pytest

anndata = pytest.importorskip("anndata")


def test_version_is_a_string():
    import spatial_anno_metrics as sam

    assert isinstance(sam.__version__, str) and sam.__version__


def test_public_names_resolve():
    import spatial_anno_metrics as sam

    assert len(sam.__all__) == len(set(sam.__all__))
    for name in sam.__all__:
        assert callable(getattr(sam, name)) or isinstance(getattr(sam, name), dict), name
    assert "from_spatialdata" in sam.__all__
