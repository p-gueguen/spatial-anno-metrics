"""Tests for purity.py - the marker co-expression / contamination purity battery.

Reference-free (crisp_purity, mecr, pmp) + reference-based (nmp, ncp). The fixture contrasts a CLEAN
section (mutually-exclusive lineage markers never co-occur) with a CONTAMINATED one (a block of cells
co-expresses both lineages, mimicking spillover / spatial doublets); every metric must rate clean as
purer than contaminated - the whole point of the module.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

MARKERS = {"T": ["CD3D", "CD3E", "TRAC"], "Mal": ["MLANA", "SOX10", "DCT"]}
GENES = ["CD3D", "CD3E", "TRAC", "MLANA", "SOX10", "DCT"] + [f"BG{i}" for i in range(20)]
_GI = {g: i for i, g in enumerate(GENES)}


def _section(contaminated_frac=0.0, n=200, seed=0):
    """Half T cells, half Mal cells. `contaminated_frac` of cells co-express BOTH lineages' markers."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(0.2, size=(n, len(GENES))).astype("float32")
    # zero the MARKER columns so lineage detection is controlled by the +8 spikes below (background BG
    # genes keep their Poisson noise); otherwise Poisson noise on marker genes cross-detects lineages.
    for g in MARKERS["T"] + MARKERS["Mal"]:
        X[:, _GI[g]] = 0
    half = n // 2
    for g in MARKERS["T"]:
        X[:half, _GI[g]] += 8
    for g in MARKERS["Mal"]:
        X[half:, _GI[g]] += 8
    n_contam = int(contaminated_frac * n)
    if n_contam:                                   # give the first n_contam cells the OTHER lineage too
        for i in range(n_contam):
            other = MARKERS["Mal"] if i < half else MARKERS["T"]
            for g in other:
                X[i, _GI[g]] += 8
    a = anndata.AnnData(X=X)
    a.var_names = GENES
    a.obs_names = [f"c{i}" for i in range(n)]
    a.obs["cell_type"] = pd.Categorical(["T"] * half + ["Mal"] * (n - half))
    return a


# --------------------------------------------------------------------------- #
# reference-free
# --------------------------------------------------------------------------- #
def test_crisp_purity_clean_vs_contaminated():
    from spatial_anno_metrics import purity as p

    clean = _section(0.0)
    dirty = _section(0.4, seed=1)
    pc = p.crisp_purity(clean, MARKERS)
    pd_ = p.crisp_purity(dirty, MARKERS)
    assert pc == pytest.approx(1.0)                       # exclusive markers never co-occur
    assert pd_ < pc                                       # contamination lowers purity
    assert "crisp_impure" in clean.obs and clean.obs["crisp_impure"].sum() == 0
    assert dirty.obs["crisp_impure"].sum() > 0


def test_mecr_clean_zero_contaminated_positive():
    from spatial_anno_metrics import purity as p

    assert p.mecr(_section(0.0), MARKERS) == pytest.approx(0.0)
    assert p.mecr(_section(0.4, seed=1), MARKERS) > 0.0


def test_mecr_no_pairs_returns_zero():
    from spatial_anno_metrics import purity as p

    assert p.mecr(_section(0.0), {"T": ["CD3D"]}) == 0.0  # only one lineage -> no disjoint pairs


def test_pmp_writes_bounded_column_and_ranks_pure_above_contaminated():
    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2)
    p.pmp(a, MARKERS, label_key="cell_type")
    v = a.obs["pmp"].to_numpy()
    assert "pmp" in a.obs
    finite = v[np.isfinite(v)]
    assert finite.size and finite.min() >= 0.0 and finite.max() <= 1.0
    # cells contaminated with the other lineage (first 60) have LOWER purity than the clean tail.
    assert np.nanmean(v[:60]) < np.nanmean(v[120:])


# --------------------------------------------------------------------------- #
# reference-based (guarded)
# --------------------------------------------------------------------------- #
def test_nmp_ncp_skip_without_reference():
    from spatial_anno_metrics import purity as p

    a = _section(0.0)
    assert p.nmp(a, reference=None)["status"] == "skipped"
    assert p.ncp(a, reference=None)["status"] == "skipped"


def test_nmp_flags_negative_marker_contamination():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=5)          # clean reference: T cells lack Mal genes and vice-versa
    dirty = _section(0.5, seed=6)        # query with heavy cross-lineage contamination
    out = p.nmp(dirty, reference=ref, label_key="cell_type", ref_label_key="cell_type")
    assert out["status"] == "ok"
    assert "nmp" in dirty.obs and 0.0 <= out["mean_nmp"] <= 1.0
    assert out["mean_nmp"] > p.nmp(_section(0.0, seed=7), reference=ref)["mean_nmp"]


def test_ncp_returns_fraction_in_unit_interval():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=8)
    out = p.ncp(_section(0.0, seed=9), reference=ref)
    assert out["status"] in ("ok", "skipped")
    if out["status"] == "ok":
        assert 0.0 <= out["ncp"] <= 1.0
