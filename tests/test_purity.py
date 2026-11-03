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


def _section(contaminated_frac=0.0, n=200, seed=0, off_target=0):
    """Half T cells, half Mal cells. `contaminated_frac` of cells co-express BOTH lineages' markers;
    the first `off_target` T cells carry ONE count of the Mal marker MLANA (ambient-level off-target)."""
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
    X[:off_target, _GI["MLANA"]] += 1
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
    rc = p.crisp_purity(clean, MARKERS)
    pc, pd_ = rc["crisp_purity"], p.crisp_purity(dirty, MARKERS)["crisp_purity"]
    assert rc["status"] == "ok" and rc["matrix"] == "X"   # no counts layer: X holds the counts
    assert pc == pytest.approx(1.0)                       # exclusive markers never co-occur
    assert pd_ < pc                                       # contamination lowers purity
    assert "crisp_impure" in clean.obs and clean.obs["crisp_impure"].sum() == 0
    assert dirty.obs["crisp_impure"].sum() > 0


def test_mecr_clean_zero_contaminated_positive():
    from spatial_anno_metrics import purity as p

    assert p.mecr(_section(0.0), MARKERS)["mecr"] == pytest.approx(0.0)
    assert p.mecr(_section(0.4, seed=1), MARKERS)["mecr"] > 0.0


def test_mecr_and_crisp_skip_instead_of_a_clean_score():
    """No lineage pair / no marker-positive cell: nothing was measured, so no 0.0 / 1.0 'perfectly clean'."""
    from spatial_anno_metrics import purity as p

    out = p.mecr(_section(0.0), {"T": ["CD3D"]})          # only one lineage -> no disjoint pairs
    assert out["status"] == "skipped" and "pair" in out["reason"]
    off_panel = {"T": ["nope1"], "Mal": ["nope2"]}
    assert p.crisp_purity(_section(0.0), off_panel)["status"] == "skipped"
    assert p.mecr(_section(0.0), off_panel)["status"] == "skipped"


def test_one_off_target_count_counts_as_detection():
    """The '>= 1 count' rule: a single ambient MLANA count makes a T cell co-detect both lineages."""
    from spatial_anno_metrics import purity as p

    a = _section(0.0, off_target=10)
    out = p.crisp_purity(a, MARKERS)
    assert out["n_impure"] == 10 and out["crisp_purity"] == pytest.approx(1 - 10 / 200)
    assert p.mecr(a, MARKERS)["mecr"] == pytest.approx(10 / 200)   # the one pair: 10 both / 200 either


SHARED = {"CD4 T": ["CD3E", "CD4"], "CD8 T": ["CD3E", "CD8A"], "B": ["MS4A1"]}


def _t_subsets(n_cd8_contam=0):
    """100 CD4 T (CD3E+CD4) + 100 CD8 T (CD3E+CD8A); the first `n_cd8_contam` CD4 T cells also detect CD8A."""
    X = np.zeros((200, 4), dtype="float32")             # CD3E, CD4, CD8A, MS4A1
    X[:, 0] = 3
    X[:100, 1] = 3
    X[100:, 2] = 3
    X[:n_cd8_contam, 2] = 3
    a = anndata.AnnData(X=X)
    a.var_names = ["CD3E", "CD4", "CD8A", "MS4A1"]
    a.obs_names = [f"c{i}" for i in range(200)]
    return a


def test_shared_marker_is_not_co_detection():
    from spatial_anno_metrics import purity as p

    a = _t_subsets()                                     # zero contamination; CD3E is in both T sets
    assert p.mecr(a, SHARED)["mecr"] == pytest.approx(0.0)
    assert p.crisp_purity(a, SHARED)["crisp_purity"] == pytest.approx(1.0)


def test_private_marker_co_detection_still_counted():
    from spatial_anno_metrics import purity as p

    a = _t_subsets(n_cd8_contam=20)
    # CD4 T vs CD8 T on CD4 vs CD8A: 20 both / 200 either = 0.1; both T-vs-B pairs 0 -> mean 0.1 / 3
    assert p.mecr(a, SHARED)["mecr"] == pytest.approx(0.1 / 3)
    assert p.crisp_purity(a, SHARED)["crisp_purity"] == pytest.approx(0.9)   # 20 impure of 200 marker-positive


def test_mecr_skips_nested_marker_sets():
    from spatial_anno_metrics import purity as p

    nested = {"T": ["CD3E"], "CD4 T": ["CD3E", "CD4"], "CD8 T": ["CD3E", "CD8A"]}
    # T has no marker of its own against either subset -> both pairs skipped, not averaged in as 0
    assert p.mecr(_t_subsets(n_cd8_contam=20), nested)["mecr"] == pytest.approx(0.1)


def test_pmp_writes_bounded_column_and_ranks_pure_above_contaminated():
    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2)
    out = p.pmp(a, MARKERS, label_key="cell_type")
    v = a.obs["pmp"].to_numpy()
    assert out["status"] == "ok" and out["mean_pmp"] == pytest.approx(np.nanmean(v))
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


def test_nmp_sees_spillover_on_a_large_panel():
    """A type's negatives must be OTHER lineages' markers. Taken as the bottom decile of its reference profile
    they were the panel's near-silent genes, and 10% fibroblast spillover scored like a clean section."""
    from scipy import sparse

    from spatial_anno_metrics import purity as p

    rng = np.random.default_rng(0)
    n_g = 2000                                           # 100 T, 100 Fib, 20 shared, 1780 near-silent genes
    genes = [f"g{i}" for i in range(n_g)]
    lam = np.full((2, n_g), 0.001)
    lam[0, :100], lam[0, 100:200] = 2.0, 0.02            # each type also carries the other's ambient RNA
    lam[1, 100:200], lam[1, :100] = 2.0, 0.02
    lam[:, 200:220] = 2.0                                # housekeeping: expressed by Fib, but T has it too
    ref = anndata.AnnData(X=sparse.csr_matrix(
        np.vstack([rng.poisson(lam[0], (200, n_g)), rng.poisson(lam[1], (200, n_g))]).astype("float32")))
    ref.var_names = genes
    ref.obs["cell_type"] = ["T"] * 200 + ["Fib"] * 200

    def t_cells(spill):
        lam_t = np.full(n_g, 0.001)
        lam_t[:100], lam_t[100:200], lam_t[200:220] = 1.0, spill, 1.0
        a = anndata.AnnData(X=rng.poisson(lam_t, (100, n_g)).astype("float32"))
        a.var_names = genes
        a.obs["cell_type"] = "T"
        return a

    assert p.nmp(t_cells(0.0), reference=ref)["mean_nmp"] == pytest.approx(0.0)
    # truth: 11 Fib counts / (100 T + 20 shared + 11 Fib + 1.8 silent) = 0.083
    assert p.nmp(t_cells(0.11), reference=ref)["mean_nmp"] == pytest.approx(0.083, abs=0.01)


def test_nmp_reports_types_it_cannot_score():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=5)
    # Mal is absent from the reference; T has no other reference type to take negatives from
    ref = ref[ref.obs["cell_type"] == "T"].copy()
    out = p.nmp(_section(0.0, seed=6), reference=ref)
    assert out["types_without_negatives"] == ["Mal", "T"]


def test_purity_refuses_log_X_without_a_counts_layer():
    """No silent fallback: a missing counts layer is replaced by X only when X holds counts."""
    import warnings

    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2)
    a.layers["counts"] = a.X.copy()
    a.X = np.log1p(a.X)
    assert p.pmp(a, MARKERS)["matrix"] == "layers['counts']"
    del a.layers["counts"]
    for fn in (p.crisp_purity, p.mecr, p.pmp):
        with pytest.raises(ValueError, match="layer="):
            fn(a, MARKERS)
    with warnings.catch_warnings(record=True) as w:               # layer=None: read X as it is, but say so
        warnings.simplefilter("always")
        assert p.pmp(a, MARKERS, layer=None)["matrix"] == "X"
    assert sum(issubclass(x.category, UserWarning) for x in w) == 1


def test_ncp_drops_the_unused_reference_label_key():
    import inspect

    from spatial_anno_metrics import purity as p

    params = list(inspect.signature(p.ncp).parameters)
    assert "ref_label_key" not in params and params[2] == "coexpr_threshold"


def test_ncp_returns_fraction_in_unit_interval():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=8)
    out = p.ncp(_section(0.0, seed=9), reference=ref)
    assert out["status"] in ("ok", "skipped")
    if out["status"] == "ok":
        assert 0.0 <= out["ncp"] <= 1.0
