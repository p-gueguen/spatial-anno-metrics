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
    assert "sam_crisp_impure" in clean.obs and clean.obs["sam_crisp_impure"].sum() == 0
    assert dirty.obs["sam_crisp_impure"].sum() > 0


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
    v = a.obs["sam_pmp"].to_numpy()
    assert out["status"] == "ok" and out["mean_pmp"] == pytest.approx(np.nanmean(v))
    assert "sam_pmp" in a.obs and "pmp" not in a.obs
    finite = v[np.isfinite(v)]
    assert finite.size and finite.min() >= 0.0 and finite.max() <= 1.0
    # cells contaminated with the other lineage (first 60) have LOWER purity than the clean tail.
    assert np.nanmean(v[:60]) < np.nanmean(v[120:])


def test_purity_outputs_key_added_and_uns_summary(tmp_path):
    """obs columns are named by key_added (default sam_<name>); the scalar summary goes to uns['sam'][fn]
    and stays h5ad-writable."""
    from spatial_anno_metrics import purity as p
    from spatial_anno_metrics.purity import _pmp_values

    a = _section(0.3, seed=2)
    out = p.pmp(a, MARKERS, key_added="my_pmp")
    assert "my_pmp" in a.obs and "sam_pmp" not in a.obs
    np.testing.assert_array_equal(a.obs["my_pmp"].to_numpy(), _pmp_values(a, MARKERS, "cell_type", "counts"))
    assert a.uns["sam"]["pmp"] == {"status": "ok", "mean_pmp": out["mean_pmp"], "n_defined": out["n_defined"],
                                   "no_dict_frac": 0.0, "matrix": "X", "key_added": "my_pmp"}
    p.crisp_purity(a, MARKERS)
    ref = _section(0.0, seed=5)
    p.nmp(a, reference=ref, score="raw")
    assert {"sam_crisp_impure", "sam_nmp"} <= set(a.obs) and {"pmp", "crisp_purity", "nmp"} <= set(a.uns["sam"])
    assert "types_without_negatives" not in a.uns["sam"]["nmp"]          # lists and frames stay in the result
    a.write_h5ad(tmp_path / "a.h5ad")
    back = anndata.read_h5ad(tmp_path / "a.h5ad")
    assert back.uns["sam"]["pmp"]["key_added"] == "my_pmp"


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
    out = p.nmp(dirty, reference=ref, label_key="cell_type", ref_label_key="cell_type", score="raw")
    assert out["status"] == "ok"
    assert "sam_nmp" in dirty.obs and 0.0 <= out["mean_nmp"] <= 1.0
    assert out["mean_nmp"] > p.nmp(_section(0.0, seed=7), reference=ref, score="raw")["mean_nmp"]


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

    assert p.nmp(t_cells(0.0), reference=ref, score="raw")["mean_nmp"] == pytest.approx(0.0)
    # truth: 11 Fib counts / (100 T + 20 shared + 11 Fib + 1.8 silent) = 0.083
    assert p.nmp(t_cells(0.11), reference=ref, score="raw")["mean_nmp"] == pytest.approx(0.083, abs=0.01)


def test_nmp_reports_types_it_cannot_score():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=5)
    # Mal is absent from the reference; T has no other reference type to take negatives from
    ref = ref[ref.obs["cell_type"] == "T"].copy()
    out = p.nmp(_section(0.0, seed=6), reference=ref, score="raw")
    assert out["status"] == "skipped"                      # no cell scorable: nothing measured, not 0.0
    assert out["types_absent_from_reference"] == ["Mal"] and out["types_without_negatives"] == ["T"]


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


def test_ncp_scores_clean_as_preserved_and_contamination_as_lost():
    from spatial_anno_metrics import purity as p

    mk = MARKERS["T"] + MARKERS["Mal"]
    ref = _section(0.0, seed=8)[:, mk].copy()    # T and Mal markers never co-occur: 9 non-coexpressed pairs
    clean = p.ncp(_section(0.0, seed=9)[:, mk].copy(), reference=ref)
    assert clean["status"] == "ok" and clean["ncp"] == pytest.approx(1.0)
    dirty = p.ncp(_section(0.5, seed=9)[:, mk].copy(), reference=ref)   # 100 of 200 cells co-detect both
    assert dirty["status"] == "ok" and dirty["ncp"] == pytest.approx(0.0)


# --------------------------------------------------------------------------- #
# nmp, score="raw": undefined is NaN, exact boundaries, >= 3 reference types
# --------------------------------------------------------------------------- #
def _nmp_world(n_t_detect, n_t=220, other=(("Fib", 100, 30),), sparse_ref=False):
    """A reference where COL1A1 is detected in exactly ``n_t_detect`` of ``n_t`` T cells and, per other type
    ``(name, n_cells, n_detect)``, in ``n_detect`` of its cells; plus 10 T section cells, each 8 CD3E +
    2 COL1A1 counts (truth: nmp 0.2 when COL1A1 is a T negative, else 0)."""
    from scipy import sparse

    rows, labels = [], []
    for name, n, k in (("T", n_t, n_t_detect), *other):
        block = np.zeros((n, 2), dtype="float32")
        block[:, 0] = 5 if name == "T" else 0                  # CD3E
        block[:k, 1] = 3                                       # COL1A1
        rows.append(block)
        labels += [name] * n
    X = np.vstack(rows)
    ref = anndata.AnnData(X=sparse.csr_matrix(X) if sparse_ref else X)
    ref.var_names = ["CD3E", "COL1A1"]
    ref.obs["cell_type"] = labels
    sec = anndata.AnnData(X=np.tile(np.array([[8, 2]], dtype="float32"), (10, 1)))
    sec.var_names = ["CD3E", "COL1A1"]
    sec.obs["cell_type"] = "T"
    return sec, ref


@pytest.mark.parametrize("sparse_ref", [False, True])
def test_nmp_inclusive_max_neg_frac_does_not_depend_on_storage(sparse_ref):
    """11/220 T cells is exactly 0.05: a negative (<= max_neg_frac), dense or CSR. scipy's sparse .mean()
    put it just above 0.05, so the CSR reference listed T as unscorable."""
    from spatial_anno_metrics import purity as p

    sec, ref = _nmp_world(11, sparse_ref=sparse_ref)
    assert p.nmp(sec, ref, score="raw")["mean_nmp"] == pytest.approx(0.2)
    sec, ref = _nmp_world(12, sparse_ref=sparse_ref)           # 12/220 > 0.05: not a negative
    out = p.nmp(sec, ref, score="raw")
    assert out["status"] == "skipped" and out["types_without_negatives"] == ["T"]


def test_nmp_inclusive_min_other_frac_and_max_over_other_types():
    from spatial_anno_metrics import purity as p

    sec, ref = _nmp_world(0, other=(("Fib", 100, 25),))        # exactly 25 % of Fib: >= min_other_frac
    assert p.nmp(sec, ref, score="raw")["mean_nmp"] == pytest.approx(0.2)
    sec, ref = _nmp_world(0, other=(("Fib", 100, 24),))
    assert p.nmp(sec, ref, score="raw")["status"] == "skipped"
    # SOME other type: 30 % of Fib qualifies although B never detects COL1A1 (a min over others would not)
    sec, ref = _nmp_world(0, other=(("Fib", 100, 30), ("B", 100, 0), ("Endo", 100, 0)))
    assert p.nmp(sec, ref, score="raw")["mean_nmp"] == pytest.approx(0.2)


def test_nmp_raw_unscorable_cells_are_nan_and_excluded_from_the_mean():
    from spatial_anno_metrics import purity as p

    sec, ref = _nmp_world(0)
    extra = anndata.AnnData(X=np.array([[0, 0]] * 3 + [[8, 2]] * 5, dtype="float32"))
    extra.var_names = sec.var_names
    extra.obs["cell_type"] = ["T"] * 3 + ["Mal"] * 5         # 3 zero-count T cells; Mal not in the reference
    both = anndata.concat([sec, extra], index_unique="-")
    out = p.nmp(both, ref, score="raw")
    v = both.obs["sam_nmp"].to_numpy()
    assert np.isnan(v[10:]).all() and v[:10] == pytest.approx(0.2)
    assert out["mean_nmp"] == pytest.approx(0.2) and out["n_defined"] == 10
    assert out["types_absent_from_reference"] == ["Mal"] and out["types_without_negatives"] == []


# --------------------------------------------------------------------------- #
# goldens: txsim / Moldia (tests/goldens/*/SOURCE.md)
# --------------------------------------------------------------------------- #
from pathlib import Path  # noqa: E402

GOLDENS = Path(__file__).parent / "goldens"


def _golden_adata(name, fname):
    with np.load(GOLDENS / name / fname, allow_pickle=False) as z:
        a = anndata.AnnData(X=z["counts"].astype("float64"))
        a.var_names = [str(g) for g in z["genes"]]
        a.obs_names = [str(c) for c in z["cells"]]
        a.obs["cell_type"] = [str(x) for x in z["labels"]]
    return a


def _expected(name):
    import json

    return json.loads((GOLDENS / name / "expected.json").read_text())["results"]


@pytest.mark.parametrize("section", ["clean", "spillover"])
def test_nmp_salas_matches_txsim_reads(section):
    from spatial_anno_metrics import purity as p

    exp = _expected("nmp_txsim")[section]
    sec = _golden_adata("nmp_txsim", f"section_{section}.npz")
    ref = _golden_adata("nmp_txsim", "reference.npz")
    ours = p._nmp_salas_reads(sec.X, sec.obs["cell_type"].to_numpy(), list(sec.var_names), ref.X,
                              ref.obs["cell_type"].to_numpy(), list(ref.var_names))
    want = exp["txsim_reads"]
    assert ours["nmp"] == pytest.approx(want["overall"], abs=1e-12)
    assert set(ours["per_type"]) == set(want["per_celltype"])
    for t, v in want["per_celltype"].items():
        assert ours["per_type"][t] == pytest.approx(v, abs=1e-12), t
    assert set(ours["per_gene"]) == set(want["per_gene"])
    for g, v in want["per_gene"].items():
        assert ours["per_gene"][g] == pytest.approx(v, abs=1e-12), g
    # the label-permutation baseline, on the 20 stored permutations
    with np.load(GOLDENS / "nmp_txsim" / "label_permutations.npz") as z:
        perms = z[section]
    labels = sec.obs["cell_type"].to_numpy()
    ref_labels = ref.obs["cell_type"].to_numpy()
    permuted = [p._nmp_salas_reads(sec.X, labels[perm], list(sec.var_names), ref.X, ref_labels,
                                   list(ref.var_names))["nmp"] for perm in perms]
    np.testing.assert_allclose(permuted, exp["txsim_reads_permuted"]["values"], rtol=0, atol=1e-12)
    # public call: the same unscaled score; scaled against its own 20 permutations (not the stored ones)
    out = p.nmp(sec, ref)
    assert out["status"] == "ok" and out["score"] == "salas" and out["n_permutations"] == 20
    assert out["nmp_unscaled"] == pytest.approx(want["overall"], abs=1e-12)
    b = out["nmp_null"]
    assert out["nmp"] == pytest.approx((out["nmp_unscaled"] - b) / (1 - b))
    assert b == pytest.approx(exp["txsim_reads_permuted"]["mean"], abs=0.02)
    assert out["nmp"] == pytest.approx(exp["txsim_reads_permuted"]["scaled_overall"], abs=0.01)
    assert "sam_nmp" not in sec.obs                        # a per-section score: no per-cell column


@pytest.mark.parametrize("section", ["clean", "spillover"])
def test_nmp_moldia_matches_xenium_benchmarking_reads(section):
    """score='moldia' reproduces Moldia Xenium_benchmarking's relative-to-mean reads NMP to 1e-12."""
    from spatial_anno_metrics import purity as p

    exp = _expected("nmp_txsim")[section]
    sec = _golden_adata("nmp_txsim", f"section_{section}.npz")
    ref = _golden_adata("nmp_txsim", "reference.npz")
    out = p.nmp(sec, ref, score="moldia")
    assert out["status"] == "ok" and out["score"] == "moldia"
    assert out["nmp_unscaled"] == pytest.approx(exp["moldia_seg_reads"]["overall"], abs=1e-12)
    assert abs(exp["txsim_reads"]["overall"] - exp["moldia_seg_reads"]["overall"]) > 1e-4


@pytest.mark.parametrize("section", ["clean", "spillover"])
def test_nmp_cells_matches_txsim_cells_golden(section):
    """score='cells' reproduces txsim / Moldia negative_marker_purity_cells and its 20 stored permutations."""
    from spatial_anno_metrics import purity as p

    exp = _expected("nmp_txsim")[section]
    sec = _golden_adata("nmp_txsim", f"section_{section}.npz")
    ref = _golden_adata("nmp_txsim", "reference.npz")
    out = p.nmp(sec, ref, score="cells")
    assert out["status"] == "ok" and out["score"] == "cells"
    assert out["nmp_unscaled"] == pytest.approx(exp["txsim_cells"]["overall"], abs=1e-12)
    assert out["nmp_unscaled"] == pytest.approx(exp["moldia_seg_cells"]["overall"], abs=1e-12)
    with np.load(GOLDENS / "nmp_txsim" / "label_permutations.npz") as z:
        perms = z[section]
    r = p._nmp_salas_reads(sec.X, sec.obs["cell_type"].to_numpy(), list(sec.var_names), ref.X,
                           ref.obs["cell_type"].to_numpy(), list(ref.var_names), perms, criterion="cells")
    np.testing.assert_allclose(r["permuted"], exp["txsim_cells_permuted"]["values"], rtol=0, atol=1e-12)
    b = float(np.mean(r["permuted"]))
    assert b == pytest.approx(exp["txsim_cells_permuted"]["mean"], abs=1e-12)
    assert (r["nmp"] - b) / (1.0 - b) == pytest.approx(exp["txsim_cells_permuted"]["scaled_overall"], abs=1e-12)


@pytest.mark.parametrize("section", ["clean", "spillover"])
def test_ncp_moldia_matches_xenium_benchmarking_golden(section):
    """ncp(..., score='moldia') reproduces Moldia negative_marker_purity_coexpression to 1e-12."""
    from scipy import sparse

    from spatial_anno_metrics import purity as p

    exp = _expected("coexpr_txsim")[section]["moldia_ncp"]
    sec = _golden_adata("coexpr_txsim", f"section_{section}.npz")
    ref = _golden_adata("coexpr_txsim", "reference.npz")
    out = p.ncp(sec, ref, score="moldia")
    assert out["status"] == "ok" and out["score"] == "moldia"
    assert out["ncp"] == pytest.approx(exp["overall"], abs=1e-12)
    assert out["n_negative_pairs"] == exp["n_negative_pairs"]
    assert set(out["per_gene_as_conditioning"]) == set(exp["per_gene_as_conditioning"])
    for g, v in exp["per_gene_as_conditioning"].items():
        assert out["per_gene_as_conditioning"][g] == pytest.approx(v, abs=1e-12), g
    assert set(out["per_gene_as_target"]) == set(exp["per_gene_as_target"])
    for g, v in exp["per_gene_as_target"].items():
        assert out["per_gene_as_target"][g] == pytest.approx(v, abs=1e-12), g
    # Sparse CSR input matches dense exactly
    sec_s, ref_s = sec.copy(), ref.copy()
    sec_s.X, ref_s.X = sparse.csr_matrix(sec_s.X), sparse.csr_matrix(ref_s.X)
    assert p.ncp(sec_s, ref_s, score="moldia")["ncp"] == pytest.approx(out["ncp"], abs=1e-12)
    with pytest.raises(ValueError, match="score"):
        p.ncp(sec, ref, score="bad")


def test_nmp_salas_skips_with_fewer_than_two_eligible_types():
    from spatial_anno_metrics import purity as p

    ref = _section(0.0, seed=5)
    out = p.nmp(_section(0.0, seed=6), ref[ref.obs["cell_type"] == "T"].copy())
    assert out["status"] == "skipped" and out["types_absent_from_reference"] == ["Mal"]


@pytest.mark.parametrize("section", ["clean", "spillover"])
@pytest.mark.parametrize("by_type", [False, True])
def test_coexpression_similarity_matches_txsim(section, by_type):
    from spatial_anno_metrics import purity as p

    exp = _expected("coexpr_txsim")[section]
    sec = _golden_adata("coexpr_txsim", f"section_{section}.npz")
    ref = _golden_adata("coexpr_txsim", "reference.npz")
    out = p.coexpression_similarity(sec, ref, label_key="cell_type" if by_type else None,
                                    ref_label_key="cell_type")
    want = exp["txsim_coexpression_similarity_by_celltype" if by_type else "txsim_coexpression_similarity"]
    assert out["status"] == "ok" and out["matrix"] == "X"
    assert out["coexpression_similarity"] == pytest.approx(want["overall"], abs=1e-12)
    assert set(out["per_gene"]) == set(want["per_gene"])
    for g, v in want["per_gene"].items():                 # null = NaN: a gene no cell type scores
        assert out["per_gene"][g] == pytest.approx(np.nan if v is None else v, abs=1e-12, nan_ok=True), g
    if by_type:
        assert set(out["per_type"]) == set(want["per_celltype"])
        for t, v in want["per_celltype"].items():
            assert out["per_type"][t] == pytest.approx(v, abs=1e-12), t


def test_coexpression_similarity_skips_without_reference_or_genes():
    from spatial_anno_metrics import purity as p

    a = _section(0.0)
    assert p.coexpression_similarity(a, None)["status"] == "skipped"
    out = p.coexpression_similarity(a, _section(0.0, seed=3), min_cells=10_000)
    assert out["status"] == "skipped" and "min_cells" in out["reason"]


# --------------------------------------------------------------------------- #
# #8e: CRISP / MECR against a detection-matched null
# --------------------------------------------------------------------------- #
def test_crisp_mecr_report_a_detection_matched_null():
    from spatial_anno_metrics import purity as p

    clean = _section(0.0, off_target=40)            # 40 T cells with one ambient MLANA count: purity 0.8
    rc, rm = p.crisp_purity(clean, MARKERS), p.mecr(clean, MARKERS)
    assert rc["crisp_purity"] == pytest.approx(0.8) and rc["n_null"] == 20
    # independent detection at the same per-gene rates co-detects far more often than these exclusive lineages
    assert rc["crisp_purity_null"] < rc["crisp_purity"] and rc["vs_null"] < 0.5
    assert rc["vs_null"] == pytest.approx((1 - rc["crisp_purity"]) / (1 - rc["crisp_purity_null"]))
    assert rm["mecr_null"] > rm["mecr"] and rm["vs_null"] == pytest.approx(rm["mecr"] / rm["mecr_null"])
    assert p.crisp_purity(clean, MARKERS) == rc                          # seeded: reproducible
    dirty = _section(0.4, seed=1, off_target=40)
    assert p.crisp_purity(dirty, MARKERS)["vs_null"] > rc["vs_null"]
    assert p.mecr(dirty, MARKERS)["vs_null"] > rm["vs_null"]
    off = p.crisp_purity(clean, MARKERS, n_null=0)
    assert off["crisp_purity_null"] is None and off["vs_null"] is None
    assert off["crisp_purity"] == rc["crisp_purity"]


def test_null_positive_control_independent_markers_have_no_excess():
    """Markers detected independently of each other (background genes used as markers) sit on the null."""
    from spatial_anno_metrics import purity as p

    a = _section(0.0, n=2000, seed=4)
    bg = {"A": ["BG0", "BG1"], "B": ["BG2", "BG3"], "C": ["BG4"]}
    rc, rm = p.crisp_purity(a, bg), p.mecr(a, bg)
    # one realisation vs the null mean: +-0.06 over seeds 0-5 at this size; exclusive lineages read 0.24
    assert rc["vs_null"] == pytest.approx(1.0, abs=0.1) and rm["vs_null"] == pytest.approx(1.0, abs=0.1)


# --------------------------------------------------------------------------- #
# #8c dictionary coverage; int labels with int-keyed marker sets
# --------------------------------------------------------------------------- #
def test_marker_metrics_report_labels_and_lineages_without_markers():
    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2)
    a.obs["cell_type"] = a.obs["cell_type"].astype(str)
    a.obs.loc[a.obs_names[:20], "cell_type"] = "NK"                     # a label with no marker set
    ms = {**MARKERS, "B": ["NOT_ON_PANEL"]}
    out = p.pmp(a, ms)
    assert out["labels_without_markers"] == ["NK"] and out["no_dict_frac"] == pytest.approx(0.1)
    for fn in (p.crisp_purity, p.mecr):
        assert fn(a, ms)["lineages_without_panel_markers"] == ["B"]


def test_int_labels_with_int_keyed_marker_sets():
    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2)
    want = p.pmp(a.copy(), MARKERS)
    b = a.copy()
    b.obs["cell_type"] = (b.obs["cell_type"] == "Mal").astype(int).to_numpy()   # T -> 0, Mal -> 1
    got = p.pmp(b, {0: MARKERS["T"], 1: MARKERS["Mal"]})
    assert got["status"] == "ok" and got["mean_pmp"] == pytest.approx(want["mean_pmp"])
    assert got["labels_without_markers"] == []
    with pytest.raises(ValueError, match="both"):
        p.pmp(b, {0: MARKERS["T"], "0": MARKERS["Mal"]})


# --------------------------------------------------------------------------- #
# #15: sparse input gives the dense values
# --------------------------------------------------------------------------- #
def test_sparse_and_dense_inputs_agree():
    from scipy import sparse

    from spatial_anno_metrics import purity as p

    a = _section(0.3, seed=2, off_target=15)
    s = a.copy()
    s.X = sparse.csr_matrix(s.X)
    for fn in (p.crisp_purity, p.mecr):
        assert fn(a, MARKERS) == fn(s, MARKERS)
    assert p.pmp(a, MARKERS)["mean_pmp"] == p.pmp(s, MARKERS)["mean_pmp"]
    ref = _section(0.0, seed=8)
    s_ref = ref.copy()
    s_ref.X = sparse.csr_matrix(s_ref.X)
    assert p.ncp(a, ref) == p.ncp(s, s_ref)
    for by_type in (False, True):
        kw = {"label_key": "cell_type", "ref_label_key": "cell_type"} if by_type else {}
        da = p.coexpression_similarity(a, ref, **kw)
        sa = p.coexpression_similarity(s, s_ref, **kw)
        assert sparse.issparse(p._lognorm_cols(s.X, [0, 1, 2]))
        assert sa["coexpression_similarity"] == pytest.approx(da["coexpression_similarity"], abs=1e-12)
        for g in da["per_gene"]:
            assert sa["per_gene"][g] == pytest.approx(da["per_gene"][g], abs=1e-12, nan_ok=True)
        if by_type:
            for t in da["per_type"]:
                assert sa["per_type"][t] == pytest.approx(da["per_type"][t], abs=1e-12, nan_ok=True)


def test_ncp_matches_a_pair_loop():
    """The vectorised pair counts equal a plain per-pair loop (positive control: a looser cut differs)."""
    from itertools import combinations

    from spatial_anno_metrics import purity as p

    rng = np.random.default_rng(1)
    ref, sec = _section(0.0, seed=8), _section(0.2, seed=9)
    sec.X = sec.X + (rng.random(sec.X.shape) < 0.05)
    pos_r, pos_s = ref.X > 0, sec.X > 0

    def rate(P, i, j):
        e = (P[:, i] | P[:, j]).sum()
        return None if e == 0 else (P[:, i] & P[:, j]).sum() / e

    def loop(thr):
        pairs = [(i, j) for i, j in combinations(range(len(GENES)), 2)
                 if rate(pos_r, i, j) is not None and rate(pos_r, i, j) <= thr]
        kept = sum(1 for i, j in pairs if rate(pos_s, i, j) is None or rate(pos_s, i, j) <= thr)
        return kept / len(pairs)

    assert p.ncp(sec, ref)["ncp"] == pytest.approx(loop(0.1), abs=0)
    assert loop(0.1) != loop(0.3)


def test_ncp_threshold_is_inclusive():
    """A pair co-detected in exactly 1 of 10 detecting cells sits ON coexpr_threshold=0.1: non-coexpressed."""
    from spatial_anno_metrics import purity as p

    X = np.zeros((10, 2), dtype="float32")
    X[:5, 0], X[4:, 1] = 1, 1                    # cell 4 detects both: 1 both / 10 either = 0.1
    ref = anndata.AnnData(X=X.copy())
    ref.var_names = ["a", "b"]
    sec = ref.copy()
    assert p.ncp(sec, ref)["ncp"] == pytest.approx(1.0)
    sec.X[3, 1] = 1                              # 2 / 10 = 0.2 in the section: no longer preserved
    assert p.ncp(sec, ref)["ncp"] == pytest.approx(0.0)


@pytest.mark.parametrize("b_count, negative", [(199, False), (249, True)])
def test_nmp_salas_share_cut_is_strict(b_count, negative):
    """Type A holds exactly 1/200 = 0.005 of gene g: not a negative (< 0.005); 1/250 = 0.004 is."""
    from spatial_anno_metrics import purity as p

    def data(a_g, b_g):
        X = np.zeros((20, 2))
        X[:10, 0], X[10:, 0], X[10:, 1] = a_g, b_g, 5       # gene h (col 1) is a clean negative of A
        return X

    labels = np.array(["A"] * 10 + ["B"] * 10)
    r = p._nmp_salas_reads(data(2, 198), labels, ["g", "h"], data(1, b_count), labels, ["g", "h"])
    assert ("g" in r["per_gene"]) is negative and "h" in r["per_gene"]
    assert r["n_negative_pairs"] == 1 + negative
