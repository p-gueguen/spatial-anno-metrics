"""Tests for isc.py - inter-sample consistency: ISCG / ISCL, regimes, over-partitioning, dropped types."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _isc_adata(consistent: bool, n_per=40, n_types=3, n_samples=5, n_genes=30, seed=0):
    """Each (sample, type) is a pseudobulk. Consistent: a type's signature block is the SAME
    across samples. Inconsistent: each sample uses a DIFFERENT block for a given type."""
    rng = np.random.default_rng(seed)
    block = n_genes // n_types
    X, types, samples = [], [], []
    for s in range(n_samples):
        for t in range(n_types):
            sig = t if consistent else rng.integers(0, n_types)  # which block is high
            base = rng.normal(0.2, 0.1, size=(n_per, n_genes))
            base[:, sig * block:(sig + 1) * block] += 5.0
            X.append(base)
            types += [f"T{t}"] * n_per
            samples += [f"S{s}"] * n_per
    a = anndata.AnnData(X=np.vstack(X).astype("float32"))
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.obs["cell_type"] = pd.Categorical(types)
    a.obs["sample"] = pd.Categorical(samples)
    return a


def test_inter_sample_consistency_high_for_reproducible_types():
    from spatial_anno_metrics import isc as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    assert good["consistency"] > bad["consistency"]


def test_inter_sample_consistency_two_axes_and_regime():
    """Garnica 2026: ISCG (2-label silhouette | pseudobulk cosine) and ISCL (silhouette | reciprocal
    classification match) both separate reproducible from non-reproducible types; product + regime."""
    from spatial_anno_metrics import isc as em

    good = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    bad = em.inter_sample_consistency(_isc_adata(consistent=False), "cell_type", "sample")
    for k in ("iscg", "iscl", "integrated"):
        assert good[k] > bad[k], k
    assert good["iscl"] > 0.9 and set(good["regime"].values()) == {"consistent"}
    assert good["dropped_types"] == {} and good["n_profiles"] == 15
    assert {k: good["thresholds"][k] for k in ("iscg", "iscl")} == {"iscg": 0.74, "iscl": 0.45}
    assert good["status"] == "ok" and good["matrix"] == "X"          # layer=None: X averaged as it is


def test_inter_sample_consistency_over_partitioning_drops_local_not_global():
    """Benchmark task 2 (the paper's positive control for ISCL): split one type into two identical
    halves in every sample. The halves cannot be reciprocally matched -> ISCL collapses; they stay
    compact -> ISCG stays high; regime = 'merge'. Untouched types stay 'consistent'."""
    from spatial_anno_metrics import isc as em

    base = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample")
    a = _isc_adata(consistent=True)
    lab = a.obs["cell_type"].astype(str).to_numpy()
    lab[np.flatnonzero(lab == "T0")[::2]] = "T0.b"  # every 2nd T0 cell, i.e. 20/20 per sample
    a.obs["cell_type"] = pd.Categorical(lab)
    split = em.inter_sample_consistency(a, "cell_type", "sample")
    for t in ("T0", "T0.b"):
        assert split["per_type_iscl"][t] < base["per_type_iscl"]["T0"] - 0.5, t
        assert split["per_type_iscg"][t] > 0.74, t
        assert split["regime"][t] == "merge", t
    assert split["regime"]["T1"] == "consistent" and split["regime"]["T2"] == "consistent"


def test_inter_sample_consistency_drops_underpowered_types_and_discloses():
    """A type present in fewer than min_samples samples is reported under dropped_types, not scored
    (the paper's 3 % 'not computable'); lowering min_samples brings it back."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True)
    lab, smp = a.obs["cell_type"].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    a = a[~((lab == "T2") & np.isin(smp, ["S0", "S1", "S2"]))].copy()  # T2 survives in 2 samples only
    out = em.inter_sample_consistency(a, "cell_type", "sample")
    assert out["dropped_types"] == {"T2": 2} and "T2" not in out["per_type_iscl"]
    relaxed = em.inter_sample_consistency(a, "cell_type", "sample", min_samples=2)
    assert relaxed["dropped_types"] == {} and "T2" in relaxed["per_type_iscl"]
    none = em.inter_sample_consistency(a, "cell_type", "sample", min_samples=6)    # only 5 samples
    assert none["status"] == "skipped" and ">= 6 samples" in none["reason"]
    assert none["iscg"] is None and none["consistency"] is None and len(none["dropped_types"]) == 3


def test_inter_sample_consistency_counts_layer_recipe():
    """layer= names a raw-count layer (was counts_layer=): pseudobulks are summed counts, log1p CP10k."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True)
    a.layers["counts"] = np.round(np.clip(a.X, 0, None) * 10)
    out = em.inter_sample_consistency(a, "cell_type", "sample", layer="counts")
    assert out["status"] == "ok" and out["matrix"] == "layers['counts']" and out["iscl"] > 0.9
    with pytest.raises(TypeError):
        em.inter_sample_consistency(a, "cell_type", "sample", counts_layer="counts")


def test_isc_regime_quadrants():
    from spatial_anno_metrics import isc as em

    assert em.isc_regime(0.9, 0.9) == "consistent"
    assert em.isc_regime(0.9, 0.1) == "merge"
    assert em.isc_regime(0.1, 0.9) == "split"
    assert em.isc_regime(0.1, 0.1) == "remove"
    assert em.isc_regime(None, 0.5) == "not_computable"
    assert isinstance(em.isc_regime(0.9, 0.9), str)                  # a helper: returns a str, not a dict


def test_headline_iscg_iscl_are_means_over_labels():
    """A label weighs once in iscg / iscl whatever its number of samples (paper); the pseudobulk-weighted mean
    (the pre-0.7.0 headline) is kept under *_profile_mean."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=False, n_samples=6, seed=3)
    lab, smp = a.obs["cell_type"].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    a = a[~((lab == "T2") & (smp == "S0"))].copy()                    # T2: 5 profiles, T0/T1: 6
    out = em.inter_sample_consistency(a, "cell_type", "sample")
    for k in ("iscg", "iscl"):
        per = out[f"per_type_{k}"]
        assert out[k] == pytest.approx(np.mean(list(per.values())), abs=1e-15)
        assert out[f"{k}_profile_mean"] == pytest.approx((6 * per["T0"] + 6 * per["T1"] + 5 * per["T2"]) / 17, abs=1e-12)
    assert out["iscg"] != pytest.approx(out["iscg_profile_mean"], abs=1e-6)  # the fixture tells the two apart


def test_iscg_reduction_pca_default_and_raw_space():
    """Default ISCG space = scTypeEval's (cosine on 30 centred PCs, capped at n_profiles - 1); reduction=None is
    cosine on the pseudobulks themselves (the pre-0.7.0 space). ISCL and the legacy silhouette ignore it."""
    from scipy.spatial.distance import cdist

    from spatial_anno_metrics import isc as em
    from spatial_anno_metrics.thresholds import THRESHOLDS

    a = _isc_adata(consistent=False, seed=1)
    raw = em.inter_sample_consistency(a, "cell_type", "sample", reduction=None)
    pca = em.inter_sample_consistency(a, "cell_type", "sample")
    assert pca["reduction"] == "pca" and pca["n_pcs"] == 14 and raw["reduction"] is None and raw["n_pcs"] is None
    P, lab, smp, _, _ = em._pseudobulk_profiles(a, "cell_type", "sample", min_cells=10, min_samples=5)
    D = cdist(P, P, "cosine")
    np.fill_diagonal(D, 0)
    g = em._two_label_silhouette(D, lab)
    assert raw["per_type_iscg"] == pytest.approx({t: g[lab == t].mean() for t in np.unique(lab)}, abs=1e-12)
    U, S, _ = np.linalg.svd(P - P.mean(0), full_matrices=False)
    E = U[:, :14] * S[:14]
    g = em._two_label_silhouette(cdist(E, E, "cosine"), lab)
    assert pca["per_type_iscg"] == pytest.approx({t: g[lab == t].mean() for t in np.unique(lab)}, abs=1e-12)
    assert abs(pca["iscg"] - raw["iscg"]) > 1e-3                      # positive control: the spaces differ here
    assert pca["per_type_iscl"] == raw["per_type_iscl"] and pca["per_type"] == raw["per_type"]
    import inspect
    sig = inspect.signature(em.inter_sample_consistency).parameters
    assert sig["n_pcs"].default == THRESHOLDS["isc"]["n_pcs"]["value"] == 30 and sig["reduction"].default == "pca"
    with pytest.raises(ValueError, match="reduction"):
        em.inter_sample_consistency(a, "cell_type", "sample", reduction="umap")


def _two_sample_profiles():
    """Sample A holds X, Y; sample B holds X, Z. Y and Z share one profile, so A_Y <-> B_Z is a mutual hit
    between two labels that are each absent from the other sample."""
    x, yz = np.r_[5.0, 4, 3, 2, 1, 0], np.r_[0.0, 1, 2, 3, 4, 5]
    P = np.vstack([x, yz, x + 0.1, yz + [0, 0, 0, 0, 0, 0.1]])
    return P, np.array(["A", "A", "B", "B"]), np.array(["X", "Y", "X", "Z"])


def test_rcm_half_rule_for_labels_absent_from_the_other_sample():
    """scTypeEval match_tidy / match_pairs: a mutual hit scores 0 only when both labels occur in both samples,
    0.5 when either label is absent from the other sample; get_consistency reads 0.5 as 1 (no match)."""
    from spatial_anno_metrics import isc as em

    P, smp, lab = _two_sample_profiles()
    D = em._rcm_dissimilarity(P, smp, lab)
    assert D[0, 2] == D[2, 0] == 0.0                                  # X <-> X: shared label, a match
    assert D[1, 3] == D[3, 1] == 0.5                                  # Y <-> Z: the 0.5 rule
    assert D[0, 3] == D[1, 2] == 1.0 and np.all(np.diag(D) == 0)


def test_rcm_half_rule_moves_iscl():
    """Y lives in S0-S2, Z in S3-S5, with the same expression: the cross-label Y <-> Z mutual hits no longer
    count as matches, so Y and Z are not penalised for a label the other samples cannot offer."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True, n_samples=6, seed=2)
    lab, smp = a.obs["cell_type"].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    lab = np.where((lab == "T2") & np.isin(smp, ["S3", "S4", "S5"]), "T2z", lab)
    a.obs["cell_type"] = pd.Categorical(lab)
    out = em.inter_sample_consistency(a, "cell_type", "sample", min_samples=3)
    P, plab, psmp, _, _ = em._pseudobulk_profiles(a, "cell_type", "sample", min_cells=10, min_samples=3)
    D = em._rcm_dissimilarity(P, psmp, plab)
    cross = np.isin(plab, ["T2", "T2z"])[:, None] & np.isin(plab, ["T2", "T2z"])[None, :] & (plab[:, None] != plab)
    assert (D[cross] == 0.5).sum() == 18                              # 3 x 3 pairs, both triangles
    assert out["per_type_iscl"]["T2"] == out["per_type_iscl"]["T2z"] == 1.0
    D[D == 0.5] = 0.0                                                 # the pre-0.7.0 rule: a mutual hit is a match
    from sklearn.metrics import silhouette_samples
    old = silhouette_samples(D, plab, metric="precomputed")[plab == "T2"].mean()
    assert old < 0.6                                                  # positive control: the rule moves ISCL


def test_rcm_skips_samples_with_one_profile():
    """scTypeEval recip_classif skips a sample pair when either sample has < 2 labels: distance stays 1."""
    from spatial_anno_metrics import isc as em

    P, smp, lab = _two_sample_profiles()
    P, smp, lab = P[[0, 1, 2]], smp[[0, 1, 2]], lab[[0, 1, 2]]       # B keeps only X
    D = em._rcm_dissimilarity(P, smp, lab)
    assert D[0, 2] == 1.0 and D[0, 1] == 1.0 and np.all(np.diag(D) == 0)


def test_rcm_genes_restricts_the_classifier():
    """rcm_genes=<gene list> runs the reciprocal classification on those genes only (e.g. HVGs); ISCG keeps every
    gene. With layer=None this equals scoring the subset AnnData's ISCL."""
    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=False, seed=4)
    a.var_names = [f"g{i}" for i in range(a.n_vars)]
    genes = list(a.var_names[::3])
    sub = em.inter_sample_consistency(a, "cell_type", "sample", rcm_genes=genes)
    ref = em.inter_sample_consistency(a[:, genes].copy(), "cell_type", "sample")
    full = em.inter_sample_consistency(a, "cell_type", "sample")
    assert sub["per_type_iscl"] == pytest.approx(ref["per_type_iscl"], abs=1e-12)
    assert sub["per_type_iscg"] == full["per_type_iscg"] and sub["rcm_genes"] == len(genes)
    assert full["rcm_genes"] == "all"
    with pytest.raises(KeyError, match="nope"):
        em.inter_sample_consistency(a, "cell_type", "sample", rcm_genes=["g0", "nope"])


def test_classic_markers_union_of_top_positive_differences():
    """SingleR classic de.n = round(500 * (2/3)^log2(N)); per ordered label pair the top de.n genes with a positive
    difference; the union is used. With N = 2 and 3 genes, de.n = 333 caps at the positive genes."""
    from spatial_anno_metrics import isc as em

    Pk = np.array([[3.0, 1.0, 0.0], [0.0, 1.0, 2.0]])
    assert sorted(em._classic_markers(Pk)) == [0, 2]                  # gene 1 never differs
    assert sorted(em._classic_markers(Pk, n=1)) == [0, 2]
    Pk = np.array([[3.0, 2.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
    assert sorted(em._classic_markers(Pk, n=1)) == [0, 3]             # top-1 per pair: 0 beats 1
    assert em._classic_marker_n(7) == round(500 * (2 / 3) ** np.log2(7)) == 160


def test_pseudobulk_sparse_equals_dense():
    """indicator.T @ X aggregation: a sparse count layer gives the same profiles as the dense one."""
    import scipy.sparse as sp

    from spatial_anno_metrics import isc as em

    a = _isc_adata(consistent=True)
    a.layers["counts"] = np.round(np.clip(a.X, 0, None) * 10)
    P1, l1, s1, d1, _ = em._pseudobulk_profiles(a, "cell_type", "sample", min_cells=10, min_samples=5, layer="counts")
    a.layers["counts"] = sp.csr_matrix(a.layers["counts"])
    P2, l2, s2, d2, _ = em._pseudobulk_profiles(a, "cell_type", "sample", min_cells=10, min_samples=5, layer="counts")
    assert np.allclose(P1, P2, rtol=0, atol=1e-12) and (l1 == l2).all() and (s1 == s2).all() and d1 == d2
    a.X = sp.csr_matrix(a.X)
    Pm, *_ = em._pseudobulk_profiles(a, "cell_type", "sample", min_cells=10, min_samples=5)
    dense = np.asarray(a.X.todense())
    lab, smp = a.obs["cell_type"].astype(str).to_numpy(), a.obs["sample"].astype(str).to_numpy()
    assert np.allclose(Pm[0], dense[(lab == "T0") & (smp == "S0")].mean(0), atol=1e-6)  # layer=None: the mean


def _mixture_adata(seed=0, n_samples=6, n_types=5, n_genes=80, n_cells=500, batch_sd=0.0, same_composition=False):
    """Raw counts with a graded per-gene baseline (as real genes have) and a per-sample cell-type composition
    drawn from a Dirichlet, so a label made of random cells differs between samples (as on real data). Each
    type adds a marker block. ``batch_sd`` > 0 multiplies every gene of a sample by its own log-normal factor;
    ``same_composition`` gives every sample the same (unequal) type counts, as the synthetic world does."""
    rng = np.random.default_rng(seed)
    base = rng.lognormal(0, 1, n_genes)
    rates = np.tile(base, (n_types, 1))
    block = n_genes // (2 * n_types)
    for t in range(n_types):
        rates[t, t * block:(t + 1) * block] += 8.0
    rates /= rates.sum(1, keepdims=True)
    X, types, samples = [], [], []
    fixed = 20 + np.round(np.linspace(0.5, 1.5, n_types) / n_types * (n_cells - 20 * n_types)).astype(int)
    for s in range(n_samples):
        n_t = fixed if same_composition else 20 + rng.multinomial(n_cells - 20 * n_types,
                                                                  rng.dirichlet(np.full(n_types, 1.5)))
        f = rng.lognormal(0, batch_sd, n_genes) if batch_sd else np.ones(n_genes)
        for t in range(n_types):
            lam = rates[t] * f / (rates[t] * f).sum() * rng.lognormal(np.log(300), 0.3, (n_t[t], 1))
            X.append(rng.poisson(lam))
            types += [f"T{t}"] * n_t[t]
            samples += [f"S{s}"] * n_t[t]
    X = np.vstack(X).astype("float32")
    a = anndata.AnnData(X=X)
    a.layers["counts"] = X.copy()
    a.obs_names = [f"c{i}" for i in range(a.n_obs)]
    a.var_names = [f"g{i}" for i in range(n_genes)]
    a.obs["cell_type"] = types
    a.obs["sample"] = samples
    return a


def _isc(a, **kw):
    from spatial_anno_metrics import isc as em

    return em.inter_sample_consistency(a, "cell_type", "sample", layer="counts", **kw)


def test_regime_status_and_threshold_echo():
    out = _isc(_mixture_adata())
    assert out["regime_status"] == "benchmarked"
    assert {"iscg", "iscl", "regime_partner_rate", "regime_partner_ratio", "regime_remove", "regime_low",
            "regime_floor_gap"} <= set(
        out["thresholds"])
    assert set(out["regime"].values()) == {"consistent"} and set(out["quadrant"].values()) == {"consistent"}
    ev = out["regime_evidence"]["T0"]
    assert set(ev) == {"iscg_calibrated", "self_match", "partner", "partner_rate"}
    assert ev["partner"] is None and ev["partner_rate"] == 0.0 and ev["self_match"] == 1.0
    assert out["iscg_floor"] < 0.74 and out["low_global_signal"] is False


def test_regime_merge_names_its_reciprocal_partner():
    """Benchmark task 2 on counts with real-like composition: two random halves of T0 in every sample are each
    other's reciprocal best hit across samples -> both 'merge', each naming the other; untouched types stay."""
    a = _mixture_adata(seed=1)
    lab, smp = a.obs["cell_type"].to_numpy().copy(), a.obs["sample"].to_numpy()
    rng = np.random.default_rng(1)
    for s in np.unique(smp):
        ids = np.flatnonzero((lab == "T0") & (smp == s))
        lab[rng.choice(ids, len(ids) // 2, replace=False)] = "T0.b"
    a.obs["cell_type"] = lab
    out = _isc(a)
    for t, u in (("T0", "T0.b"), ("T0.b", "T0")):
        assert out["regime"][t] == "merge", (t, out["regime"], out["regime_evidence"][t])
        assert out["regime_evidence"][t]["partner"] == u and out["regime_evidence"][t]["partner_rate"] >= 0.05
    assert {out["regime"][t] for t in ("T1", "T2", "T3", "T4")} == {"consistent"}


def test_regime_never_calls_split_merged_types_read_consistent():
    """Under-partitioning is not identifiable from ISC: a label that pools two distinct types is compact and
    replicable when its parts are (real data: 12 / 12 merged pairs, ISCG >= 0.88, ISCL 1.0). The call says
    'consistent', never a 'split' that ISC cannot support."""
    a = _mixture_adata(seed=2)
    lab = a.obs["cell_type"].to_numpy().copy()
    lab[np.isin(lab, ["T0", "T1"])] = "T0+T1"
    a.obs["cell_type"] = lab
    out = _isc(a)
    assert out["regime"]["T0+T1"] == "consistent"
    assert "split" not in set(out["regime"].values())


def test_regime_batch_effect_is_not_a_label_problem():
    """A per-sample batch effect lowers ISCG for every label (the paper's fixed 0.74 guide then reads 'split'),
    while the labels are right: calibrated against this dataset's own labels the calls stay 'consistent', and
    the dataset-level flag says the global signal is low."""
    out = _isc(_mixture_adata(seed=3, batch_sd=1.0))
    assert "split" in set(out["quadrant"].values())                     # positive control: the old guide misfires
    assert set(out["regime"].values()) == {"consistent"}, out["regime_evidence"]
    assert out["low_global_signal"] is True


def test_regime_random_label_reads_remove():
    """A label made of random cells of each sample (multiplets / low-quality: the paper's 'remove') is as
    compact as the per-sample average profile (the ISCG floor), not as a typical label -> 'remove'."""
    a = _mixture_adata(seed=4)
    lab, smp = a.obs["cell_type"].to_numpy().copy(), a.obs["sample"].to_numpy()
    rng = np.random.default_rng(4)
    for s in np.unique(smp):
        ids = np.flatnonzero(smp == s)
        lab[rng.choice(ids, 30, replace=False)] = "Junk"
    a.obs["cell_type"] = lab
    out = _isc(a)
    assert out["regime"]["Junk"] == "remove", out["regime_evidence"]["Junk"]
    assert out["regime_evidence"]["Junk"]["iscg_calibrated"] < 0.25
    assert {out["regime"][t] for t in ("T0", "T1", "T2", "T3", "T4")} == {"consistent"}


def test_regime_falls_back_to_iscl_when_the_floor_is_not_separable():
    """Identical composition in every sample (the synthetic world) makes a label of random cells as compact and
    replicable as a real type (floor ISCG 0.99): ISCG cannot flag it, iscg_calibrated is None and only ISCL and
    the partner evidence decide. The random label then reads 'consistent' - a documented blind spot, while the
    same label in a cohort whose composition varies reads 'remove' (test above)."""
    a = _mixture_adata(seed=0, same_composition=True)
    lab, smp = a.obs["cell_type"].to_numpy().copy(), a.obs["sample"].to_numpy()
    rng = np.random.default_rng(0)
    for s in np.unique(smp):
        lab[rng.choice(np.flatnonzero(smp == s), 30, replace=False)] = "Junk"
    a.obs["cell_type"] = lab
    out = _isc(a)
    assert out["iscg_floor"] > 0.95
    assert all(ev["iscg_calibrated"] is None for ev in out["regime_evidence"].values())
    assert set(out["regime"].values()) == {"consistent"}


def test_reciprocal_partners_counts_mutual_cross_sample_hits():
    """Per label: share of its (profile, other sample) comparisons that end in a reciprocal best hit with its
    own label (self_match) or with its most frequent other label (partner_rate). 0.5 entries are hits."""
    from spatial_anno_metrics import isc as em

    smp = np.array(["A", "A", "B", "B", "C", "C"])
    lab = np.array(["X", "Y", "X", "Y", "X", "Y"])
    D = np.ones((6, 6))
    np.fill_diagonal(D, 0)
    for i, j in ((0, 3), (1, 4), (2, 4)):                              # X_A~Y_B, Y_A~X_C, X_B~X_C (0.5 too)
        D[i, j] = D[j, i] = 0.0
    D[2, 4] = D[4, 2] = 0.5
    ev = em._reciprocal_partners(D, lab, smp)
    assert ev["X"] == {"self_match": 2 / 6, "partner": "Y", "partner_rate": 2 / 6}
    assert ev["Y"] == {"self_match": 0.0, "partner": "X", "partner_rate": 2 / 6}
    D[:] = 1
    np.fill_diagonal(D, 0)
    assert em._reciprocal_partners(D, lab, smp)["X"] == {"self_match": 0.0, "partner": None, "partner_rate": 0.0}
    # Z lives in A and B only: its comparisons are Z_A -> B and Z_B -> A (2, not 2 x 2); a hit of Z_A inside
    # C (which holds no Z) is not counted, so a rare label is not diluted by the samples it is absent from.
    smp = np.array(["A", "A", "A", "B", "B", "B", "C", "C"])
    lab = np.array(["X", "Y", "Z", "X", "Y", "Z", "X", "Y"])
    D = np.ones((8, 8))
    np.fill_diagonal(D, 0)
    D[2, 5] = D[5, 2] = 0.0                                            # Z_A ~ Z_B
    D[2, 6] = D[6, 2] = 0.0                                            # Z_A ~ X_C: outside Z's samples
    ev = em._reciprocal_partners(D, lab, smp)
    assert ev["Z"] == {"self_match": 1.0, "partner": None, "partner_rate": 0.0}


def test_isc_zero_pseudobulk_and_markers_mode_and_regime_boundaries():
    from spatial_anno_metrics import isc as em

    # Zero-norm pseudobulk does not poison iscg or regime calls to NaN
    a = _isc_adata(consistent=True)
    a.X[:40] = 0.0
    out_raw = em.inter_sample_consistency(a, "cell_type", "sample", reduction=None)
    out_pca = em.inter_sample_consistency(a, "cell_type", "sample", reduction="pca")
    assert np.isfinite(out_raw["iscg"]) and np.isfinite(out_pca["iscg"])
    assert "not_computable" not in set(out_raw["regime"].values())

    # rcm_genes="markers" end-to-end
    out_mk = em.inter_sample_consistency(_isc_adata(consistent=True), "cell_type", "sample", rcm_genes="markers")
    assert out_mk["status"] == "ok" and out_mk["rcm_genes"] == "markers" and out_mk["iscl"] > 0.9

    # Exact boundary equality for isc_regime (0.74, 0.45) and NaN handling
    assert em.isc_regime(0.74, 0.45) == "consistent"
    assert em.isc_regime(0.74, 0.4499) == "merge"
    assert em.isc_regime(0.7399, 0.45) == "split"
    assert em.isc_regime(0.7399, 0.4499) == "remove"
    assert em.isc_regime(float("nan"), 0.5) == "not_computable"


def test_isc_min_cells_exact_boundary():
    """Pin exact min_cells=10 default boundary: 9 cells per (sample, type) skips, 10 succeeds."""
    from spatial_anno_metrics import isc as em

    r9 = em.inter_sample_consistency(_isc_adata(consistent=True, n_per=9), "cell_type", "sample")
    r10 = em.inter_sample_consistency(_isc_adata(consistent=True, n_per=10), "cell_type", "sample")
    assert r9["status"] == "skipped" and r9["n_profiles"] == 0
    assert r10["status"] == "ok" and r10["n_profiles"] == 15


def test_merge_pairs_need_mutual_partners_matched_like_themselves():
    """Random halves of one type hit each other about as often as themselves (larger ratio >= 0.5 on real
    data); two labels sharing a contaminant are mutual partners too, but each matches itself 4x more often
    (ratio <= 0.25: Mast / Normal_Kerat with 40 % T cells each). Only the first pair is 'merge'."""
    from spatial_anno_metrics import isc as em

    thr = {"regime_partner_rate": 0.05, "regime_partner_ratio": 1 / 3}

    def ev(partner, rate, self_match):
        return {"self_match": self_match, "partner": partner, "partner_rate": rate}

    halves = {"A": ev("A.b", 0.24, 0.30), "A.b": ev("A", 0.24, 0.10), "C": ev(None, 0.0, 1.0)}
    assert em._merge_pairs(halves, thr) == {"A", "A.b"}
    shared = {"M": ev("N", 0.12, 0.48), "N": ev("M", 0.12, 0.52), "C": ev(None, 0.0, 1.0)}
    assert em._merge_pairs(shared, thr) == set()
    assert em._merge_pairs(shared, {**thr, "regime_partner_ratio": 0.2}) == {"M", "N"}   # positive control
    one_way = {"A": ev("B", 0.3, 0.3), "B": ev("C", 0.3, 0.3), "C": ev("B", 0.3, 0.3)}
    assert em._merge_pairs(one_way, thr) == {"B", "C"}                  # A -> B is not returned by B
    rare = {"A": ev("A.b", 0.04, 0.04), "A.b": ev("A", 0.04, 0.04)}
    assert em._merge_pairs(rare, thr) == set()                          # below regime_partner_rate
