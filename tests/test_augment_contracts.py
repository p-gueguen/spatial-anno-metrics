"""Tests for generic AnnData perturbations and behavioral metric-contract verification."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")


def _make_spatial_tissue(n: int = 480, seed: int = 0, sparse_matrix: bool = False):
    import scipy.sparse as sp
    from sklearn.decomposition import PCA

    rng = np.random.default_rng(seed)
    types = ["T_CD4", "T_CD8", "Fibro", "Epi"]
    genes = [f"G_{i}" for i in range(16)]
    marker_sets = {t: genes[i * 4 : (i + 1) * 4] for i, t in enumerate(types)}
    # Make T_CD4 and T_CD8 closer expression siblings while spatial domains touch along borders
    centers = np.full((4, 16), 0.2, dtype=np.float64)
    for i in range(4):
        centers[i, i * 4 : (i + 1) * 4] = 8.0
    centers[0, 4:6] = 1.5
    centers[1, 0:2] = 1.5

    per = n // 4
    n_cells = per * 4
    lab = np.repeat(types, per)
    samp = np.tile(["s0", "s1", "s2", "s3"], n_cells // 4)
    xy = np.zeros((n_cells, 2), dtype=np.float64)
    quad_centers = [(0.0, 0.0), (0.0, 2.2), (2.2, 0.0), (2.2, 2.2)]
    s_idx = np.array([int(s[1:]) for s in samp], dtype=np.float64)
    for i, (cx, cy) in enumerate(quad_centers):
        idx = np.flatnonzero(lab == types[i])
        xy[idx, 0] = rng.normal(cx, 0.9, size=len(idx))
        xy[idx, 1] = rng.normal(cy, 0.9, size=len(idx)) + s_idx[idx] * 20.0

    t_idx = np.repeat(np.arange(4), per)
    lam = centers[t_idx].copy()
    X = rng.poisson(lam).astype(np.float64)
    mat = sp.csr_matrix(X) if sparse_matrix else X
    obs = pd.DataFrame(
        {"cell_type": pd.Categorical(lab), "truth": pd.Categorical(lab), "sample": samp},
        index=[f"c{i}" for i in range(n_cells)],
    )
    a = anndata.AnnData(X=mat.copy(), obs=obs, var=pd.DataFrame(index=genes))
    a.layers["counts"] = mat.copy()
    a.obsm["spatial"] = xy
    a.uns["marker_sets"] = marker_sets
    tot = np.maximum(X.sum(axis=1, keepdims=True), 1.0)
    a.obsm["X_pca"] = PCA(n_components=10, random_state=0).fit_transform(np.log1p(X / tot * 1e4)).astype("float32")
    return a


def _make_ref(seed: int = 1):
    rng = np.random.default_rng(seed)
    types = ["T_CD4", "T_CD8", "Fibro", "Epi"]
    genes = [f"G_{i}" for i in range(16)]
    centers = np.full((4, 16), 0.2, dtype=np.float64)
    for i in range(4):
        centers[i, i * 4 : (i + 1) * 4] = 8.0
    centers[0, 4:6] = 1.5
    centers[1, 0:2] = 1.5
    t_idx = np.repeat(np.arange(4), 60)
    X = rng.poisson(centers[t_idx]).astype(np.float64)
    a = anndata.AnnData(
        X=X.copy(),
        obs=pd.DataFrame({"cell_type": pd.Categorical(np.repeat(types, 60))}, index=[f"r{i}" for i in range(240)]),
        var=pd.DataFrame(index=genes),
    )
    a.layers["counts"] = X.copy()
    return a


def test_perturb_anndata_all_operators_non_mutating_and_deterministic():
    import scipy.sparse as sp

    from spatial_anno_metrics import AUGMENTATIONS, perturb_anndata

    dense = _make_spatial_tissue(n=240, seed=0, sparse_matrix=False)
    sparse = _make_spatial_tissue(n=240, seed=0, sparse_matrix=True)
    orig_labels = dense.obs["cell_type"].astype(str).to_numpy().copy()

    for name, meta in AUGMENTATIONS.items():
        lv = meta["default_levels"][-1]
        p1 = perturb_anndata(dense, name, lv, seed=7, sample_key="sample", recompute_pca=True)
        p2 = perturb_anndata(dense, name, lv, seed=7, sample_key="sample", recompute_pca=True)
        assert p1.uns["sam_augmentation"]["perturbation"] == name
        assert (dense.obs["cell_type"].astype(str).to_numpy() == orig_labels).all()
        assert np.allclose(np.asarray(p1.layers["counts"]), np.asarray(p2.layers["counts"]))
        assert (p1.obs["cell_type"].astype(str).to_numpy() == p2.obs["cell_type"].astype(str).to_numpy()).all()

        ps = perturb_anndata(sparse, name, lv, seed=7, sample_key="sample")
        assert sp.issparse(ps.layers["counts"])
        assert np.allclose(np.asarray(p1.layers["counts"]), ps.layers["counts"].toarray())

    # infiltration preserves cell label-to-expression pairing while changing spatial coordinates
    inf = perturb_anndata(dense, "infiltration", 0.5, seed=3, sample_key="sample")
    assert (inf.obs["cell_type"].astype(str).to_numpy() == orig_labels).all()
    assert np.allclose(np.asarray(inf.layers["counts"]), np.asarray(dense.layers["counts"]))
    assert not np.allclose(inf.obsm["spatial"], dense.obsm["spatial"])

    # nuisance_shuffle preserves label counts while reordering rows cleanly
    ns = perturb_anndata(dense, "nuisance_shuffle", 1, seed=2, sample_key="sample")
    assert ns.obs["cell_type"].value_counts().to_dict() == dense.obs["cell_type"].value_counts().to_dict()
    assert not np.array_equal(ns.obs_names.to_numpy(), dense.obs_names.to_numpy())

    # panel_growth adds background genes
    pg = perturb_anndata(dense, "panel_growth", 40, seed=2)
    assert pg.n_vars == dense.n_vars + 40
    assert any(g.startswith("NegBg_") for g in pg.var_names)


def test_perturb_anndata_handles_float_expression():
    from spatial_anno_metrics import perturb_anndata

    a = _make_spatial_tissue(n=160, seed=4, sparse_matrix=False)
    a.layers["counts"] = (np.asarray(a.layers["counts"]) + 0.35).astype(np.float64)
    a.X = a.layers["counts"].copy()
    for op, lv in [
        ("spillover", 0.3),
        ("depth_loss", 0.5),
        ("ambient", 0.2),
        ("batch_shift", 0.5),
        ("coord_shuffle", 0.5),
        ("panel_growth", 10),
    ]:
        out = perturb_anndata(a, op, lv, seed=1, sample_key="sample")
        assert np.all(np.isfinite(out.layers["counts"]))
        assert float(out.layers["counts"].min()) >= 0.0


def test_verify_metric_contracts_positive_controls_catch_broken_metrics():
    from spatial_anno_metrics import verify_metric_contracts

    a = _make_spatial_tissue(n=200, seed=0)

    broken = {
        "blind": lambda x, r: {"rel_pw_ct": 0.75},
        "inverted": lambda x, r: {
            "ext.macro_f1": float((x.obs["cell_type"].astype(str) != x.obs["truth"].astype(str)).mean())
        },
        "order_bug": lambda x, r: {"pmp": 0.9 if x.obs_names[0] == "c0" else 0.1},
    }
    res = verify_metric_contracts(a, broken, seeds=(0, 1))
    assert res["status"] == "failed"
    verdicts = dict(zip(zip(res["summary"]["metric"], res["summary"]["perturbation"]), res["summary"]["verdict"]))
    assert verdicts[("rel_pw_ct", "label_noise")] == "fail_blind"
    assert verdicts[("rel_pw_ct", "type_swap")] == "fail_blind"
    assert verdicts[("ext.macro_f1", "label_noise")] == "fail_inverted"
    assert verdicts[("pmp", "nuisance_shuffle")] == "fail_false_alarm"


def test_verify_metric_contracts_passes_on_tissue_and_real_goldens():
    from pathlib import Path

    import spatial_anno_metrics as sam

    a = _make_spatial_tissue(n=480, seed=0)
    ref = _make_ref(seed=1)

    def _f(res_dict, key):
        val = res_dict.get(key)
        return np.nan if val is None else float(val)

    def _eval(x, r):
        out = {
            "ext.macro_f1": _f(sam.external_scores(y_true=x.obs["truth"], y_pred=x.obs["cell_type"]), "macro_f1"),
            "iv.integrated": _f(sam.internal_validity(x, "cell_type", embedding="X_pca"), "integrated"),
            "fid.mean_auc": _f(sam.marker_program_fidelity(x, "cell_type", x.uns["marker_sets"]), "mean_auc"),
            "pmp": _f(sam.pmp(x, x.uns["marker_sets"], "cell_type"), "mean_pmp"),
            "nmp": _f(sam.nmp(x, r, "cell_type", "cell_type"), "nmp"),
            "rel_pw_ct": _f(
                sam.relative_pairwise_ct_expression(x, r, label_key="cell_type"),
                "relative_pairwise_ct_expression",
            ),
            "isc.integrated": _f(sam.inter_sample_consistency(x, "cell_type", "sample"), "integrated"),
        }
        rm = sam.residual_spatial_autocorrelation(x, label_key="cell_type", n_neighs=6, graph=None)
        out["resid_moran.explained"] = _f(rm, "explained_fraction")
        out["resid_moran.residual"] = _f(rm, "residual_median_moran")
        org = sam.spatial_organization_consistency(x, label_key="cell_type", sample_key="sample", n_neighs=8)
        out["spatial_org.corr"] = _f(org, "mean_corr")
        out["spatial_org.offdiag"] = _f(org, "offdiag_mean_corr")
        return out

    res = sam.verify_metric_contracts(a, {"sam": _eval}, reference=ref, seeds=(0, 1), sample_key="sample")
    assert res["status"] == "ok", res["summary"][res["summary"]["verdict"].str.startswith("fail_")]
    assert res["n_passed"] >= 14

    # Real committed Xenium breast crop: verify residual_spatial_autocorrelation behavioral contract
    import sys

    bench_dir = str(Path(__file__).resolve().parents[1] / "benchmarks")
    if bench_dir not in sys.path:
        sys.path.insert(0, bench_dir)
    import real_perturbation as rp

    xenium = rp.load_xenium()
    res_xen = sam.verify_metric_contracts(
        xenium,
        {
            "resid": lambda x, r: {
                "resid_moran.explained": sam.residual_spatial_autocorrelation(x, label_key="label").get(
                    "explained_fraction", np.nan
                ),
                "resid_moran.residual": sam.residual_spatial_autocorrelation(x, label_key="label").get(
                    "residual_median_moran", np.nan
                ),
            }
        },
        seeds=(0,),
        label_key="label",
        perturbation_kwargs={"merge": {"pairs": [("Epithelial", "Fibroblast"), ("Tcell", "Myeloid")]}},
    )
    assert Path(__file__).name == "test_augment_contracts.py"
    assert res_xen["status"] == "ok", res_xen["summary"][res_xen["summary"]["verdict"].str.startswith("fail_")]


def test_must_not_degrade_abs_tol_and_roi_glob(tmp_path: Path) -> None:
    """Near-zero deterministic jitter respects abs_tol, while real degradation still fails."""
    import sys

    repo_root = str(Path(__file__).resolve().parents[1])
    if repo_root not in sys.path:
        sys.path.insert(0, repo_root)
    from benchmarks.augmentation.datasets import roi_inventory
    from spatial_anno_metrics.augment import _score_one

    levels = [1.0, 0.5, 0.25]
    # Jitter like resid_moran.residual on nuclear Slide-tags (0.0032 -> 0.0054) must pass with abs_tol=0.01
    s_ok = pd.Series([0.0032, 0.0041, 0.0054], index=levels)
    r_ok = _score_one(s_ok, levels, orientation=-1, role="must_not_degrade", min_rho=0.8, rel_tol=0.02, abs_tol=0.01, caveat=None)
    assert r_ok["verdict"] == "pass"
    assert r_ok["threshold"] == pytest.approx(0.01)

    # Real degradation (+0.05) still fails clearly
    s_bad = pd.Series([0.0032, 0.0200, 0.0532], index=levels)
    r_bad = _score_one(s_bad, levels, orientation=-1, role="must_not_degrade", min_rho=0.8, rel_tol=0.02, abs_tol=0.01, caveat=None)
    assert r_bad["verdict"] == "fail_false_alarm"

    # Rank 1 MERFISH human cortex glob must not double-count slidetags_human_cortex.h5ad
    (tmp_path / "merfish_human_cortex_2_rep1.h5ad").write_bytes(b"")
    (tmp_path / "slidetags_human_cortex.h5ad").write_bytes(b"")
    inv = roi_inventory(tmp_path).set_index("rank")
    assert inv.loc[1, "local_status"] == "cached (1 files)"
    assert inv.loc[2, "local_status"] == "cached (1 files)"
