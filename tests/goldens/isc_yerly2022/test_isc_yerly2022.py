"""ISC parity golden: Yerly 2022 BCC subsample vs scTypeEval 1.0.0 (Bioconductor 3.23), see SOURCE.md.

Exact where both tools follow the same recipe (pseudobulk, cosine, 2-label silhouette, silhouette on the
reciprocal-classification distances). The classifier inside RCM differs by design (ours: Spearman nearest
profile; scTypeEval: SingleR), so that gap is pinned pair by pair on the over-split label column.
"""
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from spatial_anno_metrics import inter_sample_consistency, isc as eval_metrics

HERE = Path(__file__).parent
EXP = json.loads((HERE / "expected.json").read_text())
COLUMNS = ("annotation", "annotation_oversplit")
ISCL, ISCG, ISCG0 = EXP["names"]["ISCL"], EXP["names"]["ISCG"], EXP["names"]["ISCG_no_pca"]
SIL0 = "silhouette | Pseudobulk:Cosine [reduction=FALSE]"  # our legacy `per_type`
HALVES = ["Fibroblast", "Fibroblast.Nspl2", "Tcell", "Tcell.Nspl2"]


@pytest.fixture(scope="module")
def adata():
    z = np.load(HERE / "counts.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    a = ad.AnnData(X=X.astype(np.float64), obs=pd.read_csv(HERE / "obs.csv", dtype=str).set_index("cell_id"))
    a.var_names = z["genes"].astype(str)
    a.layers["counts"] = a.X
    return a


def reference_rcm(adata, col):
    """scTypeEval's RCM distances in our profile order (0.5 read as 1, as get_consistency does), and ours."""
    P, lab, smp, _, _ = eval_metrics._pseudobulk_profiles(adata, col, "patient_bcc", min_cells=10, min_samples=5,
                                                       layer="counts")
    ids = [f"{s}_{re.sub(r'[ _+-]', '.', t)}" for t, s in zip(lab, smp)]  # scTypeEval's profile names
    D = EXP["columns"][col]["dissimilarity"]["recip_classif:Match"]
    order = [D["labels"].index(i) for i in ids]
    ref = np.array(D["matrix"], dtype=float)[np.ix_(order, order)]
    return ref, eval_metrics._rcm_dissimilarity(P, smp, lab), ids


def test_inputs_parse():
    z = np.load(HERE / "counts.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.read_csv(HERE / "obs.csv", dtype=str)
    assert X.shape == (2620, 2000) == (len(obs), len(z["genes"]))
    assert list(obs.columns) == ["cell_id", "patient_bcc", *COLUMNS] and obs["cell_id"].is_unique
    assert len(set(z["genes"])) == 2000 and X.data.min() > 0
    per = pd.crosstab(obs["annotation"], obs["patient_bcc"])
    assert per.shape == (9, 7)
    assert ((per >= 10).sum(axis=1) >= 5).sum() == 7  # Bcell (4 samples) and PlasmaCell (2) fall below
    split = obs["annotation_oversplit"] != obs["annotation"]
    assert set(obs.loc[split, "annotation"]) == {"Fibroblast", "Tcell"}
    assert (pd.crosstab(obs["annotation_oversplit"], obs["patient_bcc"]).loc[HALVES] == 25).all().all()


@pytest.mark.parametrize("col", COLUMNS)
def test_iscg_and_legacy_silhouette_match_sctypeeval(adata, col):
    res = inter_sample_consistency(adata, col, "patient_bcc", layer="counts", reduction=None)
    ref = EXP["columns"][col]["per_type"]
    assert set(res["dropped_types"]) == {t for t, v in ref.items() if v[ISCG0] is None} == {"Bcell", "PlasmaCell"}
    assert set(res["per_type_iscg"]) == {t for t, v in ref.items() if v[ISCG0] is not None}
    for t, g in res["per_type_iscg"].items():
        assert g == pytest.approx(ref[t][ISCG0], abs=1e-12)
        assert res["per_type"][t] == pytest.approx(ref[t][SIL0], abs=1e-12)


@pytest.mark.parametrize("col", COLUMNS)
def test_iscl_matches_sctypeeval_on_its_rcm_distances(adata, col, monkeypatch):
    ref, _, _ = reference_rcm(adata, col)
    monkeypatch.setattr(eval_metrics, "_rcm_dissimilarity", lambda P, samples, labels, genes=None: ref)
    res = inter_sample_consistency(adata, col, "patient_bcc", layer="counts")
    for t, v in res["per_type_iscl"].items():
        assert v == pytest.approx(EXP["columns"][col]["per_type"][t][ISCL], abs=1e-12)


def test_rcm_classifier_gap_vs_singler(adata):
    ref, ours, _ = reference_rcm(adata, "annotation")
    assert np.array_equal(ours, ref)  # same mutual hits on all 1128 profile pairs
    ref, ours, ids = reference_rcm(adata, "annotation_oversplit")
    iu = np.triu_indices(len(ids), 1)
    diff = [(ids[i], ids[j]) for i, j in zip(*iu) if ours[i, j] != ref[i, j]]
    assert len(iu[0]) == 1891 and len(diff) == 30  # SOURCE.md lists the 30 pairs

    def label(pid):
        return pid.split("_", 1)[1].removesuffix(".Nspl2")

    assert all(label(i) == label(j) in ("Fibroblast", "Tcell") for i, j in diff)  # halves of one label only


@pytest.mark.parametrize("col", COLUMNS)
def test_iscg_default_30_pcs_matches_sctypeeval(adata, col):
    """The default ISCG space is scTypeEval's (README recipe, reduction = TRUE): cosine on 30 centred PCs. Ours
    uses an exact SVD, scTypeEval irlba; SOURCE.md measured the gap at 1.4e-6 / 8.5e-7."""
    res = inter_sample_consistency(adata, col, "patient_bcc", layer="counts")
    ref = EXP["columns"][col]["per_type"]
    assert res["reduction"] == "pca" and res["n_pcs"] == 30
    for t, g in res["per_type_iscg"].items():
        assert g == pytest.approx(ref[t][ISCG], abs=2e-6)
        assert abs(ref[t][ISCG] - ref[t][ISCG0]) > 0.01  # positive control: the two spaces differ by far more
    assert res["iscg"] == pytest.approx(np.mean([ref[t][ISCG] for t in res["per_type_iscg"]]), abs=2e-6)


@pytest.mark.parametrize("col", COLUMNS)
def test_headline_is_the_mean_over_labels(adata, col):
    """iscg / iscl / integrated average over labels (paper, scTypeEval mean_over_types); the pseudobulk-weighted
    means stay under *_profile_mean (Mast has 6 profiles, every other label 7)."""
    res = inter_sample_consistency(adata, col, "patient_bcc", layer="counts", reduction=None)
    mean = EXP["columns"][col]["mean_over_types"]
    assert res["iscg"] == pytest.approx(mean[ISCG0], abs=1e-12)
    assert res["iscg_profile_mean"] == pytest.approx({"annotation": 0.855469, "annotation_oversplit": 0.825611}[col],
                                                     abs=1e-6)  # the pre-0.7.0 headline (SOURCE.md)
    assert res["iscg"] != pytest.approx(res["iscg_profile_mean"], abs=1e-4)
    assert res["iscl"] == pytest.approx(np.mean(list(res["per_type_iscl"].values())), abs=1e-15)


def test_rcm_markers_vs_singler(adata):
    """rcm_genes='markers' (SingleR classic marker union per training sample, no fine-tuning) vs SingleR's matrix:
    still identical on the authors' labels; on the over-split column the count of differing pairs is pinned."""
    for col, n_diff in (("annotation", 0), ("annotation_oversplit", 28)):
        ref, _, _ = reference_rcm(adata, col)
        P, lab, smp, _, _ = eval_metrics._pseudobulk_profiles(adata, col, "patient_bcc", min_cells=10, min_samples=5,
                                                              layer="counts")
        ours = eval_metrics._rcm_dissimilarity(P, smp, lab, genes="markers")
        assert int(np.triu(ours != ref, 1).sum()) == n_diff, col
