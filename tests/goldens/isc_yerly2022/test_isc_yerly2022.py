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
ISCL, ISCG0 = EXP["names"]["ISCL"], EXP["names"]["ISCG_no_pca"]
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
    ref[ref == 0.5] = 1.0
    return ref, eval_metrics._rcm_dissimilarity(P, smp), ids


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
    res = inter_sample_consistency(adata, col, "patient_bcc", layer="counts")
    ref = EXP["columns"][col]["per_type"]
    assert set(res["dropped_types"]) == {t for t, v in ref.items() if v[ISCG0] is None} == {"Bcell", "PlasmaCell"}
    assert set(res["per_type_iscg"]) == {t for t, v in ref.items() if v[ISCG0] is not None}
    for t, g in res["per_type_iscg"].items():
        assert g == pytest.approx(ref[t][ISCG0], abs=1e-12)
        assert res["per_type"][t] == pytest.approx(ref[t][SIL0], abs=1e-12)


@pytest.mark.parametrize("col", COLUMNS)
def test_iscl_matches_sctypeeval_on_its_rcm_distances(adata, col, monkeypatch):
    ref, _, _ = reference_rcm(adata, col)
    monkeypatch.setattr(eval_metrics, "_rcm_dissimilarity", lambda P, samples: ref)
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
