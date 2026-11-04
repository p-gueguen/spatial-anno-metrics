"""Tests for consensus.py - multi-annotator consensus with ontology-path voting (popV-style) and its entropy."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

anndata = pytest.importorskip("anndata")

H = {"T_CD4": "T", "T_CD8": "T", "Treg": "T_CD4", "B": "Lymphoid", "T": "Lymphoid"}


def _adata(rows):
    """One cell per row; columns a0, a1, ... are the annotators."""
    df = pd.DataFrame(rows, columns=[f"a{i}" for i in range(len(rows[0]))], dtype=object)
    df.index = [f"c{i}" for i in range(len(df))]
    return anndata.AnnData(X=np.zeros((len(df), 1), dtype="float32"), obs=df)


def test_consensus_entropy_hand_cases():
    from spatial_anno_metrics.consensus import consensus_entropy

    a = _adata([["B", "B", "B"],                 # full agreement
                ["T_CD4", "T_CD8", "T"],         # subtype split under one lineage
                ["T_CD4", "T_CD4", "T_CD8"],     # majority subtype
                ["Treg", "T_CD4", "T_CD8"],      # Treg votes for T_CD4 too -> T_CD4 has 2/3
                ["B", "Myeloid", "Epi"],         # three roots: no label has a majority
                ["B", "T_CD4", "Myeloid"]])      # B and T_CD4 both vote Lymphoid: 2/3 at the root
    r = consensus_entropy(a, ["a0", "a1", "a2"], hierarchy=H)
    assert r["status"] == "ok"
    assert list(r["consensus"]) == ["B", "T", "T_CD4", "T_CD4", None, "Lymphoid"]
    assert r["consensus_share"] == pytest.approx([1, 1, 2 / 3, 2 / 3, 1 / 3, 2 / 3])
    h23 = -(2 / 3 * np.log(2 / 3) + 1 / 3 * np.log(1 / 3)) / np.log(3)
    # entropy of the votes projected to the consensus depth: cell 1 -> T, T, T; cells 2-3 at depth 2 (Treg ->
    # T_CD4); cell 4 (no consensus) at the roots; cell 5 at the roots: Lymphoid, Lymphoid, Myeloid
    assert r["ontology_entropy"] == pytest.approx([0, 0, h23, h23, 1, h23])
    assert r["entropy"] == pytest.approx([0, 1, h23, 1, 1, 1])             # flat: subtype splits count fully
    s = r["summary"]
    assert s["n_annotators"] == 3 and s["pct_no_consensus"] == pytest.approx(1 / 6)
    assert s["mean_ontology_entropy"] < s["mean_entropy"]


def test_consensus_without_hierarchy_is_flat_majority():
    from spatial_anno_metrics.consensus import consensus_entropy

    a = _adata([["T_CD4", "T_CD8", "T"], ["B", "B", "T"]])
    r = consensus_entropy(a, ["a0", "a1", "a2"])
    assert list(r["consensus"]) == [None, "B"]
    assert r["ontology_entropy"] == pytest.approx(r["entropy"])
    assert r["entropy"][0] == pytest.approx(1.0)


def test_consensus_two_annotator_tie_has_no_consensus_without_a_shared_ancestor():
    """Strict majority: 1 of 2 is not a consensus; with a hierarchy the shared parent is."""
    from spatial_anno_metrics.consensus import consensus_entropy

    a = _adata([["T_CD4", "T_CD8"]])
    assert consensus_entropy(a, ["a0", "a1"])["consensus"][0] is None
    assert consensus_entropy(a, ["a0", "a1"], hierarchy=H)["consensus"][0] == "T"


def test_consensus_input_checks():
    from spatial_anno_metrics.consensus import consensus_entropy

    a = _adata([["B", None], ["B", "T"]])
    with pytest.raises(ValueError, match="missing"):
        consensus_entropy(a, ["a0", "a1"])
    assert consensus_entropy(a, ["a0"])["status"] == "skipped"
    with pytest.raises(ValueError, match="cycle"):
        consensus_entropy(_adata([["X", "Y"]]), ["a0", "a1"], hierarchy={"X": "Y", "Y": "X"})
    with pytest.raises(ValueError, match="min_share"):
        consensus_entropy(_adata([["X", "Y"]]), ["a0", "a1"], min_share=0.4)


def test_consensus_entropy_on_yerly2022_bcc_real_annotations():
    """Real multi-annotator table on Yerly et al. 2022 human BCC (2,620 cells): author annotation, the
    benchmark's over-split annotation (Fibroblast.Nspl2, Tcell.Nspl2), and a coarse lineage annotator.
    Flat majority leaves every over-split cell without consensus (3 distinct strings); ontology-path voting
    resolves them to their parent type (Tcell, Fibroblast) with zero ontology entropy at the parent depth."""
    from pathlib import Path

    from spatial_anno_metrics.consensus import consensus_entropy

    obs = pd.read_csv(Path(__file__).parent / "goldens" / "isc_yerly2022" / "obs.csv")
    obs.index = obs["cell_id"].astype(str)
    coarse_map = {
        "Tcell": "Lymphoid", "Bcell": "Lymphoid", "PlasmaCell": "Lymphoid",
        "MoMacDC": "Myeloid", "Mast": "Myeloid",
        "Fibroblast": "Stromal", "Endothelial": "Stromal",
        "Normal_Kerat": "Epithelial", "Melanocyte": "Melanocytic",
    }
    obs["annotation_coarse"] = obs["annotation"].map(coarse_map)
    hierarchy = {
        "Tcell.Nspl2": "Tcell",
        "Fibroblast.Nspl2": "Fibroblast",
        **coarse_map,
    }
    a = anndata.AnnData(X=np.zeros((len(obs), 1), dtype="float32"), obs=obs)
    cols = ["annotation", "annotation_oversplit", "annotation_coarse"]

    flat = consensus_entropy(a, cols)
    hier = consensus_entropy(a, cols, hierarchy=hierarchy)
    assert flat["status"] == "ok" and hier["status"] == "ok"

    # Under flat voting, every cell where annotation_oversplit != annotation has 3 distinct strings -> no consensus
    split_mask = (obs["annotation"] != obs["annotation_oversplit"]).to_numpy()
    assert split_mask.sum() == 350                               # 7 patients x (25 Tcell + 25 Fibroblast)
    assert flat["summary"]["pct_no_consensus"] == pytest.approx(split_mask.mean())
    assert all(c is None for c, m in zip(flat["consensus"], split_mask) if m)

    # With the hierarchy, Tcell + Tcell.Nspl2 both vote for Tcell (2/3 majority) and all 3 vote for the root
    h23 = -(2 / 3 * np.log(2 / 3) + 1 / 3 * np.log(1 / 3)) / np.log(3)
    assert hier["summary"]["pct_no_consensus"] == 0.0
    assert hier["summary"]["mean_ontology_entropy"] < hier["summary"]["mean_entropy"]
    assert hier["summary"]["mean_ontology_entropy"] == pytest.approx(h23)
    for i in np.flatnonzero(split_mask):
        assert hier["consensus"][i] == obs["annotation"].iloc[i]
        assert hier["consensus_share"][i] == pytest.approx(2 / 3)
        assert flat["entropy"][i] == pytest.approx(1.0)
        assert hier["ontology_entropy"][i] == pytest.approx(h23)
    for i in np.flatnonzero(~split_mask):
        assert hier["consensus"][i] == obs["annotation"].iloc[i]
        assert hier["ontology_entropy"][i] == pytest.approx(h23)
