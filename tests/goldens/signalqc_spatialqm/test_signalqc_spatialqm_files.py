"""Loader smoke test for the SpatialQM signal-QC parity golden: the files parse and agree with each
other. No metric is compared here; the owning work package adds that test (see SOURCE.md)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent


def _no_nan(token):
    raise ValueError(f"non-standard JSON constant {token}; NaN must be written as null")


def test_signalqc_spatialqm_golden_files_parse():
    c = np.load(HERE / "counts.npz", allow_pickle=False)
    n_cells, n_feat = (int(v) for v in c["shape"])
    assert c["indptr"][0] == 0 and len(c["indptr"]) == n_cells + 1
    assert c["indptr"][-1] == len(c["data"]) == len(c["indices"]) and c["indices"].max() < n_feat
    assert len(c["feature_id"]) == len(c["feature_name"]) == len(c["feature_type"]) == n_feat

    cells = pd.read_csv(HERE / "cells.csv")
    assert (cells["cell_id"].to_numpy() == c["cell_id"]).all()
    cs = np.concatenate([[0], np.cumsum(c["data"], dtype=np.int64)])
    rowsum = cs[c["indptr"][1:]] - cs[c["indptr"][:-1]]          # Xenium total_counts = all 541 features
    assert (rowsum == cells["total_counts"].to_numpy()).all()

    t = np.load(HERE / "transcripts.npz", allow_pickle=False)
    n_tx = len(t["feature_idx"])
    assert n_tx > 0 and all(len(t[k]) == n_tx for k in t.files)
    assert t["feature_idx"].max() < n_feat

    e = json.loads((HERE / "expected.json").read_text(), parse_constant=_no_nan)   # make_expected.R (SpatialQM)
    inp = e["meta"]["input"]                                   # the reference was computed on these files
    assert (inp["n_cells"], inp["n_features"], inp["n_transcripts"]) == (n_cells, n_feat, n_tx)
    assert inp["n_transcripts_unassigned"] == int((t["cell_id"] == -1).sum())
    m = e["getMorans"]
    assert len(m["feature"]) == len(m["type"]) == len(m["I"]) and set(m["feature"]) <= set(c["feature_name"])
