"""Loader test for the SpatialQM signal-QC parity golden: the files parse, agree with each other, and are the
files expected.json was computed from. make_inputs.py `outs` records the sha256 of the three fixtures it reads and
make_expected.R copies them into meta.input.sha256; expected.json's own sha256 is pinned here (regenerating it
means updating EXPECTED_SHA256, SOURCE.md says how). Metric values: test_signalqc_spatialqm_metrics.py."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
EXPECTED_SHA256 = "01a7b5aa1f62d50797c310816704ac16a7f2c5e623f4abd22872baa8a4bc07aa"
FIXTURES = ("counts.npz", "cells.csv", "transcripts.npz")


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


def _sha256(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def test_signalqc_spatialqm_golden_ties_expected_to_inputs():
    e = json.loads((HERE / "expected.json").read_text())
    assert e["meta"]["input"]["sha256"] == {f: _sha256(f) for f in FIXTURES}   # computed from these bytes
    assert _sha256("expected.json") == EXPECTED_SHA256                         # and not edited since
