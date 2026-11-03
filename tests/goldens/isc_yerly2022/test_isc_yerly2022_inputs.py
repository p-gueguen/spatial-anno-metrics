"""Loader smoke test for the ISC parity golden: the input files parse and match SOURCE.md.

Inputs only - the scTypeEval reference (expected.json) is not generated yet, see SOURCE.md.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).parent


def test_isc_yerly2022_inputs_parse():
    z = np.load(HERE / "counts.npz")
    X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"]))
    obs = pd.read_csv(HERE / "obs.csv", dtype=str)
    assert X.shape == (2620, 2000) == (len(obs), len(z["genes"]))
    assert list(obs.columns) == ["cell_id", "patient_bcc", "annotation"] and obs["cell_id"].is_unique
    assert len(set(z["genes"])) == 2000 and X.data.min() > 0
    per = pd.crosstab(obs["annotation"], obs["patient_bcc"])
    assert per.shape == (9, 7)
    assert ((per >= 10).sum(axis=1) >= 5).sum() == 7  # Bcell (4 samples) and PlasmaCell (2) fall below
