"""Smoke test: the reference parity goldens parse. The work packages that own each metric compare values."""
import json
from pathlib import Path

import numpy as np
import pytest

GOLDENS = Path(__file__).parent / "goldens"


def _no_nan(token):
    raise ValueError(f"{token} in a golden: NaN must be stored as null")


@pytest.mark.parametrize("name", ["nmp_txsim", "coexpr_txsim", "ecs_clusim", "deconv_openproblems"])
def test_golden_files_parse(name):
    d = GOLDENS / name
    assert (d / "SOURCE.md").is_file()
    assert "results" in json.loads((d / "expected.json").read_text(), parse_constant=_no_nan)
    for f in d.iterdir():
        assert f.stat().st_size < 5_000_000, f
        if f.suffix == ".json":
            json.loads(f.read_text(), parse_constant=_no_nan)
        elif f.suffix == ".npz":
            with np.load(f, allow_pickle=False) as z:
                assert all(z[k].size for k in z.files), f
