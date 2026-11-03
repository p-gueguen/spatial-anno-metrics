"""Parity golden for deconvolution metrics: OpenProblems task_spatial_decomposition ``r2`` and ``jsd`` vs our
``deconvolution_metrics``.

Writes ``tests/goldens/deconv_openproblems/`` (``proportions.json`` + ``expected.json``) and prints the
ours-vs-reference table for its SOURCE.md. Not run in CI. Runs the two upstream metric scripts unchanged
except for their VIASH parameter block, so it needs a clone at the pinned commit and a venv with anndata,
scikit-learn and scipy (see SOURCE.md). Run from anywhere:

    <venv>/bin/python benchmarks/parity/deconv_openproblems/make_golden.py [--upstream DIR]
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "tests" / "goldens" / "deconv_openproblems"
UPSTREAM = Path("/srv/GT/analysis/pgueguen/sam_dev/runs/wf2-prep-ref/upstream/task_spatial_decomposition")
PIN = "892ab4818b3d28656cd1a05f71b34c4d31af1de2"
N_SPOTS, N_TYPES = 40, 5


def rows(m):
    return m / m.sum(1, keepdims=True)


def make_pairs() -> dict:
    """close: small Dirichlet error; noisy: an unrelated prediction; absent_type: cell type 4 is absent from
    every true spot (an all-zero column) while the prediction gives it a little mass."""
    rng = np.random.default_rng(0)
    true = rng.dirichlet(np.ones(N_TYPES), N_SPOTS)
    absent = rows(np.column_stack([rng.dirichlet(np.ones(N_TYPES - 1), N_SPOTS), np.zeros(N_SPOTS)]))
    pairs = {"close": (true, rows(0.9 * true + 0.1 * rng.dirichlet(np.ones(N_TYPES), N_SPOTS))),
             "noisy": (true, rng.dirichlet(np.ones(N_TYPES), N_SPOTS)),
             "absent_type": (absent, rows(0.95 * absent + 0.05 * rng.dirichlet(np.ones(N_TYPES), N_SPOTS)))}
    return {k: {"true": t.tolist(), "pred": p.tolist()} for k, (t, p) in pairs.items()}


def run_metric(script: Path, true: np.ndarray, pred: np.ndarray, tmp: Path) -> float:
    """Execute an upstream metric script with its VIASH block replaced by our file paths; return its value."""
    import anndata as ad
    sol, meth, out = tmp / "solution.h5ad", tmp / "output.h5ad", tmp / f"score_{script.parent.name}.h5ad"
    ad.AnnData(np.zeros(true.shape), obsm={"proportions_true": true}).write_h5ad(sol)
    ad.AnnData(np.zeros(pred.shape), obsm={"proportions_pred": pred},
               uns={"dataset_id": "synthetic", "method_id": "synthetic"}).write_h5ad(meth)
    src = script.read_text()
    par = (f"par = {{'input_method': {str(meth)!r}, 'input_solution': {str(sol)!r}, "
           f"'output': {str(out)!r}}}\nmeta = {{'name': {script.parent.name!r}}}\n")
    src, n = re.subn(r"## VIASH START\n.*?## VIASH END\n", par, src, flags=re.S)
    assert n == 1, script
    exec(compile(src, str(script), "exec"), {"__name__": "__main__"})
    score = ad.read_h5ad(out)
    assert list(score.uns["metric_ids"]) == [script.parent.name]
    return float(score.uns["metric_values"][0])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", type=Path, default=UPSTREAM, help="task_spatial_decomposition clone")
    args = ap.parse_args()
    head = subprocess.run(["git", "-C", str(args.upstream), "rev-parse", "HEAD"], capture_output=True,
                          text=True, check=True).stdout.strip()
    if head != PIN:
        sys.exit(f"{args.upstream} is at {head}, expected {PIN}")
    sys.path.insert(0, str(REPO / "src"))
    import spatial_anno_metrics as sam
    assert Path(sam.__file__).resolve().is_relative_to(REPO / "src"), sam.__file__

    OUT.mkdir(parents=True, exist_ok=True)
    pairs = make_pairs()
    (OUT / "proportions.json").write_text(json.dumps(pairs) + "\n")
    results, ours = {}, {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, p in pairs.items():
            t, q = np.array(p["true"]), np.array(p["pred"])
            vals = {m: run_metric(args.upstream / "src" / "metrics" / m / "script.py", t, q, Path(tmp))
                    for m in ("r2", "jsd")}
            results[name] = {m: (None if np.isnan(v) else v) for m, v in vals.items()}
            ours[name] = sam.deconvolution_metrics(t, q)
    (OUT / "expected.json").write_text(json.dumps({
        "generated": date.today().isoformat(),
        "upstream": f"openproblems-bio/task_spatial_decomposition@{PIN}",
        "call": "src/metrics/{r2,jsd}/script.py, unchanged but for the VIASH par block; true -> "
                "solution.obsm['proportions_true'], pred -> method output obsm['proportions_pred']",
        "results": results}, indent=1, allow_nan=False) + "\n")

    print("| pair | OpenProblems r2 | ours r2 | OpenProblems jsd | ours jsd |")
    print("|---|---|---|---|---|")
    for name in pairs:
        r, o = results[name], ours[name]
        print(f"| {name} | {r['r2']!r} | {o['r2']!r} | {r['jsd']!r} | {o['jsd']!r} |")


if __name__ == "__main__":
    main()
