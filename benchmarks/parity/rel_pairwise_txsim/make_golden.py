"""Parity golden for relative pairwise cell-type expression similarity: txsim vs our ``relative_pairwise_ct_expression``.

Writes ``tests/goldens/rel_pairwise_txsim/`` (the same inputs as nmp_txsim + ``expected.json``) and prints
the ours-vs-reference table for its SOURCE.md. Not run in CI. Run from anywhere:

    MPLBACKEND=Agg <python_with_scanpy> benchmarks/parity/rel_pairwise_txsim/make_golden.py [--upstream DIR]
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("nmp_make_golden",
                                               HERE.parent / "nmp_txsim" / "make_golden.py")
nmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nmp)

OUT = nmp.REPO / "tests" / "goldens" / "rel_pairwise_txsim"
TARGET_SUM = 1e4


def lognorm(a):
    """``layers['lognorm']`` = log1p(1e4 * counts / total); X too, as OpenProblems does."""
    tot = a.layers["counts"].sum(1, keepdims=True)
    assert (tot > 0).all()
    a.layers["lognorm"] = np.log1p(TARGET_SUM * a.layers["counts"] / tot)
    a.X = a.layers["lognorm"].copy()
    return a


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", type=Path, default=Path("/tmp/op_audit"),
                    help="dir containing the txsim clone")
    args = ap.parse_args()
    nmp.setup_paths()
    txsim_dir = args.upstream / "txsim"
    txsim_sha = subprocess.run(
        ["git", "-C", str(txsim_dir), "rev-parse", "HEAD"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()
    mod = nmp.load_module(
        txsim_dir / "txsim" / "metrics" / "_relative_pairwise_celltype_expression.py",
        "txsim_rel_pw_ct",
    )
    import spatial_anno_metrics as sam
    assert Path(sam.__file__).resolve().is_relative_to(nmp.REPO / "src"), sam.__file__

    nmp.export_inputs(OUT)
    ref = lognorm(nmp.to_adata(OUT / "reference.npz"))
    results, ours = {}, {}
    for s in nmp.SECTIONS:
        sp, sc = nmp.shared_genes(lognorm(nmp.to_adata(OUT / f"section_{s}.npz")), ref)
        overall, per_gene, per_ct = mod.relative_pairwise_celltype_expression(
            sp.copy(), sc.copy(), key="celltype", layer="lognorm", pipeline_output=False
        )
        results[s] = {
            "section": f"section_{s}.npz",
            "reference": "reference.npz",
            "txsim_relative_pairwise_ct_expr": {
                "overall": nmp.num(overall),
                "per_celltype": nmp.series(per_ct.sort_index()),
                "per_gene": nmp.series(per_gene.sort_index()),
            },
        }
        a = nmp.to_adata(OUT / f"section_{s}.npz", key="cell_type")
        ref_a = nmp.to_adata(OUT / "reference.npz", key="cell_type")
        ours[s] = sam.relative_pairwise_ct_expression(a, ref_a)

    nmp.write_json(OUT / "expected.json", {
        "generated": date.today().isoformat(),
        "upstream": {"txsim": f"txsim@{txsim_sha}"},
        "lognorm_target_sum": TARGET_SUM,
        "results": results,
    })

    print("| section | txsim rel_pairwise_ct | ours relative_pairwise_ct_expression | max |diff| |")
    print("|---|---|---|---|")
    for s in nmp.SECTIONS:
        r = results[s]["txsim_relative_pairwise_ct_expr"]
        o = ours[s]
        d_ov = abs(r["overall"] - o["relative_pairwise_ct_expression"])
        d_ct = max(abs(r["per_celltype"][k] - o["per_type"][k]) for k in r["per_celltype"])
        d_g = max(abs(r["per_gene"][k] - o["per_gene"][k]) for k in r["per_gene"])
        print(f"| {s} | {r['overall']:.12f} | {o['relative_pairwise_ct_expression']:.12f} | "
              f"{max(d_ov, d_ct, d_g):.3e} |")


if __name__ == "__main__":
    main()
