"""Parity golden for co-expression: txsim ``coexpression_similarity`` and the Moldia negative co-expression
purity (NCP, Salas et al. 2025) vs our ``ncp``.

Writes ``tests/goldens/coexpr_txsim/`` (the same inputs as nmp_txsim + ``expected.json``) and prints the
ours-vs-reference table for its SOURCE.md. Not run in CI; same venv and clones as nmp_txsim, whose helpers
this script reuses. Run from anywhere:

    MPLBACKEND=Agg <venv>/bin/python benchmarks/parity/coexpr_txsim/make_golden.py [--upstream DIR]
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from datetime import date
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("nmp_make_golden",
                                               HERE.parent / "nmp_txsim" / "make_golden.py")
nmp = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nmp)

OUT = nmp.REPO / "tests" / "goldens" / "coexpr_txsim"
TARGET_SUM = 1e4                 # lognorm = log1p(TARGET_SUM * counts / per-cell total), float64


def lognorm(a):
    """``layers['lognorm']`` = log1p(1e4 x counts / total); X too, as OpenProblems does (txsim's gene
    filter reads X; its detection pattern is the same as the counts')."""
    tot = a.layers["counts"].sum(1, keepdims=True)
    assert (tot > 0).all()
    a.layers["lognorm"] = np.log1p(TARGET_SUM * a.layers["counts"] / tot)
    a.X = a.layers["lognorm"].copy()
    return a


def main() -> None:
    # txsim's by_celltype loops over a set of cell types: its order, and so the last bit of the averaged
    # per-gene values, follows Python's string hash seed. Pin it so a re-run gives identical bytes.
    if os.environ.get("PYTHONHASHSEED") != "0":
        os.environ["PYTHONHASHSEED"] = "0"
        os.execv(sys.executable, [sys.executable, *sys.argv])
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", type=Path, default=nmp.UPSTREAM,
                    help="dir with the txsim and Moldia clones")
    args = ap.parse_args()
    nmp.setup_paths()
    txsim_dir = nmp.pinned(args.upstream / "txsim", "txsim")
    moldia_dir = nmp.pinned(args.upstream / "Xenium_benchmarking", "Xenium_benchmarking")
    sys.path.insert(0, str(txsim_dir))
    import matplotlib.pyplot as plt
    import txsim as tx
    assert Path(tx.__file__).resolve().is_relative_to(txsim_dir.resolve()), tx.__file__
    xb = nmp.load_module(moldia_dir / "xb" / "calculating.py", "moldia_xb_calculating")
    import spatial_anno_metrics as sam
    assert Path(sam.__file__).resolve().is_relative_to(nmp.REPO / "src"), sam.__file__

    nmp.export_inputs(OUT)
    ref = lognorm(nmp.to_adata(OUT / "reference.npz"))
    kw = {"key": "celltype", "pipeline_output": False}
    results, ours = {}, {}
    for s in nmp.SECTIONS:
        sp, sc = nmp.shared_genes(lognorm(nmp.to_adata(OUT / f"section_{s}.npz")), ref)
        r = {"section": f"section_{s}.npz", "reference": "reference.npz"}
        v, per_gene, _ = tx.metrics.coexpression_similarity(sp.copy(), sc.copy(), layer="lognorm", **kw)
        r["txsim_coexpression_similarity"] = {"overall": nmp.num(v), "per_gene": nmp.series(per_gene)}
        v, per_gene, per_ct = tx.metrics.coexpression_similarity(sp.copy(), sc.copy(), layer="lognorm",
                                                                 by_celltype=True, **kw)
        r["txsim_coexpression_similarity_by_celltype"] = {
            "overall": nmp.num(v), "per_celltype": nmp.series(per_ct.sort_index()),
            "per_gene": nmp.series(per_gene)}
        v, by_col, by_row, low_sc, _, _ = xb.negative_marker_purity_coexpression(
            sp.copy(), sc.copy(), minexp=0.0, **kw)
        plt.close("all")                                                # it draws a clustermap
        r["moldia_ncp"] = {"overall": nmp.num(v), "n_negative_pairs": int(len(low_sc)),
                           "per_gene_as_conditioning": nmp.series(by_col),
                           "per_gene_as_target": nmp.series(by_row)}
        results[s] = r
        a = nmp.to_adata(OUT / f"section_{s}.npz", key="cell_type")
        ours[s] = sam.ncp(a, nmp.to_adata(OUT / "reference.npz", key="cell_type"))

    nmp.write_json(OUT / "expected.json", {
        "generated": date.today().isoformat(),
        "upstream": {k: f"{k}@{v}" for k, v in nmp.PINS.items()},
        "lognorm_target_sum": TARGET_SUM,
        "results": results})

    print("| section | txsim coexpression_similarity | txsim by cell type | Moldia NCP (pairs) | ours ncp |")
    print("|---|---|---|---|---|")
    for s in nmp.SECTIONS:
        r, m = results[s], results[s]["moldia_ncp"]
        print(f"| {s} | {r['txsim_coexpression_similarity']['overall']:.6f} | "
              f"{r['txsim_coexpression_similarity_by_celltype']['overall']:.6f} | "
              f"{m['overall']:.6f} ({m['n_negative_pairs']}) | {ours[s].get('ncp', ours[s])} |")
    for s in nmp.SECTIONS:
        by_ct = results[s]["txsim_coexpression_similarity_by_celltype"]["per_celltype"]
        print(s, "by cell type:", {k: round(v, 6) for k, v in by_ct.items()},
              "| genes kept:", len(results[s]["txsim_coexpression_similarity"]["per_gene"]))


if __name__ == "__main__":
    main()
