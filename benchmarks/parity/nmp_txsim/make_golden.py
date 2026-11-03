"""Parity golden for negative marker purity (Salas et al. 2025): txsim and the Moldia code vs our ``nmp``.

Writes ``tests/goldens/nmp_txsim/`` (inputs + ``expected.json``) and prints the ours-vs-reference tables
that go into its SOURCE.md. Not run in CI. Needs a venv with txsim's import-time dependencies and clones of
the two upstream repositories at the pinned commits (both are only run, never copied: neither has a licence).
See ``tests/goldens/nmp_txsim/SOURCE.md`` for the venv recipe. Run from anywhere:

    MPLBACKEND=Agg <venv>/bin/python benchmarks/parity/nmp_txsim/make_golden.py [--upstream DIR]

``benchmarks/parity/coexpr_txsim/make_golden.py`` reuses the helpers below (input export included).
"""
from __future__ import annotations

import argparse
import io
import json
import subprocess
import sys
import zipfile
from datetime import date
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "tests" / "goldens" / "nmp_txsim"
UPSTREAM = Path("/srv/GT/analysis/pgueguen/sam_dev/runs/wf2-prep-ref/upstream")
PINS = {"txsim": "bfdfaa9094c1f465d79b8b4ade855de8f2d2fb78",
        "Xenium_benchmarking": "50249f7d03d9cfc1453c4b81ef3b1053c4d23052"}
SECTIONS = ("clean", "spillover")
N_PERM = 20                      # Salas et al.: baseline = mean NMP over 20 random label permutations


def setup_paths() -> None:
    """This worktree's package and simulator first on sys.path (never an installed copy)."""
    for p in (REPO / "benchmarks" / "synthetic", REPO / "src"):
        sys.path.insert(0, str(p))


def pinned(clone: Path, name: str) -> Path:
    head = subprocess.run(["git", "-C", str(clone), "rev-parse", "HEAD"], capture_output=True, text=True,
                          check=True).stdout.strip()
    if head != PINS[name]:
        sys.exit(f"{clone} is at {head}, expected {PINS[name]}")
    return clone


def load_module(path: Path, name: str):
    """Import one upstream file by path (keeps the rest of the clone off sys.path)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def save_npz(path: Path, **arrays) -> None:
    """``np.savez_compressed`` with a fixed zip timestamp: identical arrays give identical bytes on re-run."""
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for key, arr in arrays.items():
            buf = io.BytesIO()
            np.lib.format.write_array(buf, np.asanyarray(arr), allow_pickle=False)
            z.writestr(zipfile.ZipInfo(f"{key}.npy", date_time=(1980, 1, 1, 0, 0, 0)), buf.getvalue(),
                       compress_type=zipfile.ZIP_DEFLATED)


def export_inputs(out: Path) -> None:
    """``reference.npz`` + ``section_{clean,spillover}.npz``: integer ``counts`` (cells x genes), ``genes``,
    ``cells``, ``labels`` (reference: ``cell_type``; sections: ``label``, which equals the truth here)."""
    import sim
    out.mkdir(parents=True, exist_ok=True)
    items = [("reference", sim.reference(seed=0, n_per_type=100), "cell_type"),
             ("section_clean", sim.simulate(seed=0, n_cells=600), "label"),
             ("section_spillover", sim.spillover(0.2, 0), "label")]
    for name, a, key in items:
        counts = np.asarray(a.layers["counts"].todense())
        assert counts.min() >= 0 and np.array_equal(counts, np.round(counts))
        save_npz(out / f"{name}.npz", counts=counts.astype(np.int32),
                 genes=np.asarray(a.var_names, dtype=str), cells=np.asarray(a.obs_names, dtype=str),
                 labels=a.obs[key].astype(str).to_numpy().astype(str))


def to_adata(npz: Path, key: str = "celltype"):
    """AnnData with float64 counts in X, ``layers['counts']`` (txsim ``layer=``) and ``layers['raw']`` (the
    Moldia functions read 'raw'); labels in ``obs[key]``."""
    import anndata as ad
    import pandas as pd
    d = np.load(npz)
    a = ad.AnnData(X=d["counts"].astype(np.float64), var=pd.DataFrame(index=d["genes"]),
                   obs=pd.DataFrame({key: d["labels"]}, index=d["cells"]))
    a.layers["counts"] = a.X.copy()
    a.layers["raw"] = a.X.copy()
    return a


def shared_genes(sp, sc):
    """Subset both to the reference's genes present in the section, as OpenProblems task_ist_preprocessing
    does before calling txsim (txsim indexes the reference by the section's genes: NegCtrl* would raise)."""
    genes = [g for g in sc.var_names if g in sp.var_names]
    return sp[:, genes].copy(), sc[:, genes].copy()


def num(x):
    x = float(x)
    return None if np.isnan(x) else x


def series(s) -> dict:
    return {str(k): num(v) for k, v in s.items()}


def write_json(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=1, allow_nan=False) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", type=Path, default=UPSTREAM, help="dir with the txsim and Moldia clones")
    args = ap.parse_args()
    setup_paths()
    txsim_dir = pinned(args.upstream / "txsim", "txsim")
    moldia_dir = pinned(args.upstream / "Xenium_benchmarking", "Xenium_benchmarking")
    sys.path.insert(0, str(txsim_dir))
    import txsim as tx
    assert Path(tx.__file__).resolve().is_relative_to(txsim_dir.resolve()), tx.__file__
    seg = load_module(moldia_dir / "notebooks" / "5_segmentation_benchmark" / "metrics.py",
                      "moldia_seg_metrics")
    import spatial_anno_metrics as sam
    assert Path(sam.__file__).resolve().is_relative_to(REPO / "src"), sam.__file__

    export_inputs(OUT)
    n_cells = {s: len(np.load(OUT / f"section_{s}.npz")["cells"]) for s in SECTIONS}
    perms = {s: np.stack([np.random.default_rng(k).permutation(n) for k in range(N_PERM)]).astype(np.int16)
             for s, n in n_cells.items()}
    save_npz(OUT / "label_permutations.npz", **perms)

    ref = to_adata(OUT / "reference.npz")
    kw = {"key": "celltype", "layer": "counts"}
    results, ours = {}, {}
    for s in SECTIONS:
        sp, sc = shared_genes(to_adata(OUT / f"section_{s}.npz"), ref)
        r = {"section": f"section_{s}.npz", "reference": "reference.npz"}
        for fl in ("reads", "cells"):
            fn = getattr(tx.metrics, f"negative_marker_purity_{fl}")
            v, per_gene, per_ct = fn(sp.copy(), sc.copy(), pipeline_output=False, **kw)
            r[f"txsim_{fl}"] = {"overall": num(v), "per_celltype": series(per_ct),
                                "per_gene": series(per_gene)}
            # the Moldia copy crashes with pipeline_output=False (`.cip` typo), so only its overall value
            r[f"moldia_seg_{fl}"] = {"overall": num(getattr(seg, f"negative_marker_purity_{fl}")(
                sp.copy(), sc.copy(), key="celltype", pipeline_output=True))}
            vals = []
            for p in perms[s]:
                spp = sp.copy()
                spp.obs["celltype"] = sp.obs["celltype"].to_numpy()[p]
                vals.append(num(fn(spp, sc.copy(), pipeline_output=True, **kw)))
            base = float(np.mean(vals))
            r[f"txsim_{fl}_permuted"] = {"values": vals, "mean": base,
                                         "scaled_overall": (r[f"txsim_{fl}"]["overall"] - base) / (1 - base)}
        results[s] = r

        a = to_adata(OUT / f"section_{s}.npz", key="cell_type")      # ours: full panel, natural call
        o = sam.nmp(a, to_adata(OUT / "reference.npz", key="cell_type"))
        ours[s] = {**o, "per_celltype": a.obs.groupby("cell_type")["nmp"].mean().to_dict()}

    write_json(OUT / "expected.json", {
        "generated": date.today().isoformat(),
        "upstream": {k: f"{k}@{v}" for k, v in PINS.items()},
        "n_permutations": N_PERM,
        "results": results})

    print("| section | txsim reads | txsim cells | Moldia reads | Moldia cells | txsim reads, scaled | "
          "txsim cells, scaled | ours mean_nmp | ours types_without_negatives |")
    print("|---|---|---|---|---|---|---|---|---|")
    for s in SECTIONS:
        r = results[s]
        print(f"| {s} | {r['txsim_reads']['overall']:.6f} | {r['txsim_cells']['overall']:.6f} | "
              f"{r['moldia_seg_reads']['overall']:.6f} | {r['moldia_seg_cells']['overall']:.6f} | "
              f"{r['txsim_reads_permuted']['scaled_overall']:.4f} | "
              f"{r['txsim_cells_permuted']['scaled_overall']:.4f} | "
              f"{ours[s]['mean_nmp']:.6f} | {ours[s]['types_without_negatives']} |")
    print("\n| section | cell type | txsim reads per-type excess | ours mean nmp |")
    print("|---|---|---|---|")
    for s in SECTIONS:
        for ct, v in results[s]["txsim_reads"]["per_celltype"].items():
            print(f"| {s} | {ct} | {v:.6f} | {ours[s]['per_celltype'][ct]:.6f} |")
    for s in SECTIONS:
        means = {fl: round(results[s][f"txsim_{fl}_permuted"]["mean"], 6) for fl in ("reads", "cells")}
        print(s, "permuted means:", means,
              "| txsim cells negatives (genes, types):", len(results[s]["txsim_cells"]["per_gene"]),
              len(results[s]["txsim_cells"]["per_celltype"]))


if __name__ == "__main__":
    main()
