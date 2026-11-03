#!/usr/bin/env python
"""Inputs for the SpatialQM signal-QC parity golden (tests/goldens/signalqc_spatialqm). Not run in CI.

Data: 10x Xenium FFPE Human Breast Cancer Rep1 (Janesick et al. 2023), files and md5 in SOURCE.md.

  python make_inputs.py section  # <data>/section.h5ad: all cells (tutorial / scale test)
  python make_inputs.py crop     # golden fixtures: counts.npz, cells.csv, transcripts.npz
  python make_inputs.py outs     # fixtures -> Xenium-style folder for make_expected.R + ours.json

`outs` reads only the committed fixtures, so the SpatialQM input is reproducible from the repo, under either venv:
cells.csv is copied byte for byte (a pandas round trip re-parses ~117 centroids one ulp differently on pandas 2.3
vs 3.1). It also writes fixtures.sha256.json, which make_expected.R copies into meta.input.sha256.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

DATA = Path("/srv/GT/analysis/pgueguen/sam_dev/data/xenium_breast")
OUTS = Path("/srv/GT/analysis/pgueguen/sam_dev/runs/wf2b-signal/crop_outs")
GOLDEN = Path(__file__).resolve().parents[3] / "tests" / "goldens" / "signalqc_spatialqm"
PREFIX = "Xenium_FFPE_Human_Breast_Cancer_Rep1_"
# 8 of the 28 'Negative Control Probe' features here are named antisense_<GENE> (SpatialQM's regexes miss them)
CONTROL_RE = r"^(NegControlProbe|NegControlCodeword|BLANK|Unassigned|antisense_)"
MARKERS = {"epithelial": ["EPCAM", "KRT7", "KRT8", "CDH1"], "stromal": ["LUM", "POSTN", "PDGFRA", "DPT"],
           "immune": ["PTPRC", "CD3E", "CD68"]}


def read_section(data: Path):
    import h5py

    with h5py.File(data / f"{PREFIX}cell_feature_matrix.h5", "r") as f:
        m = f["matrix"]
        n_feat, n_cells = (int(v) for v in m["shape"][:])
        # 10x stores features x cells as CSC: the same arrays read as cells x features CSR
        X = sp.csr_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(n_cells, n_feat))
        var = pd.DataFrame({k: m["features"][k][:].astype(str) for k in ("id", "name", "feature_type")})
        barcodes = m["barcodes"][:].astype(str)
    cells = pd.read_csv(data / f"{PREFIX}cells.csv.gz")
    if not (cells["cell_id"].astype(str).to_numpy() == barcodes).all():
        raise ValueError("cells.csv.gz and the h5 barcodes are not in the same order")
    var["control"] = (var["feature_type"] != "Gene Expression").to_numpy()
    if not (var["name"].str.match(CONTROL_RE).to_numpy() == var["control"].to_numpy()).all():
        raise ValueError("control flag from feature_type disagrees with the control-name pattern")
    return X, var, cells


def cmd_section(args):
    import anndata as ad

    X, var, cells = read_section(args.data)
    obs = cells.set_axis(cells["cell_id"].astype(str).to_numpy())
    a = ad.AnnData(X=X.astype(np.float32), obs=obs, var=var.set_axis(var.pop("name").to_numpy()))
    a.layers["counts"] = X
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()
    a.uns["source"] = {"dataset": "10x Xenium FFPE Human Breast Cancer Rep1 (Janesick et al. 2023)",
                       "url": f"https://cf.10xgenomics.com/samples/xenium/1.0.1/{PREFIX[:-1]}/",
                       "licence": "CC BY 4.0", "built_by": "benchmarks/parity/signalqc_spatialqm/make_inputs.py"}
    a.write_h5ad(args.data / "section.h5ad")
    genes = np.asarray(X[:, ~a.var["control"].to_numpy()].sum(1)).ravel()
    print(json.dumps({"n_cells": a.n_obs, "n_features": a.n_vars, "n_controls": int(a.var["control"].sum()),
                      "median_gene_counts_per_cell": float(np.median(genes))}))


def pick_window(X, names, xy, max_cells: int, max_half: float, step: float = 250.0):
    """Centre (on a `step` grid) whose `max_cells` Chebyshev-nearest cells have the most even
    epithelial / stromal / immune mix, among interior windows: half-side <= `max_half` and every
    quadrant around the centre holding >= 15 % of the cells (excludes glass and tissue edges)."""
    score = np.column_stack([np.asarray(X[:, np.isin(names, g)].sum(1)).ravel() for g in MARKERS.values()])
    cls = np.where(score.max(1) > 0, score.argmax(1), -1)
    best = (-1.0, None)
    lo, hi = xy.min(0), xy.max(0)
    for cx in np.arange(lo[0], hi[0], step):
        for cy in np.arange(lo[1], hi[1], step):
            d = np.abs(xy - (cx, cy)).max(1)
            dk = np.partition(d, [max_cells - 1, max_cells])[[max_cells - 1, max_cells]]
            h = float(dk.mean())                       # cells with d < h: exactly max_cells unless tied
            if h > max_half:
                continue
            q = 2 * (xy[d < h, 0] > cx) + (xy[d < h, 1] > cy)
            if np.bincount(q, minlength=4).min() < 0.15 * max_cells:
                continue
            frac = np.bincount(cls[(d < h) & (cls >= 0)], minlength=3) / max_cells
            ent = float(-(frac[frac > 0] * np.log(frac[frac > 0])).sum())
            if ent > best[0]:
                best = (ent, (float(cx), float(cy), h))
    return best[1]


def cmd_crop(args):
    import pyarrow as pa
    import pyarrow.parquet as pq

    X, var, cells = read_section(args.data)
    names = var["name"].to_numpy()
    xy = cells[["x_centroid", "y_centroid"]].to_numpy()
    cx, cy, h = pick_window(X, names, xy, args.max_cells, args.max_half)
    sel = np.abs(xy - (cx, cy)).max(1) < h
    t = pq.read_table(args.data / f"{PREFIX}transcripts.parquet",
                      columns=["cell_id", "overlaps_nucleus", "feature_name", "qv"],
                      filters=[("x_location", ">", cx - h), ("x_location", "<", cx + h),
                               ("y_location", ">", cy - h), ("y_location", "<", cy + h)])
    fidx = pd.Index(names).get_indexer(t["feature_name"].cast(pa.string()).to_numpy(zero_copy_only=False))
    if (fidx < 0).any():
        raise ValueError("transcript feature not in the count-matrix features")
    Xc = X[sel]
    g = args.golden
    g.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(g / "counts.npz", data=Xc.data.astype(np.int32), indices=Xc.indices.astype(np.uint16),
                        indptr=Xc.indptr.astype(np.int32), shape=np.array(Xc.shape, dtype=np.int64),
                        cell_id=cells["cell_id"].to_numpy()[sel].astype(np.int32),
                        feature_id=var["id"].to_numpy().astype(str), feature_name=names.astype(str),
                        feature_type=var["feature_type"].to_numpy().astype(str))
    cells[sel].to_csv(g / "cells.csv", index=False)
    np.savez_compressed(g / "transcripts.npz", feature_idx=fidx.astype(np.uint16),
                        cell_id=t["cell_id"].to_numpy().astype(np.int32),
                        overlaps_nucleus=t["overlaps_nucleus"].to_numpy().astype(np.uint8),
                        qv=t["qv"].to_numpy().astype(np.float32))
    print(json.dumps({"center_um": [cx, cy], "half_side_um": h, "n_cells": int(sel.sum()),
                      "n_transcripts": t.num_rows,
                      "sizes_bytes": {p.name: p.stat().st_size for p in sorted(g.iterdir()) if p.is_file()}}))


def load_golden(g: Path):
    c = np.load(g / "counts.npz", allow_pickle=False)
    X = sp.csr_matrix((c["data"], c["indices"], c["indptr"]), shape=tuple(c["shape"]))
    t = np.load(g / "transcripts.npz", allow_pickle=False)
    return X, {k: c[k] for k in ("feature_id", "feature_name", "feature_type")}, pd.read_csv(g / "cells.csv"), t


def write_outs(X, feats, cells, t, outs: Path, cells_csv: Path):
    """The crop as a Xenium-1.3+-style folder (string cell ids, 'UNASSIGNED'), the format SpatialQM parses.
    ``cells.csv.gz`` is the committed ``cells_csv`` gzipped as is (mtime 0), never re-written by pandas."""
    import scipy.io

    if outs.exists():
        shutil.rmtree(outs)
    mdir = outs / "cell_feature_matrix"
    mdir.mkdir(parents=True)
    scipy.io.mmwrite(str(mdir / "matrix.mtx"), X.T.tocoo(), field="integer")
    with open(mdir / "matrix.mtx", "rb") as src, gzip.open(mdir / "matrix.mtx.gz", "wb") as dst:
        shutil.copyfileobj(src, dst)
    (mdir / "matrix.mtx").unlink()
    pd.DataFrame(feats).to_csv(mdir / "features.tsv.gz", sep="\t", header=False, index=False)
    cells["cell_id"].to_csv(mdir / "barcodes.tsv.gz", header=False, index=False)
    with open(cells_csv, "rb") as src, gzip.GzipFile(outs / "cells.csv.gz", "wb", mtime=0) as dst:
        shutil.copyfileobj(src, dst)
    cid = t["cell_id"].astype(str).astype(object)
    cid[t["cell_id"] == -1] = "UNASSIGNED"
    pd.DataFrame({"cell_id": cid, "overlaps_nucleus": t["overlaps_nucleus"],
                  "feature_name": feats["feature_name"][t["feature_idx"]], "qv": t["qv"]}
                 ).to_csv(outs / "transcripts.csv.gz", index=False)


def _nan_to_none(o):
    if isinstance(o, dict):
        return {k: _nan_to_none(v) for k, v in o.items()}
    if isinstance(o, float) and math.isnan(o):
        return None
    return o


def run_ours(X, feats, cells) -> dict:
    """Our current signal_qc API, default arguments, on the crop (all 541 features in X)."""
    import anndata as ad
    from spatial_anno_metrics import signal_qc as sq

    ftype = feats["feature_type"]
    a = ad.AnnData(X=X.astype(np.float32),
                   var=pd.DataFrame({"feature_type": ftype, "control": ftype != "Gene Expression"},
                                    index=feats["feature_name"]))
    a.obs_names = cells["cell_id"].astype(str).to_numpy()
    a.obs["cell_area"] = cells["cell_area"].to_numpy()
    a.obs["total_counts"] = cells["total_counts"].to_numpy()
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()
    out = {"sparsity": sq.sparsity(a), "detection_entropy": sq.detection_entropy(a),
           "signal_to_noise": sq.signal_to_noise(a), "tx_per_area": sq.tx_per_area(a),
           "moran_signal": sq.moran_signal(a), "tx_per_cell": sq.tx_per_cell(a)}
    out["moran_signal"].pop("per_feature")
    for name in ("moran_signal", "signal_to_noise", "detection_entropy"):     # the SpatialQM recipes
        r = getattr(sq, name)(a.copy(), method="spatialqm")
        r.pop("per_feature", None)
        out[f"{name}_spatialqm"] = r
    mi = a.uns["moranI"]["I"].reindex(feats["feature_name"])
    by_type = {t: mi[ftype == t] for t in np.unique(ftype)}
    out["moran_I_per_feature"] = {k: float(v) for k, v in mi.items()}
    ctrl_per_cell = np.asarray(X[:, ftype != "Gene Expression"].sum(1)).ravel()
    out["explain"] = {  # numbers that locate the gaps; not part of our API
        "frac_cells_zero_control_counts": float((ctrl_per_cell == 0).mean()),
        "moran_nanmedian_by_feature_type": {t: float(np.nanmedian(v)) for t, v in by_type.items()},
        "moran_n_nan_by_feature_type": {t: int(v.isna().sum()) for t, v in by_type.items()},
    }
    return _nan_to_none(out)


def cmd_outs(args):
    X, feats, cells, t = load_golden(args.golden)
    write_outs(X, feats, cells, t, args.outs, args.golden / "cells.csv")
    sha = {f: hashlib.sha256((args.golden / f).read_bytes()).hexdigest()
           for f in ("counts.npz", "cells.csv", "transcripts.npz")}
    (args.outs / "fixtures.sha256.json").write_text(json.dumps(sha, indent=1))
    ours = run_ours(X, feats, cells)
    (args.outs / "ours.json").write_text(json.dumps(ours, indent=1, allow_nan=False))
    print(json.dumps({k: v for k, v in ours.items() if k != "moran_I_per_feature"}, indent=1))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("step", choices=["section", "crop", "outs"])
    p.add_argument("--data", type=Path, default=DATA)
    p.add_argument("--golden", type=Path, default=GOLDEN)
    p.add_argument("--outs", type=Path, default=OUTS)
    p.add_argument("--max-cells", type=int, default=2000)
    p.add_argument("--max-half", type=float, default=400.0, help="max window half-side (um)")
    a = p.parse_args()
    {"section": cmd_section, "crop": cmd_crop, "outs": cmd_outs}[a.step](a)
