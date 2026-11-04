"""Build the scDesign3 world cache: crop the real section, hand it to R, pack R's fit for ``world.py``.

    python build.py prepare   # section.h5ad + author labels -> <cache>/input/ (counts.mtx, genes.txt, obs.csv)
    Rscript generate.R <cache> <n_cores> [k_gp]           # the scDesign3 fit (SLURM: generate.sbatch)
    python build.py pack      # <cache>/fit/ -> <cache>/params.h5ad (fitted means, NB dispersion, copulas, markers)
    python build.py realism   # real crop vs scDesign3 draws (R simu_new and the Python port) vs the toy world

The cache lives outside git (``SAM_SCDESIGN3_CACHE``, default below). Inputs:

* ``section.h5ad``: 10x Xenium FFPE human breast cancer, replicate 1 (Janesick et al. 2023, CC BY), 167,780
  cells x 541 features (313 genes + negative controls), raw counts in ``layers['counts']``.
* Author cell types: GEO GSM7780153 ``Xenium_R1_Fig1-5_supervised.csv.gz`` (Janesick et al. 2023, the
  supervised annotation of Figs 1-5), one row per cell, ``Barcode`` = ``cell_id``.

The crop is x < 1300 um, y >= 4000 um (bottom-left of the section: invasive tumour, stroma and the immune
aggregate), keeping only the six author types that map onto the synthetic world's names (``TYPE_MAP``).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
DATA = Path(os.environ.get("SAM_DATA", "/srv/GT/analysis/pgueguen/sam_dev/data"))
CACHE = Path(os.environ.get("SAM_SCDESIGN3_CACHE", DATA / "sim_worlds" / "scdesign3"))
SECTION = DATA / "xenium_breast" / "section.h5ad"
LABELS = DATA / "xenium_breast" / "janesick_annotation" / "GSM7780153_Xenium_R1_Fig1-5_supervised.csv.gz"

# author label -> world type (the names benchmarks/synthetic/sim.py scenarios use)
TYPE_MAP = {"CD4+_T_Cells": "T_CD4", "CD8+_T_Cells": "T_CD8", "B_Cells": "B", "Macrophages_1": "Myeloid",
            "Stromal": "Fibro", "Invasive_Tumor": "Epi"}
CROP = {"x_max": 1300.0, "y_min": 4000.0}
MIN_DETECTION = 0.005          # a gene is fitted if detected in >= 0.5 % of crop cells


def _crop():
    import anndata as ad
    a = ad.read_h5ad(SECTION)
    lab = pd.read_csv(LABELS)
    if not (a.obs["cell_id"].astype(int).to_numpy() == lab["Barcode"].to_numpy()).all():
        raise ValueError("author labels are not in section order")
    a.obs["author_type"] = lab["Cluster"].to_numpy()
    a.obs["cell_type"] = a.obs["author_type"].map(TYPE_MAP)
    o = a.obs
    keep = ((o["x_centroid"] < CROP["x_max"]) & (o["y_centroid"] >= CROP["y_min"]) & o["cell_type"].notna()).to_numpy()
    b = a[keep].copy()
    b.layers["counts"] = sp.csr_matrix(b.layers["counts"])
    return a, b


def prepare():
    a, b = _crop()
    X = b.layers["counts"]
    ftype = b.var["feature_type"].to_numpy()
    det = np.asarray((X > 0).mean(0)).ravel()
    genes = (ftype == "Gene Expression") & (det >= MIN_DETECTION)
    ctrl = ftype == "Negative Control Probe"
    out = CACHE / "input"
    out.mkdir(parents=True, exist_ok=True)
    import scipy.io
    scipy.io.mmwrite(out / "counts.mtx", X[:, genes].T.tocoo().astype(np.int32))        # genes x cells, for R
    pd.Series(b.var_names[genes]).to_csv(out / "genes.txt", index=False, header=False)
    obs = pd.DataFrame({"cell": b.obs_names, "cell_type": b.obs["cell_type"].to_numpy(),
                        "author_type": b.obs["author_type"].to_numpy(),
                        "x": b.obs["x_centroid"].to_numpy(), "y": b.obs["y_centroid"].to_numpy(),
                        "cell_area": b.obs["cell_area"].to_numpy(),
                        "log_area": np.log(b.obs["cell_area"].to_numpy()),
                        "total_counts_genes": np.asarray(X[:, genes].sum(1)).ravel()})
    obs.to_csv(out / "obs.csv", index=False)
    # real counts of the crop (fitted genes + negative-control probes) for pack / realism
    real = b[:, genes | ctrl].copy()
    real.X = real.layers["counts"]
    del real.layers["counts"]
    real.var["control"] = ctrl[genes | ctrl]
    real.obs = real.obs[["cell_id", "x_centroid", "y_centroid", "cell_area", "author_type", "cell_type"]].copy()
    real.uns = {"source": "10x Xenium FFPE Human Breast Cancer Rep1 (Janesick 2023, CC BY); crop " + json.dumps(CROP)}
    real.write_h5ad(out / "real_crop.h5ad", compression="gzip")
    info = {"n_cells": int(b.n_obs), "n_genes_fitted": int(genes.sum()), "n_controls": int(ctrl.sum()),
            "genes_dropped_low_detection": list(b.var_names[(ftype == "Gene Expression") & ~genes]),
            "crop": CROP, "type_counts": b.obs["cell_type"].value_counts().to_dict(),
            "section_cells": int(a.n_obs)}
    (out / "prepare_info.json").write_text(json.dumps(info, indent=1))
    print(json.dumps(info, indent=1))


# --------------------------------------------------------------------------- #
# pack: R outputs -> params.h5ad
# --------------------------------------------------------------------------- #
TYPES = ["T_CD4", "T_CD8", "B", "Myeloid", "Fibro", "Epi"]


def derive_markers(type_mean: pd.DataFrame, per_type: int = 4, n_lineage: int = 2, min_mean: float = 0.1):
    """Marker sets from the FITTED type means (types x genes, mean count per cell).

    Specificity of gene g for type t = log2((mu_t + eps) / (max over the other types + eps)), eps = 0.01, among
    genes with mu_t >= ``min_mean``. As in the toy world, T_CD4 and T_CD8 share ``n_lineage`` T-lineage markers
    (specificity of min(mu_CD4, mu_CD8) against the non-T types) and take ``per_type - n_lineage`` private
    markers each (specific against every other type, the other T subtype included); every other type takes its
    ``per_type`` most specific genes. A gene goes to one set only (the lineage genes to both T sets)."""
    eps = 0.01
    used: set[str] = set()

    def top(score, n, eligible):
        s = score[eligible & ~score.index.isin(used)].sort_values(ascending=False)
        picked = list(s.index[:n])
        used.update(picked)
        return picked, s.iloc[:n]

    scores = {}
    non_t = [t for t in TYPES if not t.startswith("T_")]
    tmin = type_mean.loc[["T_CD4", "T_CD8"]].min(0)
    lin_score = np.log2((tmin + eps) / (type_mean.loc[non_t].max(0) + eps))
    lineage, scores["T_lineage"] = top(lin_score, n_lineage, tmin >= min_mean)
    ms = {}
    for t in TYPES:
        other = [u for u in TYPES if u != t]
        sc = np.log2((type_mean.loc[t] + eps) / (type_mean.loc[other].max(0) + eps))
        n = per_type - n_lineage if t.startswith("T_") else per_type
        picked, scores[t] = top(sc, n, type_mean.loc[t] >= min_mean)
        ms[t] = (lineage + picked) if t.startswith("T_") else picked
    return ms, lineage, {k: v.round(2).to_dict() for k, v in scores.items()}


def pack():
    import anndata as ad
    fit = CACHE / "fit"
    real = ad.read_h5ad(CACHE / "input" / "real_crop.h5ad")
    obs = pd.read_csv(CACHE / "input" / "obs.csv")
    mean = pd.read_csv(fit / "mean.csv.gz", index_col=0)               # cells x genes
    sigma = pd.read_csv(fit / "sigma.csv", index_col=0)["sigma"]       # NBI sigma = 1 / theta, per gene
    info = json.loads((fit / "fit_info.json").read_text())
    failed = [g for g in info.get("failed_genes", []) if g in mean.columns]
    genes = [g for g in mean.columns if g not in failed and np.isfinite(sigma.get(g, np.nan))
             and np.isfinite(mean[g]).all()]
    mean = mean[genes]
    if not (mean.index.astype(str) == obs["cell"].astype(str)).all():
        raise ValueError("fit cell order differs from obs.csv")
    cop = {}
    for t in TYPES:
        C = pd.read_csv(fit / f"copula_{t}.csv.gz", index_col=0)
        cop[t] = C.loc[genes, genes].to_numpy(np.float32)
    type_mean = mean.groupby(obs["cell_type"].to_numpy()).mean().loc[TYPES]
    ms, lineage, scores = derive_markers(type_mean)
    ctrl_names = list(real.var_names[real.var["control"].to_numpy()])
    Xr = real[:, ctrl_names].X
    ctrl_rate = np.asarray(Xr.mean(0)).ravel()
    p = ad.AnnData(X=mean.to_numpy(np.float32),
                   obs=obs.set_index(obs["cell"].astype(str)).drop(columns="cell"),
                   var=pd.DataFrame({"sigma": sigma.loc[genes].to_numpy()}, index=genes))
    p.layers["real"] = sp.csr_matrix(real[p.obs_names, genes].X.astype(np.float32))
    p.uns["copula"] = cop
    p.uns["copula_types"] = TYPES
    p.uns["control_names"] = ctrl_names
    p.uns["control_rate"] = ctrl_rate
    p.uns["marker_sets"] = ms
    p.uns["T_lineage"] = lineage
    p.uns["marker_scores_log2"] = json.dumps(scores)
    p.uns["fit_info"] = json.dumps(info)
    p.write_h5ad(CACHE / "params.h5ad", compression="gzip")
    print("genes", len(genes), "failed", failed)
    print("markers", json.dumps(ms, indent=1))
    print("scores", json.dumps(scores, indent=1))


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else ""
    if cmd == "prepare":
        prepare()
    elif cmd == "pack":
        pack()
    elif cmd == "realism":
        from realism import main as realism_main
        realism_main()
    else:
        sys.exit(__doc__)
