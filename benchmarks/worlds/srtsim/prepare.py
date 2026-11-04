"""Build the SRTsim world's cache in three steps (step 2 is ``generate.R``):

    python prepare.py crops --seeds 5 --run RUN            # 1. pick crops, write SRTsim inputs
    Rscript generate.R RUN 5 6                             # 2. SRTsim fit + count (sbatch: generate.sbatch)
    python prepare.py cache --seeds 5 --run RUN            # 3. assemble seed_<k>.h5ad + marker_sets.json

Inputs: the 10x Xenium FFPE human breast section, Rep1 (``SAM_SECTION_H5AD``: AnnData with
``layers['counts']``, ``obsm['spatial']`` in um, ``obs['cell_id']``, ``obs['cell_area']``, ``var['feature_type']``)
and the authors' supervised annotation (GEO GSM7780153, ``..._Fig1-5_supervised.csv.gz``, fetched if absent).
Output: ``SAM_SRTSIM_CACHE`` (see ``world.py``).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import world  # noqa: E402

SECTION = Path(__import__("os").environ.get(
    "SAM_SECTION_H5AD", "/srv/GT/analysis/pgueguen/sam_dev/data/xenium_breast/section.h5ad"))
ANNOT_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7780nnn/GSM7780153/suppl/"
             "GSM7780153_Xenium_R1_Fig1-5_supervised.csv.gz")
ANNOT_MD5 = "197fda5b0b9b417a4683bdc14d663afd"
STEP_UM = 25.0          # candidate window corners on a 25 um grid
MIN_FRAC, MIN_N, MIN_TOTAL = 0.04, 12, 600   # every type >= 4 % and >= 12 cells; >= 600 typed cells per crop
N_MARKERS = 4           # per type, as in the toy world (T types: 2 shared T-lineage + 2 own)


def _annotation(path: Path) -> pd.DataFrame:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(ANNOT_URL, path)
    if hashlib.md5(path.read_bytes()).hexdigest() != ANNOT_MD5:
        raise ValueError(f"{path}: md5 differs from the GEO file ({ANNOT_MD5})")
    return pd.read_csv(path)


def _section(annot: Path):
    a = ad.read_h5ad(SECTION)
    ann = _annotation(annot)
    if not (a.obs["cell_id"].astype(int).to_numpy() == ann["Barcode"].to_numpy()).all():
        raise ValueError("annotation barcodes do not match the section's cell_id order")
    a.obs["author_cluster"] = ann["Cluster"].to_numpy()
    a.obs["truth"] = ann["Cluster"].map(world.TYPEMAP).to_numpy()
    ft = a.var["feature_type"].to_numpy()
    feats = np.flatnonzero(ft == "Gene Expression").tolist() + np.flatnonzero(ft == "Negative Control Probe").tolist()
    return a[:, feats].copy()


def _candidates(xy, truth):
    """Lower-left corners (um) of every CROP_UM window, on the STEP_UM grid, that holds every type."""
    k = int(world.CROP_UM // STEP_UM)
    ex = np.arange(0, xy[:, 0].max() + STEP_UM, STEP_UM)
    ey = np.arange(0, xy[:, 1].max() + STEP_UM, STEP_UM)
    H = np.stack([np.histogram2d(xy[truth == t, 0], xy[truth == t, 1], bins=[ex, ey])[0] for t in world.TYPES])
    C = np.pad(H.cumsum(1).cumsum(2), ((0, 0), (1, 0), (1, 0)))
    S = C[:, k:, k:] - C[:, :-k, k:] - C[:, k:, :-k] + C[:, :-k, :-k]
    n = S.sum(0)
    ok = (n >= MIN_TOTAL) & (S.min(0) >= np.maximum(MIN_N, MIN_FRAC * n))
    ix, iy = np.nonzero(ok)
    return np.column_stack([ex[ix], ey[iy]])


def crops(args):
    run = Path(args.run)
    run.mkdir(parents=True, exist_ok=True)
    a = _section(Path(args.annotation))
    xy = a.obsm["spatial"]
    truth = a.obs["truth"].to_numpy().astype(object)
    cand = _candidates(xy, truth)
    print(f"{len(cand)} candidate {world.CROP_UM:g} um windows", flush=True)
    X = a.layers["counts"].tocsr() if sp.issparse(a.layers["counts"]) else sp.csr_matrix(a.layers["counts"])
    pd.DataFrame({"gene": a.var_names, "feature_type": a.var["feature_type"].to_numpy(),
                  "feature_id": a.var["id"].to_numpy() if "id" in a.var else a.var_names}).to_csv(
        run / "genes.csv", index=False)
    for seed in range(args.seeds):
        rng = np.random.default_rng(seed)
        picked = []
        for c in rng.permutation(len(cand)):
            if all(np.abs(cand[c] - cand[p]).max() >= world.CROP_UM for p in picked):
                picked.append(c)
            if len(picked) == world.N_CROPS:
                break
        rows = []
        for crop, c in enumerate(picked):
            x0, y0 = cand[c]
            inside = ((xy[:, 0] >= x0) & (xy[:, 0] < x0 + world.CROP_UM) & (xy[:, 1] >= y0)
                      & (xy[:, 1] < y0 + world.CROP_UM) & pd.notna(truth))
            for i in np.flatnonzero(inside):
                rows.append(dict(cell=f"x{a.obs['cell_id'].iloc[i]}", row=i, crop=crop, label=truth[i],
                                 author_cluster=a.obs["author_cluster"].iloc[i],
                                 x_local=xy[i, 0] - x0, y_local=xy[i, 1] - y0, x_um=xy[i, 0], y_um=xy[i, 1],
                                 cell_area=float(a.obs["cell_area"].iloc[i])))
        cells = pd.DataFrame(rows)
        cells.drop(columns="row").to_csv(run / f"seed_{seed}_cells.csv", index=False)
        scipy.io.mmwrite(run / f"seed_{seed}_counts.mtx", X[cells["row"].to_numpy()].T.tocoo(), field="integer")
        comp = cells.groupby(["crop", "label"]).size().unstack(fill_value=0)[world.TYPES]
        print(f"seed {seed}: windows {[tuple(cand[c]) for c in picked]}\n{comp}", flush=True)


def _read_params(run: Path, seed: int, genes: list[str]):
    """generate.R's long tables -> arrays (crops x types x genes) / (crops x genes); NaN = not fitted."""
    gi = {g: i for i, g in enumerate(genes)}
    out = {}
    d = pd.read_csv(run / f"seed_{seed}_params_domain.csv")
    for p in ("pi0", "theta", "mu"):
        A = np.full((world.N_CROPS, len(world.TYPES), len(genes)), np.nan)
        A[d["crop"].to_numpy(), [world.TYPES.index(t) for t in d["label"]], [gi[g] for g in d["gene"]]] = d[p]
        out[f"dom_{p}"] = A
    t = pd.read_csv(run / f"seed_{seed}_params_tissue.csv")
    for p in ("pi0", "theta", "mu"):
        A = np.full((world.N_CROPS, len(genes)), np.nan)
        A[t["crop"].to_numpy(), [gi[g] for g in t["gene"]]] = t[p]
        out[f"tis_{p}"] = A
    out["dom_mean"] = np.nan_to_num((1 - out["dom_pi0"]) * out["dom_mu"])
    models = d.groupby("model").size().to_dict()
    return out, models


def _derive_markers(ws: list[ad.AnnData]):
    """Top type-specific genes from the SIMULATED counts of every cached seed: mean depth-normalised count per
    type; specificity = log2((m_t + eps) / (max_{u != t} m_u + eps)); a gene must be detected in >= 20 % of the
    type's cells. T_CD4 / T_CD8: the 2 genes most specific to T (both T types pooled) + their 2 most specific."""
    G = ws[0].var_names[~ws[0].var["control"].to_numpy()]
    Xs, ts = [], []
    for w in ws:
        X = np.asarray(w[:, G].layers["counts"].todense(), float)
        Xs.append(X / np.maximum(X.sum(1, keepdims=True), 1) * 100)
        ts.append(w.obs["truth"].astype(str).to_numpy())
    X, t = np.vstack(Xs), np.concatenate(ts)
    det = np.vstack([(X[t == k] > 0).mean(0) for k in world.TYPES])
    M = np.vstack([X[t == k].mean(0) for k in world.TYPES])
    eps = 0.01

    def spec(row_mean, other_means, detect):
        s = np.log2((row_mean + eps) / (other_means.max(0) + eps))
        s[detect < 0.2] = -np.inf
        return s

    table, ms = [], {}
    T = [world.TYPES.index("T_CD4"), world.TYPES.index("T_CD8")]
    nonT = [i for i in range(len(world.TYPES)) if i not in T]
    tmask = np.isin(t, ["T_CD4", "T_CD8"])
    s_lin = spec(X[tmask].mean(0), M[nonT], (X[tmask] > 0).mean(0))
    lineage = list(G[np.argsort(-s_lin)[:2]])
    for i, k in enumerate(world.TYPES):
        s = spec(M[i], np.delete(M, i, 0), det[i])
        order = [g for g in G[np.argsort(-s)] if g not in lineage]
        own = order[: N_MARKERS - 2] if k in ("T_CD4", "T_CD8") else order[:N_MARKERS]
        ms[k] = (lineage + own) if k in ("T_CD4", "T_CD8") else own
        gi = {g: j for j, g in enumerate(G)}
        for g in ms[k]:
            table.append(dict(type=k, gene=g, mean_in_type=M[i, gi[g]], max_mean_other=np.delete(M[:, gi[g]], i).max(),
                              log2_specificity=(s_lin[gi[g]] if g in lineage else s[gi[g]]),
                              detected_in_type=det[i, gi[g]], role="T lineage" if g in lineage else "own"))
    return ms, lineage, pd.DataFrame(table)


def cache(args):
    run, out = Path(args.run), Path(args.cache)
    out.mkdir(parents=True, exist_ok=True)
    genes_df = pd.read_csv(run / "genes.csv")
    genes = list(genes_df["gene"])
    ctrl = (genes_df["feature_type"] == "Negative Control Probe").to_numpy()
    names = list(genes)
    nctrl = 0
    for i, c in enumerate(ctrl):                                    # NegCtrl0..27 in panel order
        if c:
            names[i] = f"NegCtrl{nctrl}"
            nctrl += 1
    ws = []
    meta = {}
    for seed in range(args.seeds):
        cells = pd.read_csv(run / f"seed_{seed}_cells.csv")
        real = scipy.io.mmread(run / f"seed_{seed}_counts.mtx").T.tocsr()
        simc = scipy.io.mmread(run / f"seed_{seed}_sim_counts.mtx").T.tocsr()
        sim_cells = pd.read_csv(run / f"seed_{seed}_sim_cells.csv")["cell"].to_numpy()
        sim_genes = pd.read_csv(run / f"seed_{seed}_sim_genes.csv")["gene"].to_numpy()
        pos = pd.Series(np.arange(len(sim_cells)), index=sim_cells)[cells["cell"]].to_numpy()
        gpos = pd.Series(np.arange(len(sim_genes)), index=sim_genes)[genes].to_numpy()
        simc = simc[pos][:, gpos]
        params, models = _read_params(run, seed, genes)
        y_off = cells["crop"].to_numpy() * (world.CROP_UM + world.GAP_UM)
        obs = pd.DataFrame({"truth": pd.Categorical(cells["label"], categories=world.TYPES),
                            "crop": cells["crop"].to_numpy(), "author_cluster": cells["author_cluster"].to_numpy(),
                            "x_local": cells["x_local"].to_numpy(), "y_local": cells["y_local"].to_numpy(),
                            "x_um": cells["x_um"].to_numpy(), "y_um": cells["y_um"].to_numpy(),
                            "cell_area": cells["cell_area"].to_numpy()},
                           index=cells["cell"].to_numpy())
        var = pd.DataFrame({"control": ctrl, "feature_id": genes_df["feature_id"].fillna("").to_numpy(),
                            "panel_name": genes}, index=names)
        w = ad.AnnData(X=None, obs=obs, var=var)
        w.layers["counts"] = simc.astype("float32")
        w.layers["real"] = real.astype("float32")
        w.obsm["spatial"] = np.column_stack([cells["x_local"].to_numpy(), cells["y_local"].to_numpy() + y_off])
        w.uns.update(params)
        lib = np.asarray(real[:, ~ctrl].sum(1)).ravel()
        w.uns["section_mean_depth"] = float(lib.mean())
        w.uns["models"] = {k: int(v) for k, v in models.items()}
        timing = run / f"seed_{seed}_timing.csv"
        if timing.exists():
            w.uns["fit_seconds"] = pd.read_csv(timing)["seconds"].to_numpy()
        w.write_h5ad(out / f"seed_{seed}.h5ad", compression="gzip")
        ws.append(w)
        meta[seed] = dict(n_cells=int(w.n_obs), models=w.uns["models"])
        print(f"seed {seed}: {w.n_obs} cells, models {models}", flush=True)
    ms, lineage, table = _derive_markers(ws)
    (out / "marker_sets.json").write_text(json.dumps({"marker_sets": ms, "t_lineage": lineage,
                                                      "method": _derive_markers.__doc__.strip()}, indent=1))
    table.to_csv(out / "marker_table.csv", index=False)
    print(json.dumps(ms, indent=1))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("step", choices=["crops", "cache"])
    p.add_argument("--seeds", type=int, default=5)
    p.add_argument("--run", default="/srv/GT/analysis/pgueguen/sam_dev/runs/sim-srtsim/srtsim_inputs")
    p.add_argument("--cache", default=str(world.CACHE))
    p.add_argument("--annotation", default=str(SECTION.parent / "janesick_annotation"
                                                / "GSM7780153_Xenium_R1_Fig1-5_supervised.csv.gz"))
    a = p.parse_args()
    {"crops": crops, "cache": cache}[a.step](a)
