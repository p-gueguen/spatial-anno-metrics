"""Build the scCube world cache once (needs scCube + torch; world.py only reads the cache).

    python benchmarks/worlds/sccube/generate.py --yerly <dir with counts.mtx, genes.txt, meta.csv> \
        [--out benchmarks/worlds/sccube/cache] [--device cuda:0] [--epochs 3000] [--seeds 5]

scCube (Qian et al., Nat Commun 2024; github.com/ZJUFanLab/scCube) is installed from its GitHub repository,
the authors' only distribution route. Each step is skipped when its output exists:

1. ``prep``: Yerly et al. 2022 BCC scRNA-seq (author ``annotation``; the T cells split by the authors' scGate
   CD4T / CD8T calls) mapped onto the six names the sweep's scenarios use (``sim.TYPES``). Cells of each type are
   split 50/50: a training half (scCube's VAE) and a reference half (real counts, ``world.reference``). Panel =
   the 10x Xenium human breast panel genes present in Yerly (306 / 313); 960 extra genes with the lowest
   between-type F statistic serve ``panel_growth``. Marker sets come from the reference half.
2. ``train``: scCube's VAE on log1p CP10k (full-transcriptome library) of the training half.
3. ``pools``: per seed, scCube generates ``POOL`` cells per type (VAE reconstructions of training cells with
   sampled latents; the training order is permuted per seed so seeds draw different templates).
4. ``layouts``: per pattern x seed x (n_cells, weights): six square tiles, one per sample, each its own scCube
   spatial pattern (``PATTERNS``). Rare-type weights and the 600-cell edge-case section are pre-generated.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "synthetic"))
import sim  # noqa: E402  (TYPES, WEIGHTS, SCENARIOS - the contract the world must honour)

TYPES = list(sim.TYPES)
POOL = 700                 # cells per type and seed: rare_type at 1 % needs 0.22 x 0.99 / 0.85 x 2400 = 615 Epi
N_EXTRA = 960              # panel_growth goes up to +960 genes
N_CTRL = 20
TRAIN_CAP = 2500           # training cells per type (Epi / T_CD4 have more)
TILES = 6                  # one tile per sample
SPATIAL_SIZE = 30          # scCube's grid (default)
GAP = 5.0                  # space between tiles, in cell spacings: no kNN edge crosses tiles
# realism targets, measured on the 10x Xenium FFPE human breast Rep1 section (167,780 cells, 313 genes):
TARGET_MEDIAN = 164.0      # median transcripts per cell on the gene panel
TARGET_SDLOG = 0.728       # SD of log(transcripts per cell)
CTRL_RATIO = 0.00105 / 0.611   # counts per control feature / counts per gene feature, per cell
AREA_MEDIAN, AREA_SDLOG, AREA_LIB_R = 140.0, 0.773, 0.562   # cell area (um^2), and its correlation with log lib

PATTERNS = {   # scCube controls: delta = GRF length scale (domain size; spatial_size 30 = tile width), lamda =
    # share of cells that keep their domain label (the rest are permuted: mixing). Types fill quantile bands of
    # one field in sim.TYPES order, so neighbouring types share borders as in the toy stripes.
    "stripes": dict(kind="stripes"),                       # scCube custom stripes: the toy geometry, real expression
    "domains": dict(kind="grf", delta=25.0, lamda=1.0),    # few large graded domains
    "fine": dict(kind="grf", delta=6.0, lamda=1.0),        # many small patches (higher spatial resolution)
    "mixed": dict(kind="grf", delta=25.0, lamda=0.7),      # the domains with 30 % of labels permuted
}


def weight_vectors():
    """Every weight vector the scenarios ask for: the default and rare_type's levels (sim.rare_type's formula)."""
    out = [np.asarray(sim.WEIGHTS, float)]
    for s in sim.SCENARIOS["rare_type"][1]:
        w = sim.WEIGHTS.copy()
        w[1] = s
        w[[0, 2, 3, 4, 5]] *= (1 - s) / w[[0, 2, 3, 4, 5]].sum()
        out.append(w)
    return out


def wkey(w) -> str:
    w = np.asarray(w, float)
    return "w" + "-".join(str(int(round(1e4 * x))) for x in w / w.sum())


# --------------------------------------------------------------------------- #
def map_types(meta: pd.DataFrame) -> pd.Series:
    a, g = meta["annotation"].astype(str), meta["scGate_multi"].astype(str)
    t = pd.Series(np.nan, index=meta.index, dtype=object)
    t[(a == "Tcell") & (g == "CD4T")] = "T_CD4"
    t[(a == "Tcell") & (g == "CD8T")] = "T_CD8"
    t[a.isin(["Bcell", "PlasmaCell"])] = "B"
    t[a == "MoMacDC"] = "Myeloid"
    t[a == "Fibroblast"] = "Fibro"
    t[a == "Cancer_cells"] = "Epi"
    return t


def lognorm(X, lib):
    X = X.multiply(1e4 / np.maximum(lib, 1)[:, None]).tocsr() if sp.issparse(X) else X * (1e4 / np.maximum(lib, 1)[:, None])
    return X.log1p() if sp.issparse(X) else np.log1p(X)


def prep(yerly: Path, out: Path):
    from scipy.io import mmread
    import anndata as ad
    t0 = time.time()
    genes = np.array(Path(yerly / "genes.txt").read_text().split())
    meta = pd.read_csv(yerly / "meta.csv")
    X = sp.csr_matrix(mmread(str(yerly / "counts.mtx")).T)          # cells x genes
    types = map_types(meta)
    keep = types.notna().to_numpy()
    X, meta, types = X[keep], meta[keep].reset_index(drop=True), types[keep].reset_index(drop=True)
    lib = np.asarray(X.sum(1)).ravel()
    rng = np.random.default_rng(0)
    half = np.zeros(len(types), dtype="U5")
    for t in TYPES:
        ids = rng.permutation(np.flatnonzero(types.to_numpy() == t))
        half[ids[: len(ids) // 2]] = "train"
        half[ids[len(ids) // 2:]] = "ref"
    panel_all = (HERE / "xenium_breast_panel.txt").read_text().split()
    gidx = {g: i for i, g in enumerate(genes)}
    panel = [g for g in panel_all if g in gidx]
    tr = half == "train"
    L = lognorm(X[tr], lib[tr])
    # extra genes: expressed (>= 10 % of training cells), not ribosomal / mitochondrial, lowest F between types
    det = np.asarray((X[tr] > 0).mean(0)).ravel()
    bad = np.array([g.startswith(("RPL", "RPS", "MT-", "MRPL", "MRPS")) for g in genes])
    cand = np.flatnonzero((det >= 0.10) & ~bad & ~np.isin(genes, panel))
    Lc = L[:, cand].toarray()
    ty = types.to_numpy()[tr]
    m = Lc.mean(0)
    between = sum((ty == t).sum() * (Lc[ty == t].mean(0) - m) ** 2 for t in TYPES) / (len(TYPES) - 1)
    within = sum(((Lc[ty == t] - Lc[ty == t].mean(0)) ** 2).sum(0) for t in TYPES) / (len(ty) - len(TYPES))
    F = between / np.maximum(within, 1e-12)
    extra = list(genes[cand[np.argsort(F)[:N_EXTRA]]])
    feats = panel + extra
    fidx = [gidx[g] for g in feats]
    # training matrix: log1p CP10k (full-transcriptome library), capped per type
    tr_ids = []
    for t in TYPES:
        ids = np.flatnonzero(tr & (types.to_numpy() == t))
        tr_ids += list(rng.choice(ids, min(TRAIN_CAP, len(ids)), replace=False))
    tr_ids = np.sort(np.array(tr_ids))
    Ltr = lognorm(X[tr_ids][:, fidx].astype(np.float64), lib[tr_ids]).toarray().astype(np.float32)
    obs = pd.DataFrame({"Cell": [f"y{i}" for i in tr_ids], "Cell_type": types.to_numpy()[tr_ids],
                        "patient": meta["patient_bcc"].to_numpy()[tr_ids]}, index=[f"y{i}" for i in tr_ids])
    a = ad.AnnData(Ltr, obs=obs, var=pd.DataFrame(index=feats))
    a.write_h5ad(out / "train.h5ad")
    # reference half: real counts on panel + extra genes
    rf = np.flatnonzero(half == "ref")
    r = ad.AnnData(X[rf][:, fidx].astype(np.float32), var=pd.DataFrame(index=feats),
                   obs=pd.DataFrame({"cell_type": types.to_numpy()[rf], "patient": meta["patient_bcc"].to_numpy()[rf],
                                     "lib_full": lib[rf]}, index=[f"y{i}" for i in rf]))
    r.write_h5ad(out / "ref_pool.h5ad")
    # marker sets from the reference half (log1p CP10k means; what a user with a reference would do)
    Lr = lognorm(X[rf][:, [gidx[g] for g in panel]].astype(np.float64), lib[rf]).toarray()
    rt = types.to_numpy()[rf]
    M = pd.DataFrame([Lr[rt == t].mean(0) for t in TYPES], index=TYPES, columns=panel)
    nonT = ["B", "Myeloid", "Fibro", "Epi"]
    lin = (M.loc[["T_CD4", "T_CD8"]].min(0) - M.loc[nonT].max(0)).sort_values(ascending=False).index[:2].tolist()
    ms = {}
    for t in TYPES:
        others = [o for o in TYPES if o != t]
        score = (M.loc[t] - M.loc[others].max(0)).drop(lin).sort_values(ascending=False)
        ms[t] = (lin + score.index[:2].tolist()) if t in ("T_CD4", "T_CD8") else score.index[:4].tolist()
    info = dict(types=TYPES, panel=panel, panel_missing=[g for g in panel_all if g not in gidx], extra=extra,
                marker_sets=ms, t_lineage=lin,
                n_cells={h: {t: int(((half == h) & (types.to_numpy() == t)).sum()) for t in TYPES}
                         for h in ("train", "ref")},
                n_train_used={t: int((obs.Cell_type == t).sum()) for t in TYPES},
                type_map="T_CD4/T_CD8 = annotation Tcell & scGate_multi CD4T/CD8T; B = Bcell + PlasmaCell; "
                         "Myeloid = MoMacDC; Fibro = Fibroblast; Epi = Cancer_cells",
                prep_seconds=round(time.time() - t0, 1))
    (out / "prep.json").write_text(json.dumps(info, indent=1))
    print("prep done", {k: info[k] for k in ("n_cells", "marker_sets", "t_lineage")}, flush=True)


def _train_adata(out):
    import anndata as ad
    a = ad.read_h5ad(out / "train.h5ad")
    a.X = np.asarray(a.X, dtype=np.float32)
    return a


def train(out: Path, device: str, epochs: int):
    import torch
    from scCube import scCube
    torch.manual_seed(0)
    np.random.seed(0)
    random.seed(0)
    a = _train_adata(out)
    t0 = time.time()
    scCube().train_vae_and_generate_cell(sc_adata=a, celltype_key="Cell_type", cell_key="Cell",
                                         target_num={t: 1 for t in TYPES}, batch_size=512, epoch_num=epochs,
                                         lr=1e-4, hidden_size=128, save_model=True, save_path=str(out),
                                         project_name="vae", used_device=device)
    (out / "train.json").write_text(json.dumps(dict(epochs=epochs, device=device, batch_size=512, lr=1e-4,
                                                    hidden_size=128, seconds=round(time.time() - t0, 1))))


def pools(out: Path, device: str, seeds: int):
    import anndata as ad
    import torch
    from scCube import scCube
    a = _train_adata(out)
    for seed in range(seeds):
        f = out / f"pool_s{seed}.h5ad"
        if f.exists():
            continue
        torch.manual_seed(seed)
        np.random.seed(seed)
        random.seed(seed)
        perm = np.random.default_rng(1000 + seed).permutation(a.n_obs)
        meta, data = scCube().load_vae_and_generate_cell(sc_adata=a[perm].copy(), celltype_key="Cell_type",
                                                         cell_key="Cell", target_num={t: POOL for t in TYPES},
                                                         hidden_size=128, load_path=str(out / "vae.pth"),
                                                         used_device=device)
        rates = np.expm1(data.T.to_numpy(dtype=np.float64)).astype(np.float32)     # CP10k-scale rates
        ct = meta["Cell_type"].to_numpy()
        order = np.argsort(pd.Categorical(ct, categories=TYPES).codes, kind="stable")
        p = ad.AnnData(rates[order], obs=pd.DataFrame({"cell_type": ct[order]}, index=[f"g{i}" for i in range(len(ct))]),
                       var=pd.DataFrame(index=list(data.index)))
        p.write_h5ad(f)
        print("pool", seed, p.shape, flush=True)


# --------------------------------------------------------------------------- #
def _counts_per_tile(w, n_tile):
    """Largest-remainder rounding of w x n_tile, at least one cell per type."""
    raw = np.asarray(w, float) / np.sum(w) * n_tile
    c = np.floor(raw).astype(int)
    c[np.argsort(-(raw - c))[: n_tile - c.sum()]] += 1
    while (c == 0).any():
        c[np.argmax(c)] -= 1
        c[np.argmin(c)] += 1
    return c


def _tile_grf(counts, seed, delta, lamda):
    from scCube.utils import SPatternGenerator
    np.random.seed(seed)
    random.seed(seed)
    types = np.repeat(TYPES, counts)             # sim.TYPES order = band order along the field
    meta = pd.DataFrame({"Cell": [f"t{i}" for i in range(len(types))], "Cell_type": types})
    res = SPatternGenerator().run(generate_meta=meta.copy(), set_seed=True, seed=seed, celltype_key="Cell_type",
                                  spatial_dim=2, spatial_size=SPATIAL_SIZE, delta=delta, lamda=lamda)
    return res["point_x"].to_numpy(float), res["point_y"].to_numpy(float), res["Cell_type"].to_numpy(str)


def _tile_stripes(w, n_tile, seed):
    import anndata as ad
    from scCube.utils import SPatternGeneratorCustom
    dummy = ad.AnnData(np.zeros((len(TYPES), 1), dtype=np.float32), obs=pd.DataFrame({"Cell_type": TYPES},
                                                                                   index=TYPES))
    g = SPatternGeneratorCustom(sc_adata=dummy, cell_num=n_tile, celltype_key="Cell_type", set_seed=True, seed=seed,
                                spatial_size=SPATIAL_SIZE, select_celltype=TYPES)
    width = np.asarray(w, float) / np.sum(w) * SPATIAL_SIZE
    top = np.cumsum(width)                           # scCube's stripe covers [y1 - width, y1]
    top[-1] += 1e-6
    sp_ = g.simulate_stripes(y1_list=list(top), y2_list=list(top), stripe_width_list=list(width),
                             stripe_purity_list=[1.0] * len(TYPES), stripe_celltype_list=TYPES,
                             infiltration_celltype_list=[[] for _ in TYPES],
                             infiltration_prop_list=[[] for _ in TYPES], background_celltype=["Epi"])
    # scCube stripes run along x; transpose so types change along x, as in the toy world
    return sp_["point_y"].to_numpy(float), sp_["point_x"].to_numpy(float), sp_["Cell_type"].to_numpy(str)


def layout(pattern, seed, n_cells, w):
    cfg = PATTERNS[pattern]
    n_tile = n_cells // TILES
    side = np.sqrt(n_tile)                           # one cell per unit area, like the toy world
    rows = []
    for t in range(TILES):
        tseed = 100_000 * list(PATTERNS).index(pattern) + 1_000 * seed + 10 * t + (1 if n_cells != 2400 else 0)
        if cfg["kind"] == "grf":
            x, y, ct = _tile_grf(_counts_per_tile(w, n_tile), tseed, cfg["delta"], cfg["lamda"])
        else:
            x, y, ct = _tile_stripes(w, n_tile, tseed)
        col, row = t % 3, t // 3
        rows.append(pd.DataFrame({"x": x / SPATIAL_SIZE * side + col * (side + GAP),
                                  "y": y / SPATIAL_SIZE * side + row * (side + GAP), "cell_type": ct, "tile": t}))
    return pd.concat(rows, ignore_index=True)


def layouts(out: Path, seeds: int):
    t0 = time.time()
    combos = [(2400, w) for w in weight_vectors()] + [(600, np.asarray(sim.WEIGHTS, float))]
    for pattern in PATTERNS:
        d = out / "layouts" / pattern
        d.mkdir(parents=True, exist_ok=True)
        for seed in range(seeds):
            for n, w in combos:
                f = d / f"s{seed}_n{n}_{wkey(w)}.csv.gz"
                if not f.exists():
                    layout(pattern, seed, n, w).to_csv(f, index=False)
        print("layouts", pattern, round(time.time() - t0, 1), "s", flush=True)
    return round(time.time() - t0, 1)


def finalize(out: Path, seeds: int, layout_seconds):
    import anndata as ad
    import scCube
    info = json.loads((out / "prep.json").read_text())
    p = ad.read_h5ad(out / "pool_s0.h5ad")
    npan = len(info["panel"])
    load = np.asarray(p.X[:, :npan]).sum(1)
    scale = TARGET_MEDIAN / np.median(load)
    size_sd = float(np.sqrt(max(TARGET_SDLOG ** 2 - np.log(load).var(), 0.05 ** 2)))
    info.update(scale=float(scale), size_sd=size_sd, ctrl_ratio=CTRL_RATIO, n_ctrl=N_CTRL, pool_per_type=POOL,
                tiles=TILES, gap=GAP, spatial_size=SPATIAL_SIZE, patterns=PATTERNS, seeds=seeds,
                weights={wkey(w): list(map(float, w)) for w in weight_vectors()},
                area=dict(median=AREA_MEDIAN, sdlog=AREA_SDLOG, r_lib=AREA_LIB_R),
                targets=dict(median_tx=TARGET_MEDIAN, sdlog_tx=TARGET_SDLOG,
                             source="10x Xenium FFPE human breast Rep1"),
                sccube=dict(version=getattr(scCube, "__version__", "2.0.0"),
                            source="git+https://github.com/ZJUFanLab/scCube@d5b423d"),
                train=json.loads((out / "train.json").read_text()), layout_seconds=layout_seconds)
    (out / "meta.json").write_text(json.dumps(info, indent=1))
    print("meta written: scale", scale, "size_sd", size_sd, flush=True)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--yerly", required=True, help="Yerly export: counts.mtx (genes x cells), genes.txt, meta.csv")
    ap.add_argument("--out", default=str(HERE / "cache"), help="cache directory (world.py: SAM_SCCUBE_CACHE)")
    ap.add_argument("--device", default="cuda:0")
    ap.add_argument("--epochs", type=int, default=3000)
    ap.add_argument("--seeds", type=int, default=5)
    ap.add_argument("--steps", nargs="*", default=["prep", "train", "pools", "layouts"])
    a = ap.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    timings = {}
    for step in a.steps:
        t0 = time.time()
        if step == "prep" and not (out / "prep.json").exists():
            prep(Path(a.yerly), out)
        elif step == "train" and not (out / "vae.pth").exists():
            train(out, a.device, a.epochs)
        elif step == "pools":
            pools(out, a.device, a.seeds)
        elif step == "layouts":
            timings["layouts"] = layouts(out, a.seeds)
        timings[step] = round(time.time() - t0, 1)
    finalize(out, a.seeds, timings.get("layouts"))
    print("timings (s)", timings, flush=True)


if __name__ == "__main__":
    main()
