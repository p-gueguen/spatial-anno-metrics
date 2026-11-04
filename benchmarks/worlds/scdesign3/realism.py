"""Realism of the scDesign3 world: the real crop vs scDesign3 draws vs the toy world (``build.py realism``).

Same 8,714 cells for the real crop, the R draw (``scDesign3::simu_new``, fit/rdraw.mtx) and the Python port
(``world.py``, seed 0); the toy world is ``sim.simulate(0)`` (2,400 cells, its own genes and markers).
Writes ``results/realism.csv`` (summary), ``results/realism_genes.csv`` (per gene) and
``results/realism.png``.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import world  # noqa: E402

sim = world.sim
OUT = HERE / "results"


def _dense(X):
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X, dtype=float)


def _moran(X, xy, k=6):
    """Per-gene Moran's I on a row-normalised kNN graph (log1p CP10k)."""
    from sklearn.neighbors import NearestNeighbors
    idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    L = _lognorm(X)
    Z = L - L.mean(0)
    num = (Z * Z[idx].mean(1)).sum(0)
    den = (Z ** 2).sum(0)
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)


def _lognorm(X):
    lib = X.sum(1, keepdims=True)
    lib[lib == 0] = 1
    return np.log1p(X / lib * 1e4)


def _auc(score, pos):
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(pos, score)


def summarise(X, types, ms, xy=None, real=None):
    """Summary statistics of one count matrix (cells x genes, columns aligned to ``real`` if given)."""
    lib = X.sum(1)
    m, v = X.mean(0), X.var(0, ddof=1)
    zf = (X == 0).mean(0)
    ok = m > 0
    slope = np.polyfit(np.log10(m[ok]), np.log10(v[ok]), 1)[0]
    L = _lognorm(X)
    tm = pd.DataFrame(X).groupby(types).mean()
    col = {g: i for i, g in enumerate(ms["__genes__"])}
    ratios, aucs = [], []
    for t, gs in ((t, gs) for t, gs in ms.items() if t != "__genes__"):
        c = [col[g] for g in gs]
        own = tm.loc[t].iloc[c].to_numpy()
        oth = tm.drop(index=t).iloc[:, c].max(0).to_numpy()
        ratios += list(own / np.maximum(oth, 1e-3))
        aucs.append(_auc(L[:, c].mean(1), types == t))
    s = {"cells": len(X), "genes": X.shape[1], "lib_median": np.median(lib),
         "lib_iqr_over_median": (np.percentile(lib, 75) - np.percentile(lib, 25)) / np.median(lib),
         "lib_cv": lib.std() / lib.mean(), "zero_fraction": (X == 0).mean(),
         "meanvar_slope_log10": slope, "median_var_over_mean": np.median(v[ok] / m[ok]),
         "marker_own_over_max_other_median": np.median(ratios), "marker_own_over_max_other_min": np.min(ratios),
         "marker_score_auc_mean": np.mean(aucs), "marker_score_auc_min": np.min(aucs)}
    per = {"mean": m, "var": v, "zero": zf}
    for t in ("Epi", "Fibro"):           # gene-gene correlation inside a type, genes with mean >= 0.2 in it
        sel = types == t
        if sel.sum() > 50:
            per[f"corr_{t}"] = np.corrcoef(L[sel].T)
    if xy is not None:
        per["moran"] = _moran(X, xy)
        s["moran_median"] = np.nanmedian(per["moran"])
    if real is not None:
        rm, rv, rz = real["mean"], real["var"], real["zero"]
        k = (rm > 0) & (m > 0)
        s["gene_mean_r_log"] = np.corrcoef(np.log10(rm[k]), np.log10(m[k]))[0, 1]
        s["gene_var_r_log"] = np.corrcoef(np.log10(rv[k]), np.log10(v[k]))[0, 1]
        s["gene_var_ratio_median"] = np.median(v[k] / rv[k])
        s["gene_zero_r"] = np.corrcoef(rz, zf)[0, 1]
        s["gene_zero_mad"] = np.mean(np.abs(rz - zf))
        s["lib_ks"] = stats.ks_2samp(real["lib"], lib).statistic
        for t in ("Epi", "Fibro"):
            if f"corr_{t}" in real and f"corr_{t}" in per:
                sel = real["type_mean"].loc[t].to_numpy() >= 0.2
                iu = np.triu_indices(sel.sum(), 1)
                a, b = real[f"corr_{t}"][np.ix_(sel, sel)][iu], per[f"corr_{t}"][np.ix_(sel, sel)][iu]
                f = np.isfinite(a) & np.isfinite(b)
                a, b = a[f], b[f]
                s[f"genecorr_{t}_r"] = np.corrcoef(a, b)[0, 1]
                s[f"genecorr_{t}_real_mean_abs"] = np.mean(np.abs(a))
                s[f"genecorr_{t}_mean_abs"] = np.nanmean(np.abs(b))
        if "moran" in per:
            k2 = np.isfinite(real["moran"]) & np.isfinite(per["moran"])
            s["moran_r"] = np.corrcoef(real["moran"][k2], per["moran"][k2])[0, 1]
    per["lib"] = lib
    per["type_mean"] = tm
    return s, per


def main():
    import anndata as ad
    P = world.params()
    pa = ad.read_h5ad(world.CACHE / "params.h5ad")
    genes = P["genes"]
    types = P["types"]
    xy = P["xy"]
    ms = {**world.marker_sets(), "__genes__": genes}
    real = _dense(pa.layers["real"])
    rd = scipy.io.mmread(world.CACHE / "fit" / "rdraw.mtx").tocsr()
    rg = [ln.strip() for ln in open(world.CACHE / "fit" / "rdraw_genes.txt")]
    rdraw = _dense(rd).T[:, [rg.index(g) for g in genes]]
    allc = np.arange(len(types))
    py = world._nb_quantile(stats.norm.cdf(world._z_pool(0, 0).astype(float)), world._mean(allc, 1.0, 0.0),
                            P["sigma"])
    rows, pers = {}, {}
    rows["real crop"], pers["real"] = summarise(real, types, ms, xy)
    pers["real"]["type_mean"].index = pers["real"]["type_mean"].index.astype(str)
    rows["scDesign3 (R simu_new)"], pers["r"] = summarise(rdraw, types, ms, xy, pers["real"])
    rows["scDesign3 (world.py, seed 0)"], pers["py"] = summarise(py, types, ms, xy, pers["real"])
    toy = sim.simulate(0)
    tg = [g for g in toy.var_names if not toy.var.loc[g, "control"]]
    Xt = _dense(toy[:, tg].layers["counts"])
    rows["toy world (sim.py, seed 0)"], _ = summarise(Xt, toy.obs["truth"].to_numpy(),
                                                     {**sim.marker_sets(), "__genes__": tg}, toy.obsm["spatial"])
    summ = pd.DataFrame(rows)
    OUT.mkdir(exist_ok=True)
    summ.to_csv(OUT / "realism.csv", float_format="%.4g")
    pd.DataFrame({"gene": genes, "mean_real": pers["real"]["mean"], "mean_R": pers["r"]["mean"],
                  "mean_py": pers["py"]["mean"], "var_real": pers["real"]["var"], "var_R": pers["r"]["var"],
                  "var_py": pers["py"]["var"], "zero_real": pers["real"]["zero"], "zero_R": pers["r"]["zero"],
                  "zero_py": pers["py"]["zero"], "moran_real": pers["real"]["moran"],
                  "moran_R": pers["r"]["moran"], "moran_py": pers["py"]["moran"]}).to_csv(
        OUT / "realism_genes.csv", index=False, float_format="%.5g")
    _plot(pers, xy, types, genes, real, py)
    with pd.option_context("display.width", 200, "display.max_columns", 10):
        print(summ.round(3))


def _plot(pers, xy, types, genes, real, py):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 3, figsize=(15, 9.5))
    r, p, q = pers["real"], pers["py"], pers["r"]
    for lab, d, c in (("real", r, "k"), ("scDesign3 (py)", p, "tab:red"), ("scDesign3 (R)", q, "tab:blue")):
        ax[0, 0].loglog(d["mean"], d["var"], ".", ms=3, c=c, alpha=0.6, label=lab)
    lim = [1e-3, 1e3]
    ax[0, 0].plot(lim, lim, "--", c="grey", lw=0.8)
    ax[0, 0].set(xlabel="gene mean", ylabel="gene variance", title="Mean-variance (304 genes)")
    ax[0, 0].legend()
    ax[0, 1].plot(r["zero"], p["zero"], ".", c="tab:red", ms=3, label="py")
    ax[0, 1].plot(r["zero"], q["zero"], ".", c="tab:blue", ms=3, label="R")
    ax[0, 1].plot([0, 1], [0, 1], "--", c="grey", lw=0.8)
    ax[0, 1].set(xlabel="zero fraction, real", ylabel="zero fraction, simulated", title="Per-gene zero fraction")
    ax[0, 1].legend()
    bins = np.logspace(0, 3.3, 50)
    for lab, d, c in (("real", r, "k"), ("py", p, "tab:red"), ("R", q, "tab:blue")):
        ax[0, 2].hist(d["lib"], bins=bins, histtype="step", color=c, label=lab, density=True)
    ax[0, 2].set(xscale="log", xlabel="counts per cell", title="Library size")
    ax[0, 2].legend()
    ax[1, 0].plot(r["moran"], p["moran"], ".", c="tab:red", ms=3)
    ax[1, 0].plot([0, 1], [0, 1], "--", c="grey", lw=0.8)
    ax[1, 0].set(xlabel="Moran's I, real", ylabel="Moran's I, scDesign3 (py)", title="Spatial autocorrelation per gene")
    g = "ACTA2" if "ACTA2" in genes else genes[0]
    j = genes.index(g)
    for k, (lab, X) in enumerate((("real", real), ("scDesign3 (py)", py))):
        axk = ax[1, 1 + k]
        axk.scatter(xy[:, 0], xy[:, 1], s=1, c=np.log1p(X[:, j]), cmap="viridis", vmax=np.log1p(np.percentile(real[:, j], 99)))
        axk.set(title=f"{g}, {lab}", aspect=1)
        axk.invert_yaxis()
        axk.set_xticks([])
        axk.set_yticks([])
    fig.tight_layout()
    fig.savefig(OUT / "realism.png", dpi=60)


if __name__ == "__main__":
    main()
