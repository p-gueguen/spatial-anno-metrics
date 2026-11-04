"""Realism checks: SRTsim counts vs the real Xenium counts of the same cells (and the toy world for scale).

    python benchmarks/worlds/srtsim/realism.py [--out benchmarks/worlds/srtsim/results]

Every cached seed's full crops (all cells, genes only unless stated). Writes ``realism.csv`` (one row per check:
real, SRTsim, toy), ``realism_markers.csv`` (per marker gene: Moran's I and specificity, real vs SRTsim) and
``realism.png``.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import ks_2samp, pearsonr, spearmanr
from sklearn.neighbors import NearestNeighbors

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import world  # noqa: E402
from world import sim  # noqa: E402


def _dense(X):
    return np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)


def _knn(xy, groups, k=6):
    """Row-standardised kNN weights, neighbours searched within each group (crop / sample)."""
    rows, cols = [], []
    for g in np.unique(groups):
        ids = np.flatnonzero(groups == g)
        nn = NearestNeighbors(n_neighbors=min(k + 1, len(ids))).fit(xy[ids]).kneighbors(xy[ids],
                                                                                      return_distance=False)[:, 1:]
        rows.append(np.repeat(ids, nn.shape[1]))
        cols.append(ids[nn.ravel()])
    r, c = np.concatenate(rows), np.concatenate(cols)
    W = sp.csr_matrix((np.ones(len(r)), (r, c)), shape=(len(xy), len(xy)))
    return sp.diags(1 / np.maximum(np.asarray(W.sum(1)).ravel(), 1)) @ W


def morans_i(X, W):
    Z = X - X.mean(0)
    den = (Z ** 2).sum(0)
    num = (Z * (W @ Z)).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)           # row-standardised: n / S0 = 1


def _corr_offdiag(X):
    C = np.corrcoef(X, rowvar=False)
    return C[np.triu_indices_from(C, 1)]


def summarise(X, types, ms, W, label):
    """One world's numbers: X = cells x genes counts (genes only)."""
    lib = X.sum(1)
    out = {"library size: median": np.median(lib), "library size: mean": lib.mean(),
           "library size: CV": lib.std() / lib.mean(), "library size: 5th pct": np.percentile(lib, 5),
           "library size: 95th pct": np.percentile(lib, 95),
           "zero fraction (cells x genes)": (X == 0).mean(),
           "genes detected per cell: median": np.median((X > 0).sum(1))}
    m, v = X.mean(0), X.var(0)
    ok = m > 0
    out["mean-variance: slope of log10 var on log10 mean"] = np.polyfit(np.log10(m[ok]), np.log10(v[ok]), 1)[0]
    out["mean-variance: median var / mean"] = np.median(v[ok] / m[ok])
    mk = [g for t in ms for g in ms[t]]
    mi = morans_i(X, W)
    out["Moran's I, all genes: median"] = np.nanmedian(mi)
    out["Moran's I, marker genes: median"] = np.nanmedian(mi[mk])
    spec = []
    for t, gs in ms.items():
        own = X[types == t][:, gs].mean()
        other = max(X[types == u][:, gs].mean() for u in ms if u != t)
        spec.append(own / max(other, 1e-9))
    out["marker specificity: median own / max other type (fold)"] = np.median(spec)
    pos = X[:, mk] > 0
    lin = {g: (t if g not in sim_tlin(ms) else "T") for t in ms for g in ms[t]}
    lin_of = np.array([lin[g] for g in mk])
    n_lin = np.array([len(set(lin_of[row])) for row in pos])
    out["cells with markers of >= 2 lineages (share of marker+ cells)"] = (n_lin >= 2).sum() / max((n_lin >= 1).sum(), 1)
    top = np.argsort(-m)[:50]
    epi = X[types == "Epi"][:, top]
    out["gene-gene |r| within Epi, top-50 genes: median"] = np.nanmedian(np.abs(_corr_offdiag(epi)))
    return {k: float(v) for k, v in out.items()}, mi


def sim_tlin(ms):
    return set(ms["T_CD4"]) & set(ms["T_CD8"])


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=str(HERE / "results"))
    a = p.parse_args(argv)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    seeds = world._seeds()
    ms = world.marker_sets()
    reals, sims, types, xys, groups, crops = [], [], [], [], [], []
    for s in seeds:
        w = world._load(s)
        g = ~w.var["control"].to_numpy()
        reals.append(_dense(w.layers["real"])[:, g])
        sims.append(_dense(w.layers["counts"])[:, g])
        types.append(w.obs["truth"].astype(str).to_numpy())
        xys.append(w.obsm["spatial"])
        groups.append(np.array([f"{s}_{c}" for c in w.obs["crop"]]))
        crops.append(w)
    genes = list(crops[0].var_names[~crops[0].var["control"].to_numpy()])
    R, S, T, XY, G = np.vstack(reals), np.vstack(sims), np.concatenate(types), np.vstack(xys), np.concatenate(groups)
    W = _knn(XY, G)
    col = {g: i for i, g in enumerate(genes)}
    msi = {t: [col[g] for g in gs] for t, gs in ms.items()}
    real, I_real = summarise(R, T, msi, W, "real")
    srt, I_sim = summarise(S, T, msi, W, "srtsim")
    toy = _toy_summary()
    rows = [dict(check=k, real=real[k], srtsim=srt[k], toy=toy.get(k, np.nan)) for k in real]
    # agreement of per-gene / per-cell statistics, real vs SRTsim (same cells, same genes)
    m_r, m_s = R.mean(0), S.mean(0)
    ok = (m_r > 0) & (m_s > 0)
    lv_r, lv_s = np.log10(R.var(0)[ok]), np.log10(S.var(0)[ok])
    z_r, z_s = (R == 0).mean(0), (S == 0).mean(0)
    lib_r, lib_s = R.sum(1), S.sum(1)
    within = []
    for t in world.TYPES:
        sel = T == t
        vr, vs = R[sel].var(0), S[sel].var(0)
        k = (vr > 0) & (vs > 0)
        within.append(pearsonr(np.log10(vr[k]), np.log10(vs[k]))[0])
    mk = [col[g] for t in ms for g in ms[t]]
    agree = {
        "per-gene log10 mean: Pearson r real vs SRTsim": pearsonr(np.log10(m_r[ok]), np.log10(m_s[ok]))[0],
        "per-gene log10 variance: Pearson r real vs SRTsim": pearsonr(lv_r, lv_s)[0],
        "per-gene log10 variance within type: mean Pearson r over types": float(np.mean(within)),
        "per-gene zero fraction: mean |real - SRTsim|": np.abs(z_r - z_s).mean(),
        "per-gene zero fraction: Pearson r": pearsonr(z_r, z_s)[0],
        "per-cell library size: Spearman rho real vs SRTsim": spearmanr(lib_r, lib_s)[0],
        "library size distribution: KS statistic": ks_2samp(lib_r, lib_s).statistic,
        "Moran's I per gene: Pearson r real vs SRTsim": pearsonr(*_finite(I_real, I_sim))[0],
        "Moran's I per marker gene: Pearson r real vs SRTsim": pearsonr(*_finite(I_real[mk], I_sim[mk]))[0],
        "Moran's I per gene: median SRTsim / real": float(np.nanmedian(I_sim / np.where(I_real > 0, I_real, np.nan))),
        "gene-gene r within Epi (top-50): Pearson r of the two correlation matrices":
            pearsonr(_corr_offdiag(R[T == "Epi"][:, np.argsort(-m_r)[:50]]),
                     _corr_offdiag(S[T == "Epi"][:, np.argsort(-m_r)[:50]]))[0],
    }
    rows += [dict(check=k, real=np.nan, srtsim=float(v), toy=np.nan) for k, v in agree.items()]
    cr = sum(_dense(w.layers["real"])[:, w.var["control"].to_numpy()].sum() for w in crops)
    cs = sum(_dense(w.layers["counts"])[:, w.var["control"].to_numpy()].sum() for w in crops)
    rows.append(dict(check="negative-control probe counts per 1,000 cells (28 probes)", real=1000 * cr / len(R),
                     srtsim=1000 * cs / len(S), toy=toy.get("ctrl_per_1000", np.nan)))
    rows.append(dict(check="cells (all cached seeds)", real=len(R), srtsim=len(S), toy=np.nan))
    df = pd.DataFrame(rows)
    df.to_csv(out / "realism.csv", index=False, float_format="%.4g")
    mrows = []
    for t, gs in ms.items():
        for g in gs:
            j = col[g]
            row = dict(type=t, gene=g, moran_real=I_real[j], moran_srtsim=I_sim[j])
            for name, X in (("real", R), ("srtsim", S)):
                own = X[T == t, j].mean()
                row[f"own_mean_{name}"] = own
                row[f"fold_vs_max_other_{name}"] = own / max(max(X[T == u, j].mean() for u in world.TYPES if u != t),
                                                             1e-9)
            mrows.append(row)
    pd.DataFrame(mrows).to_csv(out / "realism_markers.csv", index=False, float_format="%.4g")
    _plot(out / "realism.png", R, S, I_real, I_sim, mk, lib_r, lib_s)
    with pd.option_context("display.width", 200, "display.max_colwidth", 80):
        print(df.to_string(index=False))


def _finite(x, y):
    k = np.isfinite(x) & np.isfinite(y)
    return x[k], y[k]


def _toy_summary():
    """The toy world's numbers on the same checks (3 seeds, default size)."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("toy_sim", Path(sim.__file__))
    toy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(toy)                              # a fresh copy: unaffected by world.install
    Xs, Ts, XYs, Gs = [], [], [], []
    for s in range(3):
        b = toy.simulate(s)
        g = ~b.var["control"].to_numpy()
        Xs.append(_dense(b.layers["counts"])[:, g])
        Ts.append(b.obs["truth"].to_numpy())
        XYs.append(b.obsm["spatial"])
        Gs.append(np.array([f"{s}_{x}" for x in b.obs["sample"]]))
    genes = list(b.var_names[~b.var["control"].to_numpy()])
    col = {g: i for i, g in enumerate(genes)}
    ms = {t: [col[g] for g in gs] for t, gs in toy.marker_sets().items()}
    W = _knn(np.vstack(XYs), np.concatenate(Gs))
    res, _ = summarise(np.vstack(Xs), np.concatenate(Ts), ms, W, "toy")
    ctrl = b.var["control"].to_numpy()
    C = [_dense(toy.simulate(s).layers["counts"])[:, ctrl] for s in range(3)]
    res["ctrl_per_1000"] = 1000 * sum(c.sum() for c in C) / sum(len(c) for c in C)
    return res


def _plot(path, R, S, I_real, I_sim, mk, lib_r, lib_s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
    for X, c, lab in ((R, "0.3", "real"), (S, "tab:orange", "SRTsim")):
        m, v = X.mean(0), X.var(0)
        k = m > 0
        ax[0].scatter(m[k], v[k], s=5, c=c, alpha=0.6, label=lab)
    lim = [1e-4, 1e3]
    ax[0].plot(lim, lim, "k--", lw=0.8, label="Poisson")
    ax[0].set(xscale="log", yscale="log", xlabel="mean count per cell", ylabel="variance", title="Mean-variance")
    ax[0].legend(fontsize=8)
    ax[1].scatter((R == 0).mean(0), (S == 0).mean(0), s=5, c="tab:blue")
    ax[1].plot([0, 1], [0, 1], "k--", lw=0.8)
    ax[1].set(xlabel="zero fraction, real", ylabel="zero fraction, SRTsim", title="Zero fraction per gene")
    bins = np.linspace(0, np.percentile(np.r_[lib_r, lib_s], 99.5), 60)
    ax[2].hist(lib_r, bins, alpha=0.5, color="0.3", label="real")
    ax[2].hist(lib_s, bins, alpha=0.5, color="tab:orange", label="SRTsim")
    ax[2].set(xlabel="counts per cell (genes)", ylabel="cells", title="Library size")
    ax[2].legend(fontsize=8)
    ax[3].scatter(I_real, I_sim, s=5, c="0.6", label="all genes")
    ax[3].scatter(I_real[mk], I_sim[mk], s=12, c="tab:red", label="markers")
    hi = np.nanmax(np.r_[I_real, I_sim])
    ax[3].plot([0, hi], [0, hi], "k--", lw=0.8)
    ax[3].set(xlabel="Moran's I, real", ylabel="Moran's I, SRTsim", title="Spatial autocorrelation (6-NN)")
    ax[3].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)


if __name__ == "__main__":
    main()
