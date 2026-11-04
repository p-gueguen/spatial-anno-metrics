"""Realism checks: the toy world, the scCube world (each pattern), real cells, and a real Xenium section.

    PYTHONPATH=src python benchmarks/worlds/sccube/realism.py --out benchmarks/worlds/sccube/results/realism.csv

Per source and seed (clean sections, 2400 cells): transcripts per cell (median, SD of log), zero fraction,
control counts per control feature relative to per gene feature, overdispersion (median var/mean over genes;
Poisson = 1), marker specificity (median over marker genes of mean count in the marker's type / mean count in
other types; for the shared T markers 'other' = non-T), off-target marker detection (share of cells with >= 1
count of another type's markers, CRISP's rule), Moran's I (kNN-6, median over panel genes and over marker genes)
and the kNN-6 same-type share.

Sources: ``toy`` (benchmarks/synthetic/sim.py); ``sccube:<pattern>``; ``real cells`` = held-out Yerly cells
(``world.reference``) at scRNA-seq depth (no spatial columns); ``xenium`` = central windows of the 10x Xenium
FFPE human breast Rep1 section (313 genes, real controls; no cell types in the file, so no marker columns).
Each source is reported at its native depth and binomially thinned to a common depth (``--depth``, 60 counts;
cells below it dropped), so expression statistics compare like for like: zero fraction, overdispersion and
marker detection all depend on depth.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2] / "src"))
sys.path.insert(0, str(HERE.parents[1] / "synthetic"))
sys.path.insert(0, str(HERE))
import sim  # noqa: E402
import world  # noqa: E402

toy_simulate = sim.simulate


def _dense(X):
    return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)


def _moran(L, xy, k=6):
    from sklearn.neighbors import NearestNeighbors
    idx = NearestNeighbors(n_neighbors=k + 1).fit(xy).kneighbors(xy, return_distance=False)[:, 1:]
    Z = L - L.mean(0)
    den = (Z ** 2).sum(0)
    num = (Z * Z[idx].mean(1)).sum(0)
    return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan), idx


def stats(X, genes, ctrl, types=None, xy=None, ms=None, lineage=()):
    X = _dense(X).astype(float)
    g = ~ctrl
    Xg = X[:, g]
    lib = Xg.sum(1)
    out = dict(n_cells=len(X), n_genes=int(g.sum()), median_tx=np.median(lib), sdlog_tx=np.log(lib[lib > 0]).std(),
               zero_fraction=(Xg == 0).mean())
    out["ctrl_to_gene"] = (X[:, ctrl].mean() / Xg.mean()) if ctrl.any() else np.nan
    mu, var = Xg.mean(0), Xg.var(0)
    out["overdispersion"] = np.median(var[mu > 0.05] / mu[mu > 0.05])
    if types is not None:                         # within each type, averaged over types: biology minus type mix
        od = []
        for t in np.unique(types):
            Xt = Xg[types == t]
            m_, v_ = Xt.mean(0), Xt.var(0)
            if (m_ > 0.05).any():
                od.append(np.median(v_[m_ > 0.05] / m_[m_ > 0.05]))
        out["overdispersion_within_type"] = np.mean(od)
    gi = {x: i for i, x in enumerate(np.asarray(genes)[g])}
    if types is not None and ms:
        ratios = []
        for t, mk in ms.items():
            for m_ in mk:
                if m_ not in gi:
                    continue
                if m_ in lineage:
                    on, off = np.isin(types, ["T_CD4", "T_CD8"]), ~np.isin(types, ["T_CD4", "T_CD8"])
                else:
                    on, off = types == t, types != t
                ratios.append(Xg[on, gi[m_]].mean() / max(Xg[off, gi[m_]].mean(), 1e-9))
        out["marker_specificity"] = np.median(ratios)
        hits = []
        for i, t in enumerate(types):
            other = [gi[m_] for u, mk in ms.items() if u != t for m_ in mk
                     if m_ in gi and m_ not in ms.get(t, []) and m_ not in lineage]
            hits.append((Xg[i, other] > 0).any())
        out["offtarget_marker_detect"] = np.mean(hits)
    if xy is not None:
        L = np.log1p(Xg / np.maximum(lib, 1)[:, None] * 100)
        mi, idx = _moran(L, xy)
        out["moran_median_all"] = np.nanmedian(mi)
        mk = [gi[m_] for v in (ms or {}).values() for m_ in v if m_ in gi]
        out["moran_median_markers"] = np.nanmedian(mi[mk]) if mk else np.nan
        if types is not None:
            out["knn_same_type"] = (types[idx] == types[:, None]).mean()
    return out


def thin_to(X, ctrl, depth, rng):
    """Cells with >= ``depth`` gene counts, binomially thinned to ``depth`` on average (controls thinned alike)."""
    X = _dense(X).astype(int)
    lib = X[:, ~ctrl].sum(1)
    keep = lib >= depth
    return rng.binomial(X[keep], (depth / lib[keep])[:, None]), keep


def plot_layouts(path, seed=0):
    """Seed-0 clean sections of the four patterns, coloured by cell type (the README figure)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = dict(zip(sim.TYPES, ["#1f77b4", "#17becf", "#9467bd", "#ff7f0e", "#2ca02c", "#d62728"]))
    fig, axes = plt.subplots(1, len(world.PATTERNS), figsize=(4.2 * len(world.PATTERNS), 3.4))
    for ax, p in zip(axes, world.PATTERNS):
        s = world.simulate(seed, pattern=p)
        xy, t = s.obsm["spatial"], s.obs["truth"].to_numpy()
        for ty in sim.TYPES:
            m = t == ty
            ax.scatter(xy[m, 0], xy[m, 1], s=2, c=col[ty], label=ty, linewidths=0)
        ax.set_title(p)
        ax.set_aspect("equal")
        ax.axis("off")
    axes[-1].legend(markerscale=5, fontsize=7, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=110)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot", default=None, help="also write the layout figure to this path")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--xenium", default=None, help="10x Xenium breast Rep1 section as h5ad (layers['counts'], var['control'], obsm['spatial']); omitted = no Xenium rows")
    ap.add_argument("--depth", type=int, default=60, help="common depth for the depth-matched block")
    ap.add_argument("--out", default=str(HERE / "results" / "realism.csv"))
    a = ap.parse_args(argv)
    items = []          # (source, seed, X, genes, ctrl, types, xy, marker sets, lineage)
    for seed in range(a.seeds):
        t = toy_simulate(seed)
        items.append(("toy", seed, t.layers["counts"], t.var_names, t.var["control"].to_numpy(),
                      t.obs["truth"].to_numpy(), t.obsm["spatial"], t.uns["marker_sets"], sim.T_LINEAGE))
        for p in world.PATTERNS:
            s = world.simulate(seed, pattern=p)
            items.append((f"sccube:{p}", seed, s.layers["counts"], s.var_names, s.var["control"].to_numpy(),
                          s.obs["truth"].to_numpy(), s.obsm["spatial"], s.uns["marker_sets"], world.t_lineage()))
        r = world.reference(seed, n_per_type=400)     # real held-out cells at scRNA-seq depth
        items.append(("real cells (Yerly ref)", seed, r.layers["counts"], r.var_names, np.zeros(r.n_vars, bool),
                      r.obs["cell_type"].to_numpy(), None, world.marker_sets(), world.t_lineage()))
    if a.xenium:
        import anndata as ad
        x = ad.read_h5ad(a.xenium)
        xy = x.obsm["spatial"] if "spatial" in x.obsm else x.obs[["x_centroid", "y_centroid"]].to_numpy()
        c = np.median(xy, 0)
        for i, half in enumerate([400, 600]):                  # centred windows, ~2,500 and ~5,000 cells
            keep = (np.abs(xy - c) < half).all(1)
            items.append((f"xenium {2 * half} um window", i, x.layers["counts"][keep], x.var_names,
                          x.var["control"].to_numpy(), None, xy[keep], None, ()))
    rows = []
    for src, seed, X, genes, ctrl, types, xy_, ms, lin in items:
        rows.append(dict(block="native depth", source=src, seed=seed, **stats(X, genes, ctrl, types, xy_, ms, lin)))
        Xt, keep = thin_to(X, ctrl, a.depth, np.random.default_rng(seed))
        rows.append(dict(block=f"thinned to {a.depth}", source=src, seed=seed,
                         **stats(Xt, genes, ctrl, None if types is None else types[keep], None, ms, lin)))
    df = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False, float_format="%.5g")
    summ = df.drop(columns="seed").groupby(["block", "source"], sort=False).mean()
    print(summ.round(3).to_string())
    summ.to_csv(Path(a.out).with_name("realism_summary.csv"), float_format="%.4g")
    if a.plot:
        plot_layouts(a.plot)


if __name__ == "__main__":
    main()
