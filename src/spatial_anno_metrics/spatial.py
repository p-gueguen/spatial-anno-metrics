"""Spatial label-quality metrics: does the cell-type label make SPATIAL sense.

The reference-free spatial regime the package was named for but only covered on the *signal* side
(``signal_qc.moran_signal`` scores the data, not the labels). These score the LABELS against the
tissue geometry:

* :func:`spatial_coherence` - per-cell fraction of spatial neighbors sharing the label, the dataset PAS
  (Percentage of Abnormal Spots) and CHAOS, both as defined by Shang & Zhou 2022.
* :func:`neighborhood_sanity` - flags cell types whose strongest neighborhood-enrichment partner is
  a DIFFERENT type (a systematic-mislabeling smell), via a squidpy permutation null.
* :func:`neighbor_contamination` - per-cell share of counts on markers of the ADJACENT cell types
  (spillover from neighbors) + a directed source x target matrix (SegTraQ ``contamination_strength``),
  with a baseline from same-type cells away from the source and a permutation p-value.

Per-cell values go to ``obs[key_added]`` (default ``sam_<name>``) and the scalar summary to
``uns['sam'][<function name>]``. A spatial graph these functions build is cached in
``obsp['spatial_connectivities']`` when none is there (see ``graph=`` in each function).

**Caveat (load-bearing):** coherence, PAS and CHAOS detect *random* label noise only. They are blind to a
compact mislabelled patch, they *rise* when two neighbouring types are merged into one label, and they fall
for correctly labelled infiltrating cells. Use them to DOWN-WEIGHT / flag, never to hard-delete. Source:
SpatialScribe QC Layer-6; grounded in Salas et al. Nat Methods 2025. Depends on: numpy, scipy; squidpy
(builds the spatial graph). See ``docs/cell-annotation-quality-metrics.md`` §3e.
"""
from __future__ import annotations

import numpy as np

from ._utils import _squidpy, _store
from ._validate import check_labels, dictionary_coverage, get_counts, skip
from .thresholds import value as _thr

__all__ = ["spatial_coherence", "neighborhood_sanity", "neighbor_contamination"]
_GRAPH_TAG = "sam_graph"          # uns['spatial_neighbors'][_GRAPH_TAG]: which graph this module cached


def _spatial_graph(adata, n_neighs: int, graph: str | None = None):
    """``(CSR connectivities, description)``; row i lists cell i's neighbours.

    ``graph=None``: a graph in ``obsp['spatial_connectivities']`` that was put there by hand (no
    ``uns['spatial_neighbors']``) is used as it is; one this module cached with the same settings is reused;
    otherwise a squidpy kNN graph of ``n_neighs`` is built. ``"delaunay"``: a squidpy Delaunay graph.
    ``"obsp"``: ``obsp['spatial_connectivities']`` as it is. A built graph is cached on ``adata`` only when
    ``obsp`` holds none, so a graph built elsewhere (another ``n_neighs``, ``moran_signal``) is never reused
    by mistake nor overwritten.
    """
    import anndata as ad
    import pandas as pd
    from scipy.sparse import csr_matrix

    if graph not in (None, "delaunay", "obsp"):
        raise ValueError(f"graph must be None, 'delaunay' or 'obsp', got {graph!r}")
    have = "spatial_connectivities" in adata.obsp
    meta = adata.uns.get("spatial_neighbors", {}) if have else {}
    if graph == "obsp" or (have and graph is None and not meta):
        if not have:
            raise ValueError("graph='obsp' needs adata.obsp['spatial_connectivities']")
        return csr_matrix(adata.obsp["spatial_connectivities"]), "obsp"
    want = "delaunay" if graph == "delaunay" else f"knn(n_neighs={n_neighs})"
    if have and meta.get(_GRAPH_TAG) == want:
        return csr_matrix(adata.obsp["spatial_connectivities"]), want
    if "spatial" not in adata.obsm:
        return None, "no obsm['spatial'] or reusable obsp['spatial_connectivities']"
    kw = {"delaunay": True} if graph == "delaunay" else {"n_neighs": n_neighs}
    target = (ad.AnnData(obs=pd.DataFrame(index=adata.obs_names), obsm={"spatial": adata.obsm["spatial"]})
              if have else adata)
    _squidpy().gr.spatial_neighbors(target, coord_type="generic", **kw)
    if target is adata:
        adata.uns["spatial_neighbors"][_GRAPH_TAG] = want
    return csr_matrix(target.obsp["spatial_connectivities"]), want


def _chaos(adata, labels: np.ndarray) -> float:
    """Shang & Zhou 2022 CHAOS (Methods; SpatialPCA ``fx_CHAOS``).

    Coordinates scaled per axis (mean 0, sd 1), then the mean over all cells of the distance to the nearest
    cell with the same label; a label held by one cell adds 0 but counts in N. Re-implemented from the
    paper's definition (no code copied). NaN without ``obsm['spatial']``.
    """
    from scipy.spatial import cKDTree

    if "spatial" not in adata.obsm:
        return float("nan")
    xy = np.asarray(adata.obsm["spatial"], dtype=float)
    sd = xy.std(axis=0, ddof=1)
    z = (xy - xy.mean(axis=0)) / np.where(sd > 0, sd, 1.0)
    tot = 0.0
    for t in np.unique(labels):
        p = z[labels == t]
        if len(p) > 1:
            tot += float(cKDTree(p).query(p, k=2)[0][:, 1].sum())
    return tot / len(labels)


def spatial_coherence(adata, label_key: str = "cell_type", n_neighs: int = _thr("spatial", "coherence_n_neighs"),
                      pas_threshold: float = _thr("spatial", "pas_threshold"), *, graph: str | None = None,
                      key_added: str = "sam_spatial_coherence") -> dict:
    """Per-cell fraction of same-label spatial neighbours, plus the dataset PAS and CHAOS.

    PAS (Percentage of Abnormal Spots; Shang & Zhou 2022, doi 10.1038/s41467-022-34879-1) is the fraction of
    cells whose same-label neighbour fraction is ``<= pas_threshold``; the defaults (10 nearest neighbours,
    0.4) are the paper's "label differs from at least six of its ten neighbours". A cell with no neighbour
    counts as normal. CHAOS (same paper) is the mean distance from each cell to its nearest same-label cell
    on per-axis standardised coordinates. Lower PAS and CHAOS = spatially smoother labels.

    These detect *random* label noise only: they are blind to a compact mislabelled patch, they rise
    (look better) when two neighbouring types are merged into one label, and they fall for correctly
    labelled infiltrating cells. Use them to down-weight or flag, never to drop cells or to rank
    annotations by correctness.

    Parameters
    ----------
    adata : anndata.AnnData
        Cells with ``obs[label_key]`` and ``obsm['spatial']`` (or a graph in ``obsp``, see ``graph``).
    label_key : str
        Label column in ``adata.obs``.
    n_neighs : int
        Size of the kNN graph built when one is needed (Shang & Zhou: 10; 0.6.x: 15).
    pas_threshold : float
        Same-label fraction at or below which a cell is abnormal (inclusive; 0.6.x used ``< 0.2``).
    graph : {None, "delaunay", "obsp"}
        ``None``: a hand-made ``obsp['spatial_connectivities']`` (no ``uns['spatial_neighbors']``) is used as
        it is and ``n_neighs`` is then ignored; otherwise a kNN graph of ``n_neighs`` (a cached graph with
        other settings is not reused). ``"delaunay"``: squidpy Delaunay. ``"obsp"``: the obsp graph as it is.
    key_added : str
        ``obs`` column for the per-cell fraction (NaN for a cell without neighbours).

    Returns
    -------
    dict
        ``{status, pas, chaos, mean_coherence, n_neighs, pas_threshold, graph}``; ``chaos`` is NaN without
        ``obsm['spatial']``; ``graph`` names the graph used. Fewer than 2 labels -> a skip.
        Writes ``obs[key_added]`` and ``uns['sam']['spatial_coherence']``, and caches the graph it builds when
        ``obsp`` holds none.
    """
    from scipy.sparse import csr_matrix

    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    A, how = _spatial_graph(adata, n_neighs, graph)
    if A is None:
        return skip(how)

    cats = {c: i for i, c in enumerate(np.unique(labels))}
    lab_idx = np.array([cats[lb] for lb in labels])
    onehot = csr_matrix((np.ones(adata.n_obs), (np.arange(adata.n_obs), lab_idx)),
                        shape=(adata.n_obs, len(cats)))
    same = np.asarray((A @ onehot)[np.arange(adata.n_obs), lab_idx]).ravel()
    deg = np.asarray(A.sum(1)).ravel()
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.where(deg > 0, same / deg, np.nan)
    adata.obs[key_added] = frac

    pas = float(np.mean(np.nan_to_num(frac, nan=1.0) <= pas_threshold))
    out = {"status": "ok", "pas": pas, "chaos": _chaos(adata, labels), "mean_coherence": float(np.nanmean(frac)),
           "n_neighs": n_neighs, "pas_threshold": pas_threshold, "graph": how}
    _store(adata, "spatial_coherence", out, key_added)
    return out


def neighborhood_sanity(adata, label_key: str = "cell_type", n_neighs: int = _thr("spatial", "sanity_n_neighs"),
                        seed: int = 0, *, graph: str | None = None) -> dict:
    """Flag cell types whose strongest spatial neighbourhood-enrichment partner is another type.

    A healthy annotation self-associates spatially. Runs squidpy's neighbourhood enrichment (permutation
    z-score) and flags any type whose STRONGEST partner is a different type - a mild systematic-mislabeling
    smell (two labels that are really one, a swap, or two types forced together). Reference-free (labels +
    coords). Runs on a copy of the label column with unused levels dropped, so the caller's column is left as
    it is.

    Parameters
    ----------
    adata : anndata.AnnData
        Cells with ``obs[label_key]`` and ``obsm['spatial']`` (or a graph in ``obsp``).
    label_key : str
        Label column in ``adata.obs``.
    n_neighs : int
        Size of the kNN graph built when one is needed.
    seed : int
        Seed of squidpy's permutations.
    graph : {None, "delaunay", "obsp"}
        Which spatial graph to use, as in :func:`spatial_coherence`.

    Returns
    -------
    dict
        ``{status, categories, self_enrichment (diagonal z), suspicious [{type, stronger_with, self_z,
        cross_z}], graph}``. Fewer than 2 labels -> a skip. Writes nothing but the graph cache
        (``obsp['spatial_connectivities']``, ``obsp['spatial_distances']``, ``uns['spatial_neighbors']``) when
        it builds a graph and ``obsp`` holds none.
    """
    import anndata as ad
    import pandas as pd

    labels = check_labels(adata, label_key)
    if len(np.unique(labels)) < 2:
        return skip(f"need >= 2 labels in adata.obs[{label_key!r}], found {len(np.unique(labels))}")
    sq = _squidpy()
    A, how = _spatial_graph(adata, n_neighs, graph)
    if A is None:
        return skip(how)
    col = adata.obs[label_key]
    col = (col.cat.remove_unused_categories() if isinstance(col.dtype, pd.CategoricalDtype)
           else col.astype("category"))
    tmp = ad.AnnData(obs=pd.DataFrame({label_key: col}))
    tmp.obsp["spatial_connectivities"] = A
    z = np.asarray(sq.gr.nhood_enrichment(tmp, cluster_key=label_key, seed=seed, show_progress_bar=False,
                                          copy=True)[0])
    cats = list(col.cat.categories)
    diag = np.diag(z) if z.ndim == 2 and z.shape[0] == z.shape[1] else np.array([])

    suspicious = []
    if z.ndim == 2 and z.shape[0] > 1:
        for i, c in enumerate(cats):
            j = int(np.argmax(z[i]))
            if j != i:
                suspicious.append({"type": str(c), "stronger_with": str(cats[j]),
                                   "self_z": float(z[i][i]), "cross_z": float(z[i][j])})
    return {"status": "ok", "categories": [str(c) for c in cats],
            "self_enrichment": [float(x) for x in diag], "suspicious": suspicious, "graph": how}


def _perm_p(v1: np.ndarray, v0: np.ndarray, n_perm: int, rng) -> float:
    """One-sided permutation p of ``mean(v1) - mean(v0) > 0``.

    The group labels are shuffled over the pooled values; ``(1 + #null >= observed) / (1 + n_perm)``. NaN
    when ``n_perm`` is 0 or a group is empty.
    """
    if n_perm <= 0 or not len(v1) or not len(v0):
        return float("nan")
    pool = np.concatenate([v1, v0])
    n1, n0, n_all, tot = len(v1), len(v0), len(pool), float(pool.sum())
    obs = v1.mean() - v0.mean()
    idx = rng.permuted(np.broadcast_to(np.arange(n_all), (n_perm, n_all)), axis=1)[:, :n1]
    s1 = pool[idx].sum(axis=1)
    null = s1 / n1 - (tot - s1) / n0
    return float((1 + np.sum(null >= obs - 1e-12)) / (1 + n_perm))


def _marker_sets(marker_sets: dict, panel) -> tuple[dict, dict]:
    """``(positive, negative)`` on-panel gene sets keyed by ``str(label)``.

    ``negative`` holds only the types given in SegTraQ's ``{"positive": [...], "negative": [...]}`` form with
    a ``"negative"`` list.
    """
    pos, neg = {}, {}
    for ct, mk in marker_sets.items():
        if isinstance(mk, dict):
            pos[str(ct)] = {g for g in mk.get("positive", []) if g in panel}
            if "negative" in mk:
                neg[str(ct)] = {g for g in mk["negative"] if g in panel}
        else:
            pos[str(ct)] = {g for g in mk if g in panel}
    return pos, neg


def neighbor_contamination(adata, marker_sets: dict, label_key: str = "cell_type",
                           n_neighs: int = _thr("spatial", "contamination_n_neighs"),
                           require_neighbor_expression: bool = True, *, graph: str | None = None,
                           n_perm: int = _thr("spatial", "contamination_n_perm"), seed: int = 0,
                           layer: str | None = "counts", key_added: str = "sam_neighbor_contamination") -> dict:
    """Spillover from ADJACENT cells, per cell and per directed (source, target) type pair.

    For a cell of type t, a gene is a local contaminant when it is a marker of a type s touching the cell and
    a negative of t; with ``require_neighbor_expression`` a touching s cell must also detect it. A gene two
    touching sources share counts once. The per-cell value is the cell's counts on its local contaminants /
    its total counts. Re-implements SegTraQ's ``contamination_strength`` (github.com/LazDaria/SegTraQ),
    vectorised and sparse over cells; row i of the graph lists cell i's neighbours, as in SegTraQ.

    t's negatives are, by default, every marker of another type that is not a marker of t. A nested pair
    (one type's set inside the other's, e.g. ``T`` in ``CD4 T``) is hierarchical, not spillover, and is
    skipped, as ``purity.mecr`` skips it. SegTraQ's explicit form ``{type: {"positive": [...], "negative":
    [...]}}`` sets t's negatives to its list (contaminant = a listed negative of t that s lists as positive).

    The raw value cannot separate spillover from real co-expression: a gene t truly expresses scores the
    same next to s as anywhere. ``baseline`` is the same share in t cells with no s neighbour, ``excess`` the
    raw value minus that baseline, and ``p_value`` a one-sided permutation test of the excess (the "touches
    s" label shuffled over the t cells). The baseline counts all of s's contaminant genes while the raw value
    counts only those a neighbour detects, so the excess errs low; spatially clustered border cells make the
    p-value optimistic when expression drifts smoothly across the tissue.

    Parameters
    ----------
    adata : anndata.AnnData
        Cells with counts, ``obs[label_key]`` and ``obsm['spatial']`` (or a graph in ``obsp``).
    marker_sets : dict
        ``{label: [genes]}`` or ``{label: {"positive": [...], "negative": [...]}}`` (``"negative"`` optional;
        both forms may be mixed). Keys are matched to labels as strings, so int keys work for int labels.
        Off-panel genes are dropped.
    label_key : str
        Label column in ``adata.obs``.
    n_neighs : int
        Size of the kNN graph built when one is needed.
    require_neighbor_expression : bool
        Count a gene only when a touching source cell detects it (SegTraQ's default).
    graph : {None, "delaunay", "obsp"}
        Which spatial graph to use, as in :func:`spatial_coherence`; ``"delaunay"`` is SegTraQ's default.
    n_perm : int
        Permutations for ``p_value`` (0 = no test, NaN p-values).
    seed : int
        Seed of the permutations.
    layer : str or None
        Count matrix (see ``_validate.get_counts``).
    key_added : str
        ``obs`` column for the per-cell value (NaN when no other marker-set type that could spill into the
        cell touches it, when its type has no marker set, or when it has no counts).

    Returns
    -------
    dict
        ``{status, mean_contamination, n_evaluable, source_target, n_source_target, baseline, n_baseline,
        excess, p_value, split_slope, split_corr, distance_corr, distance_curve, no_dict_frac,
        labels_without_markers, graph, n_perm, matrix}``. The frames are source (rows) x target (columns):
        ``source_target`` the mean per-cell contamination over the ``n_source_target`` evaluable target cells
        (NaN, not 0, where no target cell touches the source, as in SegTraQ), ``baseline`` over the
        ``n_baseline`` target cells touching no source cell. ``split_slope`` and ``split_corr`` (after Bilous
        et al. *Nat Methods* 2026, SPLIT) regress/correlate each target cell's ``(s, t)`` contaminant marker
        share against ``s``'s share of the cell's spatial neighbours (> 0 under spillover, ~0 under uniform
        co-expression). When ``obsm['spatial']`` is present, ``distance_corr`` (Spearman r) and
        ``distance_curve`` (4 quantile distance bins) track contaminant marker share vs Euclidean distance to
        the nearest ``s`` cell (CellSPA readout: < 0 and decaying with distance under spillover).
        ``no_dict_frac`` / ``labels_without_markers`` disclose the cells whose label has no on-panel marker
        set: they are never scored. Fewer than 2 labels with an on-panel marker set, or no evaluable cell,
        -> a skip (nothing written; the coverage keys are added).
    """
    import pandas as pd
    import scipy.sparse as sps
    from scipy.spatial import cKDTree
    from scipy.stats import spearmanr

    labels = check_labels(adata, label_key)
    mat, matrix = get_counts(adata, layer)
    panel = {g: i for i, g in enumerate(adata.var_names)}
    pos, neg = _marker_sets(marker_sets, panel)
    cov = dictionary_coverage(labels, pos)
    types = [t for t in np.unique(labels) if pos.get(t)]
    if len(types) < 2:
        return {**skip(f"need >= 2 labels in adata.obs[{label_key!r}] with an on-panel marker set, "
                       f"found {len(types)}"), **cov}
    A, how = _spatial_graph(adata, n_neighs, graph)
    if A is None:
        return {**skip(how), **cov}

    def pair_genes(s, t):
        if t in neg:
            return pos[s] & neg[t]
        return pos[s] - pos[t] if pos[t] - pos[s] else set()     # nested (t inside s): hierarchical, skip

    pairs = {(s, t): sorted(pair_genes(s, t)) for s in types for t in types if s != t}
    pairs = {k: g for k, g in pairs.items() if g}
    genes = sorted(set().union(*pairs.values())) if pairs else []
    col = {g: j for j, g in enumerate(genes)}
    Xm = sps.csr_matrix(mat[:, [panel[g] for g in genes]])          # cells x contaminant genes, sparse
    total = np.asarray(mat.sum(axis=1)).ravel().astype(float)
    Ac = A.tocsc()
    A_bin = (A > 0).astype(float).tocsc()
    deg = np.asarray(A_bin.sum(axis=1)).ravel()
    xy = np.asarray(adata.obsm["spatial"], dtype=float) if "spatial" in adata.obsm else None
    if xy is not None and not np.all(np.isfinite(xy)):
        xy = None
    rng = np.random.default_rng(seed)

    frames = {k: pd.DataFrame(np.nan, index=types, columns=types)
              for k in ("source_target", "n_source_target", "baseline", "n_baseline", "excess", "p_value",
                        "split_slope", "split_corr", "distance_corr")}
    curve_rows = []
    rows, cols = [], []                                # (cell, gene) local contaminants, OR-ed over sources
    for s in types:
        src = np.flatnonzero(labels == s)
        As = Ac[:, src]
        adj = np.asarray(As.sum(axis=1)).ravel() > 0   # touches an s cell
        s_nb_all = np.asarray(A_bin[:, src].sum(axis=1)).ravel()
        tree_s = cKDTree(xy[src]) if xy is not None and len(src) > 0 else None
        nb = (As @ (Xm[src] > 0).astype(float)).tocsr() if require_neighbor_expression else None
        for t in types:
            g = [col[x] for x in pairs.get((s, t), [])]
            if not g:
                continue
            tgt = np.flatnonzero((labels == t) & (total > 0))
            if not len(tgt):
                continue
            xs = Xm[tgt][:, g]
            y_all = np.asarray(xs.sum(axis=1)).ravel() / total[tgt]
            deg_t = deg[tgt]
            with np.errstate(divide="ignore", invalid="ignore"):
                share = np.where(deg_t > 0, s_nb_all[tgt] / deg_t, np.nan)
            m_sh = np.isfinite(share) & np.isfinite(y_all)
            if m_sh.sum() >= 3 and share[m_sh].std() > 0:
                sh_c = share[m_sh] - share[m_sh].mean()
                y_c = y_all[m_sh] - y_all[m_sh].mean()
                denom_sh = float(np.dot(sh_c, sh_c))
                frames["split_slope"].loc[s, t] = float(np.dot(sh_c, y_c) / denom_sh)
                y_norm = float(np.sqrt(np.dot(y_c, y_c)))
                frames["split_corr"].loc[s, t] = (float(np.dot(sh_c, y_c) / (np.sqrt(denom_sh) * y_norm))
                                                  if y_norm > 0 else 0.0)
            if tree_s is not None:
                d_st = tree_s.query(xy[tgt], k=1)[0]
                m_d = np.isfinite(d_st) & np.isfinite(y_all)
                if m_d.sum() >= 3 and d_st[m_d].std() > 0:
                    frames["distance_corr"].loc[s, t] = (
                        float(spearmanr(d_st[m_d], y_all[m_d]).statistic) if y_all[m_d].std() > 0 else 0.0)
                if m_d.sum() >= 4 and d_st[m_d].max() > d_st[m_d].min():
                    edges = np.unique(np.quantile(d_st[m_d], np.linspace(0.0, 1.0, 5)))
                    if len(edges) >= 2:
                        b_idx = np.clip(np.searchsorted(edges[1:-1], d_st[m_d], side="left"), 0, len(edges) - 2)
                        for b in range(len(edges) - 1):
                            in_b = b_idx == b
                            if in_b.any():
                                curve_rows.append({
                                    "source": s, "target": t, "bin": int(b),
                                    "r_min": float(edges[b]), "r_max": float(edges[b + 1]),
                                    "mean_distance": float(d_st[m_d][in_b].mean()),
                                    "mean_contamination": float(y_all[m_d][in_b].mean()),
                                    "n_cells": int(in_b.sum()),
                                })
            if require_neighbor_expression:           # an s neighbour detects the gene
                hit = nb[tgt][:, g] > 0
            else:                                      # any s neighbour
                hit = sps.csr_matrix(np.repeat(adj[tgt, None], len(g), axis=1))
            ok = hit.getnnz(axis=1) > 0
            if not ok.any():
                continue
            v1 = np.asarray(xs.multiply(hit).sum(axis=1)).ravel()[ok] / total[tgt][ok]
            far = ~adj[tgt]
            v0 = np.asarray(xs[far].sum(axis=1)).ravel() / total[tgt][far]
            base = float(v0.mean()) if len(v0) else float("nan")
            for k, v in (("source_target", v1.mean()), ("n_source_target", ok.sum()), ("baseline", base),
                         ("n_baseline", len(v0)), ("excess", v1.mean() - base),
                         ("p_value", _perm_p(v1, v0, n_perm, rng))):
                frames[k].loc[s, t] = v
            r, c = hit.nonzero()
            rows.append(tgt[r])
            cols.append(np.asarray(g)[c])

    if not rows:
        return {**skip("no cell touches another marker-set type that could spill into it"), **cov}
    distance_curve = pd.DataFrame(
        curve_rows,
        columns=["source", "target", "bin", "r_min", "r_max", "mean_distance", "mean_contamination", "n_cells"],
    )
    r, c = np.concatenate(rows), np.concatenate(cols)
    local = sps.csr_matrix((np.ones(len(r)), (r, c)), shape=Xm.shape) > 0   # duplicates (two sources) -> once
    evaluable = local.getnnz(axis=1) > 0
    frac = np.full(adata.n_obs, np.nan)
    frac[evaluable] = np.asarray(Xm.multiply(local).sum(axis=1)).ravel()[evaluable] / total[evaluable]
    adata.obs[key_added] = frac
    out = {"status": "ok", "mean_contamination": float(np.nanmean(frac)), "n_evaluable": int(evaluable.sum()),
           **frames, "distance_curve": distance_curve, **cov, "graph": how, "n_perm": int(n_perm),
           "matrix": matrix}
    _store(adata, "neighbor_contamination", out, key_added)
    return out
