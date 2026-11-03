"""Marker co-expression / contamination purity - the annotation-quality regime that asks
"are these cells transcriptionally PURE, or contaminated by spillover / doublets?"

Two evidence regimes (like the rest of the package):

* **Reference-free** (panel markers only, no ground truth): :func:`crisp_purity`, :func:`mecr`,
  :func:`pmp`. Quantify spillover as how often *mutually exclusive* lineage markers co-occur in a
  cell - the direct read-out of transcript diffusion / spatial doublets - using only a curated
  ``marker_sets`` dict.
* **Reference-based** (need a matched scRNA/snRNA reference): :func:`nmp`, :func:`ncp`. Use the
  reference's per-type profile / gene-pair co-expression structure to flag contamination the
  panel markers alone can't see.

All are cheap (cell x gene only) and panel-arm-independent. Source: the SpatialScribe QC Layer-3
contamination battery (github.com/p-gueguen/spatial-scribe), grounded in the spatial-QC literature
(Salas et al. Nat Methods 2025; SPLIT). See ``docs/cell-annotation-quality-metrics.md``.

Every function reads the count matrix chosen by ``layer=`` (default ``layers['counts']``, or ``X`` when it
holds raw counts; see ``_validate.get_counts``) and names it under ``"matrix"`` in its result.

Depends on: numpy, anndata. Heavy imports stay inside functions.
"""
from __future__ import annotations

from itertools import combinations

from ._validate import check_labels, get_counts, skip

__all__ = ["crisp_purity", "mecr", "pmp", "nmp", "ncp"]


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _private_detection(adata, marker_sets: dict[str, list[str]], layer):
    """``(marker_pos, pairs, matrix)``: per cell, does it detect any on-panel marker; and, lazily for each
    lineage pair ``(a, b)``, ``(da, db)`` = does it detect >=1 marker of ``a`` that is NOT a marker of ``b``
    (and vice versa). A gene listed for both lineages is evidence for neither against the other, so it never
    counts as co-detection; a nested pair (one set has no marker of its own against the other, e.g. ``T`` vs
    ``CD4 T``) is hierarchical, not exclusive, and is skipped. ``matrix`` names the matrix read."""
    import numpy as np

    mat, matrix = get_counts(adata, layer)
    panel = set(adata.var_names)
    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    sets = {ct: {g for g in mk if g in panel} for ct, mk in marker_sets.items()}
    sets = {ct: s for ct, s in sets.items() if s}
    genes = sorted(set().union(*sets.values()))
    col = {g: j for j, g in enumerate(genes)}
    sub = mat[:, [var_idx[g] for g in genes]]
    det = (sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)) > 0   # cells x marker genes

    def any_of(gs):
        return det[:, [col[g] for g in gs]].any(axis=1)

    pairs = ((any_of(sets[a] - sets[b]), any_of(sets[b] - sets[a]))
             for a, b in combinations(sets, 2) if sets[a] - sets[b] and sets[b] - sets[a])
    return det.any(axis=1), pairs, matrix


# --------------------------------------------------------------------------- #
# reference-free (panel markers only)
# --------------------------------------------------------------------------- #
def crisp_purity(adata, marker_sets: dict[str, list[str]], *, layer: str | None = "counts") -> dict:
    """CRISP purity (reference-free): a marker-positive cell is IMPURE if it detects markers from
    ``>= 2`` lineages, judging each pair on the markers the two sets do not share (a shared gene such as
    CD3E in both CD4 T and CD8 T sets is not co-detection). A marker counts as detected at >= 1 count.
    Writes per-cell ``adata.obs['crisp_impure']`` (bool) and returns ``{status, crisp_purity, n_impure,
    n_marker_positive, matrix}`` with ``crisp_purity = 1 - n_impure / n_marker_positive`` (1.0 = perfectly
    pure). No marker-positive cell -> a skip (nothing measured), not 1.0."""
    import numpy as np

    marker_pos, pairs, matrix = _private_detection(adata, marker_sets, layer)
    impure = np.zeros(adata.n_obs, dtype=bool)
    for da, db in pairs:
        impure |= da & db
    n_pos = int(marker_pos.sum())
    if n_pos == 0:
        return skip("no cell detects an on-panel marker of any lineage")
    adata.obs["crisp_impure"] = impure
    return {"status": "ok", "crisp_purity": float(1.0 - impure.sum() / n_pos), "n_impure": int(impure.sum()),
            "n_marker_positive": n_pos, "matrix": matrix}


def mecr(adata, marker_sets: dict[str, list[str]], *, layer: str | None = "counts") -> dict:
    """MECR (Mutually Exclusive Co-expression Rate, reference-free): mean over lineage pairs of
    ``#(both detected) / #(either detected)``, each lineage detected on the markers the pair does not share;
    nested pairs (``T`` vs ``CD4 T``) are skipped. 0 = mutually-exclusive markers never co-occur (clean);
    higher = more spillover / doublets. Uses the curated ``marker_sets`` pairs (contrast :func:`ncp`,
    which derives its pairs from a reference). Returns ``{status, mecr, n_pairs, matrix}`` (``n_pairs`` =
    pairs with >= 1 detecting cell, the ones averaged); no such pair -> a skip, not 0.0."""
    import numpy as np

    _, pairs, matrix = _private_detection(adata, marker_sets, layer)
    n_pairs, rates = 0, []
    for da, db in pairs:
        n_pairs += 1
        either = (da | db).sum()
        if either:
            rates.append((da & db).sum() / either)
    if not n_pairs:
        return skip("no pair of lineages with private on-panel markers")
    if not rates:
        return skip(f"none of the {n_pairs} lineage pairs is detected in any cell")
    return {"status": "ok", "mecr": float(np.mean(rates)), "n_pairs": len(rates), "matrix": matrix}


def pmp(adata, marker_sets: dict[str, list[str]], label_key: str = "cell_type", *,
        layer: str | None = "counts") -> dict:
    """Per-cell Marker Purity (reference-free, needs an assigned label): of a cell's *lineage-marker*
    transcripts, the fraction from its ASSIGNED type's on-panel markers. Writes ``obs['pmp']`` in [0,1].

    **Panel-size invariant.** The denominator is the cell's counts over the UNION of all lineage
    markers (its "marker transcriptome"), NOT its total panel counts - so a 2-7 gene marker set is not
    structurally ~0 on a 5K / WTA panel. Low PMP => the marker signal is genuinely dominated by OTHER
    lineages (impure / doublet / mislabeled). ``NaN`` when the cell has no lineage-marker transcripts
    at all (purity undefined; callers should treat NaN as neutral). Returns ``{status, mean_pmp, n_defined,
    matrix}`` (``mean_pmp`` over the defined cells); no cell with a defined PMP -> a skip."""
    import numpy as np

    labels = check_labels(adata, label_key)
    panel = set(adata.var_names)
    mat, matrix = get_counts(adata, layer)
    var_idx = {g: i for i, g in enumerate(adata.var_names)}

    def _colsum(idx, rowmask=None):
        m = mat if rowmask is None else mat[rowmask]
        sub = m[:, idx]
        arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        return np.asarray(arr).sum(1).ravel().astype(float)

    all_idx = sorted({var_idx[g] for genes in marker_sets.values() for g in genes if g in panel})
    denom = _colsum(all_idx) if all_idx else np.zeros(adata.n_obs)

    out = np.full(adata.n_obs, np.nan, dtype=float)   # NaN = purity undefined (no marker signal)
    for ct in np.unique(labels):
        gi = [var_idx[g] for g in marker_sets.get(ct, []) if g in panel]
        mask = labels == ct
        if not gi or not mask.any():
            continue
        num = _colsum(gi, rowmask=mask)
        d = denom[mask]
        with np.errstate(divide="ignore", invalid="ignore"):
            out[mask] = np.where(d > 0, np.clip(num / d, 0, 1), np.nan)
    n_defined = int(np.isfinite(out).sum())
    if not n_defined:
        return skip("PMP is undefined for every cell: no label has an on-panel marker set, or no cell of such "
                    "a label detects a lineage marker")
    adata.obs["pmp"] = out
    return {"status": "ok", "mean_pmp": float(np.nanmean(out)), "n_defined": n_defined, "matrix": matrix}


# --------------------------------------------------------------------------- #
# reference-based (need a matched scRNA/snRNA reference)
# --------------------------------------------------------------------------- #
def nmp(adata, reference=None, label_key: str = "cell_type", ref_label_key: str = "cell_type",
        max_neg_frac: float = 0.05, min_other_frac: float = 0.25, *, layer: str | None = "counts") -> dict:
    """Negative-marker proportion (reference-based): mean fraction of each cell's counts on genes its
    assigned type should NOT express. A gene shared with the reference is a negative for type t when the
    reference detects it in ``<= max_neg_frac`` of t's cells and ``>= min_other_frac`` of some other type's:
    another lineage's marker that t lacks, i.e. where that lineage's spillover lands (the rarity rule of
    SegTraQ's ``markers_from_reference``, without its DE step). High nmp flags contamination / mislabeling.
    The reference-based counterpart to :func:`pmp`. Reference-guarded (no-ops -> ``{'status':'skipped'}``
    without a reference). Writes ``obs['nmp']``; returns ``{status, mean_nmp, types_without_negatives, matrix}``
    - cells of the listed types (absent from the reference, or no gene passes the rule) score 0 by
    construction, not because they are clean. The reference's own ``X`` is read for detection only."""
    if reference is None:
        return {"status": "skipped", "reason": "no reference"}
    import numpy as np

    shared = [g for g in adata.var_names if g in set(reference.var_names)]
    if not shared:
        return {"status": "skipped", "reason": "no shared genes with reference"}

    labels = check_labels(adata, label_key)
    ref_labels = check_labels(reference, ref_label_key)

    ref_idx = {g: i for i, g in enumerate(reference.var_names)}
    ref_det = reference.X[:, [ref_idx[g] for g in shared]] > 0
    det = {ct: np.asarray(ref_det[ref_labels == ct].mean(axis=0)).ravel() for ct in np.unique(ref_labels)}

    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    mat, matrix = get_counts(adata, layer)
    sub = mat[:, [var_idx[g] for g in shared]]
    total = np.asarray(sub.sum(axis=1)).ravel().astype(float)

    out = np.zeros(adata.n_obs, dtype=float)
    no_neg = []
    for ct in np.unique(labels):
        others = [d for o, d in det.items() if o != ct]
        neg = (np.flatnonzero((det[ct] <= max_neg_frac) & (np.max(others, axis=0) >= min_other_frac))
               if ct in det and others else [])
        if not len(neg):
            no_neg.append(str(ct))  # no negative-set evidence -> leave at 0, disclosed
            continue
        cell_mask = labels == ct
        neg_counts = np.asarray(sub[cell_mask][:, neg].sum(axis=1)).ravel().astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[cell_mask] = np.where(total[cell_mask] > 0, neg_counts / total[cell_mask], 0.0)

    adata.obs["nmp"] = out
    return {"status": "ok", "mean_nmp": float(np.mean(out)), "types_without_negatives": no_neg, "matrix": matrix}


def ncp(adata, reference=None, coexpr_threshold: float = 0.1, max_genes: int = 400, max_pairs: int = 2000, *,
        layer: str | None = "counts") -> dict:
    """Non-coexpression preservation (reference-based, MECR-style): fraction of reference-defined
    non-coexpressed gene pairs that STAY non-coexpressed in the section - a spillover/diffusion sanity
    check complementary to :func:`mecr` (which uses curated lineage pairs instead of reference-derived
    ones). A shared-gene pair is "non-coexpressed" when its co-detection rate ``#(both>0)/#(either>0)``
    is ``<= coexpr_threshold``; ``ncp`` is the fraction of reference-flagged pairs still below it in the
    section. Reference-guarded. To stay tractable on WTA/large panels, the shared gene set is capped to
    the top ``max_genes`` by reference mean expression and the candidate pairs to ``max_pairs`` (stable
    order, no RNG). Label-free (the reference's labels are not used). Returns ``{status, ncp, matrix}``."""
    if reference is None:
        return {"status": "skipped", "reason": "no reference"}
    import numpy as np

    shared = [g for g in adata.var_names if g in set(reference.var_names)]
    if len(shared) < 2:
        return {"status": "skipped", "reason": "fewer than 2 shared genes with reference"}

    ref_idx = {g: i for i, g in enumerate(reference.var_names)}
    ref_gi = [ref_idx[g] for g in shared]
    ref_sub = reference.X[:, ref_gi]
    ref_arr = ref_sub.toarray() if hasattr(ref_sub, "toarray") else np.asarray(ref_sub)

    if len(shared) > max_genes:
        top = np.sort(np.argsort(ref_arr.mean(axis=0))[::-1][:max_genes])
        shared = [shared[i] for i in top]
        ref_arr = ref_arr[:, top]
    ref_pos = ref_arr > 0

    n_genes = len(shared)
    pairs = []
    for i, j in combinations(range(n_genes), 2):
        either = (ref_pos[:, i] | ref_pos[:, j]).sum()
        if either == 0:
            continue
        rate = (ref_pos[:, i] & ref_pos[:, j]).sum() / either
        if rate <= coexpr_threshold:
            pairs.append((i, j))
    pairs = pairs[:max_pairs]  # deterministic cap
    if not pairs:
        return {"status": "skipped", "reason": "no non-coexpressed reference gene pairs found"}

    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    gi = [var_idx[g] for g in shared]
    mat, matrix = get_counts(adata, layer)
    sub = mat[:, gi]
    arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    pos = arr > 0

    preserved = 0
    for i, j in pairs:
        either = (pos[:, i] | pos[:, j]).sum()
        if either == 0:
            preserved += 1  # neither gene ever detected -> trivially preserved
            continue
        rate = (pos[:, i] & pos[:, j]).sum() / either
        if rate <= coexpr_threshold:
            preserved += 1
    return {"status": "ok", "ncp": float(preserved / len(pairs)), "matrix": matrix}
