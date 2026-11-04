"""Marker co-expression / contamination purity.

The annotation-quality regime that asks "are these cells transcriptionally PURE, or contaminated by
spillover / doublets?". Two evidence regimes, like the rest of the package:

* **Reference-free** (panel markers only, no ground truth): :func:`crisp_purity`, :func:`mecr`,
  :func:`pmp`. Quantify spillover as how often *mutually exclusive* lineage markers co-occur in a
  cell - the direct read-out of transcript diffusion / spatial doublets - using only a curated
  ``marker_sets`` dict. CRISP and MECR are reported beside a detection-matched null (markers detected
  independently at their own rates). None of them can tell spillover from genuine co-expression.
* **Reference-based** (need a matched scRNA/snRNA reference): :func:`nmp`, :func:`ncp`,
  :func:`coexpression_similarity`. The reference carries the real co-expression, so these tolerate it:
  NMP is the metric to trust where cell types genuinely share markers.

Every function reads the count matrix chosen by ``layer=`` (default ``layers['counts']``, or ``X`` when it
holds raw counts; see ``_validate.get_counts``) and names it under ``"matrix"`` in its result. Per-cell values
go to ``obs[key_added]`` (default ``sam_<name>``) and the scalar summary to ``uns['sam'][<function name>]``.
``marker_sets`` keys are compared as str, like the labels, so integer cluster ids work.

Sources: Salas et al. *Nat Methods* 2025 (NMP, NCP), txsim (``coexpression_similarity``), Hartman &
Satija 2024 (MECR), Rose et al. *bioRxiv* 2026 (CRISP); see ``docs/cell-annotation-quality-metrics.md``.
"""
from __future__ import annotations

from collections import Counter
from itertools import combinations

import numpy as np

from ._utils import _store
from ._validate import check_labels, dictionary_coverage, get_counts, matrix_name, skip
from .thresholds import value as _thr

__all__ = ["crisp_purity", "mecr", "pmp", "nmp", "ncp", "coexpression_similarity"]
_DETECT = _thr("purity", "detection_threshold")         # detected = count > _DETECT, i.e. one count
_SALAS_MAX_SHARE = _thr("purity", "nmp_salas_max_share")
_SALAS_MIN_CELLS = _thr("purity", "nmp_salas_min_cells")
_TARGET_SUM = _thr("purity", "coexpr_target_sum")


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _str_keys(marker_sets: dict) -> dict[str, list[str]]:
    """Return ``marker_sets`` with str keys and str genes, as labels come out of ``check_labels``.

    Without this, integer cluster labels with an int-keyed dict matched nothing and scored NaN silently.
    Two keys that become the same string (``0`` and ``"0"``) raise ``ValueError``.
    """
    out: dict[str, list[str]] = {}
    for k, genes in marker_sets.items():
        if str(k) in out:
            raise ValueError(f"marker_sets has both {k!r} and {str(k)!r} as keys; keys are compared as str")
        out[str(k)] = [str(g) for g in genes]
    return out


def _on_panel(adata, marker_sets: dict[str, list[str]]) -> dict[str, list[str]]:
    """Restrict each marker set to the genes in ``adata.var_names``."""
    panel = set(map(str, adata.var_names))
    return {t: [g for g in gs if g in panel] for t, gs in marker_sets.items()}


def _private_detection(adata, marker_sets: dict, layer):
    """Detect markers and list the exclusive lineage pairs.

    Returns ``(det, pair_cols, off_panel, matrix)``: ``det`` is a CSC bool matrix, cells x
    on-panel marker genes, count > ``_DETECT``. ``pair_cols`` holds, per lineage pair ``(a, b)``, the
    columns of ``a``'s markers that are not ``b``'s and vice versa. A gene listed for both lineages is
    evidence for neither against the other; a nested pair (one set has no marker of its own against the
    other, e.g. ``T`` vs ``CD4 T``) is hierarchical, not exclusive, and is left out. ``off_panel`` lists
    the lineages none of whose markers is on the panel; ``matrix`` names the matrix read.
    """
    from scipy import sparse

    mat, matrix = get_counts(adata, layer)
    sets = {ct: set(gs) for ct, gs in _on_panel(adata, _str_keys(marker_sets)).items()}
    off_panel = sorted(ct for ct, s in sets.items() if not s)
    sets = {ct: s for ct, s in sets.items() if s}
    genes = sorted(set().union(*sets.values()))
    var_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}
    sub = mat[:, [var_idx[g] for g in genes]]
    det = sparse.csc_matrix(sub > _DETECT if sparse.issparse(sub) else np.asarray(sub) > _DETECT, dtype=bool)
    det.eliminate_zeros()
    col = {g: j for j, g in enumerate(genes)}
    cols = {ct: {col[g] for g in s} for ct, s in sets.items()}
    pair_cols = [(sorted(cols[a] - cols[b]), sorted(cols[b] - cols[a])) for a, b in combinations(cols, 2)
                 if cols[a] - cols[b] and cols[b] - cols[a]]
    return det, pair_cols, off_panel, matrix


def _codetection(det, pair_cols):
    """Yield ``(da, db)`` per pair: does each cell detect a private marker of ``a`` / of ``b``."""
    for ca, cb in pair_cols:
        yield (np.asarray(det[:, ca].sum(axis=1)).ravel() > 0, np.asarray(det[:, cb].sum(axis=1)).ravel() > 0)


def _crisp(det, pair_cols):
    """Return ``(impure, n_marker_positive)`` for one detection matrix."""
    impure = np.zeros(det.shape[0], dtype=bool)
    for da, db in _codetection(det, pair_cols):
        impure |= da & db
    return impure, int((np.diff(det.tocsr().indptr) > 0).sum())


def _mecr_rates(det, pair_cols):
    """Return ``#(both) / #(either)`` for each pair detected in >= 1 cell."""
    rates = []
    for da, db in _codetection(det, pair_cols):
        either = (da | db).sum()
        if either:
            rates.append((da & db).sum() / either)
    return rates


def _shuffle_columns(det, rng):
    """Reassign each column's detections to random cells: per-gene counts kept, co-occurrence lost."""
    from scipy import sparse

    n = det.shape[0]
    rows = [rng.choice(n, k, replace=False) for k in np.diff(det.indptr)]
    idx = np.concatenate(rows) if rows else np.zeros(0, dtype=np.int64)
    out = sparse.csc_matrix((np.ones(idx.size, dtype=bool), idx, det.indptr.copy()), shape=det.shape)
    out.sort_indices()
    return out


def _ratio(a: float, b: float) -> float:
    """Return ``a / b``, NaN when ``b`` is 0."""
    return float(a / b) if b > 0 else float("nan")


def _null_mean(det, n_null: int, seed: int, stat):
    """Mean of ``stat(det)`` over ``n_null`` column-shuffled copies of ``det``; None when ``n_null`` is 0."""
    if n_null < 1:
        return None
    # ponytail: O(cells x marker genes) random draws per copy; draw only the marker columns if it gets slow
    rng = np.random.default_rng(seed)
    return float(np.mean([stat(_shuffle_columns(det, rng)) for _ in range(n_null)]))


# --------------------------------------------------------------------------- #
# reference-free (panel markers only)
# --------------------------------------------------------------------------- #
def crisp_purity(adata, marker_sets: dict[str, list[str]], *, layer: str | None = "counts",
                 n_null: int = _thr("purity", "null_draws"), seed: int = 0,
                 key_added: str = "sam_crisp_impure") -> dict:
    """CRISP purity (reference-free): the share of marker-positive cells that detect markers of ONE lineage.

    A marker-positive cell is impure if it detects markers from ``>= 2`` lineages, judging each pair on the
    markers the two sets do not share (CD3E in both a CD4 T and a CD8 T set is not co-detection); nested
    pairs (``T`` vs ``CD4 T``) are skipped. A marker counts as detected at >= 1 count. Rule and formula from
    CRISP (Rose et al. *bioRxiv* 2026); CRISP's default detection threshold is > 2 counts and it expects
    disjoint marker lists, so the values are not comparable to published CRISP scores.

    The null answers "how pure would these markers look if each were detected independently at its own
    rate?": each marker gene's detections are reassigned to random cells (``n_null`` draws). ``vs_null`` =
    observed impurity / null impurity: 0 = perfectly exclusive lineages, 1 = no more exclusive than
    independently detected markers. It carries the section's depth and detection rates, which the raw
    purity does not (a clean synthetic section reads 0.80 purity at >= 1 count, ``vs_null`` 0.22).
    Spillover, doublets and genuine co-expression all raise it: use :func:`nmp` with a reference to tell real
    co-expression from contamination.

    Parameters
    ----------
    adata
        Section with raw counts.
    marker_sets
        ``{lineage: [genes]}``; genes off the panel are ignored.
    layer
        Count matrix (see ``_validate.get_counts``).
    n_null
        Draws of the detection-matched null; 0 skips it (null keys None).
    seed
        Seed of the null draws.
    key_added
        ``obs`` column for the per-cell impure flag (bool).

    Returns
    -------
    dict
        ``{status, crisp_purity, crisp_purity_null, vs_null, n_null, n_impure, n_marker_positive,
        lineages_without_panel_markers, matrix}`` with ``crisp_purity = 1 - n_impure / n_marker_positive``
        (1.0 = perfectly pure), ``crisp_purity_null`` its mean over the null draws and ``vs_null = (1 -
        crisp_purity) / (1 - crisp_purity_null)`` (NaN when the null is perfectly pure). No marker-positive
        cell -> a skip (nothing measured), not 1.0.
    """
    det, pair_cols, off_panel, matrix = _private_detection(adata, marker_sets, layer)
    impure, n_pos = _crisp(det, pair_cols)
    if n_pos == 0:
        return {**skip("no cell detects an on-panel marker of any lineage"),
                "lineages_without_panel_markers": off_panel}

    def stat(d):
        imp, n = _crisp(d, pair_cols)
        return 1.0 - imp.sum() / n

    purity = float(1.0 - impure.sum() / n_pos)
    null = _null_mean(det, n_null, seed, stat)
    adata.obs[key_added] = impure
    out = {"status": "ok", "crisp_purity": purity, "crisp_purity_null": null,
           "vs_null": None if null is None else _ratio(1.0 - purity, 1.0 - null), "n_null": int(n_null),
           "n_impure": int(impure.sum()), "n_marker_positive": n_pos,
           "lineages_without_panel_markers": off_panel, "matrix": matrix}
    _store(adata, "crisp_purity", out, key_added)
    return out


def mecr(adata, marker_sets: dict[str, list[str]], *, layer: str | None = "counts",
         n_null: int = _thr("purity", "null_draws"), seed: int = 0) -> dict:
    """MECR (Mutually Exclusive Co-expression Rate, reference-free; Hartman & Satija 2024).

    Mean over lineage pairs of ``#(both detected) / #(either detected)``, each lineage detected on the
    markers the pair does not share; nested pairs (``T`` vs ``CD4 T``) are skipped. 0 = mutually exclusive
    markers never co-occur (clean); higher = more spillover / doublets / co-expression. Uses the curated
    ``marker_sets`` pairs (contrast :func:`ncp`, which derives its pairs from a reference). The null is the
    one of :func:`crisp_purity`: each marker gene's detections reassigned to random cells; ``vs_null`` =
    ``mecr / mecr_null`` (0 = never co-detected, 1 = as often as independently detected markers).

    Parameters
    ----------
    adata
        Section with raw counts.
    marker_sets
        ``{lineage: [genes]}``; genes off the panel are ignored.
    layer
        Count matrix (see ``_validate.get_counts``).
    n_null
        Draws of the detection-matched null; 0 skips it (null keys None).
    seed
        Seed of the null draws.

    Returns
    -------
    dict
        ``{status, mecr, mecr_null, vs_null, n_null, n_pairs, lineages_without_panel_markers, matrix}``:
        ``n_pairs`` = pairs with >= 1 detecting cell, the ones averaged. No such pair -> a skip, not 0.0.
    """
    det, pair_cols, off_panel, matrix = _private_detection(adata, marker_sets, layer)
    rates = _mecr_rates(det, pair_cols)
    if not pair_cols:
        return {**skip("no pair of lineages with private on-panel markers"),
                "lineages_without_panel_markers": off_panel}
    if not rates:
        return {**skip(f"none of the {len(pair_cols)} lineage pairs is detected in any cell"),
                "lineages_without_panel_markers": off_panel}
    value = float(np.mean(rates))
    null = _null_mean(det, n_null, seed, lambda d: np.mean(_mecr_rates(d, pair_cols)))
    return {"status": "ok", "mecr": value, "mecr_null": null,
            "vs_null": None if null is None else _ratio(value, null), "n_null": int(n_null), "n_pairs": len(rates),
            "lineages_without_panel_markers": off_panel, "matrix": matrix}


def _pmp_values(adata, marker_sets: dict[str, list[str]], label_key: str, layer: str | None):
    """Per-cell PMP as a float ndarray (NaN = undefined).

    The computation behind :func:`pmp`, without writing to ``adata``. The AQI's purity term reads it too.
    """
    labels = check_labels(adata, label_key)
    ms = _on_panel(adata, _str_keys(marker_sets))
    mat, _ = get_counts(adata, layer)
    var_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}

    def _colsum(idx, rowmask=None):
        m = mat if rowmask is None else mat[rowmask]
        sub = m[:, idx]
        arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        return np.asarray(arr).sum(1).ravel().astype(float)

    all_idx = sorted({var_idx[g] for genes in ms.values() for g in genes})
    denom = _colsum(all_idx) if all_idx else np.zeros(adata.n_obs)

    out = np.full(adata.n_obs, np.nan, dtype=float)   # NaN = purity undefined (no marker signal)
    for ct in np.unique(labels):
        gi = [var_idx[g] for g in ms.get(ct, [])]
        mask = labels == ct
        if not gi or not mask.any():
            continue
        num = _colsum(gi, rowmask=mask)
        d = denom[mask]
        with np.errstate(divide="ignore", invalid="ignore"):
            out[mask] = np.where(d > 0, np.clip(num / d, 0, 1), np.nan)
    return out


def pmp(adata, marker_sets: dict[str, list[str]], label_key: str = "cell_type", *,
        layer: str | None = "counts", key_added: str = "sam_pmp") -> dict:
    """Per-cell Marker Purity (reference-free, needs an assigned label).

    Of a cell's *lineage-marker* transcripts, the fraction from its ASSIGNED type's on-panel markers.

    **Panel-size invariant.** The denominator is the cell's counts over the UNION of all lineage markers
    (its "marker transcriptome"), NOT its total panel counts - so a 2-7 gene marker set is not structurally
    ~0 on a 5K / WTA panel. Low PMP => the marker signal is genuinely dominated by OTHER lineages (impure /
    doublet / mislabeled). ``NaN`` when the cell has no lineage-marker transcripts, or its label has no
    on-panel marker set (purity undefined; callers should treat NaN as neutral); those labels are listed.

    Parameters
    ----------
    adata
        Section with raw counts.
    marker_sets
        ``{label: [genes]}``; keys compared as str.
    label_key
        ``obs`` column with the assigned labels.
    layer
        Count matrix (see ``_validate.get_counts``).
    key_added
        ``obs`` column for the per-cell PMP in [0, 1].

    Returns
    -------
    dict
        ``{status, mean_pmp, n_defined, no_dict_frac, labels_without_markers, matrix}``: ``mean_pmp`` over the
        defined cells; ``labels_without_markers`` = labels with no on-panel marker set (their cells are NaN)
        and ``no_dict_frac`` their share of cells. No cell with a defined PMP -> a skip.
    """
    values = _pmp_values(adata, marker_sets, label_key, layer)
    cov = dictionary_coverage(check_labels(adata, label_key), _on_panel(adata, _str_keys(marker_sets)))
    n_defined = int(np.isfinite(values).sum())
    if not n_defined:
        return {**skip("PMP is undefined for every cell: no label has an on-panel marker set, or no cell of "
                       "such a label detects a lineage marker"), **cov}
    adata.obs[key_added] = values
    out = {"status": "ok", "mean_pmp": float(np.nanmean(values)), "n_defined": n_defined, **cov,
           "matrix": matrix_name(adata, layer)}
    _store(adata, "pmp", out, key_added)
    return out


# --------------------------------------------------------------------------- #
# reference-based (need a matched scRNA/snRNA reference)
# --------------------------------------------------------------------------- #
def _type_means(M, labels, types):
    """Return the types x genes matrix of mean values of ``M`` (float64) over each type's cells."""
    from scipy import sparse

    rows = {t: i for i, t in enumerate(types)}
    keep = np.flatnonzero(np.isin(labels, types))
    ind = sparse.csr_matrix((np.ones(keep.size), ([rows[t] for t in labels[keep]], keep)),
                            shape=(len(types), M.shape[0]))
    sums = ind @ (M.astype(np.float64) if sparse.issparse(M) else np.asarray(M, dtype=np.float64))
    sums = sums.toarray() if sparse.issparse(sums) else np.asarray(sums)
    return sums / np.asarray(ind.sum(axis=1)).reshape(-1, 1)


def _nmp_salas_reads(X, labels, genes, ref_X, ref_labels, ref_genes, permutations=(), *, criterion: str = "txsim"):
    """Salas et al. 2025 negative marker purity, as txsim ``negative_marker_purity_reads`` computes it.

    Re-implemented from the definitions in ``tests/goldens/nmp_txsim/SOURCE.md`` (txsim has no licence: run
    for the golden, never copied). Types in both datasets with >= ``_SALAS_MIN_CELLS`` cells in each are
    used. ``x[g, t]`` = the mean count of gene ``g`` over type ``t``'s cells divided by its sum over the types
    (a type-balanced share). ``(g, t)`` is a negative pair when the reference share is < ``_SALAS_MAX_SHARE``
    (``criterion="txsim"``) or when the reference mean is < ``_SALAS_MAX_SHARE`` of the mean over types, i.e.
    share < ``_SALAS_MAX_SHARE / len(types)`` (``criterion="moldia"``, Moldia's older reads criterion).
    NMP = ``1 - max(0, mean x_section - mean x_reference)`` over the negative pairs (NaN shares skipped);
    ``per_type`` / ``per_gene`` = mean over a type's (gene's) negative pairs of ``max(0, x_sp - x_sc)``.
    ``permutations`` (index arrays) relabel the section and give ``permuted`` NMPs.
    """
    ref_col = {g: i for i, g in enumerate(ref_genes)}
    keep = [i for i, g in enumerate(genes) if g in ref_col]
    gnames = [genes[i] for i in keep]
    Xs, Xr = X[:, keep], ref_X[:, [ref_col[g] for g in gnames]]
    labels, ref_labels = np.asarray(labels), np.asarray(ref_labels)
    n_sp, n_sc = Counter(labels.tolist()), Counter(ref_labels.tolist())
    both = set(n_sp) & set(n_sc)
    types = sorted(t for t in both if n_sp[t] >= _SALAS_MIN_CELLS and n_sc[t] >= _SALAS_MIN_CELLS)
    info = {"types_absent_from_reference": sorted(set(n_sp) - set(n_sc)),
            "types_below_min_cells": sorted(both - set(types)), "types_scored": types}
    if len(types) < 2:
        return {**info, "nmp": np.nan, "reason": f"fewer than 2 cell types with >= {_SALAS_MIN_CELLS} cells "
                                                 "in both the section and the reference"}

    def shares(M, lab):
        from scipy import sparse

        if criterion == "cells":
            det_m = M > _DETECT if sparse.issparse(M) else np.asarray(M) > _DETECT
            return _type_means(det_m, lab, types)
        m = _type_means(M, lab, types)
        with np.errstate(divide="ignore", invalid="ignore"):
            return m / m.sum(axis=0)

    x_sc = shares(Xr, ref_labels)
    cut = _SALAS_MAX_SHARE / len(types) if criterion == "moldia" else _SALAS_MAX_SHARE
    neg = np.zeros(x_sc.shape, dtype=bool)
    np.less(x_sc, cut, out=neg, where=np.isfinite(x_sc))
    if not neg.any():
        return {**info, "nmp": np.nan, "reason": "no negative (gene, type) pair in the reference"}
    base = x_sc[neg].mean()

    def overall(lab):
        x_sp = shares(Xs, lab)
        with np.errstate(invalid="ignore"):
            m = np.nanmean(x_sp[neg]) if np.isfinite(x_sp[neg]).any() else np.nan
        return (1.0 - (m - base) if m > base else 1.0), x_sp

    value, x_sp = overall(labels)
    excess = np.where(neg, np.clip(x_sp - x_sc, 0, None), np.nan)
    fin = np.isfinite(excess)
    per_type = {t: float(np.nanmean(excess[i])) for i, t in enumerate(types) if fin[i].any()}
    per_gene = {g: float(np.nanmean(excess[:, j])) for j, g in enumerate(gnames) if fin[:, j].any()}
    return {**info, "nmp": float(value), "per_type": per_type, "per_gene": per_gene,
            "n_negative_pairs": int(neg.sum()),
            "permuted": [float(overall(labels[p])[0]) for p in permutations]}


def nmp(adata, reference=None, label_key: str = "cell_type", ref_label_key: str = "cell_type",
        max_neg_frac: float = _thr("purity", "nmp_max_neg_frac"),
        min_other_frac: float = _thr("purity", "nmp_min_other_frac"), *, score: str = "salas",
        n_permutations: int = _thr("purity", "nmp_n_permutations"), seed: int = 0,
        layer: str | None = "counts", key_added: str = "sam_nmp") -> dict:
    """Negative marker purity (reference-based): reads on genes a cell's type should not express.

    The reference carries the genuine co-expression, so a gene two types really share is never a negative:
    this is the purity metric that tolerates real co-expression (the reference-free :func:`crisp_purity`,
    :func:`mecr` and :func:`pmp` cannot).

    ``score="salas"`` (default): Salas et al. *Nat Methods* 2025, as txsim / OpenProblems compute it
    (``negative_marker_purity_reads``; golden-tested). Per type, the type-balanced share of each gene's mean
    count; a (gene, type) pair is negative when the type holds < 0.5 % of the gene's reference share; NMP =
    1 - the section's mean share on the negative pairs in excess of the reference's. Scaled as in the paper:
    ``nmp = (NMP - b) / (1 - b)``, ``b`` = mean NMP over ``n_permutations`` random permutations of the section
    labels, so 0 = random labels and 1 = no more negative-marker reads than the reference. Only types with
    >= 10 cells in both datasets are scored. A per-section score: nothing is written to ``adata``. The
    reference matrix is read with the same ``layer`` rule as the section. ``score="moldia"`` runs the same
    Salas et al. reads formula with Moldia ``Xenium_benchmarking``'s older relative-to-mean negative-pair
    criterion (reference share < ``0.005 / n_types``; golden-tested). ``score="cells"`` runs the cell-fraction
    variant ``negative_marker_purity_cells`` (txsim / Moldia ``Xenium_benchmarking``; golden-tested), where
    ``f[g, t]`` is the fraction of type-``t`` cells detecting gene ``g`` and negative pairs have
    ``f_reference[g, t] < 0.005``.

    ``score="raw"`` (the v0.6 score): per cell, the fraction of its counts (over the genes shared with the
    reference) on its type's negatives, lower = cleaner. A gene is a negative for type ``t`` when the
    reference detects it in ``<= max_neg_frac`` of ``t``'s cells (SegTraQ's ``markers_from_reference``
    rarity cut) and ``>= min_other_frac`` of SOME other type's (a package choice; SegTraQ uses 0.1 plus a DE
    step, a 10-cell minimum and drops MT-/RPL/RPS genes, none of which is done here). Detection fractions
    are detecting cells / type cells, exact for dense and sparse references. Cells of a type absent from the
    reference, of a type with no negative gene, or with no counts are NaN (undefined, not clean). Writes
    ``obs[key_added]``; the reference's own ``X`` is read for detection only.

    Parameters
    ----------
    adata
        Section with raw counts and labels.
    reference
        Matched scRNA/snRNA reference with labels; None -> a skip.
    label_key, ref_label_key
        Label columns of the section and the reference.
    max_neg_frac, min_other_frac
        ``score="raw"`` negative rule (inclusive cuts).
    score
        ``"salas"``, ``"moldia"``, ``"cells"``, or ``"raw"``.
    n_permutations
        ``score="salas"`` / ``"moldia"`` / ``"cells"``: label permutations of the baseline ``b`` (>= 1).
    seed
        ``score="salas"`` / ``"moldia"`` / ``"cells"``: seed of the permutations.
    layer
        Count matrix (see ``_validate.get_counts``).
    key_added
        ``score="raw"``: ``obs`` column for the per-cell value.

    Returns
    -------
    dict
        ``score="salas"`` / ``"moldia"`` / ``"cells"``: ``{status, score, nmp, nmp_unscaled, nmp_null,
        n_permutations, per_type, per_gene, n_negative_pairs, types_scored, types_absent_from_reference,
        types_below_min_cells, matrix}``; ``per_type`` / ``per_gene`` are mean excess shares (0 = clean).
        ``nmp`` is NaN when ``b >= 1``.
        ``score="raw"``: ``{status, score, mean_nmp, n_defined, types_absent_from_reference,
        types_without_negatives, matrix}``; ``mean_nmp`` over the defined cells. Nothing scorable -> a skip
        carrying the type lists.
    """
    if score not in ("salas", "moldia", "cells", "raw"):
        raise ValueError(f"score must be 'salas', 'moldia', 'cells' or 'raw', got {score!r}")
    if reference is None:
        return skip("no reference")
    labels = check_labels(adata, label_key)
    ref_labels = check_labels(reference, ref_label_key)
    ref_names = set(map(str, reference.var_names))
    shared = [g for g in map(str, adata.var_names) if g in ref_names]
    if not shared:
        return skip("no shared genes with reference")
    mat, matrix = get_counts(adata, layer)

    if score in ("salas", "moldia", "cells"):
        if n_permutations < 1:
            raise ValueError(f"n_permutations must be >= 1, got {n_permutations}")
        ref_mat, _ = get_counts(reference, layer)
        rng = np.random.default_rng(seed)
        perms = [rng.permutation(adata.n_obs) for _ in range(n_permutations)]
        crit = "txsim" if score == "salas" else score
        r = _nmp_salas_reads(mat, labels, list(map(str, adata.var_names)), ref_mat, ref_labels,
                             list(map(str, reference.var_names)), perms, criterion=crit)
        lists = {k: r[k] for k in ("types_scored", "types_absent_from_reference", "types_below_min_cells")}
        if not np.isfinite(r["nmp"]):
            return {**skip(r["reason"]), **lists}
        b = float(np.mean(r["permuted"]))
        scaled = (r["nmp"] - b) / (1.0 - b) if b < 1 else float("nan")
        return {"status": "ok", "score": score, "nmp": float(scaled), "nmp_unscaled": r["nmp"],
                "nmp_null": b, "n_permutations": int(n_permutations), "per_type": r["per_type"],
                "per_gene": r["per_gene"], "n_negative_pairs": r["n_negative_pairs"], **lists,
                "matrix": matrix}

    ref_idx = {g: i for i, g in enumerate(map(str, reference.var_names))}
    ref_det = reference.X[:, [ref_idx[g] for g in shared]] > _DETECT
    det = {}
    for ct in np.unique(ref_labels):
        m = ref_labels == ct
        det[ct] = np.asarray(ref_det[m].sum(axis=0)).ravel() / int(m.sum())   # exact: count / n, any storage

    var_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}
    sub = mat[:, [var_idx[g] for g in shared]]
    total = np.asarray(sub.sum(axis=1)).ravel().astype(float)

    out = np.full(adata.n_obs, np.nan, dtype=float)
    absent, no_neg = [], []
    for ct in np.unique(labels):
        if ct not in det:
            absent.append(str(ct))
            continue
        others = [d for o, d in det.items() if o != ct]
        neg = (np.flatnonzero((det[ct] <= max_neg_frac) & (np.max(others, axis=0) >= min_other_frac))
               if others else [])
        if not len(neg):
            no_neg.append(str(ct))
            continue
        cell_mask = labels == ct
        neg_counts = np.asarray(sub[cell_mask][:, neg].sum(axis=1)).ravel().astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[cell_mask] = np.where(total[cell_mask] > 0, neg_counts / total[cell_mask], np.nan)

    lists = {"types_absent_from_reference": absent, "types_without_negatives": no_neg}
    n_defined = int(np.isfinite(out).sum())
    if not n_defined:
        return {**skip("no cell is scorable: every type is absent from the reference or has no negative "
                       "gene, or its cells have no counts on the shared genes"), **lists}
    adata.obs[key_added] = out
    res = {"status": "ok", "score": "raw", "mean_nmp": float(np.nanmean(out)), "n_defined": n_defined,
           **lists, "matrix": matrix}
    _store(adata, "nmp", res, key_added)
    return res


def _cooccurrence(M):
    """Return ``(both, either)``: gene x gene counts of cells detecting both / either gene of a pair."""
    from scipy import sparse

    P = sparse.csc_matrix(M > _DETECT if sparse.issparse(M) else np.asarray(M) > _DETECT, dtype=np.int64)
    both = (P.T @ P).toarray()
    n = np.diag(both)
    return both, n[:, None] + n[None, :] - both


def ncp(adata, reference=None, coexpr_threshold: float = _thr("purity", "ncp_coexpr_threshold"),
        max_genes: int = _thr("purity", "ncp_max_genes"), max_pairs: int = _thr("purity", "ncp_max_pairs"), *,
        score: str = "pairs", moldia_threshold: float = _thr("purity", "ncp_moldia_threshold"),
        layer: str | None = "counts") -> dict:
    """Non-coexpression preservation (reference-based, MECR-style; after Salas et al. 2025's NCP).

    ``score="pairs"`` (default): fraction of reference-defined non-coexpressed gene pairs that STAY
    non-coexpressed in the section - a spillover/diffusion sanity check complementary to :func:`mecr` (which
    uses curated lineage pairs instead of reference-derived ones). A shared-gene pair is "non-coexpressed"
    when its co-detection rate ``#(both>0)/#(either>0)`` is ``<= coexpr_threshold``; ``ncp`` is the fraction
    of reference-flagged pairs still below it in the section (the paper's wording). To stay tractable on
    WTA/large panels, the shared gene set is capped to the top ``max_genes`` by reference mean expression and
    the candidate pairs to ``max_pairs`` (stable order, no RNG). Label-free.

    ``score="moldia"``: Moldia ``Xenium_benchmarking``'s ``negative_marker_purity_coexpression``
    (``tests/goldens/coexpr_txsim``; golden-tested). For each ordered pair of shared genes ``(g, h)`` (in
    reference order, capped to ``max_genes``), ``K[g, h] = #(g>0 & h>0) / #(h>0)`` (1.0 when ``#(h>0) == 0``);
    negative pairs have ``K_reference[g, h] < moldia_threshold`` (default 0.05), and
    ``ncp = 1 - mean(max(0, K_section - K_reference))`` over those negative pairs, with per-gene excess means
    ``per_gene_as_conditioning`` (by column ``h``) and ``per_gene_as_target`` (by row ``g``).

    Parameters
    ----------
    adata
        Section with raw counts.
    reference
        Matched reference (its ``X`` is read for detection under ``score="pairs"``, or its ``layer`` under
        ``score="moldia"``); None -> a skip.
    coexpr_threshold
        ``score="pairs"``: co-detection rate at or below which a pair is non-coexpressed.
    max_genes, max_pairs
        Tractability caps.
    score
        ``"pairs"`` or ``"moldia"``.
    moldia_threshold
        ``score="moldia"``: conditional co-detection cut ``K_reference[g, h] < moldia_threshold`` (default 0.05).
    layer
        Count matrix of the section (see ``_validate.get_counts``).

    Returns
    -------
    dict
        ``score="pairs"``: ``{status, ncp, matrix}``.
        ``score="moldia"``: ``{status, score, ncp, n_negative_pairs, per_gene_as_conditioning,
        per_gene_as_target, matrix}``.
    """
    if score not in ("pairs", "moldia"):
        raise ValueError(f"score must be 'pairs' or 'moldia', got {score!r}")
    if reference is None:
        return skip("no reference")
    ref_names = set(map(str, reference.var_names))
    shared = [g for g in map(str, adata.var_names) if g in ref_names]
    if len(shared) < 2:
        return skip("fewer than 2 shared genes with reference")

    ref_idx = {g: i for i, g in enumerate(map(str, reference.var_names))}
    var_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}
    if score == "moldia":
        genes = [g for g in map(str, reference.var_names) if g in var_idx]
        ref_mat, _ = get_counts(reference, layer)
        ref_sub = ref_mat[:, [ref_idx[g] for g in genes]]
        if len(genes) > max_genes:
            top = np.sort(np.argsort(np.asarray(ref_sub.mean(axis=0)).ravel())[::-1][:max_genes])
            genes = [genes[i] for i in top]
            ref_sub = ref_sub[:, top]
        mat, matrix = get_counts(adata, layer)
        sec_sub = mat[:, [var_idx[g] for g in genes]]

        def _cond_coocc(M):
            both_m, _ = _cooccurrence(M)
            denom = np.diag(both_m).astype(np.float64)[None, :]
            with np.errstate(divide="ignore", invalid="ignore"):
                k = both_m.astype(np.float64) / denom
            return np.where(np.isnan(k), 1.0, k)

        K_sc, K_sp = _cond_coocc(ref_sub), _cond_coocc(sec_sub)
        neg = K_sc < moldia_threshold
        if not neg.any():
            return skip("no non-coexpressed reference gene pairs found")
        val = float(1.0 - np.mean(np.clip(K_sp[neg] - K_sc[neg], 0.0, None)))
        D = np.where(neg, K_sp - K_sc, np.nan)
        fin_col, fin_row = np.isfinite(D).any(axis=0), np.isfinite(D).any(axis=1)
        by_cond = {g: float(np.nanmean(D[:, j])) for j, g in enumerate(genes) if fin_col[j]}
        by_targ = {g: float(np.nanmean(D[i])) for i, g in enumerate(genes) if fin_row[i]}
        return {"status": "ok", "score": "moldia", "ncp": val, "n_negative_pairs": int(neg.sum()),
                "per_gene_as_conditioning": by_cond, "per_gene_as_target": by_targ, "matrix": matrix}

    ref_sub = reference.X[:, [ref_idx[g] for g in shared]]
    if len(shared) > max_genes:
        top = np.sort(np.argsort(np.asarray(ref_sub.mean(axis=0)).ravel())[::-1][:max_genes])
        shared = [shared[i] for i in top]
        ref_sub = ref_sub[:, top]

    iu = np.triu_indices(len(shared), k=1)                  # the itertools.combinations order
    both, either = (m[iu] for m in _cooccurrence(ref_sub))
    with np.errstate(divide="ignore", invalid="ignore"):
        keep = np.flatnonzero((either > 0) & (both / either <= coexpr_threshold))[:max_pairs]  # stable cap
    if not keep.size:
        return skip("no non-coexpressed reference gene pairs found")

    mat, matrix = get_counts(adata, layer)
    both, either = (m[iu][keep] for m in _cooccurrence(mat[:, [var_idx[g] for g in shared]]))
    with np.errstate(divide="ignore", invalid="ignore"):
        preserved = (either == 0) | (both / either <= coexpr_threshold)   # never detected: preserved
    return {"status": "ok", "ncp": float(preserved.sum() / keep.size), "matrix": matrix}


def _lognorm_cols(M, cols):
    """``log1p(target_sum * counts / total)`` of columns ``cols`` (preserves CSR sparsity), totals over all genes."""
    from scipy import sparse

    tot = np.asarray(M.sum(axis=1), dtype=np.float64).ravel()
    scale = np.divide(_TARGET_SUM, tot, out=np.zeros_like(tot), where=tot > 0)
    sub = M[:, cols]
    if sparse.issparse(sub):
        sub = sub.tocsr().astype(np.float64, copy=True)
        sub.data *= np.repeat(scale, np.diff(sub.indptr))
        np.log1p(sub.data, out=sub.data)
        sub.eliminate_zeros()
        return sub
    return np.log1p(np.asarray(sub, dtype=np.float64) * scale[:, None])


def _corrcoef_cols(L):
    """Gene-by-gene Pearson correlation over rows of ``L`` without densifying sparse ``cells x genes``."""
    from scipy import sparse

    if not sparse.issparse(L):
        return np.atleast_2d(np.corrcoef(L, rowvar=False))
    n, g = L.shape
    if n < 2:
        return np.full((g, g), np.nan)
    mu = np.asarray(L.mean(axis=0), dtype=np.float64).ravel()
    gram = (L.T @ L).toarray().astype(np.float64, copy=False)
    cov = gram - n * np.outer(mu, mu)
    sd = np.sqrt(np.maximum(np.diag(cov), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        C = cov / np.outer(sd, sd)
    np.clip(C, -1.0, 1.0, out=C, where=np.isfinite(C))
    np.fill_diagonal(C, np.where(sd > 0, 1.0, np.nan))
    return C


def _similarity(L_sp, L_sc):
    """Return ``1 - |rho_sc - rho_sp| / 2`` (Pearson over cells), diagonal NaN."""
    with np.errstate(divide="ignore", invalid="ignore"):
        S = 1.0 - np.abs(_corrcoef_cols(L_sc) - _corrcoef_cols(L_sp)) / 2.0
    S = np.atleast_2d(S)
    np.fill_diagonal(S, np.nan)
    return S


def coexpression_similarity(adata, reference=None, *, label_key: str | None = None,
                            ref_label_key: str = "cell_type",
                            min_cells: int = _thr("purity", "coexpr_min_cells"),
                            layer: str | None = "counts") -> dict:
    """Co-expression similarity (reference-based; txsim ``coexpression_similarity``, golden-tested).

    Do gene-gene correlations in the section match the reference's? Both datasets are log-normalised as
    ``log1p(1e4 * counts / total)`` (totals over each dataset's own genes); the genes kept are the shared ones
    detected in ``>= min_cells`` cells of each. ``rho`` = Pearson correlation over cells per dataset;
    ``S = 1 - |rho_reference - rho_section| / 2`` per gene pair (1 = same correlation, 0 = opposite sign at
    full strength). Spillover and diffusion add correlations the reference does not have. Complements
    :func:`ncp`, which only looks at the pairs the reference never co-detects.

    With ``label_key``, the score is computed within each cell type present in both (genes re-filtered within
    the type) and averaged over the types, as txsim's ``by_celltype=True``.

    Parameters
    ----------
    adata
        Section with raw counts.
    reference
        Matched scRNA/snRNA reference with raw counts (same ``layer`` rule); None -> a skip.
    label_key
        Section label column; None scores all cells together.
    ref_label_key
        Reference label column (used with ``label_key``).
    min_cells
        A gene must be detected in at least this many cells of each dataset (of each type, by type).
    layer
        Count matrix (see ``_validate.get_counts``).

    Returns
    -------
    dict
        ``{status, coexpression_similarity, per_gene, n_genes, matrix}``, plus ``per_type`` with ``label_key``
        (NaN for a type with < 2 genes passing). ``per_gene`` = mean over a gene's pairs (by type: of the
        type-averaged ``S``; NaN for a gene no type scores). No gene passing the filter -> a skip.
    """
    if reference is None:
        return skip("no reference")
    if label_key is not None:
        labels = check_labels(adata, label_key)
        ref_labels = check_labels(reference, ref_label_key)
    ref_idx = {g: i for i, g in enumerate(map(str, reference.var_names))}
    var_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}
    shared = [g for g in var_idx if g in ref_idx]
    if len(shared) < 2:
        return skip("fewer than 2 shared genes with reference")
    mat, matrix = get_counts(adata, layer)
    ref_mat, _ = get_counts(reference, layer)
    L_sp = _lognorm_cols(mat, [var_idx[g] for g in shared])
    L_sc = _lognorm_cols(ref_mat, [ref_idx[g] for g in shared])

    def expressed(L_a, L_b):
        na = np.asarray((L_a > 0).sum(axis=0)).ravel()
        nb = np.asarray((L_b > 0).sum(axis=0)).ravel()
        return (na >= min_cells) & (nb >= min_cells)

    keep = expressed(L_sp, L_sc)
    if not keep.any():
        return skip(f"no shared gene is detected in >= {min_cells} cells of both datasets (min_cells)")
    genes = [g for g, k in zip(shared, keep) if k]
    L_sp, L_sc = L_sp[:, keep], L_sc[:, keep]
    iu = np.triu_indices(len(genes), k=1)

    def upper_mean(S):
        u = S[iu]
        return float(np.nanmean(u)) if np.isfinite(u).any() else float("nan")

    out = {"status": "ok"}
    if label_key is None:
        S = _similarity(L_sp, L_sc)
        out["coexpression_similarity"] = upper_mean(S)
    else:
        total, count, per_type = np.zeros((len(genes),) * 2), np.zeros((len(genes),) * 2), {}
        for t in sorted(set(labels) & set(ref_labels)):
            a, b = labels == t, ref_labels == t
            k = expressed(L_sp[a], L_sc[b])
            S = np.full((len(genes),) * 2, np.nan)
            if k.any():
                S[np.ix_(k, k)] = _similarity(L_sp[a][:, k], L_sc[b][:, k])
            per_type[t] = upper_mean(S)
            fin = np.isfinite(S)
            total[fin] += S[fin]
            count += fin
        with np.errstate(divide="ignore", invalid="ignore"):
            S = total / count                                     # NaN where no type scored the pair
        vals = [v for v in per_type.values() if np.isfinite(v)]
        out["coexpression_similarity"] = float(np.mean(vals)) if vals else float("nan")
        out["per_type"] = per_type
    fin = np.isfinite(S).any(axis=1)
    out["per_gene"] = {g: float(np.nanmean(S[i])) if fin[i] else float("nan") for i, g in enumerate(genes)}
    out.update({"n_genes": len(genes), "matrix": matrix})
    return out
