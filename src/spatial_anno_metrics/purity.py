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

Depends on: numpy, anndata (X or the ``counts`` layer). Heavy imports stay inside functions.
"""
from __future__ import annotations

from itertools import combinations


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _lineage_detection(adata, marker_sets: dict[str, list[str]]) -> dict:
    """{lineage: bool array over cells} = 'cell detects >=1 on-panel marker of this lineage'."""
    import numpy as np

    panel = set(adata.var_names)
    detect: dict[str, "np.ndarray"] = {}
    for ct, mk in marker_sets.items():
        genes = [g for g in mk if g in panel]
        if not genes:
            continue
        sub = adata[:, genes].X
        arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        detect[ct] = (arr > 0).any(axis=1)
    return detect


# --------------------------------------------------------------------------- #
# reference-free (panel markers only)
# --------------------------------------------------------------------------- #
def crisp_purity(adata, marker_sets: dict[str, list[str]]) -> float:
    """CRISP purity (reference-free): a marker-positive cell is IMPURE if it detects markers from
    ``>= 2`` disjoint lineages. Writes per-cell ``adata.obs['crisp_impure']`` (bool) and returns the
    dataset purity ``1 - N_impure / N_marker_positive`` (1.0 = perfectly pure)."""
    import numpy as np

    detect = _lineage_detection(adata, marker_sets)
    if not detect:
        adata.obs["crisp_impure"] = False
        return 1.0
    D = np.column_stack(list(detect.values()))       # cells x lineages (bool)
    n_lin = D.sum(axis=1)
    impure = n_lin >= 2
    marker_pos = n_lin >= 1
    adata.obs["crisp_impure"] = impure
    return float(1.0 - impure.sum() / max(1, marker_pos.sum()))


def mecr(adata, marker_sets: dict[str, list[str]]) -> float:
    """MECR (Mutually Exclusive Co-expression Rate, reference-free): mean over disjoint lineage pairs
    of ``#(both detected) / #(either detected)``. 0 = mutually-exclusive markers never co-occur (clean);
    higher = more spillover / doublets. Uses the curated ``marker_sets`` pairs (contrast :func:`ncp`,
    which derives its pairs from a reference)."""
    import numpy as np

    detect = _lineage_detection(adata, marker_sets)
    rates = []
    for a, b in combinations(detect, 2):
        da, db = detect[a], detect[b]
        either = (da | db).sum()
        if either:
            rates.append((da & db).sum() / either)
    return float(np.mean(rates)) if rates else 0.0


def pmp(adata, marker_sets: dict[str, list[str]], label_key: str = "cell_type") -> None:
    """Per-cell Marker Purity (reference-free, needs an assigned label): of a cell's *lineage-marker*
    transcripts, the fraction from its ASSIGNED type's on-panel markers. Writes ``obs['pmp']`` in [0,1].

    **Panel-size invariant.** The denominator is the cell's counts over the UNION of all lineage
    markers (its "marker transcriptome"), NOT its total panel counts - so a 2-7 gene marker set is not
    structurally ~0 on a 5K / WTA panel. Low PMP => the marker signal is genuinely dominated by OTHER
    lineages (impure / doublet / mislabeled). ``NaN`` when the cell has no lineage-marker transcripts
    at all (purity undefined; callers should treat NaN as neutral). Uses the raw ``counts`` layer when
    present, else ``X``."""
    import numpy as np

    panel = set(adata.var_names)
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    var_idx = {g: i for i, g in enumerate(adata.var_names)}

    def _colsum(idx, rowmask=None):
        m = mat if rowmask is None else mat[rowmask]
        sub = m[:, idx]
        arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
        return np.asarray(arr).sum(1).ravel().astype(float)

    all_idx = sorted({var_idx[g] for genes in marker_sets.values() for g in genes if g in panel})
    denom = _colsum(all_idx) if all_idx else np.zeros(adata.n_obs)

    labels = adata.obs[label_key].astype(str).to_numpy()
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
    adata.obs["pmp"] = out


# --------------------------------------------------------------------------- #
# reference-based (need a matched scRNA/snRNA reference)
# --------------------------------------------------------------------------- #
def nmp(adata, reference=None, label_key: str = "cell_type", ref_label_key: str = "cell_type") -> dict:
    """Negative-marker proportion (reference-based): mean fraction of each cell's counts on genes its
    assigned type should NOT express (the bottom decile of that type's mean profile in the reference,
    over genes shared with the section panel). High nmp flags contamination / mislabeling. The
    reference-based counterpart to :func:`pmp`. Reference-guarded (no-ops -> ``{'status':'skipped'}``
    without a reference). Writes ``obs['nmp']``; returns ``{status, mean_nmp}``."""
    if reference is None:
        return {"status": "skipped", "reason": "no reference"}
    import numpy as np

    shared = [g for g in adata.var_names if g in set(reference.var_names)]
    if not shared:
        return {"status": "skipped", "reason": "no shared genes with reference"}

    labels = adata.obs[label_key].astype(str).to_numpy()
    ref_labels = reference.obs[ref_label_key].astype(str).to_numpy()

    ref_idx = {g: i for i, g in enumerate(reference.var_names)}
    ref_gi = [ref_idx[g] for g in shared]
    ref_sub = reference.X[:, ref_gi]
    ref_arr = ref_sub.toarray() if hasattr(ref_sub, "toarray") else np.asarray(ref_sub)

    var_idx = {g: i for i, g in enumerate(adata.var_names)}
    gi = [var_idx[g] for g in shared]
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
    sub = mat[:, gi]
    arr = sub.toarray() if hasattr(sub, "toarray") else np.asarray(sub)
    total = arr.sum(axis=1).astype(float)

    out = np.zeros(adata.n_obs, dtype=float)
    for ct in np.unique(labels):
        cell_mask = labels == ct
        if not cell_mask.any():
            continue
        ref_mask = ref_labels == ct
        if not ref_mask.any():
            continue  # type absent from reference -> no negative-set evidence, leave at 0
        profile = ref_arr[ref_mask].mean(axis=0)
        n_neg = max(1, int(np.ceil(profile.size * 0.1)))  # bottom decile = lowest-expressed genes
        neg_idx = np.argsort(profile)[:n_neg]
        neg_counts = arr[np.ix_(cell_mask, neg_idx)].sum(axis=1).astype(float)
        with np.errstate(divide="ignore", invalid="ignore"):
            out[cell_mask] = np.where(total[cell_mask] > 0, neg_counts / total[cell_mask], 0.0)

    out = np.nan_to_num(out, nan=0.0)
    adata.obs["nmp"] = out
    return {"status": "ok", "mean_nmp": float(np.mean(out))}


def ncp(adata, reference=None, ref_label_key: str = "cell_type",
        coexpr_threshold: float = 0.1, max_genes: int = 400, max_pairs: int = 2000) -> dict:
    """Non-coexpression preservation (reference-based, MECR-style): fraction of reference-defined
    non-coexpressed gene pairs that STAY non-coexpressed in the section - a spillover/diffusion sanity
    check complementary to :func:`mecr` (which uses curated lineage pairs instead of reference-derived
    ones). A shared-gene pair is "non-coexpressed" when its co-detection rate ``#(both>0)/#(either>0)``
    is ``<= coexpr_threshold``; ``ncp`` is the fraction of reference-flagged pairs still below it in the
    section. Reference-guarded. To stay tractable on WTA/large panels, the shared gene set is capped to
    the top ``max_genes`` by reference mean expression and the candidate pairs to ``max_pairs`` (stable
    order, no RNG). Returns ``{status, ncp}``."""
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
    mat = adata.layers["counts"] if "counts" in adata.layers else adata.X
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
    return {"status": "ok", "ncp": float(preserved / len(pairs))}
