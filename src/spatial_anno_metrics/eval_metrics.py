"""Cell-annotation quality metrics - the computable battery behind
``docs/research/cell-annotation-quality-metrics.md``.

What it does
------------
Scores how good a cell-type labeling is, across the three regimes of the catalog:

  * ``internal_validity``          - reference-free cluster validity (the scTypeEval family:
                                     silhouette, 2-label silhouette, neighborhood purity,
                                     orbital-medoid, Ward-PropMatch, average similarity, and a
                                     geometric-mean ``integrated`` score).
  * ``inter_sample_consistency``   - reference-free: is a cell type reproducible across
                                     biological replicates (scTypeEval ISC).
  * ``marker_program_fidelity``    - reference-free given marker sets: per-type AUC-ROC and
                                     Cohen's d of the marker-enrichment score (§3g).
  * ``external_scores``            - needs ground truth: balanced accuracy, macro/weighted/
                                     per-class F1, ARI/AMI/NMI, Cohen's kappa, and ECS.
  * ``element_centric_similarity`` - resolution-agnostic partition agreement (Gates & Ahn 2019).
  * ``hierarchical_accuracy``      - partial credit along a cell-type lineage hierarchy.
  * ``composition_accuracy``       - predicted vs true cell-type proportions (L1 / Pearson / JSD).
  * ``annotation_quality``         - orchestrates the reference-free battery into one headline;
                                     returns one headline dict.

How to use it
-------------
>>> from spatial_anno_metrics import eval_metrics as em
>>> em.internal_validity(adata, label_key="cell_type")          # reference-free
>>> em.external_scores(pred_labels, truth_labels)               # benchmark with truth

Depends on: numpy, pandas, scikit-learn, scipy (all in the main env). Heavy imports are kept
inside the functions so the package imports cheaply.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _dense(X):
    return np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)


def _scale01(x: float) -> float:
    """Map a silhouette-like value in [-1, 1] to [0, 1]."""
    return float((np.clip(x, -1.0, 1.0) + 1.0) / 2.0)


def _embedding(adata, embedding: str | None, seed: int = 0) -> np.ndarray:
    """Return the cell embedding: ``obsm[embedding]`` if present, else a quick PCA of ``X``."""
    if embedding and embedding in adata.obsm:
        return np.asarray(adata.obsm[embedding], dtype=float)
    from sklearn.decomposition import PCA

    X = _dense(adata.X).astype(float)
    n_comp = int(min(50, X.shape[1], max(2, X.shape[0] - 1)))
    return PCA(n_components=n_comp, random_state=seed).fit_transform(X)


# --------------------------------------------------------------------------- #
# Internal cluster-validity (scTypeEval family) - reference-free
# --------------------------------------------------------------------------- #
def internal_validity(adata, label_key: str = "cell_type", embedding: str | None = "X_pca",
                      k: int = 30, subsample: int = 2000, seed: int = 0) -> dict:
    """Reference-free cluster-validity of a labeling on an embedding (scTypeEval family).

    Returns ``{silhouette, silhouette_2label, neighborhood_purity, orbital_medoid,
    ward_propmatch, avg_similarity, integrated, n_used}``. ``integrated`` is the geometric mean
    of all six (each scaled to [0,1]) - a single "is this labeling internally coherent?" number.
    O(N^2) parts (silhouettes) run on a ``subsample`` of cells for scalability; the rest use all.
    """
    from sklearn.cluster import AgglomerativeClustering
    from sklearn.metrics import pairwise_distances, silhouette_score
    from sklearn.neighbors import NearestNeighbors

    E = _embedding(adata, embedding, seed=seed)
    y = adata.obs[label_key].astype(str).to_numpy()
    types = np.unique(y)
    n = len(y)
    if len(types) < 2:
        return {"status": "skipped: <2 cell types"}

    # neighborhood purity - fraction of k neighbors (excluding self) sharing the label; all cells.
    kk = int(min(k, n - 1))
    nn = NearestNeighbors(n_neighbors=kk + 1).fit(E)
    idx = nn.kneighbors(E, return_distance=False)[:, 1:]           # drop self
    purity = float(np.mean([(y[nbrs] == y[i]).mean() for i, nbrs in enumerate(idx)]))

    # orbital-medoid - fraction of cells whose nearest type-medoid is their own type.
    medoids = {}
    for t in types:
        pts = E[y == t]
        c = pts.mean(0)
        medoids[t] = pts[np.argmin(((pts - c) ** 2).sum(1))]      # actual cell nearest the centroid
    Mnames = list(types)
    M = np.vstack([medoids[t] for t in Mnames])
    d2 = ((E[:, None, :] - M[None, :, :]) ** 2).sum(2) if n * len(Mnames) < 5_000_000 else None
    if d2 is None:
        nearest = np.array([Mnames[np.argmin(((E[i] - M) ** 2).sum(1))] for i in range(n)])
    else:
        nearest = np.array(Mnames)[d2.argmin(1)]
    orbital = float((nearest == y).mean())

    # Ward-PropMatch - agreement between labels and an unsupervised Ward clustering.
    ward = AgglomerativeClustering(n_clusters=len(types), linkage="ward").fit_predict(E)
    prop = []
    for t in types:
        w = ward[y == t]
        prop.append(np.bincount(w, minlength=len(types)).max() / max(1, len(w)))
    ward_prop = float(np.mean(prop))

    # average similarity - mean cosine of each cell to its own-type centroid (cohesion).
    cents = {t: E[y == t].mean(0) for t in types}
    def _cos(u, v):
        nu, nv = np.linalg.norm(u), np.linalg.norm(v)
        return float(u @ v / (nu * nv)) if nu > 0 and nv > 0 else 0.0
    avg_sim = float(np.mean([_cos(E[i], cents[y[i]]) for i in range(n)]))

    # silhouettes - subsample the O(N^2) distance work.
    rng = np.random.default_rng(seed)
    sel = np.arange(n) if n <= subsample else rng.choice(n, subsample, replace=False)
    Es, ys = E[sel], y[sel]
    sil = float(silhouette_score(Es, ys)) if len(np.unique(ys)) > 1 else 0.0
    D = pairwise_distances(Es)
    np.fill_diagonal(D, np.nan)
    s2 = []
    for t in np.unique(ys):
        own = ys == t
        if own.sum() < 2 or (~own).sum() < 1:
            continue
        a_in = np.nanmean(D[np.ix_(own, own)], axis=1)             # mean within-type distance
        b_out = D[np.ix_(own, ~own)].mean(1)                       # mean distance to all others
        s2.extend((b_out - a_in) / np.maximum(a_in, b_out))
    sil_2label = float(np.mean(s2)) if s2 else 0.0

    parts = [_scale01(sil), _scale01(sil_2label), purity, orbital, ward_prop, _scale01(avg_sim)]
    integrated = float(np.exp(np.mean(np.log(np.clip(parts, 1e-9, 1.0)))))
    return {
        "silhouette": sil, "silhouette_2label": sil_2label, "neighborhood_purity": purity,
        "orbital_medoid": orbital, "ward_propmatch": ward_prop, "avg_similarity": avg_sim,
        "integrated": integrated, "n_used": int(len(sel)),
    }


# --------------------------------------------------------------------------- #
# Inter-sample consistency (ISC) - reference-free, needs biological replicates
# --------------------------------------------------------------------------- #
def inter_sample_consistency(adata, label_key: str, sample_key: str,
                             metric: str = "cosine") -> dict:
    """Is each cell type reproducible across biological samples (scTypeEval ISC)?

    Aggregate each (sample, type) to a pseudobulk profile, then take the **silhouette of those
    pseudobulks labeled by type** (in ``metric`` space): high = same-type profiles cluster tightly
    across samples and separate from other types. Returns ``{consistency, per_type, n_profiles}``.
    """
    from sklearn.metrics import silhouette_samples, silhouette_score

    X = _dense(adata.X).astype(float)
    y = adata.obs[label_key].astype(str).to_numpy()
    s = adata.obs[sample_key].astype(str).to_numpy()
    profiles, plabels = [], []
    for t in np.unique(y):
        for samp in np.unique(s):
            m = (y == t) & (s == samp)
            if m.sum() == 0:
                continue
            profiles.append(X[m].mean(0))
            plabels.append(t)
    P = np.vstack(profiles)
    plabels = np.asarray(plabels)
    if len(np.unique(plabels)) < 2 or len(plabels) <= len(np.unique(plabels)):
        return {"status": "skipped: need >=2 types with >=2 samples each"}
    overall = float(silhouette_score(P, plabels, metric=metric))
    svals = silhouette_samples(P, plabels, metric=metric)
    per_type = {t: float(svals[plabels == t].mean()) for t in np.unique(plabels)}
    return {"consistency": overall, "per_type": per_type, "n_profiles": int(len(plabels))}


# --------------------------------------------------------------------------- #
# Marker-program fidelity (reference-free, needs marker sets) - catalog §3g
# --------------------------------------------------------------------------- #
def marker_program_fidelity(adata, label_key: str, marker_sets: dict[str, list[str]]) -> dict:
    """Do the labels track the marker programs? Per type, AUC-ROC + Cohen's d of that type's
    marker-load score at ranking its own labeled cells above the rest (reference-free given the
    marker sets). Returns ``{per_type: {auc, cohens_d, n}, mean_auc, mean_cohens_d}``."""
    from sklearn.metrics import roc_auc_score

    var_index = {g: i for i, g in enumerate(map(str, adata.var_names))}
    labels = adata.obs[label_key].astype(str).to_numpy()
    X = _dense(adata.X).astype(float)
    per, aucs, ds = {}, [], []
    for t, genes in marker_sets.items():
        idx = [var_index[g] for g in genes if g in var_index]
        if not idx:
            continue
        score = X[:, idx].sum(1)                                   # marker load per cell
        y = (labels == t).astype(int)
        if y.sum() == 0 or y.sum() == len(y):
            per[t] = {"auc": None, "cohens_d": None, "n": int(y.sum())}
            continue
        try:
            auc = float(roc_auc_score(y, score))
        except Exception:
            auc = None
        pos, neg = score[y == 1], score[y == 0]
        num = (len(pos) - 1) * pos.var(ddof=1) + (len(neg) - 1) * neg.var(ddof=1)
        sd = np.sqrt(num / max(1, len(pos) + len(neg) - 2))
        d = float((pos.mean() - neg.mean()) / sd) if sd > 0 else 0.0
        per[t] = {"auc": auc, "cohens_d": d, "n": int(y.sum())}
        if auc is not None:
            aucs.append(auc)
        ds.append(d)
    return {
        "per_type": per,
        "mean_auc": float(np.mean(aucs)) if aucs else None,
        "mean_cohens_d": float(np.mean(ds)) if ds else None,
    }


# --------------------------------------------------------------------------- #
# Panel adequacy - can a target panel resolve a reference's cell types? (§4 battery A #1)
# --------------------------------------------------------------------------- #
def panel_resolvability(reference, label_key: str, panel_genes, target_depth: float = 50.0, *,
                        max_cells_per_type: int = 250, cv: int = 3, seed: int = 0) -> dict:
    """Can a target panel resolve this reference's cell types AT THE PANEL'S OWN DEPTH?

    Restrict the reference to the genes it shares with ``panel_genes``, **binomially thin each cell
    to ``target_depth`` counts** (the section's real per-cell depth), cross-validate a classifier, and
    score per-type with F1 (via :func:`external_scores`). Low per-type F1 + a dominant confuser means
    the panel cannot separate that type from a look-alike -> the reference's labels are too granular
    for the probe set (merge them).

    The thinning is load-bearing: at full scRNA depth every type separates (one-vs-rest AUC > 0.98
    even for 25 kidney types on a breast panel) - so this **supersedes**
    ``panel_check.identifiability_auc`` (one-vs-rest AUC at full depth, which inflates to ~1.0). Uses
    multiclass per-class F1, not one-vs-rest AUC, because AUC hides sibling confusion.

    IMPORTANT: this measures *separability given the panel*, not tissue relevance - it is only
    meaningful for a **tissue-matched** reference x panel pair (a wrong-tissue reference can score
    high here and still annotate nonsense).

    Returns ``{status, n_shared_genes, n_types, macro_f1, balanced_accuracy, ecs, frac_resolvable,
    frac_not_resolvable, per_type: {ct: {f1, recall, precision, confused_with, confused_frac, tier}}}``.
    """
    import scipy.sparse as sp
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import precision_recall_fscore_support
    from sklearn.model_selection import StratifiedKFold, cross_val_predict

    f1_ok, f1_weak, conf_hi = 0.50, 0.30, 0.25
    panel = {str(g) for g in panel_genes}
    shared = [g for g in reference.var_names if str(g) in panel]
    if len(shared) < 25:
        return {"status": "insufficient_overlap", "n_shared_genes": len(shared), "per_type": {}}
    sub = reference[:, shared].copy()
    xs = sub.X[:200].toarray() if sp.issparse(sub.X) else np.asarray(sub.X[:200])
    raw = bool(np.allclose(xs, np.round(xs)) and xs.max() > 1)
    if not raw and "counts" in sub.layers:
        sub.X = sub.layers["counts"]
        raw = True
    rng = np.random.default_rng(seed)
    y = sub.obs[label_key].astype(str).to_numpy()
    idx: list = []
    for ct in np.unique(y):
        pos = np.where(y == ct)[0]
        idx.extend(pos if len(pos) <= max_cells_per_type
                   else rng.choice(pos, max_cells_per_type, replace=False))
    sub = sub[np.array(sorted(idx))].copy()
    y = sub.obs[label_key].astype(str).to_numpy()
    keep = {ct for ct, n in zip(*np.unique(y, return_counts=True)) if n >= 10}
    m = np.isin(y, list(keep))
    sub, y = sub[m].copy(), y[m]
    if len(np.unique(y)) < 2:
        return {"status": "too_few_types", "n_shared_genes": len(shared), "per_type": {}}
    if raw:                                                        # thin to the panel's real depth
        A = (sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)).astype(np.int64)
        tot = A.sum(1, keepdims=True)
        p = np.clip(target_depth / np.maximum(tot, 1), 0.0, 1.0)
        sub.X = sp.csr_matrix(rng.binomial(A, np.broadcast_to(p, A.shape)).astype("float32"))
    Xn = _dense(sub.X).astype(float)                             # library-size normalize + log1p (no scanpy dep)
    lib = Xn.sum(1, keepdims=True)
    lib[lib == 0] = 1.0
    X = np.log1p(Xn / lib * 1e4)
    classes = np.unique(y)
    proba = cross_val_predict(LogisticRegression(max_iter=300, C=0.3), X, y,
                              cv=StratifiedKFold(cv, shuffle=True, random_state=seed),
                              method="predict_proba")
    pred = classes[np.argmax(proba, 1)]
    ext = external_scores(pred, y)
    prec, rec, f1, _ = precision_recall_fscore_support(y, pred, labels=classes, zero_division=0)
    per = {}
    for i, ct in enumerate(classes):
        mask = y == ct
        wrong = pred[mask][pred[mask] != ct]
        cw, cf = (None, 0.0)
        if len(wrong):
            v, c = np.unique(wrong, return_counts=True)
            j = int(c.argmax())
            cw, cf = str(v[j]), float(c[j] / mask.sum())
        tier = ("not_resolvable" if f1[i] < f1_weak or cf >= 0.30
                else "weak" if f1[i] < f1_ok or cf >= conf_hi else "resolvable")
        per[str(ct)] = {"f1": float(f1[i]), "recall": float(rec[i]), "precision": float(prec[i]),
                        "confused_with": cw, "confused_frac": round(cf, 2), "tier": tier}
    tiers = [d["tier"] for d in per.values()]
    return {"status": "ok", "n_shared_genes": len(shared), "n_types": int(len(classes)),
            "macro_f1": ext["macro_f1"], "balanced_accuracy": ext["balanced_accuracy"], "ecs": ext["ecs"],
            "frac_resolvable": float(np.mean([t == "resolvable" for t in tiers])),
            "frac_not_resolvable": float(np.mean([t == "not_resolvable" for t in tiers])),
            "per_type": per}


# --------------------------------------------------------------------------- #
# External harness (needs ground truth)
# --------------------------------------------------------------------------- #
def external_scores(pred, truth) -> dict:
    """Reference-based scores of ``pred`` vs ground-truth ``truth`` (benchmarks only).

    Returns accuracy, balanced accuracy, macro/weighted/per-class F1, macro precision/recall,
    ARI, AMI, NMI, Cohen's kappa, and ECS. Thin wrapper over ``sklearn.metrics`` + ``ecs``.
    """
    from sklearn.metrics import (accuracy_score, adjusted_mutual_info_score,
                                  adjusted_rand_score, balanced_accuracy_score, cohen_kappa_score,
                                  f1_score, normalized_mutual_info_score, precision_score,
                                  recall_score)

    pred = np.asarray(pred).astype(str)
    truth = np.asarray(truth).astype(str)
    labels = sorted(set(truth) | set(pred))
    per_class = f1_score(truth, pred, average=None, labels=labels, zero_division=0)
    return {
        "accuracy": float(accuracy_score(truth, pred)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, pred)),
        "macro_f1": float(f1_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "weighted_f1": float(f1_score(truth, pred, average="weighted", labels=labels, zero_division=0)),
        "precision_macro": float(precision_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "recall_macro": float(recall_score(truth, pred, average="macro", labels=labels, zero_division=0)),
        "per_class_f1": {c: float(v) for c, v in zip(labels, per_class)},
        "ari": float(adjusted_rand_score(truth, pred)),
        "ami": float(adjusted_mutual_info_score(truth, pred)),
        "nmi": float(normalized_mutual_info_score(truth, pred)),
        "kappa": float(cohen_kappa_score(truth, pred)),
        "ecs": element_centric_similarity(truth, pred),
    }


# --------------------------------------------------------------------------- #
# Element-Centric Similarity (Gates & Ahn 2019) - resolution-agnostic partition agreement
# --------------------------------------------------------------------------- #
def element_centric_similarity(a, b, alpha: float = 0.9) -> float:
    """ECS between two hard labelings ``a`` and ``b`` over the same N elements.

    Exact, vectorized closed form for hard partitions (the paper's affinity-matrix definition, with
    the self-loop terms cancelling): for element k in a-cluster A (size nA) and b-cluster B (size
    nB) with overlap nAB, the L1 affinity distance is
    ``nAB*|α/nA − α/nB| + (nA−nAB)*α/nA + (nB−nAB)*α/nB`` and ``ECS_k = 1 − dist/(2α)``. Returns the
    mean over elements, in [0, 1]; 1.0 iff the partitions are identical. Resolution-agnostic.
    """
    import pandas as pd

    a = pd.factorize(np.asarray(a))[0]
    b = pd.factorize(np.asarray(b))[0]
    _, inva, ca = np.unique(a, return_inverse=True, return_counts=True)
    _, invb, cb = np.unique(b, return_inverse=True, return_counts=True)
    nA = ca[inva].astype(float)
    nB = cb[invb].astype(float)
    pair = inva.astype(np.int64) * (invb.max() + 1) + invb
    _, invp, cp = np.unique(pair, return_inverse=True, return_counts=True)
    nAB = cp[invp].astype(float)
    dist = nAB * np.abs(alpha / nA - alpha / nB) + (nA - nAB) * (alpha / nA) + (nB - nAB) * (alpha / nB)
    return float(np.mean(1.0 - dist / (2.0 * alpha)))


# --------------------------------------------------------------------------- #
# Hierarchical & composition accuracy (external)
# --------------------------------------------------------------------------- #
def hierarchical_accuracy(pred, truth, hierarchy: dict[str, str], partial: float = 0.5) -> dict:
    """Accuracy with partial credit along a cell-type lineage tree. ``hierarchy`` maps
    subtype -> parent lineage (labels absent from it are their own parent). A correct subtype
    scores 1, a right-lineage/wrong-subtype scores ``partial``, a wrong lineage scores 0. Returns
    ``{hierarchical_accuracy, subtype_accuracy, lineage_accuracy}``."""
    pred = np.asarray(pred).astype(str)
    truth = np.asarray(truth).astype(str)
    pp = np.array([hierarchy.get(x, x) for x in pred])
    pt = np.array([hierarchy.get(x, x) for x in truth])
    exact = pred == truth
    lineage = pp == pt
    score = np.where(exact, 1.0, np.where(lineage, partial, 0.0))
    return {
        "hierarchical_accuracy": float(score.mean()),
        "subtype_accuracy": float(exact.mean()),
        "lineage_accuracy": float(lineage.mean()),
    }


def composition_accuracy(pred, truth) -> dict:
    """Predicted vs true cell-type **proportions**: L1 (total-variation-ish), Pearson r, and
    Jensen-Shannon divergence. Returns ``{l1, pearson, jsd, predicted, true}``."""
    import pandas as pd
    from scipy.spatial.distance import jensenshannon

    pred = pd.Series(np.asarray(pred).astype(str))
    truth = pd.Series(np.asarray(truth).astype(str))
    types = sorted(set(truth) | set(pred))
    p = pred.value_counts(normalize=True).reindex(types).fillna(0).to_numpy()
    q = truth.value_counts(normalize=True).reindex(types).fillna(0).to_numpy()
    l1 = float(np.abs(p - q).sum())
    pear = float(np.corrcoef(p, q)[0, 1]) if len(types) > 1 and p.std() > 0 and q.std() > 0 else 1.0
    jsd = float(jensenshannon(p, q) ** 2)                         # jensenshannon returns the distance (sqrt of div)
    return {"l1": l1, "pearson": pear, "jsd": jsd,
            "predicted": dict(zip(types, p.tolist())), "true": dict(zip(types, q.tolist()))}


def deconvolution_metrics(true_prop, pred_prop) -> dict:
    """Score predicted vs ground-truth cell-type PROPORTION matrices (spots × cell-types) - the
    OpenProblems ``task_spatial_decomposition`` metrics (+ RMSE / Pearson). This is the SOFT
    per-spot sibling of ``composition_accuracy`` (which compares hard-label GLOBAL proportions):
    use it to score a deconvolver's per-spot/per-cell proportion output (RCTD/TACCO/cell2location)
    against known mixtures (pseudospots). The two matrices must share cell-type columns + spot rows.

    Returns ``{r2, jsd, jsd_per_spot, rmse, pearson}``:
      * ``r2``          - ``sklearn.r2_score(true, pred, multioutput='uniform_average')`` (per cell
                          type across spots, averaged; higher better) - the OpenProblems ``r2``.
      * ``jsd``         - mean over cell types of ``scipy.jensenshannon(true[:,k], pred[:,k])``
                          (axis=0; lower better) - the OpenProblems ``jsd`` (JS distance).
      * ``jsd_per_spot``- mean over spots of the row-wise JS distance - the more standard per-spot
                          deconvolution JSD (each spot's proportion vector is a distribution).
      * ``rmse``        - RMSE over the whole matrix (lower better).
      * ``pearson``     - mean per-cell-type Pearson r across spots (higher better).
    """
    from scipy.spatial.distance import jensenshannon
    from sklearn.metrics import r2_score

    true = np.asarray(true_prop, dtype=float)
    pred = np.asarray(pred_prop, dtype=float)
    if true.shape != pred.shape:
        raise ValueError(f"proportion matrices differ in shape: {true.shape} vs {pred.shape}")
    n_ct = true.shape[1]

    r2 = float(r2_score(true, pred, multioutput="uniform_average"))
    jsd_ct = float(np.nanmean([jensenshannon(true[:, k], pred[:, k]) for k in range(n_ct)]))
    jsd_spot = float(np.nanmean([jensenshannon(true[i], pred[i]) for i in range(true.shape[0])]))
    rmse = float(np.sqrt(np.mean((true - pred) ** 2)))
    pears = [np.corrcoef(true[:, k], pred[:, k])[0, 1]
             for k in range(n_ct) if true[:, k].std() > 0 and pred[:, k].std() > 0]
    pearson = float(np.mean(pears)) if pears else float("nan")
    return {"r2": r2, "jsd": jsd_ct, "jsd_per_spot": jsd_spot, "rmse": rmse, "pearson": pearson}


# --------------------------------------------------------------------------- #
# Orchestrator (reference-free battery)
# --------------------------------------------------------------------------- #
def annotation_quality(adata, label_key: str = "cell_type",
                       marker_sets: dict[str, list[str]] | None = None,
                       sample_key: str | None = None, embedding: str | None = "X_pca",
                       k: int = 30, subsample: int = 5000, seed: int = 0) -> dict:
    """Run the reference-free quality battery and return a compact headline: internal validity,
    marker-program fidelity (if ``marker_sets``), and inter-sample consistency (if ``sample_key``
    present). Each sub-metric is guarded - a failure returns ``{'status': 'skipped: ...'}`` rather
    than raising, so this is safe to call from the QC funnel on any labeled section."""
    out: dict = {}
    try:
        out["internal_validity"] = internal_validity(
            adata, label_key=label_key, embedding=embedding, k=k, subsample=subsample, seed=seed)
    except Exception as e:  # pragma: no cover - defensive
        out["internal_validity"] = {"status": f"skipped: {e}"}
    if marker_sets:
        try:
            out["marker_fidelity"] = marker_program_fidelity(adata, label_key, marker_sets)
        except Exception as e:  # pragma: no cover
            out["marker_fidelity"] = {"status": f"skipped: {e}"}
    if sample_key and sample_key in adata.obs.columns:
        try:
            out["inter_sample_consistency"] = inter_sample_consistency(adata, label_key, sample_key)
        except Exception as e:  # pragma: no cover
            out["inter_sample_consistency"] = {"status": f"skipped: {e}"}
    return out
