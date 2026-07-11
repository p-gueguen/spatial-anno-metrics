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

    # Davies-Bouldin (lower = better separated) + Calinski-Harabasz (higher = better): standard sklearn
    # internal indices, reported RAW (not folded into `integrated`) - the catalog battery-C over-cluster check.
    from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
    yenc = np.unique(ys, return_inverse=True)[1]
    db = float(davies_bouldin_score(Es, yenc)) if len(np.unique(ys)) > 1 else float("nan")
    ch = float(calinski_harabasz_score(Es, yenc)) if len(np.unique(ys)) > 1 else float("nan")

    parts = [_scale01(sil), _scale01(sil_2label), purity, orbital, ward_prop, _scale01(avg_sim)]
    integrated = float(np.exp(np.mean(np.log(np.clip(parts, 1e-9, 1.0)))))
    return {
        "silhouette": sil, "silhouette_2label": sil_2label, "neighborhood_purity": purity,
        "orbital_medoid": orbital, "ward_propmatch": ward_prop, "avg_similarity": avg_sim,
        "davies_bouldin": db, "calinski_harabasz": ch,
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


def marker_gene_overlap(adata, label_key: str, marker_sets: dict[str, list[str]], n_top: int = 25) -> dict:
    """Do the labels reproduce known markers? Per type, the fraction of its curated markers that land in
    the top-``n_top`` data-derived DE genes for that label (standardized in-vs-out mean difference).
    Reference-free given marker sets; the DE-overlap complement to :func:`marker_program_fidelity`'s
    ranking view. Returns ``{per_type: {overlap, n_curated}, mean_overlap}``."""
    genes = np.array(list(map(str, adata.var_names)))
    gene_set = set(genes)
    labels = adata.obs[label_key].astype(str).to_numpy()
    X = _dense(adata.X).astype(float)
    per, fracs = {}, []
    for t, curated in marker_sets.items():
        cur = [g for g in map(str, curated) if g in gene_set]
        y = labels == t
        if not cur or y.sum() < 2 or (~y).sum() < 2:
            per[t] = {"overlap": None, "n_curated": len(cur)}
            continue
        score = (X[y].mean(0) - X[~y].mean(0)) / np.sqrt(X[~y].var(0) + 1e-9)
        top = set(genes[np.argsort(-score)[:n_top]])
        ov = len(top & set(cur)) / len(cur)
        per[t] = {"overlap": float(ov), "n_curated": len(cur)}
        fracs.append(ov)
    return {"per_type": per, "mean_overlap": float(np.mean(fracs)) if fracs else None}


def avg_bio(adata, label_key: str, truth_key: str, embedding: str | None = "X_pca") -> dict:
    """AvgBIO (Zhu 2026 / scIB): mean(ARI, NMI, scaled-ASW) - biological-structure conservation of a
    labeling vs a ground truth. ARI/NMI compare ``label_key`` to ``truth_key``; ASW is the silhouette
    of the TRUE labels in ``embedding`` (how well real types separate), scaled to [0,1]. External (needs
    ground truth). Returns ``{ari, nmi, asw, avg_bio}``."""
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score

    pred = adata.obs[label_key].astype(str).to_numpy()
    truth = adata.obs[truth_key].astype(str).to_numpy()
    ari = float(adjusted_rand_score(truth, pred))
    nmi = float(normalized_mutual_info_score(truth, pred))
    E = _embedding(adata, embedding)
    asw = float(silhouette_score(E, truth)) if len(np.unique(truth)) > 1 else 0.0
    return {"ari": ari, "nmi": nmi, "asw": asw,
            "avg_bio": float(np.mean([ari, nmi, (asw + 1) / 2]))}


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
def annotation_quality_index(adata, label_key: str, marker_sets: dict[str, list[str]] | None,
                             internal_validity_result: dict | None = None, *,
                             reference=None, ref_key: str = "cell_type", panel_genes=None,
                             median_depth: float | None = None, method_label_cols=None,
                             abstention_labels=None, provenance: dict | None = None,
                             min_type_n: int = 50, p_softmin: float = -4.0) -> dict:
    """AQI - one Annotation Quality Index in [0, 1].

    ``AQI = w_coh * min(A, soft_min_{p=-4}({C, M}))`` - the soft-min of the ACHIEVED-quality terms,
    capped by the panel-adequacy ceiling: **C** contamination/purity (macro median PMP x retention;
    panel-size invariant), **M** marker-program fidelity (macro one-vs-rest AUC), capped by **A** panel x
    depth adequacy (base-rate-normalized depth-matched macro F1 - a resolvable CEILING: you cannot type
    better than the panel resolves, but a rich panel does not by itself mean the labels are right).
    Coherence **H** is a bounded <=15% multiplier (self-consistency can never manufacture score);
    completeness is a side flag (honest abstention is not low quality). Every component is macro over
    predicted types with ``n >= min_type_n`` and any UNDEFINED term is DROPPED (never set to 0).

    Cross-method AGREEMENT **G** (>=3 voters) is NOT a cross-section term: its absolute level tracks
    REFERENCE quality, not accuracy, and does not transfer across sections (2026-07 validation). It is
    reported under ``abstention`` as a WITHIN-section per-cell signal (which cells to grey out), and
    ``regime`` (``with_agreement`` vs ``index_only``) says whether that signal is available.

    It is an INDEX (section-comparable, monotone-in-quality), **not** ``P(correct)`` - the no-transfer
    finding forbids a universal reference-free accuracy curve. Ground-truth accuracy is a caller concern
    (:func:`external_scores`), not this index.

    Returns ``{aqi, aqi_core, regime, components:{A,C,G,M,H,completeness}, abstention:{signal,n_voters,
    available}, w_coh, active_set, argmin, n_ge50_types, flags, provenance}``. Pure/guarded; QC-safe.
    """
    from .purity import pmp as _pmp

    labels = adata.obs[label_key].astype(str).to_numpy()
    abst = {str(x) for x in (abstention_labels or [])}
    vc = adata.obs[label_key].astype(str).value_counts()
    big = [str(t) for t, n in vc.items() if int(n) >= min_type_n and str(t) not in abst]
    big_mask = np.isin(labels, big)

    def _c01(x):
        return float(min(max(float(x), 0.0), 1.0))

    # A - panel x depth adequacy: base-rate-null-normalized macro F1 at the section's own depth.
    A = None
    frac_resolvable = None
    if reference is not None and panel_genes is not None:
        try:
            pr = panel_resolvability(reference, ref_key, panel_genes,
                                     target_depth=float(median_depth) if median_depth else 50.0)
            mf, nt = pr.get("macro_f1"), pr.get("n_types")
            frac_resolvable = pr.get("frac_resolvable")
            if mf is not None and nt and int(nt) >= 2:
                null = 1.0 / int(nt)
                A = _c01((float(mf) - null) / (1.0 - null))
        except Exception:  # pragma: no cover - defensive
            A = None

    # C - contamination/purity: min(macro median PMP, retention). Deleting transcripts -> NaN PMP ->
    # retention falls, so it cannot be gamed by dropping counts. Panel-size invariant (lineage-marker denom).
    C = None
    if marker_sets:
        try:
            _pmp(adata, marker_sets, label_key)                     # writes obs['pmp']
            pv = np.asarray(adata.obs["pmp"], dtype=float)
            meds, defined = [], 0
            for t in big:
                col = pv[labels == t]
                ok = col[~np.isnan(col)]
                if ok.size:
                    meds.append(float(np.median(ok)))
                defined += int((~np.isnan(col)).sum())
            if meds and big_mask.any():
                C = _c01(min(float(np.mean(meds)), defined / int(big_mask.sum())))
        except Exception:  # pragma: no cover - defensive
            C = None

    # G - cross-method agreement (the only reference-free correctness term); needs >=3 diverse voters.
    G = None
    cols = [c for c in (method_label_cols or []) if c in adata.obs.columns]
    if len(cols) >= 3:
        try:
            votes = np.vstack([adata.obs[c].astype(str).to_numpy() for c in cols])   # voters x cells
            agree = (votes == labels[None, :]).mean(0)              # frac of voters backing the label
            per = [float(agree[labels == t].mean()) for t in big if (labels == t).any()]
            if per:
                G = _c01(float(np.mean(per)))
        except Exception:  # pragma: no cover - defensive
            G = None

    # M - marker-program fidelity: macro one-vs-rest AUC over n>=50 types, rescaled 2*(auc-0.5).
    M = None
    if marker_sets:
        try:
            mpf = marker_program_fidelity(adata, label_key, marker_sets)
            aucs = [v["auc"] for t, v in mpf.get("per_type", {}).items()
                    if t in big and v.get("auc") is not None and int(v.get("n", 0)) >= min_type_n]
            if aucs:
                M = _c01(2.0 * (float(np.mean(aucs)) - 0.5))
        except Exception:  # pragma: no cover - defensive
            M = None

    # H - coherence: bounded <=15% multiplier only (coherence != correctness).
    H = internal_validity_result.get("integrated") if isinstance(internal_validity_result, dict) else None
    w_coh = float(np.clip(0.85 + 0.15 * H, 0.85, 1.0)) if isinstance(H, (int, float)) else 1.0

    # completeness - side usability flag, NOT a term (honest abstention on a shallow panel is correct).
    typed_frac = float(np.mean(~np.isin(labels, list(abst)))) if abst else 1.0
    completeness = _c01(typed_frac / frac_resolvable) if frac_resolvable else _c01(typed_frac)

    # combine: soft-min over the ACHIEVED-quality terms {C, M}, CAPPED by the panel x depth ceiling A,
    # then the bounded coherence multiplier. This is the 2026-07 VALIDATION correction (validate_aqi.py,
    # 3 GT sections, common 6-lineage axis). Measured section-level Spearman(component, true balanced
    # accuracy): M = +1.0 (marker fidelity orders sections perfectly), C = +0.5, but A = -1.0 and
    # G = +0.5. A ANTI-orders accuracy because it is a resolvable-CEILING (high on any deep panel
    # regardless of the labels actually assigned) - so it must CAP, never drive. G's absolute level
    # tracks REFERENCE quality, not accuracy, and does not transfer across sections - so it is a WITHIN-
    # section per-cell abstention signal (reported under `abstention`), never a cross-section term. The
    # old soft_min(A,C,G,M) mixed a ceiling and a non-transferring term into the bottleneck (rho +0.5,
    # wrong section order); min(A, soft_min(C,M)) restores rho +1.0.
    core = {"C": C, "M": M}
    active = {k: v for k, v in core.items() if v is not None}
    names = {"C": "contamination", "M": "marker fidelity", "A": "panel/depth ceiling"}
    aqi_core = argmin = aqi = None
    if active:
        clipped = {k: min(max(v, 0.01), 1.0) for k, v in active.items()}
        aqi_core = _c01(float(np.mean([c ** p_softmin for c in clipped.values()])) ** (1.0 / p_softmin))
        capped = aqi_core if A is None else min(aqi_core, A)          # adequacy is a CEILING, not a driver
        argmin = ("panel/depth ceiling" if (A is not None and A < aqi_core)
                  else names[min(clipped, key=clipped.get)])
        aqi = _c01(w_coh * capped)

    return {
        "aqi": aqi, "aqi_core": aqi_core,
        "regime": "with_agreement" if G is not None else "index_only",
        "components": {"A": A, "C": C, "G": G, "M": M, "H": H, "completeness": completeness},
        # cross-method agreement is a WITHIN-section per-cell abstention signal, not a cross-section term.
        "abstention": {"signal": G, "n_voters": len(cols), "available": G is not None},
        "w_coh": w_coh, "active_set": list(active.keys()) + (["A"] if A is not None else []),
        "argmin": argmin, "n_ge50_types": len(big),
        "flags": {"adequacy_unknown": A is None, "low_support": len(big) < 2,
                  "abstention_available": G is not None},
        "provenance": dict(provenance or {}, n_voters=len(cols)),
    }


def annotation_quality(adata, label_key: str = "cell_type",
                       marker_sets: dict[str, list[str]] | None = None,
                       sample_key: str | None = None, embedding: str | None = "X_pca",
                       k: int = 30, subsample: int = 5000, seed: int = 0, *,
                       reference=None, ref_key: str = "cell_type", panel_genes=None,
                       median_depth: float | None = None, method_label_cols=None,
                       abstention_labels=None, platform=None, n_panel_genes=None,
                       normalization=None, marker_set_name=None, compute_aqi: bool = True) -> dict:
    """Run the reference-free quality battery and return a compact headline: internal validity,
    marker-program fidelity (if ``marker_sets``), inter-sample consistency (if ``sample_key``), and
    the single **AQI** headline (:func:`annotation_quality_index`) in ``out['aqi']``. Each sub-metric is
    guarded - a failure returns ``{'status': 'skipped: ...'}`` rather than raising, so this is safe to
    call from the QC funnel on any labeled section.

    The AQI extra inputs are all optional and degrade gracefully: pass ``reference``+``panel_genes``
    (+``median_depth``) for the adequacy term, ``method_label_cols`` (>=3 obs columns) for the
    agreement term, and ``marker_sets`` for purity + fidelity. With none of them AQI is a pure
    coherence index. ``platform``/``n_panel_genes``/``normalization``/``marker_set_name`` are stamped
    into ``provenance`` (every number must name its normalization + marker set)."""
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
    if compute_aqi:
        try:
            out["aqi"] = annotation_quality_index(
                adata, label_key, marker_sets, out.get("internal_validity"),
                reference=reference, ref_key=ref_key, panel_genes=panel_genes,
                median_depth=median_depth, method_label_cols=method_label_cols,
                abstention_labels=abstention_labels,
                provenance={"platform": platform, "n_panel_genes": n_panel_genes,
                            "median_depth": median_depth, "normalization": normalization,
                            "marker_set": marker_set_name, "reference": reference is not None})
        except Exception as e:  # pragma: no cover - defensive
            out["aqi"] = {"status": f"skipped: {e}"}
    return out


# --------------------------------------------------------------------------- #
# Conformal prediction sets (catalog s3d) - distribution-free typing uncertainty
# --------------------------------------------------------------------------- #
def conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha: float = 0.1,
                              class_conditional: bool = False,
                              hierarchy: dict | None = None, y_query=None) -> dict:
    """Distribution-free conformal prediction SETS for cell-type annotation (catalog s3d
    "conformal set size" - the cleanest reference-free formalization of abstention).

    Model-agnostic: bring any annotator's per-cell class-probability matrix (RCTD weights,
    SingleR/CellTypist/scANVI posteriors, marker-score softmax). Split-conformal LAC: the
    calibration nonconformity is ``s = 1 - p(true class)``; the (1-alpha) quantile ``qhat`` sets
    the inclusion threshold and the query set is ``C(x) = {y : p(x)_y >= 1 - qhat}``. Set size
    **1 = confident, >1 = ambiguous, 0 = novel/OOD**, with the finite-sample guarantee
    ``P(true in C(x)) >= 1 - alpha`` **provided calibration and query are exchangeable**. On
    native spatial data calibrated against a dissociated reference, exchangeability is broken by
    platform shift, so the nominal coverage is approximate (recalibrate on a small platform-matched
    labeled subset, e.g. protein-confirmed cells) - the set SIZE stays a valid relative uncertainty.

    Parameters
    ----------
    proba_query : (n_q, K) query class probabilities. classes : (K,) labels matching the columns.
    proba_cal, y_cal : (n_c, K) calibration probabilities + their TRUE labels.
    class_conditional : per-true-class (Mondrian) quantile - **recommended**; marginal conformal
        silently under-covers rare cell types while class-conditional restores per-type coverage.
    hierarchy : optional ``{leaf: parent}`` map - also returns lineage-collapsed sets (scConform's
        ontology idea: an ambiguous SUBTYPE set becomes a confident LINEAGE call).
    y_query : optional TRUE query labels -> also returns empirical coverage + per-type coverage.

    Returns ``{sets, set_size, summary, [collapsed_sets, pct_ambiguous_one_parent],
    [coverage, per_type_coverage]}``. Depends on numpy only.
    """
    proba_query = np.asarray(proba_query, dtype=float)
    proba_cal = np.asarray(proba_cal, dtype=float)
    classes = np.asarray([str(c) for c in classes])
    y_cal = np.asarray(y_cal).astype(str)
    cls_idx = {c: j for j, c in enumerate(classes)}

    true_col = np.array([cls_idx[c] for c in y_cal])
    s_cal = 1.0 - proba_cal[np.arange(len(y_cal)), true_col]

    def _qhat(scores):
        n = len(scores)
        if n == 0:
            return 1.0
        level = min(1.0, np.ceil((n + 1) * (1 - alpha)) / n)
        return float(np.quantile(scores, level, method="higher"))

    if class_conditional:
        thr = {c: (_qhat(s_cal[y_cal == c]) if (y_cal == c).sum() >= 10 else _qhat(s_cal))
               for c in classes}
        incl_thresh = np.array([1.0 - thr[c] for c in classes])      # per-class inclusion threshold
        include = proba_query >= incl_thresh[None, :]
    else:
        q = _qhat(s_cal)
        include = proba_query >= (1.0 - q)

    sets = [[classes[j] for j in np.where(include[i])[0]] for i in range(include.shape[0])]
    size = include.sum(axis=1)
    out = {
        "sets": sets, "set_size": size,
        "summary": {
            "mean_set_size": float(size.mean()) if size.size else float("nan"),
            "pct_singleton": float((size == 1).mean()) if size.size else float("nan"),
            "pct_ambiguous": float((size > 1).mean()) if size.size else float("nan"),
            "pct_empty_novel": float((size == 0).mean()) if size.size else float("nan"),
            "alpha": alpha, "class_conditional": class_conditional,
        },
    }
    if hierarchy is not None:
        collapsed = [sorted({hierarchy.get(x, x) for x in s}) for s in sets]
        out["collapsed_sets"] = collapsed
        amb = size > 1
        out["pct_ambiguous_one_parent"] = float(
            sum(1 for i in range(len(sets)) if amb[i] and len(collapsed[i]) == 1) / max(1, amb.sum()))
    if y_query is not None:
        y_query = np.asarray(y_query).astype(str)
        covered = np.array([y_query[i] in sets[i] for i in range(len(sets))])
        out["coverage"] = float(covered.mean())
        out["per_type_coverage"] = {
            c: float(np.mean([y_query[i] in sets[i] for i in np.where(y_query == c)[0]]))
            for c in np.unique(y_query)}
    return out
