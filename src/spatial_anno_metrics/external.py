"""External harness (needs ground truth): AvgBIO, accuracy / F1 / ARI / kappa / ECS, hierarchical and
composition accuracy, and the OpenProblems deconvolution proportion metrics."""
from __future__ import annotations

import numpy as np

from ._utils import _embedding

__all__ = ["avg_bio", "external_scores", "element_centric_similarity", "hierarchical_accuracy",
           "composition_accuracy", "deconvolution_metrics"]


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
