"""Ontology-aware annotation metrics: whether labels sit in the embedding the way the Cell Ontology relates them."""
from __future__ import annotations

import numpy as np

from ._utils import _embedding
from ._validate import check_labels, skip
from .thresholds import value as _thr

__all__ = ["ontology_consistency"]
_MIN_LABELS = _thr("ontology", "min_labels")


def _read_obo(path) -> dict[str, list[str]]:
    """``{term id: [is_a parent ids]}`` of the non-obsolete ``[Term]`` stanzas of an OBO 1.2/1.4 file.

    Parameters
    ----------
    path : str or path-like
        OBO file, e.g. ``cl-basic.obo`` from the Cell Ontology releases.

    Returns
    -------
    dict
        Term id -> its ``is_a`` parents, in file order. ``relationship:`` lines (part_of, develops_from, ...),
        ``[Typedef]`` stanzas and ``is_obsolete: true`` terms are ignored.
    """
    stanzas, cur = [], None
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith("["):
                cur = {"kind": line, "is_a": []}
                stanzas.append(cur)
            elif cur is not None and ": " in line:
                key, val = line.split(": ", 1)
                if key == "is_a":
                    cur["is_a"].append(val.split()[0])            # drops "{qualifiers}" and "! name"
                elif key in ("id", "is_obsolete"):
                    cur[key] = val.split()[0]
    return {s["id"]: s["is_a"] for s in stanzas
            if s["kind"] == "[Term]" and "id" in s and s.get("is_obsolete") != "true"}


def _ontology_distances(ids, parents: dict) -> np.ndarray:
    """Pairwise is_a path length between terms: the shortest walk up from both to a common ancestor.

    Parameters
    ----------
    ids : list of str
        Term ids, all keys of ``parents``.
    parents : dict
        ``{term: [is_a parents]}`` from :func:`_read_obo`.

    Returns
    -------
    numpy.ndarray
        Symmetric ``len(ids) x len(ids)`` float matrix, 0 on the diagonal. Two terms with no common ancestor get
        one more than the largest finite distance (only ranks are used downstream).
    """
    up = []
    for t in ids:                                                  # BFS up the DAG: shortest hop count
        dist, frontier = {t: 0}, [t]
        while frontier:
            nxt = []
            for c in frontier:
                for p in parents.get(c, []):
                    if p not in dist:
                        dist[p] = dist[c] + 1
                        nxt.append(p)
            frontier = nxt
        up.append(dist)
    k = len(ids)
    D = np.zeros((k, k))
    for i in range(k):
        for j in range(i + 1, k):
            common = up[i].keys() & up[j].keys()
            D[i, j] = D[j, i] = min((up[i][c] + up[j][c] for c in common), default=np.inf)
    if np.isinf(D).any():
        D[np.isinf(D)] = D[np.isfinite(D)].max() + 1
    return D


def _rank_corr(x, y) -> float:
    """Spearman correlation (Pearson of average ranks); NaN when either side is constant."""
    from scipy.stats import rankdata

    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def ontology_consistency(adata, label_key: str, obo_path, label_to_cl: dict, *, embedding: str | None = "X_pca",
                         n_perm: int = _thr("ontology", "n_perm"), min_cells: int = _thr("ontology", "min_cells"),
                         seed: int = 0) -> dict:
    """Score whether label centroids sit in the embedding the way the Cell Ontology relates the labels.

    Reference-free and truth-free: a "macrophage" whose centroid lies inside the epithelium is far from its
    ontology neighbours (monocyte) and close to its ontology strangers. For every pair of labels, the Euclidean
    distance between their embedding centroids is compared with the ``is_a`` path length between their Cell
    Ontology terms (shortest walk up to a common ancestor); ``rho`` is the Spearman correlation over the pairs,
    and ``p_value`` a one-sided Mantel test that permutes which label holds which ontology term
    (``(1 + #{rho_perm >= rho}) / (1 + n_perm)``; with ``K`` labels there are only ``K!`` distinct permutations,
    so with 4 labels the p cannot go below ~1/24). ``loo_delta`` is how much ``rho`` rises when that label is
    left out; the label with the largest positive one is ``most_inconsistent``.

    The ontology is read with the standard library (``is_a`` only; ``part_of`` and other relationships are
    ignored). Siblings far apart in a good embedding (two T-cell states) lower ``rho`` as well: the score says
    how ontology-like the layout is, not that a given label is wrong.

    Parameters
    ----------
    adata : AnnData
        Annotated cells.
    label_key : str
        ``adata.obs`` column with the labels.
    obo_path : str or path-like
        OBO file of the ontology, e.g. ``cl-basic.obo`` from github.com/obophenotype/cell-ontology releases.
    label_to_cl : dict
        Label -> ontology term id (``"CL:0000235"``). Labels missing here or whose id is not in the file are
        left out and listed in ``excluded``.
    embedding : str or None
        ``adata.obsm`` key; absent or None -> a PCA of ``X``.
    n_perm : int
        Label permutations for the p-value.
    min_cells : int
        Labels with fewer cells are left out (their centroid is noise).
    seed : int
        Seed for the permutations and the fallback PCA.

    Returns
    -------
    dict
        ``{status, rho, p_value, n_perm, n_labels, most_inconsistent, per_label: {label: {cl_id, loo_delta,
        nearest_embedding, nearest_ontology}}, excluded: {label: reason}}``. ``nearest_ontology`` lists every
        label tied at the smallest ontology distance. Fewer than 4 usable labels, or ontology distances that are
        all equal, -> a skip (with ``excluded``). ``rho`` is NaN when all centroid distances are equal.
    """
    y = check_labels(adata, label_key)
    parents = _read_obo(obo_path)
    labels, counts = np.unique(y, return_counts=True)
    excluded, keep = {}, []
    for lab, n in zip(labels.tolist(), counts.tolist()):
        cl = label_to_cl.get(lab)
        if cl is None:
            excluded[lab] = "no CL id in label_to_cl"
        elif cl not in parents:
            excluded[lab] = f"{cl} not in the ontology file"
        elif n < min_cells:
            excluded[lab] = f"{n} cells < min_cells={min_cells}"
        else:
            keep.append(str(lab))
    if len(keep) < _MIN_LABELS:
        return {**skip(f"{len(keep)} labels with a CL id and >= {min_cells} cells, need >= {_MIN_LABELS}"),
                "excluded": excluded}
    ids = [label_to_cl[lab] for lab in keep]
    Do = _ontology_distances(ids, parents)
    E = _embedding(adata, embedding, seed=seed)
    C = np.vstack([E[y == lab].mean(0) for lab in keep])
    De = np.sqrt(((C[:, None, :] - C[None, :, :]) ** 2).sum(-1))
    k = len(keep)
    iu = np.triu_indices(k, 1)
    if np.unique(Do[iu]).size < 2:
        return {**skip("all ontology distances between the kept labels are equal"), "excluded": excluded}
    rho = _rank_corr(De[iu], Do[iu])
    rng = np.random.default_rng(seed)
    perm = [_rank_corr(De[iu], Do[np.ix_(p, p)][iu]) for p in (rng.permutation(k) for _ in range(n_perm))]
    p_value = (1 + int(np.sum(np.asarray(perm) >= rho))) / (1 + n_perm) if np.isfinite(rho) else float("nan")
    per = {}
    for i, lab in enumerate(keep):
        rest = [j for j in range(k) if j != i]
        sub = np.ix_(rest, rest)
        iu_r = np.triu_indices(k - 1, 1)
        others = np.array(rest)
        o_row = Do[i, others]
        per[lab] = {"cl_id": ids[i], "loo_delta": _rank_corr(De[sub][iu_r], Do[sub][iu_r]) - rho,
                    "nearest_embedding": keep[int(others[np.argmin(De[i, others])])],
                    "nearest_ontology": [keep[j] for j in others[o_row == o_row.min()]]}
    deltas = {lab: d["loo_delta"] for lab, d in per.items() if np.isfinite(d["loo_delta"])}
    worst = max(deltas, key=deltas.get) if deltas else None
    return {"status": "ok", "rho": rho, "p_value": p_value, "n_perm": int(n_perm), "n_labels": k,
            "most_inconsistent": worst if worst is not None and deltas[worst] > 0 else None,
            "per_label": per, "excluded": excluded}
