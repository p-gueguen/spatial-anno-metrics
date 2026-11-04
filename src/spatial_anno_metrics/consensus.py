"""Multi-annotator consensus with ontology-path voting, and its entropy.

An agreement signal between several annotations of the same cells that does not penalise one annotator calling a
subtype and another its parent.
"""
from __future__ import annotations

from collections import Counter

import numpy as np

from ._hierarchy import ancestors, validate_hierarchy
from ._validate import check_labels, skip
from .thresholds import value as _thr

__all__ = ["consensus_entropy"]


def _entropy(labels) -> float:
    p = np.array(list(Counter(labels).values()), float) / len(labels)
    return float(-(p * np.log(p)).sum())


def consensus_entropy(adata, label_keys, *, hierarchy: dict | None = None,
                      min_share: float = _thr("consensus", "min_share")) -> dict:
    """Per-cell consensus label and agreement entropy across annotators, with ontology-path voting.

    Every annotator's label votes for itself and for all its ancestors in ``hierarchy`` (popV's ontology
    voting, Ergen et al. 2024, Nat Genet, doi 10.1038/s41588-024-01993-3; popV propagates only OnClass's
    votes, here every annotator's do). The consensus is the deepest label whose vote share exceeds
    ``min_share`` (strict majority by default; in a tree it is unique), or None when no label has one - a
    plurality rule would always pick the root once votes reach ancestors. Two entropies, both divided by
    ``log(n_annotators)`` so 0 = all agree and 1 = all differ:

    - ``entropy``: flat Shannon entropy of the raw votes (a CD4 / CD8 / T split counts as full disagreement);
    - ``ontology_entropy``: entropy of the votes projected to the depth of the consensus label (to the roots
      when there is none), so disagreement *below* the consensus does not count; votes shallower than that
      depth stay as they are. Without a hierarchy the two are equal.

    Parameters
    ----------
    adata : AnnData
        Cells, one annotation per ``obs`` column.
    label_keys : list of str
        The ``obs`` columns to compare, one per annotator (>= 2). Missing labels raise ``ValueError``.
    hierarchy : dict, optional
        ``{child: parent}`` map (a label absent from it is a root). A cycle raises ``ValueError``.
    min_share : float
        Vote share a label must exceed to be the consensus, in [0.5, 1).

    Returns
    -------
    dict
        ``{status, consensus, consensus_share, entropy, ontology_entropy, summary: {n_annotators, mean_entropy,
        mean_ontology_entropy, mean_consensus_share, pct_full_agreement, pct_no_consensus}}``; arrays are per
        cell. ``consensus_share`` is the consensus label's vote share, or the top share of any label when
        there is no consensus. Fewer than 2 annotators is a skip. Writes nothing to ``adata``.
    """
    if not 0.5 <= min_share < 1:
        raise ValueError(f"min_share must be in [0.5, 1) so the consensus is unique, got {min_share}")
    label_keys = list(label_keys)
    if len(label_keys) < 2:
        return skip(f"need >= 2 annotations to compare, got {len(label_keys)}")
    L = np.column_stack([check_labels(adata, k) for k in label_keys]).astype(str)
    n_a = L.shape[1]
    hierarchy = hierarchy or {}
    validate_hierarchy(hierarchy)
    paths: dict = {}

    def path(label):                                      # root-first: [root, ..., label]
        if label not in paths:
            paths[label] = [*reversed(ancestors(label, hierarchy)), label]
        return paths[label]

    uniq, inv = np.unique(L, axis=0, return_inverse=True)
    res = []
    for row in uniq:
        ps = [path(x) for x in row]
        votes = Counter(x for p in ps for x in p)
        best = max((x for x, v in votes.items() if v / n_a > min_share), key=lambda x: len(path(x)), default=None)
        share = votes[best] / n_a if best is not None else max(votes.values()) / n_a
        d = len(path(best)) - 1 if best is not None else 0
        res.append((best, share, _entropy(list(row)) / np.log(n_a),
                    _entropy([p[min(d, len(p) - 1)] for p in ps]) / np.log(n_a)))
    inv = np.asarray(inv).ravel()
    cons = np.array([r[0] for r in res], dtype=object)[inv]
    share, ent, onto = (np.array([r[i] for r in res], float)[inv] for i in (1, 2, 3))
    return {"status": "ok", "consensus": cons, "consensus_share": share, "entropy": ent, "ontology_entropy": onto,
            "summary": {"n_annotators": n_a, "mean_entropy": float(ent.mean()),
                        "mean_ontology_entropy": float(onto.mean()), "mean_consensus_share": float(share.mean()),
                        "pct_full_agreement": float((ent == 0).mean()),
                        "pct_no_consensus": float(np.mean([c is None for c in cons]))}}
