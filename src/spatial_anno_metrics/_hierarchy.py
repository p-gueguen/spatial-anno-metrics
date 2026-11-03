"""Label-hierarchy helpers on a ``{child: parent}`` map, for hierarchy-aware metrics.

A label absent from the map, or mapped to itself, is a root (the convention of ``hierarchical_accuracy``).
"""
from __future__ import annotations

__all__ = []


def ancestors(label, parent_map: dict) -> list:
    """Return ``label``'s ancestors, nearest first: ``[parent, grandparent, ..., root]``; ``[]`` for a root.

    A cycle in ``parent_map`` raises ``ValueError``.
    """
    out, seen, cur = [], {label}, label
    while cur in parent_map and parent_map[cur] != cur:
        cur = parent_map[cur]
        if cur in seen:
            raise ValueError(f"cycle in parent_map through {cur!r}")
        seen.add(cur)
        out.append(cur)
    return out


def lca_depth(a, b, parent_map: dict) -> int:
    """Return the depth of the lowest common ancestor of ``a`` and ``b``, counting a root as depth 0.

    ``a == b`` gives the depth of ``a``; an ancestor of the other label is the LCA. ``-1`` when they share no
    ancestor, i.e. sit under different roots.
    """
    on_b = {b, *ancestors(b, parent_map)}
    lca = next((x for x in [a, *ancestors(a, parent_map)] if x in on_b), None)
    return -1 if lca is None else len(ancestors(lca, parent_map))
