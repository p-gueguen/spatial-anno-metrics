"""Parity golden for element-centric similarity: clusim ``element_sim`` vs our ``element_centric_similarity``.

Writes ``tests/goldens/ecs_clusim/`` (``pairs.json`` + ``expected.json``) and prints the ours-vs-reference
table for its SOURCE.md. Not run in CI. Needs a venv with clusim 0.4 (PyPI, MIT) and pandas; run from
anywhere:

    <venv>/bin/python benchmarks/parity/ecs_clusim/make_golden.py
"""
from __future__ import annotations

import json
import sys
from datetime import date
from importlib.metadata import version
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "tests" / "goldens" / "ecs_clusim"
CLUSIM = "0.4"
ALPHAS = (0.9, 0.5)              # 0.9 = both defaults; for a partition the score should not depend on alpha
SIZES = [100, 70, 50, 40, 25, 15]                                     # 300 elements, uneven clusters


def make_pairs() -> dict:
    rng = np.random.default_rng(0)
    a = rng.permutation(np.repeat(np.arange(len(SIZES)), SIZES))
    split = a.copy()
    big = np.flatnonzero(a == 0)
    split[big[:40]] = len(SIZES)                                        # cluster 0 (100) -> 60 + 40
    merged = np.where(a == 5, 4, a)                                     # clusters 4 (25) + 5 (15) -> one
    relabel = np.array([3, 5, 0, 1, 4, 2])[a]                           # same partition, other names
    pairs = {"identical": a.copy(), "relabelled": relabel, "split": split, "merged": merged,
             "random": rng.integers(0, len(SIZES), a.size)}
    return {k: {"a": a.tolist(), "b": b.tolist()} for k, b in pairs.items()}


def main() -> None:
    sys.path.insert(0, str(REPO / "src"))
    from clusim import sim as csim
    from clusim.clustering import Clustering

    import spatial_anno_metrics as sam
    assert Path(sam.__file__).resolve().is_relative_to(REPO / "src"), sam.__file__
    assert version("clusim") == CLUSIM, version("clusim")

    OUT.mkdir(parents=True, exist_ok=True)
    pairs = make_pairs()
    (OUT / "pairs.json").write_text(json.dumps(pairs) + "\n")

    def clu(labels):
        c = Clustering()
        c.from_membership_list(labels)
        return c

    results = {}
    for alpha in ALPHAS:
        res = {}
        for name, p in pairs.items():
            c1, c2 = clu(p["a"]), clu(p["b"])
            scores, idx = csim.element_sim_elscore(c1, c2, alpha=alpha)
            r = {"ecs": float(csim.element_sim(c1, c2, alpha=alpha))}
            if alpha == ALPHAS[0]:
                r["per_element"] = [float(scores[idx[e]]) for e in range(len(p["a"]))]
            res[name] = r
        results[f"alpha={alpha}"] = res
    (OUT / "expected.json").write_text(json.dumps({
        "generated": date.today().isoformat(), "upstream": f"clusim=={CLUSIM} (PyPI)",
        "call": "clusim.sim.element_sim(Clustering().from_membership_list(a), ...(b), alpha=alpha); "
                "per_element from element_sim_elscore, in element order",
        "results": results}, indent=1, allow_nan=False) + "\n")

    print("| pair | clusim alpha=0.9 | clusim alpha=0.5 | ours | ours - clusim |")
    print("|---|---|---|---|---|")
    for name, p in pairs.items():
        ref = results["alpha=0.9"][name]["ecs"]
        ours = sam.element_centric_similarity(p["a"], p["b"])
        half = results["alpha=0.5"][name]["ecs"]
        print(f"| {name} | {ref:.17g} | {half:.17g} | {ours:.17g} | {ours - ref:.1e} |")


if __name__ == "__main__":
    main()
