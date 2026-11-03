# ecs_clusim - element-centric similarity (clusim), reference outputs

Parity golden for `element_centric_similarity` against clusim, the reference implementation by the method's
authors. Generated 2026-11-03. Ours matches to <= 1.1e-16 on all five pairs.

## Files

| File | Content |
|---|---|
| `pairs.json` | `{pair: {"a": [300 ints], "b": [300 ints]}}`, element `i` is position `i` in both lists |
| `expected.json` | `results["alpha=0.9"][pair]` = `{"ecs", "per_element"}` (300 scores, element order); `results["alpha=0.5"][pair]` = `{"ecs"}` |

## Data

Synthetic label vectors made by `benchmarks/parity/ecs_clusim/make_golden.py` (this repository, MIT); no
external data, no DOI. `a` = 6 clusters of sizes 100/70/50/40/25/15, shuffled with
`np.random.default_rng(0).permutation`; `b` per pair:

| pair | `b` | expected by hand |
|---|---|---|
| `identical` | `a` | 1 |
| `relabelled` | `a` with names mapped 0->3, 1->5, 2->0, 3->1, 4->4, 5->2 (same partition) | 1 |
| `split` | the first 40 (by position) of cluster 0's 100 elements moved to a new cluster 6 | (60 x 0.6 + 40 x 0.4 + 200) / 300 = 0.84 |
| `merged` | cluster 5 (15) merged into cluster 4 (25) | (25 x 0.625 + 15 x 0.375 + 260) / 300 = 0.9375 |
| `random` | 300 draws from `integers(0, 6)` of the same generator, after the shuffle | - |

## Upstream

- clusim 0.4 from PyPI (`clusim-0.4-py3-none-any.whl`, uploaded 2020-05-13, sha256
  `f2fc70ec9238610408072604a235b1598f02f3cca28b4bc4a3a6ec25024affa2`), MIT,
  https://github.com/Hoosier-Clusters/clusim.
- Paper: Gates, Wood, Hetrick & Ahn, "Element-centric clustering comparison unifies overlaps and hierarchy",
  *Sci Rep* 2019, [10.1038/s41598-019-44892-y](https://doi.org/10.1038/s41598-019-44892-y).

## Exact call

```python
from clusim.clustering import Clustering
import clusim.sim as sim
c1 = Clustering(); c1.from_membership_list(a)
c2 = Clustering(); c2.from_membership_list(b)
sim.element_sim(c1, c2, alpha=alpha)                       # ecs, alpha in (0.9, 0.5); 0.9 is the default
scores, index = sim.element_sim_elscore(c1, c2, alpha=0.9)
per_element = [scores[index[i]] for i in range(300)]       # back to element order
```

Defaults otherwise (`r=1.`, `rescale_path_type='max'`, `ppr_implementation='prpack'`; the partition path does
not use them).

## Definition, as implemented

For a partition, clusim builds each element's personalized-PageRank affinity in closed form,
`p_i(j) = alpha / |C_i|` for `j` in `i`'s cluster plus `(1 - alpha)` when `j = i`; the element score is
`1 - ||p_i(a) - p_i(b)||_1 / (2 alpha)` and `ecs` is its mean over elements. `alpha` cancels for partitions
(the alpha = 0.5 and 0.9 values agree to 1e-16), which is the closed form our function uses.

## Ours vs reference

Ours = `spatial_anno_metrics.element_centric_similarity(a, b)` at `f71512f` (v0.6.0, default `alpha=0.9`).

| pair | clusim alpha=0.9 | clusim alpha=0.5 | ours | ours - clusim |
|---|---|---|---|---|
| identical | 1 | 1 | 1 | 0.0e+00 |
| relabelled | 1 | 1 | 1 | 0.0e+00 |
| split | 0.83999999999999986 | 0.83999999999999997 | 0.83999999999999997 | 1.1e-16 |
| merged | 0.9375 | 0.9375 | 0.9375 | 0.0e+00 |
| random | 0.16213666399788138 | 0.16213666399788146 | 0.16213666399788143 | 5.6e-17 |

## Regenerate

```bash
uv venv <venv> --python 3.12
uv pip install --python <venv>/bin/python clusim==0.4 pandas   # pandas: for our function only
<venv>/bin/python benchmarks/parity/ecs_clusim/make_golden.py
```

Used here: Python 3.12.13, clusim 0.4, numpy 2.5.4, scipy 1.18.1, python-igraph 1.0.0, networkx 3.7,
mpmath 1.4.1, pandas 3.1.0. The script rewrites this directory and prints the table above.
