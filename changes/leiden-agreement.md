# leiden-agreement (branch feat/leiden-agreement)

## CHANGELOG

### Added
- `leiden_label_agreement(adata, label_key, embedding="X_pca", *, n_neighbors=15, resolutions=(0.2, ..., 2.0),
  seed=0, embedding_label_trained=False)`: scib's Leiden NMI / ARI (Luecken et al. 2022, Nat Methods,
  10.1038/s41592-021-01336-8), the agreement between an unsupervised Leiden clustering of the embedding and the
  labels, as in the SPLIT paper's cell-type separation panel. Mirrors scib-metrics 0.5.6
  `nmi_ari_cluster_labels_leiden` (BSD-3): UMAP connectivities of the 15-NN graph (self included; umap-learn's
  `fuzzy_simplicial_set` ported to numpy, no umap dependency), igraph modularity Leiden at resolutions
  0.2, 0.4, ..., 2.0 with a fresh `random.Random(seed)` per run, and the first resolution with the best NMI;
  `ari` is that clustering's ARI. Returns `{status, nmi, ari, resolution, n_clusters, n_labels, n_neighbors,
  clusters, sweep, n_used}`; `clusters` is the kept clustering, `sweep` the per-resolution values. Given the same
  kNN and seed it equals scib-metrics to the digit (clustering, NMI, ARI and the whole NMI sweep) on 3 fixtures;
  the kNN here is exact (scikit-learn), scib-metrics' benchmark uses approximate pynndescent. Needs igraph, the
  new `leiden` extra; without it the result is a skip naming the extra. igraph's global random generator is
  restored to the `random` module after each run. Every cell is used: 1e5 cells (30-d) in 23 s at 1.2 GB peak
  RSS, most of it the 10 Leiden runs.
- `kmeans_label_agreement(adata, label_key, embedding="X_pca", *, n_init=1, seed=0,
  embedding_label_trained=False)`: scib-metrics' k-means NMI / ARI (`nmi_ari_cluster_labels_kmeans`: k-means++
  with k = the number of labels, `n_init=1`), with scikit-learn's k-means instead of scib-metrics' JAX one, so
  comparable in kind, not to the digit (equal on 1 fixture of 3; NMI 0.776 vs 0.820 and 0.744 vs 0.735 on the
  other two). No optional dependency; 1e5 cells in 2.3 s.
- `[leiden]` optional extra: `igraph>=0.10.3` (0.10.3 renamed `community_leiden`'s `resolution_parameter` to
  `resolution`; 0.10.0-0.10.2 fail).
- `THRESHOLDS["internal"]`: `leiden_n_neighbors` (15), `leiden_resolutions` (0.2..2.0) and `kmeans_n_init` (1),
  all from scib-metrics 0.5.6 (not provisional).
- CI installs the `leiden` extra, and the lowest-direct job gains a `dev,leiden` leg (igraph 0.10.3 with the
  core floors).
- Tests: clean blobs give NMI = ARI = 1 (both functions), shuffled labels NMI < 0.1 and |ARI| < 0.05; the
  igraph-missing skip; determinism under a seed and igraph's generator handed back; NMI / ARI equal
  scikit-learn's on the returned clustering, with the ARI taken at the best-NMI resolution; scib-metrics parity
  values (`benchmarks/parity/leiden_scib/make_expected.py`); the UMAP-connectivity port's invariants.

Not added to `annotation_quality`: like `clisi` and `isolated_label_scores`, the scib-family scores are called
on their own, and the AQI does not read them.

## README rows

Reference-free table, after the `clisi` row:

| `leiden_label_agreement(adata, label_key)` / `kmeans_label_agreement(adata, label_key)` | NMI and ARI between the labels and an unsupervised clustering of the embedding: Leiden on the 15-NN UMAP graph over resolutions 0.2-2.0, best NMI kept (needs the `leiden` extra), or k-means with k = n labels | scib / scib-metrics (Luecken et al. 2022) |

Install section, after the `[spatial]` line:

```bash
pip install "spatial-anno-metrics[leiden] @ git+https://github.com/p-gueguen/spatial-anno-metrics"   # + igraph (Leiden NMI / ARI)
```

## Catalog (`docs/cell-annotation-quality-metrics.md`)

Metric table, after the cLISI row:

| **Leiden / k-means NMI and ARI** | agreement between the labels and an unsupervised clustering of the embedding: Leiden (15-NN UMAP graph, resolution with the best NMI over 0.2-2.0) or k-means (k = n labels) | 0–1, higher | **scib** (Luecken et al. 2022); `leiden_label_agreement` (igraph, `leiden` extra), `kmeans_label_agreement` |

Function list, after `clisi(...)`:

```python
leiden_label_agreement(adata, label_key, *, n_neighbors=15) -> dict  # scib Leiden NMI / ARI (igraph)
kmeans_label_agreement(adata, label_key) -> dict                     # scib k-means NMI / ARI
```

Pitfalls, after the circularity bullet:

- **Leiden NMI / ARI picks its resolution with the labels** (scib's rule): the best of 10 clusterings, so it
  flatters a labeling a little, and a labeling coarser or finer than any resolution reaches is capped below 1.
  Read the `sweep`, and compare labelings at the same settings, never against a number from another dataset.
