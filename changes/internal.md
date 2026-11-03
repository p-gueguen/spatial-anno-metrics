# internal (package wp/internal)

## CHANGELOG

### Changed
- `internal_validity`: Ward-PropMatch runs on the `subsample` (as the silhouettes, Davies-Bouldin and
  Calinski-Harabasz already did), not on every cell. Ward on every cell is O(N^2) in memory and
  time; it now runs on 1.64e5 Xenium cells in 20 s at 0.74 GB peak RSS (was
  2,951 s / 211 GB). **Breaking** when `n_obs > subsample`: `ward_propmatch` and `integrated` move (synthetic sweep,
  2,400 cells vs `subsample=2000`: `ward_propmatch` by up to 0.094, `integrated` by up to 0.014; the other 36 sweep
  metrics are unchanged). old -> new: Ward on every cell -> pass `subsample >= n_obs`, which reproduces every
  0.6.x / wp-base output to 4e-16 (relative).
- `internal_validity`: neighbourhood purity and average similarity are vectorised; orbital-medoid uses an (n, K)
  `cdist` instead of an (n, K, D) broadcast. Values equal to 4e-16.
- `internal_validity`: when the subsample has fewer than 2 labels or only singletons, `silhouette`,
  `ward_propmatch`, `davies_bouldin` and `calinski_harabasz` are NaN, and `silhouette_2label` is NaN when no label
  has 2 cells. **Breaking**: the silhouettes were a made-up 0.0 there (or sklearn raised `ValueError` for
  singletons). The new `nan_reasons` dict names each NaN metric and why; `integrated` is the geometric mean of the
  defined parts. old -> new: `silhouette == 0.0` -> `isnan(silhouette)` with a `nan_reasons["silhouette"]`.
- `label_stability`: each rep hands `annotate_fn` a new AnnData (thinned counts in `X` and `layers[layer]`,
  copies of `obs`/`var`, shallow copies of `uns`/`obsm`/`obsp`) instead of `adata.copy()`. Other layers and
  `raw` are no longer carried: they were the memory cost, and a labeller that read them saw unthinned data.
  Flip rates are unchanged (0 difference on 3 synthetic sections and in the sweep).
- `_utils._embedding` (the PCA fallback of `internal_validity`, `avg_bio`, `isolated_label_scores`, `clisi`):
  a sparse `X` gets PCA by implicit mean-centering (arpack, scikit-learn >= 1.4) instead of a dense copy. Same
  components as the dense PCA: on 6,000 Xenium cells (313 genes) pairwise distances in the 50-PC embedding differ
  by <= 8.1e-14 relative and `internal_validity` outputs by <= 7.4e-16. scikit-learn < 1.4, or
  `n_components >= min(X.shape)`, densifies as before.

### Added
- `isolated_label_scores(adata, label_key, *, batch_key=None, embedding="X_pca", iso_threshold=None,
  subsample=2000, seed=0)`: scib isolated-label ASW and F1 (Luecken et al. 2022, Nat Methods,
  10.1038/s41592-021-01336-8). A label is isolated when it occurs in the fewest batches; without `batch_key`
  every label is, giving macro averages over labels (a rare label weighs as much as an abundant one). ASW is a
  port of scib-metrics 0.5.6 `isolated_labels` (BSD-3) and matches it to 1e-5 (float32 reference) on 3 fixtures.
  F1 follows scib's best-cluster rule, but sweeps Ward cuts into 2 .. 2 x n_labels clusters instead of Leiden
  resolutions (no new dependency), so F1 is comparable in kind, not to the digit. Above `subsample` cells each
  label keeps at most `subsample // n_labels` cells.
- `clisi(adata, label_key, *, embedding="X_pca", n_neighbors=90, perplexity=None)`: scib cLISI,
  `(n_labels - median LISI) / (n_labels - 1)`, port of scib-metrics 0.5.6 `clisi_knn` (JAX -> numpy, BSD-3);
  matches it to 1e-5 on a fixture. 1.64e5 cells: 26 s, 1.4 GB peak RSS.
- `internal_validity(..., embedding_label_trained=False)`: True adds a `warning` to the result - internal
  metrics on an embedding trained with the labels are circular (Metric Mirages, 10.1101/2024.04.02.587824).
- `THRESHOLDS["internal"]["clisi_n_neighbors"]` (90, scib-metrics) and `["isolated_f1_cuts_per_label"]` (2).
- Tests: `avg_similarity` asserted; Ward sees only the subsample; NaN silhouettes with reasons; the
  label-trained flag; sparse PCA without densifying; a planted rare label that `isolated_label_scores` catches and
  cLISI misses (isolated ASW 0.53, F1 0.32; cLISI 1.00); scib-metrics parity values; `label_stability` negative
  control (random labeller: flip rate 0.4-0.6) and the slim per-rep object.

### Integrator (files not mine)
- `tests/test_validate.py::_NAN_LABEL_CALLS`: add
  `"isolated_label_scores": lambda sam, a: sam.isolated_label_scores(a, "cell_type")` and
  `"clisi": lambda sam, a: sam.clisi(a, "cell_type")`.
- `tests/test_thresholds.py::DEFAULTS`: add `("internal", "clisi", "n_neighbors"): ("internal", "clisi_n_neighbors")`
  and `("internal", "isolated_label_scores", "subsample"): ("internal", "subsample")`.
- `sam.__all__` grows from 33 to 35 names (the contract text says 33).

## README rows

| `internal_validity(adata, label_key, embedding="X_pca")` | silhouette, 2label-silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch, avg similarity, **Davies-Bouldin + Calinski-Harabasz**, `integrated` (geometric mean of the defined parts); the O(N^2) parts on a `subsample` (runs at 1e5+ cells); NaN + `nan_reasons` when undefined; `embedding_label_trained=True` flags circularity | scTypeEval + sklearn |
| `isolated_label_scores(adata, label_key, batch_key=None)` | rare-label separation: isolated-label ASW (per-label mean silhouette) and F1 (best cluster of a Ward-cut sweep); without `batch_key`, macro averages over labels | scib / scib-metrics (Luecken et al. 2022) |
| `clisi(adata, label_key)` | cLISI: label purity of perplexity-weighted kNN neighbourhoods, scaled to [0, 1] (median-based, blind to small labels) | scib-metrics / Harmony LISI |
| `label_stability(adata, annotate_fn, drop_frac=0.2)` | label-flip rate under transcript subsampling; annotator-agnostic (pass any labeller); the labeller gets a slim AnnData per rep (thinned counts, no other layers) | SpatialScribe QC L5 |

## Catalog

§3a (Cluster-validity / separation), add two rows after "Calinski–Harabasz":

| **Isolated-label ASW / F1** | for labels seen in the fewest batches (every label without batches): the label's mean silhouette (rescaled to 0-1), and its best F1 against any cluster of a clustering sweep. Macro over labels, so a small label counts fully | 0–1, higher | **scib** (Luecken et al. 2022); `isolated_label_scores` (F1 sweeps Ward cuts, not Leiden) |
| **cLISI** | inverse Simpson index of labels in each cell's perplexity-30 kNN neighbourhood; (n_labels − median) / (n_labels − 1) | 0–1, higher | Harmony LISI; **scib**; `clisi` |

§3a, after the "Interpreting on sparse panels" paragraph: "`internal_validity` computes neighbourhood purity,
orbital-medoid and average similarity on every cell and the O(N^2) parts (silhouettes, Ward-PropMatch,
Davies-Bouldin, Calinski-Harabasz) on a `subsample` (default 2,000; 5,000 inside `annotation_quality`). An
undefined silhouette (one label or only singletons in the subsample) is NaN with a reason, never 0."

§8 (Cross-cutting rules), new bullet after "Internal metrics are consistency, not correctness":
"- **Internal metrics on a label-trained embedding are circular** (scANVI, a fine-tuned foundation model): the
  embedding was optimised to separate these labels, so silhouette / kNN purity / cLISI / ISC score the training,
  not the annotation. Compute them on a label-agnostic embedding (PCA of the counts). `internal_validity(...,
  embedding_label_trained=True)` flags it. (Metric Mirages, doi 10.1101/2024.04.02.587824, preprint.)"

§8, new bullet: "- **cLISI hides a small label, and the silhouette does not say which label is wrong**: an 8-cell
  label planted inside another type (188 cells) leaves cLISI at 1.00 (a median over cells); the silhouette drops
  (0.86 -> 0.53) without naming the label, while that label's isolated-label ASW is 0.53 (0.98 when it stands
  apart) and its F1 0.32 (1.00). Read `isolated_label_scores` when rare types matter."

## Roadmap

- #3 (`internal_validity` at >= 1e5 cells) -> **done**. Ward on the subsample, purity / avg-similarity vectorised,
  orbital-medoid via `cdist` (`src/spatial_anno_metrics/internal.py`, `internal_validity`); test
  `tests/test_internal.py::test_internal_validity_runs_ward_on_the_subsample_only`. Xenium breast section
  (Janesick 2023, 50-PC PCA, 12 KMeans labels), SLURM, 8 CPUs, same nodes before/after:
  1e5 cells 771 s / 78.7 GB -> 8.7 s / 0.70 GB; 1.64e5 cells 2,951 s / 211 GB -> 20 s /
  0.74 GB (peak RSS of the whole process, h5ad load ~0.5-0.7 GB included). At both sizes every other output is
  identical; `ward_propmatch` moves 0.703 -> 0.746 (1e5) and 0.718 -> 0.732 (1.64e5), `integrated` by <= 0.008. `subsample >= n_obs` reproduces wp-base
  to 4e-16. (The `ncp` pair loops named in the item belong to purity, not this package.)
- #8 (single-control test gaps) -> **partial**: `avg_similarity` asserted
  (`test_internal_validity_reports_avg_similarity_clean_above_shuffled`), `label_stability` negative control
  (`test_label_stability_negative_control_random_labeller_flips_often`). The other gaps belong to other modules.
- #15 (dense hotspots) -> **partial**: `_embedding` sparse PCA (`test_embedding_pca_of_sparse_x_matches_dense_without_densifying`)
  and no `adata.copy()` in `label_stability` (`test_label_stability_hands_the_labeller_a_slim_object_not_a_full_copy`)
  done. Still open in this package: `label_stability` thins a dense float32 copy of the counts per rep
  (1e6 x 5k = 20 GB); keeping it sparse changes the RNG stream and what `annotate_fn` receives, so it was left.
  Others (fidelity, pseudobulk, purity helpers) belong to other packages. Deviation: the item says
  `TruncatedSVD`; that does not centre, so it would give different geometry from the dense PCA path for the same
  data - arpack PCA with implicit centering gives the same embedding (8e-14).
- #23 (scib `isolated_labels` / cLISI) -> **done** with one deviation: `isolated_label_scores` (ASW = scib-metrics
  to 1e-5; F1 by Ward-cut sweep instead of Leiden) and `clisi` (= scib-metrics to 1e-5), tests
  `test_isolated_label_scores_*`, `test_clisi_separates_clean_from_shuffled_and_matches_scib_metrics`.
  scib-metrics is BSD-3 (Copyright 2022 Adam Gayoso); ported code carries that attribution.
- Catalog §8 rule "internal metrics on a label-trained embedding are circular" -> **done** in code
  (module docstring, `internal_validity(embedding_label_trained=)`, test
  `test_internal_validity_flags_a_label_trained_embedding`); catalog text above.
- Foundation leftover "silhouette fabricated 0.0" -> **done** (`test_internal_validity_silhouettes_are_nan_with_a_reason_when_undefined`).
