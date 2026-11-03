# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/) (before 1.0, a minor version can change the API).

## [Unreleased] - 0.7.0

0.7.0 makes the API consistent before 1.0: one result shape, one name per concept, explicit count matrices.
Calls written for 0.6.0 need the changes listed under "Changed (breaking)", each written `old -> new`.

### Changed (breaking)

**Result shapes**

- Every metric returns a dict with `"status": "ok"`. A metric that cannot score its input returns
  `{"status": "skipped", "reason": "..."}`, sometimes with context keys (counts, dropped types). This replaces
  four skip shapes: `{"status": "skipped: <why>"}`, `{"status": "insufficient_overlap"}`,
  `{"status": "too_few_types"}` and `{"status": "skipped", "reason": ...}`. Two exceptions: `isc_regime` still
  returns a str, and `annotation_quality_index` keeps its result shape.
- `crisp_purity(...)`: float -> `{status, crisp_purity, n_impure, n_marker_positive, matrix}`.
- `mecr(...)`: float -> `{status, mecr, n_pairs, matrix}`.
- `pmp(...)`: `None` -> `{status, mean_pmp, n_defined, matrix}`.
- `element_centric_similarity(a, b)`: float -> `{status, ecs}`.
- `panel_resolvability`: `status` `"insufficient_overlap"` / `"too_few_types"` -> `"skipped"` with a `reason`
  (`n_shared_genes` and `per_type` are still returned).
- `neighbor_contamination`: result keys `matrix`, `n_matrix` (the source x target frames) ->
  `source_target`, `n_source_target`. `matrix` now names the count matrix read, as in every count-based metric.
- `spatial_coherence`: result key `k` -> `n_neighs`.

**Arguments**

- `external_scores(pred, truth)` -> `external_scores(*, y_true, y_pred)` (sklearn's order, keyword-only).
- `hierarchical_accuracy(pred, truth, hierarchy, partial=0.5)` ->
  `hierarchical_accuracy(*, y_true, y_pred, hierarchy, partial=0.5)`.
- `composition_accuracy(pred, truth)` -> `composition_accuracy(*, y_true, y_pred)`.
- `deconvolution_metrics(true_prop, pred_prop)` -> `deconvolution_metrics(*, true_prop, pred_prop)`.
- Embedding kNN size is `n_neighbors`: `internal_validity(k=30)` -> `internal_validity(n_neighbors=30)`, and
  `annotation_quality(k=30)` -> `annotation_quality(n_neighbors=30)`.
- Spatial-graph size is `n_neighs`: `spatial_coherence(k=15)` -> `spatial_coherence(n_neighs=15)`.
- `inter_sample_consistency(counts_layer=None)` -> `inter_sample_consistency(layer=None)`.
- `ncp(adata, reference, ref_label_key, coexpr_threshold, ...)` -> `ncp(adata, reference, coexpr_threshold, ...)`:
  the unused `ref_label_key` is gone, so a positional third argument now sets `coexpr_threshold`.

**Count matrices**

- Count-based metrics take `layer=` and name the matrix they read in the result (`"matrix": "layers['counts']"`
  or `"X"`): `crisp_purity`, `mecr`, `pmp`, `nmp`, `ncp`, `neighbor_contamination`, `per_cell_confidence`,
  `label_stability`, `marker_program_fidelity`, `marker_gene_overlap`, `moran_signal`, `signal_to_noise`,
  `sparsity`, `detection_entropy`, `tx_per_area`, `run_signal_qc` and `inter_sample_consistency`.
- The default is `layer="counts"`. A missing layer falls back to `X` only when `X` holds raw counts (non-negative
  integers); otherwise the metric raises `ValueError` asking for `layer=`. Before, a missing or misspelled layer
  silently fell back to `X`, so log-normalised data was scored as counts. `layer=None` reads `X` as it is, with
  one `UserWarning` when `X` is not raw counts.
- Metrics that read `X` before now read `layers['counts']` when it exists: `crisp_purity`, `mecr`,
  `marker_program_fidelity`, `marker_gene_overlap` and the signal-QC metrics. `per_cell_confidence(layer=None)`
  -> `per_cell_confidence(layer="counts")`.
- `inter_sample_consistency` keeps `layer=None` as its default (`X` averaged as it is, for log-normalised HVGs).
  A named layer selects the summed-counts recipe, as `counts_layer` did.

**Outputs written to the AnnData**

- Every function that writes per-cell values to `obs` takes `key_added=`. The default column names carry a
  `sam_` prefix, so they no longer overwrite generic names such as `pmp` or `sparsity` that you may already use:
  - `pmp`: `obs['pmp']` -> `obs['sam_pmp']`;
  - `crisp_purity`: `crisp_impure` -> `sam_crisp_impure`;
  - `nmp`: `nmp` -> `sam_nmp`;
  - `spatial_coherence`: `spatial_coherence` -> `sam_spatial_coherence`;
  - `neighbor_contamination`: `neighbor_contamination` -> `sam_neighbor_contamination`;
  - `label_stability`: `label_stability` -> `sam_label_stability`;
  - `sparsity`, `detection_entropy`, `tx_per_area`: -> `sam_sparsity`, `sam_detection_entropy`, `sam_tx_per_area`.
  - Two functions write two columns, and their `key_added` is a pair: `per_cell_confidence`
    (`marker_margin`, `marker_entropy` -> `("sam_marker_margin", "sam_marker_entropy")`) and `signal_to_noise`
    (`mean_signal_ratio`, `max_signal_ratio` -> `("sam_mean_signal_ratio", "sam_max_signal_ratio")`).
- `annotation_quality` and `annotation_quality_index` no longer write `obs['pmp']`. The AQI reads the PMP values
  directly; its results are unchanged.
- `neighborhood_sanity` no longer converts your label column to `category` and no longer writes
  `uns['<label_key>_nhood_enrichment']`. It works on a copy of the column, without unused levels; before, an
  unused level changed its result.

**Labels and undefined scores**

- A NaN / None label raises `ValueError` naming the column and the count. Before, pandas >= 3 raised a
  `TypeError` deep inside a sort, and pandas < 3 silently scored the cells as a cell type called "nan".
- Skips instead of made-up values:
  - `crisp_purity` with no marker-positive cell (was 1.0);
  - `mecr` with no lineage pair (was 0.0);
  - `conformal_prediction_sets` with an empty calibration set (was a threshold of 1.0, i.e. every class in every
    set);
  - `per_cell_confidence` with fewer than 2 marker sets on the panel (was margin 0 and entropy 1, written to
    `obs`);
  - `avg_bio` with fewer than 2 truth labels (was an ASW of 0.0);
  - with fewer than 2 labels: `internal_validity`, `marker_program_fidelity`, `marker_gene_overlap`,
    `spatial_coherence` and `neighborhood_sanity` (`neighborhood_sanity` raised).

### Added

- Each function that writes `obs` also keeps its scalar results in `adata.uns["sam"][<function name>]`, with
  the column name(s) under `key_added`. Frames, per-type dicts and lists stay in the returned dict only, so the
  AnnData still saves to h5ad.

- `from_spatialdata(sdata, table_key="table")` (new module `io`): the AnnData table of a SpatialData object.
  Needs the optional `spatialdata` package.
- `THRESHOLDS` (new module `thresholds`): every hard-coded constant (ISC 0.74 / 0.45, panel F1 tiers, kNN sizes,
  NMP / NCP rules, the PAS threshold, the SNR floor, ...) as `{"value", "source", "provisional"}`. The functions
  take their defaults from it; the values are unchanged.
- `spatial_anno_metrics._hierarchy`, for hierarchy-aware metrics: `ancestors(label, parent_map)` and
  `lca_depth(a, b, parent_map)`.
- `spatial_anno_metrics._validate`, shared by the metrics: `check_labels`, `get_counts`, `skip` and
  `dictionary_coverage` (which labels have no marker set, and the fraction of cells they cover).

### Fixed

- `label_stability` floored non-integer values before thinning (`astype(int)`: 0.6 -> 0). It now rounds them,
  and `get_counts` warns that the matrix is not raw counts.
- A squidpy import that fails with an error other than `ImportError` (for example a dask `NotImplementedError`)
  is raised as `ImportError("squidpy failed to import: ...")`.
- `annotation_quality` docstring: without `marker_sets` the AQI is `None`, not a coherence index.
- `__version__` comes from the installed package metadata ("unknown" in an uninstalled source tree).

### Deprecated

- `spatial_anno_metrics.eval_metrics` was split into the modules `internal`, `isc`, `fidelity`, `external`,
  `panel`, `conformal`, `aqi` and `battery`. Import from the package or from those modules. The old module
  still re-exports every name it defined, and warns once on import.
