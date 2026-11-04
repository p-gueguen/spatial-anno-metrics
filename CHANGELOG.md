# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/) (before 1.0, a minor version can change the API).

## [Unreleased] - 0.7.0

0.7.0 makes the API consistent before 1.0: one result shape, one name per concept, explicit count matrices.
It also aligns several metrics with their published definitions (NMP, PAS, ISCG, hierarchical credit, panel
tiers) and adds 12 metrics.
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

- Count-based metrics take a keyword-only `layer=` and name the matrix they read in the result
  (`"matrix": "layers['counts']"` or `"X"`): `crisp_purity`, `mecr`, `pmp`, `nmp`, `ncp`, `neighbor_contamination`,
  `per_cell_confidence`, `label_stability`, `marker_program_fidelity`, `marker_gene_overlap`, `moran_signal`,
  `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, `run_signal_qc`, `inter_sample_consistency`
  and `panel_resolvability`. `per_cell_confidence(adata, markers, "counts")` and a positional `layer` in
  `label_stability` -> `layer="counts"`.
- `panel_resolvability` always thins the counts chosen by `layer` to `target_depth`. Before, it guessed from the
  first 200 rows whether `X` held counts and, when `X` was log-normalised and no `counts` layer existed, scored
  it at full depth without a warning: the inflated macro F1 its docstring warns about (on a simulated reference
  at `target_depth=3`, 1.00 instead of 0.65). That input now raises `ValueError` asking for `layer=`.
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
  - `marker_program_fidelity` and `marker_gene_overlap` when no type has a marker on the panel (was status
    "ok" with `None` means and no reason);
  - `neighbor_contamination` with fewer than 2 labels that have an on-panel marker set, or with no cell next to
    another marker-set type (was status "ok", a NaN mean and an all-NaN `obs` column);
  - with fewer than 2 labels: `internal_validity`, `marker_program_fidelity`, `marker_gene_overlap`,
    `spatial_coherence` and `neighborhood_sanity` (`neighborhood_sanity` raised).
- `marker_program_fidelity`: Cohen's d of a type whose marker load is constant is `None` (was 0.0) and is left
  out of `mean_cohens_d`.
- `composition_accuracy`: `pearson` is NaN when undefined, i.e. one type or equal shares on either side (was 1.0).

**Metric definitions and defaults** (each with the call that restores the 0.6.x value)

- `nmp`: the default `score` is now `"salas"`, Salas et al. 2025 negative marker purity as txsim
  `negative_marker_purity_reads` computes it (`nmp`, higher = cleaner, scaled so 0 = the mean of 20 label
  permutations; writes nothing to `adata`). `score="moldia"` computes the older Moldia `Xenium_benchmarking`
  reads criterion (`x_sc < max_share / n_types`, relative to the mean across types; matched to 1e-12 on the
  golden), and `score="cells"` computes the cell-based variant (`negative_marker_purity_cells` in txsim and
  Moldia `Xenium_benchmarking`; matched to 1e-12 on the golden). The 0.6 score is
  `score="raw"` (`mean_nmp`, lower = cleaner, writes
  `obs['sam_nmp']`); on the golden it reproduces 0.001030 / 0.004628 exactly.
- `nmp(score="raw")`: cells of a type absent from the reference, of a type with no negative gene, or with no counts
  on the shared genes are NaN, not 0.0; `mean_nmp` is the mean over the defined cells (new `n_defined`), and a
  section with no scorable cell is a skip. `types_without_negatives` -> two lists, `types_absent_from_reference`
  and `types_without_negatives` (in the reference, no gene passes). Over-split / merged labels no longer pull
  `mean_nmp` towards 0 (synthetic `over_split` level 3: 0.000565 -> 0.001349, clean 0.0014).
- `spatial_coherence`: PAS follows Shang & Zhou 2022 (SpatialPCA `fx_PAS`): abnormal = the label differs from
  >= 6 of the 10 nearest neighbours. `n_neighs=15 -> 10`, `pas_threshold=0.2 -> 0.4`, and the comparison is
  inclusive (`frac <= pas_threshold`). `spatial_coherence(adata, n_neighs=15, pas_threshold=0.2)` reproduces 0.6.x
  except for cells with exactly 20 % same-label neighbours, which now count as abnormal (Xenium breast, 167,780
  cells: PAS 0.153952 -> 0.205150, exactly the 8,590 cells at 0.2). The new PAS reads higher and responds more
  (synthetic `label_noise` 0.4: 0.344 -> 0.553; `coord_shuffle` 1.0: 0.512 -> 0.976).
- `inter_sample_consistency`: ISCG is computed in scTypeEval's default space, cosine on the first 30 centred PCs of
  the pseudobulks (`reduction="pca"`, `n_pcs=30`, capped at `min(n_profiles, n_genes) - 1`; exact SVD). It matches
  scTypeEval 1.0.0's ISCG within 1.4e-6 on the Yerly 2022 golden. Values move by +0.012 to +0.073 per label on the
  golden and by -0.012 to +0.025 on the synthetic sweep. Raw-space ISCG: `reduction=None` (0.6.x to 2.3e-14).
- `inter_sample_consistency`: headline `iscg` / `iscl` are means over labels (paper, scTypeEval
  `mean_over_types`), like `integrated`; they move when labels have unequal sample counts (golden ISCG 0.855469 ->
  0.854361). The 0.6 pseudobulk means are `iscg_profile_mean` / `iscl_profile_mean`.
- `inter_sample_consistency`: the reciprocal-classification match follows scTypeEval's `match_pairs` rule: a
  mutual best hit between two labels counts only when both labels occur in both samples, and a sample with fewer
  than 2 scored profiles matches nothing. This moves ISCL only when a scored label is missing from some samples and
  is part of a cross-label mutual hit (0 values moved on the golden and the synthetic sweep).
- `internal_validity`: Ward-PropMatch runs on the `subsample` (as the silhouettes, Davies-Bouldin and
  Calinski-Harabasz already did), not on every cell: O(N^2) memory made it unusable at 1e5 cells (1.64e5 Xenium
  cells: 2,951 s / 211 GB -> 20 s / 0.74 GB). When `n_obs > subsample`, `ward_propmatch` and `integrated` move
  (synthetic sweep, 2,400 cells: by up to 0.094 and 0.014). `subsample >= n_obs` reproduces 0.6.x to 4e-16.
- `internal_validity`: when the subsample has fewer than 2 labels or only singletons, `silhouette`,
  `ward_propmatch`, `davies_bouldin` and `calinski_harabasz` are NaN (and `silhouette_2label` when no label has 2
  cells); they were a made-up 0.0, or sklearn raised `ValueError` for singletons. A new `nan_reasons` dict names
  each NaN metric and why; `integrated` is the geometric mean of the defined parts.
- `hierarchical_accuracy`: the flat partial credit becomes credit by LCA depth. The default `credit="lca"` is the
  per-cell hierarchical F-measure (Kiritchenko et al. 2006): an honest coarse call (the parent) scores 2/3, a wrong
  sibling subtype 1/2 (both scored 0.5 before), a wrong lineage 0. On a two-level map nothing else changes.
  `credit="lineage", partial=0.5` reproduces the old scores exactly. `hierarchy` may be any depth.
- `panel_resolvability`: tiers use the chance-corrected F1, `(F1 - prevalence) / (1 - prevalence)`: `resolvable`
  at >= 0.8 and `not_resolvable` below 0.5 (was raw F1 >= 0.5 / < 0.3; the confusion cuts 0.25 / 0.30 are
  unchanged). With no private markers at depth 50 (macro-F1 0.71) `frac_resolvable` falls from 1.00 to 0.11.
  `f1_resolvable=0.5, f1_weak=0.3, null_adjusted=False` reproduces the old tiers exactly; `macro_f1` and the
  per-type `f1` are unchanged.
- `signal_to_noise` (default `method="cell"`) no longer returns the eps floor: per-cell ratios are defined only for
  cells with a control count (the others are NaN, not log2(x / 1e-9)), the medians and
  `frac_cells_signal_below_control` are taken over those cells, and the metric skips when no cell has a control
  count. On the Xenium crop the median moves 27.956 -> 6.491 (224 of 2,000 cells have a control count). `eps=1e-9`
  reproduces 0.6 exactly.
- Controls are flagged by feature type, never by name: `var['control']` wins when present, otherwise every feature
  whose `var['feature_type']` (or `feature_types`) is not `"Gene Expression"` is a control. Before, without
  `var['control']` there were no controls.
- `neighbor_contamination`: `n_source_target` is NaN, not 0, for a pair with no evaluable target cell (as in
  SegTraQ), and a nested marker-set pair (`T` inside `CD4 T`) is skipped as hierarchical, as `mecr` skips it (0.6.x
  counted the subtype's private markers as spillover).
- `spatial_coherence`, `neighborhood_sanity`, `neighbor_contamination`: a spatial graph cached by another call with
  other settings (another `n_neighs`, `moran_signal`'s graph, a user's squidpy graph) is no longer reused silently;
  the requested kNN graph is built (and not cached over the existing one). A graph put in
  `obsp['spatial_connectivities']` by hand (no `uns['spatial_neighbors']`) is still used as it is.
- `label_stability`: each rep hands `annotate_fn` a new AnnData (thinned counts in `X` and `layers[layer]`, copies
  of `obs`/`var`, shallow copies of `uns`/`obsm`/`obsp`) instead of `adata.copy()`. Other layers and `raw` are no
  longer carried (a labeller that read them saw unthinned data). Flip rates are unchanged.
- `conformal_prediction_sets`: `class_conditional` accepts False, True or `"clustered"`; any other value raises
  `ValueError`.

### Added

- Each function that writes `obs` also keeps its scalar results in `adata.uns["sam"][<function name>]`, with
  the column name(s) under `key_added`. Frames, per-type dicts and lists stay in the returned dict only, so the
  AnnData still saves to h5ad.

- `from_spatialdata(sdata, table_key="table")` and `segtraq_covariates(adata_or_sdata, table_key="table")`
  (new module `io`): the AnnData table of a SpatialData object and SegTraQ's transcript-level admixture covariates
  from `obs`. Needs the optional `spatialdata` package for `SpatialData` inputs.
- `THRESHOLDS` and `thresholds_for(module, **overrides)` (new module `thresholds`): every hard-coded constant
  (ISC 0.74 / 0.45, panel F1 tiers, kNN sizes, NMP / NCP rules, the PAS threshold, the SNR floor, ...) as
  `{"value", "source", "provisional"}`. The functions take their defaults from it; the values are unchanged.
- `spatial_anno_metrics._hierarchy`, for hierarchy-aware metrics: `ancestors(label, parent_map)` and
  `lca_depth(a, b, parent_map)`.
- `spatial_anno_metrics._validate`, shared by the metrics: `check_labels`, `get_counts`, `skip` and
  `dictionary_coverage` (which labels have no marker set, and the fraction of cells they cover).

**New metrics**

- `isolated_label_scores(adata, label_key, *, batch_key=None, embedding="X_pca", iso_threshold=None,
  subsample=2000, seed=0)`: scib isolated-label ASW and F1 (Luecken et al. 2022, Nat Methods,
  10.1038/s41592-021-01336-8). A label is isolated when it occurs in the fewest batches; without `batch_key` every
  label is, giving macro averages over labels. ASW is a port of scib-metrics 0.5.6 `isolated_labels` (BSD-3) and
  matches it to 1e-5; F1 sweeps Ward cuts into 2 .. 2 x n_labels clusters instead of Leiden resolutions (no new
  dependency), so it is comparable in kind, not to the digit.
- `clisi(adata, label_key, *, embedding="X_pca", n_neighbors=90, perplexity=None)`: scib cLISI,
  `(n_labels - median LISI) / (n_labels - 1)`, port of scib-metrics 0.5.6 `clisi_knn` (BSD-3), matches it to 1e-5.
  A planted 8-cell label inside another type leaves cLISI at 1.00; its isolated-label ASW is 0.53 and F1 0.32.
- `coexpression_similarity(adata, reference, *, label_key=None, ref_label_key="cell_type", min_cells=20,
  layer="counts")`: txsim's co-expression similarity (`1 - |rho_ref - rho_section| / 2` over gene pairs, on
  `log1p(1e4 * counts / total)`); with `label_key`, per cell type. Matches txsim to 6.7e-16 on a golden.
- `spatial_ari(*, y_true, y_pred, coords, alpha=0.8, max_cells=5000, seed=0)`: spatially aware Rand index and
  adjusted Rand index (Yan, Feng & Luo, *Biometrics* 2025, 10.1093/biomtc/ujaf127): a disagreeing pair gets partial
  credit by distance. Written from the paper (the authors' R package is GPL and was not copied); reproduces the
  vignette values to all 7 printed digits. `alpha=0` gives the ARI exactly; above `max_cells` a seeded subset is
  scored (O(n^2) pairs).
- `deconvolution_rank_score(results, metrics=("rmse", "jsd_per_spot", "aupr"))`: ranks several deconvolution
  methods on one dataset per metric (1 = best, ties share the mean rank) and averages the ranks (Spotless).
- `decomposition_purity(adata, weights_key, *, key_added="sam_decomposition_purity")`: 1 - the largest normalised
  deconvolution weight per cell, as in Spatial Touchstone; ported from SegTraQ `sp.rctd_contamination` (MIT,
  attributed). Rows with a non-finite or negative weight, or a non-positive total, are NaN.
- `tx_per_cell(adata, *, layer="counts", key_added="sam_tx_per_cell")`: transcripts per cell, mean (= SpatialQM
  `getTxPerCell`) and median; `run_signal_qc` adds a `tx_per_cell` entry.
- `transcript_qc(transcripts, negative_features, *, feature_key="feature_name", cell_key="cell_id")`, from a
  transcripts DataFrame (no spatialdata dependency): `frac_tx_in_cells` (= `getCellTxFraction`; Xenium <= 1.2 `-1`
  and >= 1.3 `"UNASSIGNED"` ids) and `fdr` (negative transcripts / all transcripts x genes / negatives detected;
  SpatialQM's `getGlobalFDR` is `fdr / 100`). Bytes feature names and cell ids are decoded as UTF-8.
- `calibration_error(proba, classes, y_true, *, n_bins=15)`: top-label expected calibration error (Guo et al. 2017,
  15 equal-width bins) and class-wise ECE; class-wise ECE matches values built from sklearn's `calibration_curve`
  to 1e-12.
- `conformal_novelty(X_query, X_ref, *, method="pca", fdr=0.1, n_components=20, frac_train=0.5, seed=0)`:
  split-conformal novelty p-values in the reference's feature space, BH-adjusted (Bates et al. 2023, *Ann Stat*,
  10.1214/22-AOS2244); the score is the PCA reconstruction error or `IsolationForest` (`method="iforest"`).
  Exchangeability matters: thin the reference to the section's per-cell depths first (an un-thinned 400-count
  reference flags 85-100 % of known 150-count cells in the synthetic check).
- `consensus_entropy(adata, label_keys, *, hierarchy=None, min_share=0.5)` (new module `consensus`): per-cell
  consensus across annotators with ontology-path voting (popV, Ergen et al. 2024, *Nat Genet*,
  10.1038/s41588-024-01993-3): every label votes for itself and its ancestors, the consensus is the deepest label
  with more than `min_share` of the votes; flat `entropy` and `ontology_entropy` (votes projected to the consensus
  depth, so a CD4 / CD8 / T split scores 0), both divided by log(number of annotators). Writes nothing to `adata`.
- `ontology_consistency(adata, label_key, obo_path, label_to_cl, *, embedding="X_pca", distance="path",
  rwr_restart=0.5, n_perm=999, min_cells=10, seed=0)`: do the label centroids sit in the embedding the way the
  Cell Ontology relates the labels? Spearman `rho` between centroid distances and `is_a` path lengths (`distance="path"`)
  or Random Walk with Restart stationary-distribution Jensen-Shannon distances (`distance="rwr"`, `scGraph-OntoRWR`,
  Wang et al. *Genome Biol* 2025, 10.1186/s13059-025-03781-6) over label pairs, a one-sided Mantel permutation
  `p_value`, per-label leave-one-out `loo_delta` and `most_inconsistent`. The OBO file is read with the standard
  library.

**New options and result keys**

- `internal_validity(..., embedding_label_trained=False)`: True adds a `warning` (internal metrics on an embedding
  trained with the labels are circular; Metric Mirages, 10.1101/2024.04.02.587824); result key `nan_reasons`.
- `inter_sample_consistency(..., reduction="pca" | None, n_pcs=30, rcm_genes=None | "markers" | <gene list>)`:
  `rcm_genes="markers"` runs the reciprocal classifier on SingleR's classic marker union of each training sample;
  a gene list (e.g. HVGs) restricts it. Result keys `iscg_profile_mean`, `iscl_profile_mean`, `regime_status`
  (`isc_regime` / `regime` are marked experimental), `reduction`, `n_pcs`, `rcm_genes`.
- `crisp_purity` and `mecr` report a detection-matched null (each marker's detections reassigned to random cells,
  `n_null=20` draws, `seed=0`): `crisp_purity_null` / `mecr_null` and `vs_null` (0 = exclusive lineages, 1 = no
  more exclusive than independent markers). `n_null=0` skips it.
- Dictionary coverage: `pmp`, `marker_program_fidelity`, `marker_gene_overlap` and `neighbor_contamination` report
  `labels_without_markers` and `no_dict_frac`; `crisp_purity` and `mecr` report `lineages_without_panel_markers`.
- `spatial_coherence`: `chaos` (Shang & Zhou 2022 CHAOS: mean distance to the nearest same-label cell on per-axis
  standardised coordinates; NaN without `obsm['spatial']`) and `graph`.
- `graph=` (`None` | `"delaunay"` | `"obsp"`) on `spatial_coherence`, `neighborhood_sanity` and
  `neighbor_contamination`; `"delaunay"` is SegTraQ's default graph.
- `neighbor_contamination`: `baseline`, `n_baseline`, `excess` and `p_value` source x target frames (the same
  share in target cells with no source neighbour, the raw value minus it, and a one-sided permutation p; `n_perm=199`,
  `seed=0`), which separate spillover from genuine co-expression; `split_slope` and `split_corr` frames (Bilous et al.
  *Nat Methods* 2026 SPLIT readout: regression slope and Pearson correlation of target cells' contaminant share vs
  the source type's share of their spatial neighbours); `distance_corr` and `distance_curve` (CellSPA readout:
  Spearman correlation and 4-quantile distance curve of contaminant share vs Euclidean distance to the nearest
  source cell); and accepts SegTraQ's marker format `{type: {"positive": [...], "negative": [...]}}` (mixable with
  plain lists).
- `spatial_coherence`, `neighborhood_sanity`, `neighbor_contamination` and `moran_signal` return `status="skipped"`
  with a reason instead of raising `KeyError` when `obsm['spatial']` (and any reusable `obsp['spatial_connectivities']`)
  is missing.
- `get_counts` and `conformal_novelty` reject matrices containing `NaN` or infinite values with a `ValueError`.
- `moran_signal(..., method="spatialqm")`: SpatialQM `getMorans` (per group, cells with a count, log-normalised,
  20 nearest neighbours weighted 1/d), reproducing all 341 per-feature values of the Xenium golden to 5.3e-15;
  needs no squidpy and writes nothing. Result keys `n_genes_excluded`, `n_controls_excluded`, `per_feature`,
  `method`.
- `signal_to_noise(..., method="spatialqm")` (`getMeanSignalRatio` / `getMaxRatio`), keyword-only `eps=`, result
  keys `pooled_mean_signal_ratio`, `n_cells_with_controls`, `frac_cells_with_controls`;
  `detection_entropy(..., method="spatialqm")` (`getEntropy`); `tx_per_area` result key `mean_tx_per_area`.
- `deconvolution_metrics(..., aupr_method="sklearn" | "precrec")` result keys `aupr` (presence/absence average
  precision, or Davis & Goadrich 2006 non-linear PR interpolation when `aupr_method="precrec"`), `aupr_precrec`
  (always reported: Spotless / R `precrec::evalmod` curve area, Sang-aram et al. *eLife* 2024) and `ssim` /
  `ssim_per_type` (per-type SSIM after max-scaling, Li et al. *Nat Methods* 2022; ported from
  QuKunLab/SpatialBenchmarking, BSD-2-Clause).
- `ncp(..., score="pairs" | "moldia", moldia_threshold=0.05)`: `score="moldia"` computes Moldia
  `Xenium_benchmarking`'s `non_coexpression_purity` (`1 - mean(clip(p_B|A - p_B|A^sc, 0, None))` over reference
  negative pairs `p_B|A^sc < moldia_threshold`) and returns `per_gene_as_conditioning` and `per_gene_as_target`,
  reproducing `tests/goldens/coexpr_txsim/expected.json` to 1e-12.
- `hierarchical_accuracy` result keys `credit`, `hierarchical_precision`, `hierarchical_recall`,
  `coarse_call_frac` and `coarse_call_accuracy`.
- `conformal_prediction_sets(..., class_conditional="clustered", *, n_clusters=None, frac_clustering=None, seed=0)`:
  clustered Mondrian conformal (Ding et al. 2023, arXiv 2306.09335), which pools classes with similar score
  distributions; result keys `clusters` and, with `y_query`, `covgap` (Ding et al.'s CovGap as a fraction).
- `panel_resolvability` per-type keys `f1_null` (the type's share of the scored cells) and `f1_adjusted`, and
  keyword-only `f1_resolvable`, `f1_weak`, `null_adjusted`.
- `THRESHOLDS` entries for the new defaults: `internal` (`clisi_n_neighbors`, `isolated_f1_cuts_per_label`), `isc`
  (`n_pcs`, `rcm_marker_n`), `purity` (`null_draws`, `nmp_salas_max_share`, `nmp_salas_min_cells`,
  `nmp_n_permutations`, `coexpr_min_cells`, `coexpr_target_sum`, `ncp_moldia_threshold`), `spatial`
  (`contamination_n_perm`), `external` (`spari_alpha`, `spari_max_cells`), `signal_qc` (`spatialqm_moran_n_neighs`,
  `spatialqm_pseudocount`), `conformal` (`ece_n_bins`, `clustered_quantiles`, `novelty_fdr`, `novelty_frac_train`,
  `novelty_n_components`), and the new blocks `consensus` (`min_share`) and `ontology` (`min_labels`, `min_cells`,
  `n_perm`, `rwr_restart`); `panel` gains `f1_adjusted_resolvable` and `f1_adjusted_weak`.
- `sam.__all__` gains the 12 new metrics and helpers above (41 names).

**Goldens and tests**

- Value tests against reference implementations: `element_centric_similarity` vs clusim 0.4 (1e-12),
  `deconvolution_metrics` vs OpenProblems `task_spatial_decomposition`, Spotless `precrec` AUPR and Li et al. 2022
  SSIM (1e-12), `nmp` (`score="salas" | "moldia" | "cells"`), `ncp` (`score="moldia"`) and
  `coexpression_similarity` vs txsim / Moldia (1e-12), the 30-PC ISCG and the label-mean headline vs scTypeEval
  1.0.0, and every SpatialQM scalar offered vs SpatialQM. The SpatialQM golden records the sha256 of its three
  input fixtures and pins its own. `.gitattributes` locks `tests/goldens/**` and `tests/data/*.obo` to `-text`.
- `tests/data/cl_tiny.obo`: a 14-term Cell Ontology slice for the ontology tests.
- Synthetic and real-data perturbation benchmarks: `benchmarks/synthetic/sweep.py` (5 seeds across all scenarios,
  recording `onto.rwr_rho`, `ncp.moldia`, and `aupr_precrec`) and `benchmarks/real_perturbation.py` (controlled
  perturbations on the committed Yerly 2022 human BCC cohort and 10x Xenium human breast cancer crop).
- Realistic simulation worlds for the synthetic sweep (`benchmarks/worlds/`), each behind the `simulate` /
  `reference` contract of `benchmarks/synthetic/sim.py`, so the 15 scenarios and the scoring are unchanged:
  scDesign3 fitted to the 10x Xenium breast section (one crop, as bands or as six Dirichlet-composition pieces),
  SRTsim fitted to random crops of the same section (one crop set per seed), and scCube (VAE on Yerly 2022 BCC
  scRNA-seq, Xenium breast panel) in four spatial patterns. `sweep.py --world {toy,scdesign3,srtsim,sccube}` (with
  `--pattern`, `--sample-mode`, `--alpha0`, `--fixed-crops`) installs a world through one helper
  (`benchmarks/worlds/__init__.py`); the per-world `run.py` scripts are thin wrappers around it. The default
  `--world toy` writes the same files byte for byte. `benchmarks/worlds/README.md` compares the verdicts across
  worlds (`compare_worlds.py`), and the catalog's caveats on spatial metrics, purity under co-expression and
  spillover, conformal novelty and sibling confusion now cite them.

### Changed

- Minimum versions: `scipy>=1.11` (was 1.10; `scipy.stats.false_discovery_control`) and `scikit-learn>=1.4`
  (was 1.3; sparse-input PCA). Both proven by the lowest-direct CI job.

- `benchmarks/synthetic/sweep.py` uses the 0.7.0 API, with the same metric names in `values.csv`. Warnings are
  recorded instead of silenced. `values.csv` gains an `n_user_warnings` column: the UserWarnings raised by the
  adapter that produced the value. A new `warnings.csv` counts every warning message. `errors.csv` is always
  written and now holds every adapter error, edge cases included, with its full traceback (before, only a
  one-line message).

- Large sparse inputs are no longer densified (same values):
  - `internal_validity`: neighbourhood purity and average similarity are vectorised and orbital-medoid uses an
    (n, K) `cdist` (values equal to 4e-16);
  - `_utils._embedding` (the PCA fallback of `internal_validity`, `avg_bio`, `isolated_label_scores`, `clisi`): a
    sparse `X` gets PCA by implicit centring (arpack, scikit-learn >= 1.4) instead of a dense copy; pairwise
    distances in the 50-PC embedding differ from the dense PCA by <= 8.1e-14 (6,000 Xenium cells). scikit-learn <
    1.4, or `n_components >= min(X.shape)`, densifies as before;
  - `inter_sample_consistency`: pseudobulks are `indicator.T @ X`, sparse in;
  - `marker_program_fidelity` / `marker_gene_overlap` slice marker columns and take column moments on sparse
    matrices; CRISP / MECR keep the marker detections sparse (167,780 cells x 541 genes: fidelity peak memory
    1,132 -> 57 MB, overlap 2,184 -> 650 MB and 12.5 -> 3.3 s);
  - `ncp` pair counts are one sparse `P.T @ P` instead of a Python loop over pairs;
  - `neighbor_contamination` is sparse throughout (400,000 cells x 300 markers: 1.88 -> 0.14 GB,
    10.2 -> 3.5 s with `n_perm=0`);
  - `coexpression_similarity` keeps `log1p(1e4 * counts / total)` sparse in CSR and computes the gene x gene Pearson
    correlation matrix via the sparse Gram matrix `(L.T @ L).toarray() - n * outer(mu, mu)` without densifying
    cells x genes;
  - `conformal_novelty` fits `PCA(svd_solver="arpack")` or `IsolationForest` directly on CSR inputs and computes
    squared PCA reconstruction error via the sparse projection identity `||x - mu||^2 - ||(x - mu) V^T||^2`;
  - `per_cell_confidence` averages marker columns on sparse CSR slices directly, and `label_stability` thins CSR
    `.data` in place (`eliminate_zeros()`) without densifying the count matrix.
- Signal QC: `moran_signal`, `transcript_qc` and every `method="spatialqm"` call return their result and write
  nothing to `uns['sam']`. Docstrings no longer claim SpatialQM ports for different quantities; the default recipes
  say how they differ.
- Docstrings of the spatial metrics say what they detect: spatial coherence / PAS / CHAOS detect random label noise
  only (blind to a compact mislabelled patch, rise when types merge, fall for correctly labelled infiltrating
  cells).
- `benchmarks/synthetic/sweep.py`: the `nmp` adapter reads the Salas score (higher = cleaner).
### Fixed

- `label_stability` floored non-integer values before thinning (`astype(int)`: 0.6 -> 0). It now rounds them,
  and `get_counts` warns that the matrix is not raw counts.
- A squidpy import that fails with an error other than `ImportError` (for example a dask `NotImplementedError`)
  is raised as `ImportError("squidpy failed to import: ...")`.
- `annotation_quality` docstring: without `marker_sets` the AQI is `None`, not a coherence index.
- `__version__` comes from the installed package metadata ("unknown" in an uninstalled source tree).

- Integer marker-set keys: `crisp_purity`, `mecr`, `pmp`, `marker_program_fidelity`, `marker_gene_overlap` and
  `neighbor_contamination` compare `marker_sets` keys as str (like the labels), so integer cluster ids with an
  int-keyed dict score instead of giving NaN (or an all-NaN 0x0 matrix) silently; `0` and `"0"` both present
  raises `ValueError`.
- `neighbor_contamination` and `spatial_coherence` crashed (`AttributeError: ... 'tocsc'`) on a dense ndarray in
  `obsp['spatial_connectivities']`.
- `moran_signal` summaries were NaN whenever a feature had no count (on the crop, 139 features made the medians
  and `separation` NaN and `frac_genes_above_control_p95` 0.0); undefined features are now left out and counted.
- `nmp(score="raw")` detection fractions are detecting cells / type cells, so the inclusive cuts no longer depend on
  storage (with COL1A1 in exactly 11/220 T reference cells a CSR reference gave 0.000 and a dense one 0.200; both
  give 0.200 now).
- Attribution: only the `<= 5 %` cut of the raw NMP rule is SegTraQ's (`markers_from_reference`); the `>= 25 %`
  other-type floor is a package choice (`thresholds.py` source text and `provisional` flag corrected). CRISP is Rose
  et al. 2026 (Santangelo lab), and the segmentation-admixture paper is Mitchel et al. *Nat Genet* 2026.
- `hierarchical_accuracy` scored an honest coarse call ("T" for a CD4 T cell) the same as a confident wrong
  sibling subtype.
- Goldens: `benchmarks/parity/isc_yerly2022/make_golden.py` ran against the removed 0.6 API (re-run:
  `counts.npz` / `obs.csv` byte-identical); `run_sctypeeval.R` calls `scTypeEval:::purge_label` instead of copying
  its regex (`expected.json` byte-identical); the SpatialQM `make_inputs.py outs` gzips the committed `cells.csv`
  byte for byte, so the regenerated `expected.json` is identical on Python 3.10 and 3.12.
### Deprecated

- `spatial_anno_metrics.eval_metrics` was split into the modules `internal`, `isc`, `fidelity`, `external`,
  `panel`, `conformal`, `aqi` and `battery`. Import from the package or from those modules. The old module
  still re-exports every name it defined, and warns once on import.
