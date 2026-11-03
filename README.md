# spatial-anno-metrics

**Cell-type annotation quality metrics — spatial + scRNA-seq.** A small, dependency-light
(`numpy` / `pandas` / `scikit-learn` / `scipy` / `anndata`; `squidpy` optional) battery for the
question *"is this cell-type annotation any good?"* — organised by **what evidence you have**,
because in real projects there is almost never a ground truth.

It consolidates metrics from several sources into one place, with attribution (see
[**the catalog**](docs/cell-annotation-quality-metrics.md) for the full reference, decision guide,
default metric batteries, and per-metric definitions).

## Install

Not on PyPI yet - install from GitHub:

```bash
pip install "git+https://github.com/p-gueguen/spatial-anno-metrics"   # core (no spatial)
pip install "spatial-anno-metrics[spatial] @ git+https://github.com/p-gueguen/spatial-anno-metrics"   # + squidpy (Moran's I)
```

## What's in it

### Reference-free (the default regime — no ground truth)

| Function | Metrics | From |
|---|---|---|
| `internal_validity(adata, label_key, embedding="X_pca")` | silhouette, 2label-silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch, avg similarity, **Davies-Bouldin + Calinski-Harabasz**, `integrated` (geometric mean) | scTypeEval + sklearn |
| `inter_sample_consistency(adata, label_key, sample_key)` | **ISC on two axes** (Garnica et al. 2026): `iscg` (2-label silhouette on pseudobulk cosine = cross-sample compactness), `iscl` (silhouette on reciprocal-classification match = cross-sample replicability), `integrated` = product, per-label `regime` via `isc_regime` (consistent / merge / split / remove), `dropped_types` below 5 samples × 10 cells | scTypeEval / Garnica et al. 2026 |
| `marker_program_fidelity(adata, label_key, marker_sets)` | per-type AUC-ROC + Cohen's d of the marker-enrichment score | Zhu et al. 2026 |
| `marker_gene_overlap(adata, label_key, marker_sets, n_top=25)` | fraction of a type's curated markers landing in its top-N data-derived DE genes (the fidelity DE-overlap view) | Zhu et al. 2026 |
| `panel_resolvability(reference, label_key, panel_genes)` | can a target panel resolve a reference's types **at the panel's own depth**? (depth-thinned CV classifier → per-type F1 + confuser) | — |
| `crisp_purity(adata, marker_sets)`, `mecr(adata, marker_sets)`, `pmp(adata, marker_sets, label_key)` | marker co-expression **contamination purity**: CRISP per-cell impurity + dataset purity, MECR (mutually-exclusive co-detection rate over lineage pairs; a marker two sets share is not co-detection, nested sets are skipped), PMP (per-cell marker purity, panel-size-invariant) | SpatialScribe QC L3 |
| `annotation_quality(adata, label_key, marker_sets=, sample_key=)` | orchestrates the reference-free battery into one headline | — |

### Reference-based purity (needs a matched scRNA/snRNA reference)

| Function | Metric | From |
|---|---|---|
| `nmp(adata, reference, label_key)` | negative-marker proportion: fraction of a cell's counts on genes its assigned type should NOT express (detected in ≤5% of its reference cells but ≥25% of another type's) — the reference-based counterpart to `pmp`; `types_without_negatives` lists the types it cannot score | SpatialScribe QC L3; negative rule from SegTraQ |
| `ncp(adata, reference)` | non-coexpression preservation: fraction of reference-defined non-coexpressed gene pairs that stay non-coexpressed in the section (MECR-style, reference-derived pairs) | SpatialScribe QC L3 |

### Spatial label quality (does the label make spatial sense? — needs coords)

| Function | Metric | From |
|---|---|---|
| `spatial_coherence(adata, label_key, k=15)` | per-cell fraction of spatial neighbors sharing the label + dataset **PAS** (Proportion of Abnormal Spots). Down-weight only — rare infiltrating cells are legitimately incoherent | SpatialScribe QC L6 |
| `neighborhood_sanity(adata, label_key)` | flags types whose strongest neighborhood-enrichment partner is a *different* type (a mislabeling smell) via a squidpy permutation null | SpatialScribe QC L6 |
| `neighbor_contamination(adata, marker_sets, label_key)` | per-cell share of counts on markers of the **adjacent** cell types (not its own) + a directed source → target matrix: spillover from neighbors, which `pmp` / `nmp` cannot tell from a mislabel. It cannot tell spillover from real co-expression either: compare with same-type cells that have no such neighbor | SegTraQ `contamination_strength` |

### Per-cell confidence

| Function | Metric | From |
|---|---|---|
| `per_cell_confidence(adata, marker_sets)` | marker-score **margin** (top1−top2) + softmax **entropy** — the margin, not the absolute score, is the confidence | SpatialScribe QC L5 |
| `label_stability(adata, annotate_fn, drop_frac=0.2)` | label-flip rate under transcript subsampling; annotator-agnostic (pass any labeller) | SpatialScribe QC L5 |
| `conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha=0.1)` | split-conformal **prediction sets** from any annotator's class probabilities: set size 1 = confident, > 1 = ambiguous, 0 = novel. Pass `class_conditional=True` (Mondrian; default False) to keep coverage on rare types; an optional `hierarchy` collapses an ambiguous subtype set to a lineage call. Coverage holds only when calibration and query cells are exchangeable (a dissociated reference vs a spatial section is not) | split conformal (LAC) |

### Signal QC (annotation-independent — score the *data*, not the labels)

| Function | Metric | From |
|---|---|---|
| `moran_signal(adata)` | **Moran's I gene-vs-control** (real genes spatially autocorrelated, controls I≈0; the separation is a signal-quality score) — via `squidpy.gr.spatial_autocorr` | SpatialQM (CSO) |
| `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, `run_signal_qc` | SNR vs controls, zero fraction, per-cell detection entropy, transcript density | SpatialQM (CSO) |

### External (needs ground truth — benchmarks only)

| Function | Metrics | From |
|---|---|---|
| `external_scores(pred, truth)` | balanced accuracy, macro/weighted/per-class F1, ARI/AMI/NMI, Cohen's kappa, ECS | classic + Gates & Ahn 2019 |
| `element_centric_similarity(a, b)` | exact vectorized ECS for hard partitions (resolution-agnostic) | Gates & Ahn 2019 |
| `hierarchical_accuracy(pred, truth, hierarchy)` | partial credit along the lineage tree | Zhu et al. 2026 |
| `avg_bio(adata, label_key, truth_key, embedding)` | AvgBIO = mean(ARI, NMI, scaled-ASW) — biological-structure conservation | Zhu et al. 2026 / scIB |
| `composition_accuracy(pred, truth)` | hard-label global proportions: L1 / Pearson / JSD | — |
| `deconvolution_metrics(true_prop, pred_prop)` | soft per-spot proportion matrices: **R² / JSD** (+ RMSE, per-spot JSD) | OpenProblems `task_spatial_decomposition` |

## Quick start

```python
import spatial_anno_metrics as sam

# Reference-free: you have labels + an embedding (the common case)
sam.internal_validity(adata, label_key="cell_type", embedding="X_pca")
# -> {'silhouette': .., 'neighborhood_purity': .., 'integrated': .., ...}

# Reference-free fidelity: you have marker sets but no truth
sam.marker_program_fidelity(adata, "cell_type", {"T cell": ["CD3D", "CD3E", "TRAC"], ...})

# Contamination purity: reference-free from curated markers (spillover / spatial doublets)
sam.crisp_purity(adata, marker_sets)   # dataset purity + obs['crisp_impure']; sam.mecr(adata, marker_sets)

# Signal QC: is the section trustworthy at all? (needs coords + control probes)
sam.moran_signal(adata)          # median gene Moran's I >> control -> real signal

# Benchmark with ground truth
sam.external_scores(pred_labels, truth_labels)     # macro_f1, ARI, kappa, ECS, ...
sam.deconvolution_metrics(true_prop, pred_prop)    # R2, JSD (OpenProblems)
```

## Choosing metrics

Don't compute everything — pick the battery for your situation. The
[decision guide](docs/cell-annotation-quality-metrics.md#1-decision-guide) maps
*what you have* (ground truth? replicates? markers? coords?) to *what to compute*. Two rules that
save you: **margin > absolute score** for confidence, and keep **≥1 method-independent metric**
(marker-only) so you're not circular with a reference.

With ≥5 biological replicates, read ISC on **two axes** rather than one number: low ISCL alone → merge
(over-partitioned), low ISCG alone → split (too coarse), both low → remove; the product ranks annotators
or reference atlases without ground truth (R² > 0.98 with macro-F1 in the 2026 benchmark). Labels with
fewer than 5 samples × 10 cells are reported as not computable, and malignant / patient-private
populations score low by construction - leave them out of the verdict. Details and the benchmark's
design lessons: catalog §3b and §7.

## Sources & attribution

The metrics are consolidated, with attribution, from:

- **scTypeEval / ISC** — Garnica, Andreatta & Carmona, "Evaluating cell type annotations in single-cell omics in the absence of ground truth," bioRxiv 2026, [10.64898/2026.06.10.731285](https://doi.org/10.64898/2026.06.10.731285); package [github.com/carmonalab/scTypeEval](https://github.com/carmonalab/scTypeEval) (Bioconductor); benchmark code [carmonalab/ISC_benchmark_reproducibility](https://github.com/carmonalab/ISC_benchmark_reproducibility). Source of ISCG / ISCL / RCM, the 0.74 / 0.45 guide thresholds and the four refinement regimes; this package is a pseudobulk-level Python re-implementation (Spearman nearest-profile classifier in place of SingleR). `[preprint]`
- **Zhu et al. 2026** — "Benchmarking cell type annotation in spatial transcriptomics," bioRxiv [10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full) (ECS, AvgBIO, hierarchical accuracy, marker-program fidelity; accuracy ≠ fidelity).
- **SpatialQM** — Center for Spatial OMICs, [github.com/Center-for-Spatial-OMICs/SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM) (Moran's I gene-vs-control + signal-quality metrics); published with the Spatial Touchstone study, Plummer et al., *Nat Biotechnol* 2025, [10.1038/s41587-025-02811-9](https://doi.org/10.1038/s41587-025-02811-9).
- **OpenProblems** — [task_spatial_decomposition](https://github.com/openproblems-bio/task_spatial_decomposition) (deconvolution R² / JSD).
- **SegTraQ** — Lazic et al. (EMBL / UZH), segmentation and transcript-assignment QC, [github.com/LazDaria/SegTraQ](https://github.com/LazDaria/SegTraQ) (MIT). Source of `neighbor_contamination` (its `contamination_strength`, re-implemented vectorized on an AnnData), the `nmp` negative-marker rule (`markers_from_reference`, without its DE step) and the shared-marker exclusion in MECR / CRISP. Its segmentation modules (border admixture, transcript positions, z-layer overlap) need transcripts + polygons: run SegTraQ itself on the SpatialData. `[tool]`
- Classic metric definitions: ARI (Hubert & Arabie 1985), AMI/NMI (Vinh et al. 2010), silhouette (Rousseeuw 1987), ECS (Gates & Ahn 2019).

This package re-implements/consolidates metrics from those tools in Python; please **cite the
original sources** in publications, not just this package.

## Status

`v0.6.0`, experimental. Several thresholds/metrics in the catalog are flagged provisional —
version accordingly. Tests: `pip install -e ".[spatial,dev]" && pytest`.

## License

MIT.
