# spatial-anno-metrics

**Cell-type annotation quality metrics — spatial + scRNA-seq.** A small, dependency-light
(`numpy` / `pandas` / `scikit-learn` / `scipy` / `anndata`; `squidpy` optional) battery for the
question *"is this cell-type annotation any good?"* — organised by **what evidence you have**,
because in real projects there is almost never a ground truth.

It consolidates metrics from several sources into one place, with attribution (see
[**the catalog**](docs/cell-annotation-quality-metrics.md) for the full reference, decision guide,
default metric batteries, and per-metric definitions).

## Install

```bash
pip install spatial-anno-metrics                 # core (no spatial)
pip install "spatial-anno-metrics[spatial]"      # + squidpy, for the Moran's I signal metric
```

## What's in it

### Reference-free (the default regime — no ground truth)

| Function | Metrics | From |
|---|---|---|
| `internal_validity(adata, label_key, embedding="X_pca")` | silhouette, 2label-silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch, avg similarity, `integrated` (geometric mean) | scTypeEval |
| `inter_sample_consistency(adata, label_key, sample_key)` | ISC — a type reproducible across biological replicates | scTypeEval |
| `marker_program_fidelity(adata, label_key, marker_sets)` | per-type AUC-ROC + Cohen's d of the marker-enrichment score | Zhu et al. 2026 |
| `panel_resolvability(reference, label_key, panel_genes)` | can a target panel resolve a reference's types **at the panel's own depth**? (depth-thinned CV classifier → per-type F1 + confuser) | — |
| `crisp_purity(adata, marker_sets)`, `mecr(adata, marker_sets)`, `pmp(adata, marker_sets, label_key)` | marker co-expression **contamination purity**: CRISP per-cell impurity + dataset purity, MECR (mutually-exclusive co-detection rate over disjoint lineage pairs), PMP (per-cell marker purity, panel-size-invariant) | SpatialScribe QC L3 |
| `annotation_quality(adata, label_key, marker_sets=, sample_key=)` | orchestrates the reference-free battery into one headline | — |

### Reference-based purity (needs a matched scRNA/snRNA reference)

| Function | Metric | From |
|---|---|---|
| `nmp(adata, reference, label_key)` | negative-marker proportion: fraction of a cell's counts on genes its assigned type should NOT express (the reference's per-type bottom decile) — the reference-based counterpart to `pmp` | SpatialScribe QC L3 |
| `ncp(adata, reference)` | non-coexpression preservation: fraction of reference-defined non-coexpressed gene pairs that stay non-coexpressed in the section (MECR-style, reference-derived pairs) | SpatialScribe QC L3 |

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

## Sources & attribution

The metrics are consolidated, with attribution, from:

- **scTypeEval** — Carmona lab, [github.com/carmonalab/scTypeEval](https://github.com/carmonalab/scTypeEval) (reference-free internal validity + inter-sample consistency; a Python re-implementation of its metric family). `[tool / self-benchmark]`
- **Zhu et al. 2026** — "Benchmarking cell type annotation in spatial transcriptomics," bioRxiv [10.64898/2026.06.16.732716](https://www.biorxiv.org/content/10.64898/2026.06.16.732716v1.full) (ECS, AvgBIO, hierarchical accuracy, marker-program fidelity; accuracy ≠ fidelity).
- **SpatialQM** — Center for Spatial OMICs, [github.com/Center-for-Spatial-OMICs/SpatialQM](https://github.com/Center-for-Spatial-OMICs/SpatialQM) (Moran's I gene-vs-control + signal-quality metrics). ⚠ NOT the identically-named Plummer/Spatial-Touchstone tool.
- **OpenProblems** — [task_spatial_decomposition](https://github.com/openproblems-bio/task_spatial_decomposition) (deconvolution R² / JSD).
- Classic metric definitions: ARI (Hubert & Arabie 1985), AMI/NMI (Vinh et al. 2010), silhouette (Rousseeuw 1987), ECS (Gates & Ahn 2019).

This package re-implements/consolidates metrics from those tools in Python; please **cite the
original sources** in publications, not just this package.

## Status

`v0.1.0`, experimental. Several thresholds/metrics in the catalog are flagged provisional —
version accordingly. Tests: `pip install -e ".[spatial,dev]" && pytest`.

## License

MIT.
