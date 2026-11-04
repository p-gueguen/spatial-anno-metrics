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
| `internal_validity(adata, label_key, embedding="X_pca")` | silhouette, 2label-silhouette, neighborhood purity, orbital-medoid, Ward-PropMatch, avg similarity, **Davies-Bouldin + Calinski-Harabasz**, `integrated` (geometric mean of the defined parts); the O(N^2) parts on a `subsample` (runs at 1e5+ cells); NaN + `nan_reasons` when undefined; `embedding_label_trained=True` flags circularity | scTypeEval + sklearn |
| `isolated_label_scores(adata, label_key, batch_key=None)` | rare-label separation: isolated-label ASW (per-label mean silhouette) and F1 (best cluster of a Ward-cut sweep); without `batch_key`, macro averages over labels | scib / scib-metrics (Luecken et al. 2022) |
| `clisi(adata, label_key)` | cLISI: label purity of perplexity-weighted kNN neighbourhoods, scaled to [0, 1] (median-based, blind to small labels) | scib-metrics / Harmony LISI |
| `inter_sample_consistency(adata, label_key, sample_key, *, layer=None, reduction="pca", n_pcs=30, rcm_genes=None)` | **ISC on two axes** (Garnica et al. 2026), checked against scTypeEval 1.0.0: `iscg` (2-label silhouette on pseudobulk cosine in 30 centred PCs = cross-sample compactness; `reduction=None` for the raw space), `iscl` (silhouette on reciprocal-classification match = cross-sample replicability; `rcm_genes="markers"` or a gene list restricts the classifier), `integrated` = mean over labels of ISCG × ISCL; headline values are means over labels; per-label `regime` benchmarked on real cohorts (consistent / merge / remove / inconclusive; never split, which ISC cannot see) with the paper's fixed-threshold `quadrant` alongside, `dropped_types` below 5 samples × 10 cells | scTypeEval / Garnica et al. 2026 |
| `marker_program_fidelity(adata, label_key, marker_sets)` | per-type AUC-ROC + Cohen's d of the marker-enrichment score; reports `labels_without_markers` | Zhu et al. 2026 |
| `marker_gene_overlap(adata, label_key, marker_sets, n_top=25)` | fraction of a type's curated markers landing in its top-N data-derived DE genes (the fidelity DE-overlap view); reports `labels_without_markers` | Zhu et al. 2026 |
| `panel_resolvability(reference, label_key, panel_genes)` | can a target panel resolve a reference's types **at the panel's own depth**? (depth-thinned CV classifier → per-type F1 + confuser; tiers on the chance-corrected F1, `(F1 − prevalence) / (1 − prevalence)` ≥ 0.8 resolvable, < 0.5 not) | — |
| `ontology_consistency(adata, label_key, obo_path, label_to_cl, *, distance="path")` | do label centroids in the embedding follow **Cell Ontology** distances (`distance="path"` `is_a` shortest path or `distance="rwr"` Random Walk with Restart, `scGraph-OntoRWR`)? Spearman over label pairs + Mantel permutation p; names the label whose removal raises it most (a "macrophage" inside the epithelium) | Cell Ontology (`is_a`); Wang et al. 2025 (`scGraph-OntoRWR`) |
| `crisp_purity(adata, marker_sets)`, `mecr(adata, marker_sets)`, `pmp(adata, marker_sets, label_key)` | marker co-expression **contamination purity**: CRISP per-cell impurity + dataset purity (rule and formula from CRISP, Rose et al. 2026, detected at >= 1 count), MECR (mutually-exclusive co-detection rate over lineage pairs; a marker two sets share is not co-detection, nested sets are skipped), both against a detection-matched null (`vs_null`: 0 = exclusive, 1 = as random as independent markers); PMP (per-cell marker purity, panel-size-invariant). None of them tells spillover from genuine co-expression: use `nmp` | Rose et al. 2026; Hartman & Satija 2024; SpatialScribe QC L3 |
| `annotation_quality(adata, label_key, marker_sets=, sample_key=)` | orchestrates the reference-free battery into one headline | — |

### Reference-based purity (needs a matched scRNA/snRNA reference)

| Function | Metric | From |
|---|---|---|
| `nmp(adata, reference, label_key)` | negative marker purity, Salas et al. 2025 as txsim computes it (default `score="salas"`): reads on (gene, type) pairs the reference does not express, scaled so 0 = random labels and 1 = clean; tolerates co-expression the reference has. `score="moldia"`: older Moldia `Xenium_benchmarking` reads criterion (`x_sc < max_share / n_types`). `score="cells"`: txsim and Moldia cell-based NMP (`negative_marker_purity_cells`). `score="raw"`: per-cell fraction of counts on genes the reference detects in <= 5 % of the type's cells and >= 25 % of another type's (NaN when a type cannot be scored) | Salas et al. 2025; txsim; Moldia `Xenium_benchmarking`; raw `<= 5 %` cut from SegTraQ |
| `ncp(adata, reference)` | non-coexpression preservation: fraction of reference-defined non-coexpressed gene pairs that stay non-coexpressed in the section (`score="pairs"`, MECR-style), or Moldia `Xenium_benchmarking` conditional clipped-excess NCP (`score="moldia"`) | Salas et al. 2025; Moldia `Xenium_benchmarking` |
| `coexpression_similarity(adata, reference)` | gene-gene correlation agreement with the reference, `1 - \|rho_ref - rho_section\| / 2` averaged over pairs; `label_key=` per cell type | txsim |

### Spatial label quality (does the label make spatial sense? — needs coords)

| Function | Metric | From |
|---|---|---|
| `spatial_coherence(adata, label_key, n_neighs=10, pas_threshold=0.4)` | per-cell fraction of spatial neighbors sharing the label + dataset **PAS** (Percentage of Abnormal Spots: label differs from >= 6 of 10 neighbours) and **CHAOS** (mean distance to the nearest same-label cell), both as in Shang & Zhou 2022. Detects *random* label noise only: blind to a compact mislabelled patch, better when types are merged, worse for correctly labelled infiltrating cells. Down-weight only | Shang & Zhou 2022 (SpatialPCA); SpatialScribe QC L6 |
| `neighborhood_sanity(adata, label_key)` | flags types whose strongest neighborhood-enrichment partner is a *different* type (a mislabeling smell, e.g. two types forced together) via a squidpy permutation null | SpatialScribe QC L6 |
| `neighbor_contamination(adata, marker_sets, label_key)` | per-cell share of counts on markers of the **adjacent** cell types (not its own) + directed source → target frames: raw spillover, a **baseline** from same-type cells with no such neighbour, the **excess** over it, a permutation **p-value**, the SPLIT readout (`split_slope`, `split_corr` vs neighbour share), and CellSPA's negative-marker-vs-distance readout (`distance_corr`, `distance_curve`). Accepts SegTraQ's `{"positive", "negative"}` marker format; `graph="delaunay"` for SegTraQ's graph; reports labels without a marker set | SegTraQ `contamination_strength`; SPLIT (Bilous et al. 2026); CellSPA |

### Per-cell confidence

| Function | Metric | From |
|---|---|---|
| `per_cell_confidence(adata, marker_sets)` | marker-score **margin** (top1−top2) + softmax **entropy** — the margin, not the absolute score, is the confidence | SpatialScribe QC L5 |
| `label_stability(adata, annotate_fn, drop_frac=0.2)` | label-flip rate under transcript subsampling; annotator-agnostic (pass any labeller); the labeller gets a slim AnnData per rep (thinned counts, no other layers) | SpatialScribe QC L5 |
| `conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha=0.1)` | split-conformal **prediction sets** from any annotator's class probabilities: set size 1 = confident, > 1 = ambiguous, 0 = novel. `class_conditional=True` (Mondrian) keeps coverage on rare types, and `class_conditional="clustered"` pools classes with similar scores when many classes have few calibration cells. An optional `hierarchy` collapses an ambiguous subtype set to a lineage call. With `y_query`, it reports coverage, per-type coverage and CovGap. Coverage holds only when calibration and query cells are exchangeable (a dissociated reference vs a spatial section is not) | split conformal (LAC); clustered conformal, Ding et al. 2023 |
| `calibration_error(proba, classes, y_true)` | **ECE** and **class-wise ECE** of an annotator's probabilities: is the confidence a probability? 0 = calibrated | Guo et al. 2017 |
| `conformal_novelty(X_query, X_ref, method="pca")` | conformal **novelty p-values** in the reference's feature space (PCA reconstruction error or isolation forest), BH-adjusted to a false discovery rate. Valid only when known cells look like the reference: thin the reference to the section's per-cell depths first | Bates et al. 2023 |
| `consensus_entropy(adata, label_keys, hierarchy=None)` | per-cell **consensus** across several annotations with ontology-path voting (a label also votes for its ancestors), plus flat and ontology-aware **entropy**: a subtype split under an agreed parent does not count as disagreement | popV (Ergen et al. 2024) |

### Signal QC (annotation-independent — score the *data*, not the labels)

| Function | Metric | From |
|---|---|---|
| `moran_signal(adata, method=)` | **Moran's I gene-vs-control**: separation of gene and control medians, NaN-safe (undefined features counted). `method="squidpy"` (default): raw counts, squidpy kNN graph; `"spatialqm"`: SpatialQM `getMorans` (log-norm, 20-NN 1/d, negative-control probes) | SpatialQM (CSO) |
| `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, `tx_per_cell`, `run_signal_qc` | SNR vs controls (per cell where a control count exists + pooled; `method="spatialqm"` = getMeanSignalRatio / getMaxRatio), zero fraction, detection entropy (`method="spatialqm"` = getEntropy), transcripts per area and per cell (mean and median) | SpatialQM (CSO) |
| `transcript_qc(transcripts, negative_features)` | from a transcripts table: fraction of transcripts in cells (FTC) and the negative-codeword false discovery rate (FDR) | SpatialQM (CSO) |

### External (needs ground truth — benchmarks only)

| Function | Metrics | From |
|---|---|---|
| `external_scores(y_true=, y_pred=)` | balanced accuracy, macro/weighted/per-class F1, ARI/AMI/NMI, Cohen's kappa, ECS | classic + Gates & Ahn 2019 |
| `element_centric_similarity(a, b)` | exact vectorized ECS for hard partitions (resolution-agnostic); matches clusim to 1e-12 | Gates et al. 2019 |
| `spatial_ari(y_true=, y_pred=, coords=)` | **spRI / spARI**: Rand index / ARI where a disagreeing pair gets partial credit by spatial distance (wrong merge of near cells, wrong split of far cells cost less); `alpha=0` gives ARI | Yan, Feng & Luo, *Biometrics* 2025 |
| `hierarchical_accuracy(y_true=, y_pred=, hierarchy=)` | credit by LCA depth on a label tree of any depth (hierarchical F, Kiritchenko 2006): parent call 2/3, wrong sibling 1/2; micro hP/hR; accuracy of coarse (non-leaf) calls. `credit="lineage"` keeps the old flat partial credit | Zhu et al. 2026; Kiritchenko et al. 2006 |
| `avg_bio(adata, label_key, truth_key, embedding)` | AvgBIO = mean(ARI, NMI, scaled-ASW) — biological-structure conservation | Zhu et al. 2026 / scIB |
| `composition_accuracy(y_true=, y_pred=)` | hard-label global proportions: L1 / Pearson / JSD | — |

### Deconvolution (soft proportions per spot or cell)

| Function | Metrics | From |
|---|---|---|
| `deconvolution_metrics(true_prop=, pred_prop=)` | soft per-spot proportion matrices: **R² / JSD** (OpenProblems), RMSE, per-spot JSD, Pearson, **presence/absence AUPR** (`aupr` and Spotless `aupr_precrec`), **per-type SSIM** (Li et al. 2022) | OpenProblems `task_spatial_decomposition`; Spotless 2024; Li et al. 2022 |
| `deconvolution_rank_score(results)` | rank-aggregated score of several methods on one dataset: mean rank over RMSE / per-spot JSD / AUPR (Spotless) or any chosen metrics | Spotless 2024; Li et al. 2022 |
| `decomposition_purity(adata, weights_key)` | 1 - largest normalised deconvolution weight per cell (higher = more mixed): segmentation / reference-gap readout | Spatial Touchstone (Plummer 2025); SegTraQ `rctd_contamination` |

## Quick start

```python
import spatial_anno_metrics as sam

# Reference-free: you have labels + an embedding (the common case)
sam.internal_validity(adata, label_key="cell_type", embedding="X_pca")
# -> {'silhouette': .., 'neighborhood_purity': .., 'integrated': .., ...}

# Reference-free fidelity: you have marker sets but no truth
sam.marker_program_fidelity(adata, "cell_type", {"T cell": ["CD3D", "CD3E", "TRAC"], ...})

# Contamination purity: reference-free from curated markers (spillover / spatial doublets)
sam.crisp_purity(adata, marker_sets)   # dataset purity + obs['sam_crisp_impure']; sam.mecr(adata, marker_sets)

# Signal QC: is the section trustworthy at all? (needs coords + control probes)
sam.moran_signal(adata)          # median gene Moran's I >> control -> real signal

# Benchmark with ground truth
sam.external_scores(y_true=truth_labels, y_pred=pred_labels)     # macro_f1, ARI, kappa, ECS, ...
sam.deconvolution_metrics(true_prop=true_prop, pred_prop=pred_prop)    # R2, JSD (OpenProblems), AUPR, SSIM
```

## Choosing metrics

Don't compute everything — pick the battery for your situation. The
[decision guide](docs/cell-annotation-quality-metrics.md#1-decision-guide--which-metrics-to-compute) maps
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
- **SegTraQ** — Lazic et al. (EMBL / UZH), segmentation and transcript-assignment QC, [github.com/LazDaria/SegTraQ](https://github.com/LazDaria/SegTraQ) (MIT). Source of `neighbor_contamination` (its `contamination_strength`, re-implemented vectorized on an AnnData), the `<= 5 %` cut of `nmp(score="raw")` (`markers_from_reference`; its other-type floor is 0.1 plus a DE step, ours 0.25) and the shared-marker exclusion in MECR / CRISP. Its segmentation modules (border admixture, transcript positions, z-layer overlap) need transcripts + polygons: run SegTraQ itself on the SpatialData. `[tool]`
- **CRISP** — Rose et al., "CRISP enables comparisons of image-based spatial transcriptomic segmentation quality across ten organs," bioRxiv 2026, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947); code [github.com/santangelo-lab/CRISP](https://github.com/santangelo-lab/CRISP) (MIT). Source of the CRISP purity rule and formula (`crisp_purity` detects at > 0 counts, so its values are not comparable to published CRISP scores). `[preprint]`
- **txsim** — [github.com/theislab/txsim](https://github.com/theislab/txsim): reference implementation of Salas et al.'s NMP (`negative_marker_purity_reads`) and of co-expression similarity. `nmp` (default `score="salas"`) and `coexpression_similarity` are re-implemented from the definitions (txsim has no licence; it was run only to build the goldens). `[tool]`
- Classic metric definitions: ARI (Hubert & Arabie 1985), AMI/NMI (Vinh et al. 2010), silhouette (Rousseeuw 1987), ECS (Gates & Ahn 2019).

This package re-implements/consolidates metrics from those tools in Python; please **cite the
original sources** in publications, not just this package.

## Status

`v0.7.0`, experimental. Several thresholds/metrics in the catalog are flagged provisional —
version accordingly. Tests: `pip install -e ".[spatial,dev]" && pytest`.

## License

MIT.
