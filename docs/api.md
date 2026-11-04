# API

Every name below is importable from the top-level package, e.g.
`from spatial_anno_metrics import internal_validity`. Sections follow the metric catalog order.

```{eval-rst}
.. currentmodule:: spatial_anno_metrics
```

## `internal`

Reference-free internal validity of a labeling on an embedding.

```{eval-rst}
.. autosummary::
   :toctree: generated

   internal_validity
   isolated_label_scores
   clisi
```

## `isc`

Inter-sample consistency: is each cell type reproducible across biological samples.

```{eval-rst}
.. autosummary::
   :toctree: generated

   inter_sample_consistency
   isc_regime
```

## `fidelity`

Marker-program fidelity (reference-free, needs marker sets).

```{eval-rst}
.. autosummary::
   :toctree: generated

   marker_program_fidelity
   marker_gene_overlap
```

## `external`

External harness (needs ground truth).

```{eval-rst}
.. autosummary::
   :toctree: generated

   avg_bio
   external_scores
   element_centric_similarity
   spatial_ari
   hierarchical_accuracy
   composition_accuracy
   deconvolution_metrics
   deconvolution_rank_score
```

## `panel`

Panel adequacy: can a target gene panel resolve a reference's cell types at the panel's own depth.

```{eval-rst}
.. autosummary::
   :toctree: generated

   panel_resolvability
```

## `battery`

The reference-free quality battery in one guarded call.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annotation_quality
```

## `aqi`

AQI: the single Annotation Quality Index headline, with its advisory lineage-coverage helpers.

```{eval-rst}
.. autosummary::
   :toctree: generated

   annotation_quality_index
```

## `conformal`

Conformal prediction sets, calibration error and conformal novelty p-values for cell-type annotation.

```{eval-rst}
.. autosummary::
   :toctree: generated

   conformal_prediction_sets
   calibration_error
   conformal_novelty
```

## `consensus`

Multi-annotator consensus with ontology-path voting, and its entropy.

```{eval-rst}
.. autosummary::
   :toctree: generated

   consensus_entropy
```

## `purity`

Marker co-expression / contamination purity.

```{eval-rst}
.. autosummary::
   :toctree: generated

   crisp_purity
   mecr
   pmp
   nmp
   ncp
   coexpression_similarity
```

## `spatial`

Spatial label-quality metrics: does the cell-type label make SPATIAL sense.

```{eval-rst}
.. autosummary::
   :toctree: generated

   spatial_coherence
   neighborhood_sanity
   neighbor_contamination
```

## `confidence`

Per-cell confidence / uncertainty (reference-free, given marker sets).

```{eval-rst}
.. autosummary::
   :toctree: generated

   per_cell_confidence
   label_stability
```

## `signal_qc`

Annotation-independent signal / QC metrics, squidpy-compatible.

```{eval-rst}
.. autosummary::
   :toctree: generated

   moran_signal
   signal_to_noise
   sparsity
   detection_entropy
   tx_per_area
   tx_per_cell
   transcript_qc
   run_signal_qc
```

## `ontology`

Ontology-aware annotation metrics: whether labels sit in the embedding the way the Cell Ontology relates them.

```{eval-rst}
.. autosummary::
   :toctree: generated

   ontology_consistency
```

## `decomposition`

Decomposition metrics: scores read from per-cell cell-type weights of a deconvolution (e.g. RCTD full mode).

```{eval-rst}
.. autosummary::
   :toctree: generated

   decomposition_purity
```

## `io`

Input adapters: get the AnnData the metrics read out of other containers.

```{eval-rst}
.. autosummary::
   :toctree: generated

   from_spatialdata
   segtraq_covariates
```

## `thresholds`

Every hard-coded constant the metrics use, with where it comes from.

```{eval-rst}
.. autosummary::
   :toctree: generated

   thresholds.THRESHOLDS
   thresholds.thresholds_for
```

