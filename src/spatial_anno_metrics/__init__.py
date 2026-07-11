"""spatial_anno_metrics - cell-type annotation quality metrics (spatial + scRNA-seq).

A small, dependency-light battery for answering "is this cell-type annotation any good?",
organised by what evidence you have:

  * reference-free internal validity (scTypeEval family) + inter-sample consistency;
  * marker-program fidelity (reference-free, given marker sets): AUC-ROC / Cohen's d;
  * marker co-expression / contamination purity: reference-free CRISP / MECR / PMP + reference-based NMP / NCP;
  * spatial label quality: spatial coherence + PAS, neighborhood-enrichment sanity (labels vs geometry);
  * per-cell confidence (reference-free, marker sets): marker-score margin + entropy, label stability;
  * external harness (needs ground truth): F1 / ARI / kappa / ECS / hierarchical / composition / AvgBIO;
  * deconvolution proportion metrics (OpenProblems): R2 / JSD / RMSE;
  * annotation-independent signal QC (SpatialQM): Moran's I gene-vs-control, SNR, sparsity, entropy.

See ``docs/cell-annotation-quality-metrics.md`` for the catalog + the sources each metric is from
(scTypeEval, Zhu et al. 2026, SpatialQM/Center-for-Spatial-OMICs, OpenProblems, SpatialScribe QC L3).
"""
from .confidence import (
    label_stability,
    per_cell_confidence,
)
from .eval_metrics import (
    annotation_quality,
    annotation_quality_index,
    avg_bio,
    composition_accuracy,
    conformal_prediction_sets,
    deconvolution_metrics,
    element_centric_similarity,
    external_scores,
    hierarchical_accuracy,
    inter_sample_consistency,
    internal_validity,
    marker_gene_overlap,
    marker_program_fidelity,
    panel_resolvability,
)
from .purity import (
    crisp_purity,
    mecr,
    ncp,
    nmp,
    pmp,
)
from .signal_qc import (
    detection_entropy,
    moran_signal,
    run_signal_qc,
    signal_to_noise,
    sparsity,
    tx_per_area,
)
from .spatial import (
    neighborhood_sanity,
    spatial_coherence,
)

__version__ = "0.4.1"

__all__ = [
    "internal_validity", "inter_sample_consistency", "marker_program_fidelity",
    "marker_gene_overlap", "avg_bio",
    "external_scores", "element_centric_similarity", "hierarchical_accuracy",
    "composition_accuracy", "deconvolution_metrics", "panel_resolvability",
    "annotation_quality", "annotation_quality_index", "conformal_prediction_sets",
    "crisp_purity", "mecr", "pmp", "nmp", "ncp",
    "spatial_coherence", "neighborhood_sanity",
    "per_cell_confidence", "label_stability",
    "moran_signal", "signal_to_noise", "sparsity", "detection_entropy",
    "tx_per_area", "run_signal_qc",
]
