"""Deprecated: ``spatial_anno_metrics.eval_metrics`` was split into modules in 0.7.0.

Import from the package (``import spatial_anno_metrics as sam``) or from the new modules: ``internal``,
``isc``, ``fidelity``, ``external``, ``panel``, ``conformal``, ``aqi`` and ``battery``. This shim re-exports
every name the old module defined (public and private) and warns once, on import.
"""
import warnings

from ._utils import _dense, _embedding, _scale01  # noqa: F401
from .aqi import _lineage_coverage, _reference_coverage  # noqa: F401
from .aqi import annotation_quality_index
from .battery import annotation_quality
from .conformal import conformal_prediction_sets
from .external import (avg_bio, composition_accuracy, deconvolution_metrics, element_centric_similarity,
                       external_scores, hierarchical_accuracy)
from .fidelity import marker_gene_overlap, marker_program_fidelity
from .internal import internal_validity
from .isc import _pseudobulk_profiles, _rcm_dissimilarity, _two_label_silhouette  # noqa: F401
from .isc import ISC_THRESHOLDS, inter_sample_consistency, isc_regime
from .panel import panel_resolvability

__all__ = [
    "ISC_THRESHOLDS", "annotation_quality", "annotation_quality_index", "avg_bio", "composition_accuracy",
    "conformal_prediction_sets", "deconvolution_metrics", "element_centric_similarity", "external_scores",
    "hierarchical_accuracy", "internal_validity", "inter_sample_consistency", "isc_regime",
    "marker_gene_overlap", "marker_program_fidelity", "panel_resolvability",
]

_MESSAGE = ("spatial_anno_metrics.eval_metrics is deprecated since 0.7.0 and will be removed: import from "
            "spatial_anno_metrics, or from its modules internal / isc / fidelity / external / panel / conformal / "
            "aqi / battery")
warnings.warn(_MESSAGE, DeprecationWarning, stacklevel=2)
