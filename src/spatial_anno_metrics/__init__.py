"""spatial_anno_metrics - cell-type annotation quality metrics (spatial + scRNA-seq).

A small, dependency-light battery for answering "is this cell-type annotation any good?",
organised by what evidence you have:

  * reference-free internal validity (scTypeEval family) + inter-sample consistency;
  * marker-program fidelity (reference-free, given marker sets): AUC-ROC / Cohen's d;
  * marker co-expression / contamination purity: reference-free CRISP / MECR / PMP + reference-based NMP / NCP;
  * spatial label quality: spatial coherence + PAS, neighborhood-enrichment sanity (labels vs geometry),
    neighbor contamination (spillover from adjacent cell types);
  * per-cell confidence (reference-free, marker sets): marker-score margin + entropy, label stability;
  * external harness (needs ground truth): F1 / ARI / kappa / ECS / hierarchical / composition / AvgBIO;
  * deconvolution proportion metrics (OpenProblems): R2 / JSD / RMSE;
  * annotation-independent signal QC (SpatialQM): Moran's I gene-vs-control, SNR, sparsity, entropy.

See ``docs/cell-annotation-quality-metrics.md`` for the catalog + the sources each metric is from
(scTypeEval, Zhu et al. 2026, SpatialQM/Center-for-Spatial-OMICs, OpenProblems, SpatialScribe QC L3).
"""
from . import (aqi, battery, confidence, conformal, decomposition, external, fidelity, internal, isc, ontology,
               panel, purity, signal_qc, spatial)
from .aqi import *  # noqa: F403
from .battery import *  # noqa: F403
from .confidence import *  # noqa: F403
from .conformal import *  # noqa: F403
from .decomposition import *  # noqa: F403
from .external import *  # noqa: F403
from .fidelity import *  # noqa: F403
from .internal import *  # noqa: F403
from .isc import *  # noqa: F403
from .ontology import *  # noqa: F403
from .panel import *  # noqa: F403
from .purity import *  # noqa: F403
from .signal_qc import *  # noqa: F403
from .spatial import *  # noqa: F403

__version__ = "0.6.0"

# The public API is each module's __all__, in catalog order. eval_metrics (deprecated shim) is not imported here.
__all__ = [name for mod in (internal, isc, fidelity, external, panel, battery, aqi, conformal, purity, spatial,
                            confidence, signal_qc, ontology, decomposition) for name in mod.__all__]


def __getattr__(name):
    """``sam.eval_metrics`` used to be bound by this file: load the deprecated shim lazily. The import system
    also lands here for ``from spatial_anno_metrics import eval_metrics``, so the warning is raised here, at the
    caller's line (raised inside the shim it would point at importlib, which default filters hide)."""
    if name != "eval_metrics":
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        shim = importlib.import_module(f"{__name__}.eval_metrics")
    warnings.warn(shim._MESSAGE, DeprecationWarning, stacklevel=2)
    return shim
