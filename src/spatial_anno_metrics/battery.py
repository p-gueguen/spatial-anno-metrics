"""The reference-free quality battery in one guarded call: internal validity, marker-program fidelity,
inter-sample consistency and the AQI headline."""
from __future__ import annotations

from ._validate import skip
from .aqi import annotation_quality_index
from .fidelity import marker_program_fidelity
from .internal import internal_validity
from .isc import inter_sample_consistency

__all__ = ["annotation_quality"]


# --------------------------------------------------------------------------- #
# Orchestrator (reference-free battery)
# --------------------------------------------------------------------------- #
def annotation_quality(adata, label_key: str = "cell_type",
                       marker_sets: dict[str, list[str]] | None = None,
                       sample_key: str | None = None, embedding: str | None = "X_pca",
                       n_neighbors: int = 30, subsample: int = 5000, seed: int = 0, *,
                       reference=None, ref_key: str = "cell_type", panel_genes=None,
                       median_depth: float | None = None, method_label_cols=None,
                       abstention_labels=None, platform=None, n_panel_genes=None,
                       normalization=None, marker_set_name=None, compute_aqi: bool = True) -> dict:
    """Run the reference-free quality battery and return a compact headline: internal validity,
    marker-program fidelity (if ``marker_sets``), inter-sample consistency (if ``sample_key``), and
    the single **AQI** headline (:func:`annotation_quality_index`) in ``out['aqi']``. Each sub-metric is
    guarded - a failure returns ``{'status': 'skipped', 'reason': ...}`` rather than raising, so this is
    safe to call from the QC funnel on any labeled section. ``n_neighbors`` is the embedding kNN size of
    :func:`internal_validity`; marker fidelity reads ``layers['counts']`` (or ``X`` when it holds counts).

    The AQI extra inputs are all optional and degrade gracefully: pass ``reference``+``panel_genes``
    (+``median_depth``) for the adequacy term, ``method_label_cols`` (>=3 obs columns) for the
    agreement term, and ``marker_sets`` for purity + fidelity. Without ``marker_sets`` the AQI has no
    term to score (its core is purity C and fidelity M; coherence only scales it), so ``out['aqi']['aqi']``
    is ``None``, not a coherence index. ``platform``/``n_panel_genes``/``normalization``/``marker_set_name``
    are stamped into ``provenance`` (every number must name its normalization + marker set)."""
    out: dict = {"status": "ok"}
    try:
        out["internal_validity"] = internal_validity(
            adata, label_key=label_key, embedding=embedding, n_neighbors=n_neighbors, subsample=subsample,
            seed=seed)
    except Exception as e:  # pragma: no cover - defensive
        out["internal_validity"] = skip(f"{type(e).__name__}: {e}")
    if marker_sets:
        try:
            out["marker_fidelity"] = marker_program_fidelity(adata, label_key, marker_sets)
        except Exception as e:
            out["marker_fidelity"] = skip(f"{type(e).__name__}: {e}")
    if sample_key and sample_key in adata.obs.columns:
        try:
            out["inter_sample_consistency"] = inter_sample_consistency(adata, label_key, sample_key)
        except Exception as e:  # pragma: no cover
            out["inter_sample_consistency"] = skip(f"{type(e).__name__}: {e}")
    if compute_aqi:
        try:
            out["aqi"] = annotation_quality_index(
                adata, label_key, marker_sets, out.get("internal_validity"),
                reference=reference, ref_key=ref_key, panel_genes=panel_genes,
                median_depth=median_depth, method_label_cols=method_label_cols,
                abstention_labels=abstention_labels, cohort_result=out.get("inter_sample_consistency"),
                provenance={"platform": platform, "n_panel_genes": n_panel_genes,
                            "median_depth": median_depth, "normalization": normalization,
                            "marker_set": marker_set_name, "reference": reference is not None})
        except Exception as e:  # pragma: no cover - defensive
            out["aqi"] = skip(f"{type(e).__name__}: {e}")
    return out
