"""Every hard-coded constant the metrics use, with where it comes from.

``THRESHOLDS[module][name] = {"value": ..., "source": "...", "provisional": bool}``. ``provisional`` is True for a
package choice that no study has validated for this use, and False for a value taken as is from the cited
method. The metric modules read these values once, at import, as their defaults: editing this dict does not
change a function's behaviour (pass the parameter instead, where there is one). Results that report the
thresholds they used (``inter_sample_consistency``'s ``thresholds``) echo these values.
"""
from __future__ import annotations

__all__ = ["THRESHOLDS"]

_PKG = "package default"


def _e(value, source, provisional=True):
    return {"value": value, "source": source, "provisional": provisional}


THRESHOLDS: dict[str, dict[str, dict]] = {
    "internal": {
        "n_neighbors": _e(30, f"{_PKG}: embedding kNN size for neighborhood purity"),
        "subsample": _e(2000, f"{_PKG}: cells used for the O(N^2) silhouettes"),
        "pca_n_components": _e(50, f"{_PKG}: PCA of X when obsm[embedding] is absent (also used by avg_bio)"),
    },
    "battery": {
        "n_neighbors": _e(30, f"{_PKG}: internal_validity kNN size inside annotation_quality"),
        "subsample": _e(5000, f"{_PKG}: internal_validity subsample inside annotation_quality"),
    },
    "isc": {
        "iscg": _e(0.74, "Garnica, Andreatta & Carmona 2026 (bioRxiv 10.64898/2026.06.10.731285), Fig. 5c: "
                         "median - 1 IQR over a 9-atlas meta-study; a guide, not a cut-off"),
        "iscl": _e(0.45, "Garnica, Andreatta & Carmona 2026, Fig. 5c (as iscg)"),
        "min_cells": _e(10, "Garnica et al. 2026 default: cells per (sample, type) pseudobulk", False),
        "min_samples": _e(5, "Garnica et al. 2026 default: pseudobulks a type needs to be scored", False),
        "pseudobulk_scale": _e(1e4, "Garnica et al. 2026 counts recipe: log1p of counts per 10k", False),
    },
    "fidelity": {
        "overlap_n_top": _e(25, f"{_PKG}: top DE genes per label in marker_gene_overlap"),
    },
    "external": {
        "ecs_alpha": _e(0.9, "Gates & Ahn 2019 element-centric similarity; clusim's default", False),
        "hierarchical_partial": _e(0.5, f"{_PKG}: credit for the right lineage, wrong subtype "
                                        "(docs/ROADMAP.md #9 proposes credit by LCA depth)"),
    },
    "panel": {
        "target_depth": _e(50.0, f"{_PKG}: counts per cell to thin the reference to"),
        "max_cells_per_type": _e(250, f"{_PKG}: reference cells sampled per type"),
        "cv_folds": _e(3, f"{_PKG}: cross-validation folds"),
        "min_shared_genes": _e(25, f"{_PKG}: fewer genes shared with the panel -> skip"),
        "min_cells_per_type": _e(10, f"{_PKG}: reference types with fewer cells are dropped"),
        "f1_resolvable": _e(0.5, f"{_PKG}: per-type F1 for 'resolvable' (docs/ROADMAP.md 8h: lenient)"),
        "f1_weak": _e(0.3, f"{_PKG}: per-type F1 below which a type is 'not_resolvable'"),
        "confused_weak": _e(0.25, f"{_PKG}: share of a type called one other type that makes it 'weak'"),
        "confused_not_resolvable": _e(0.30, f"{_PKG}: share of a type called one other type that makes it "
                                            "'not_resolvable'"),
        "logreg_C": _e(0.3, f"{_PKG}: inverse regularisation of the logistic-regression classifier"),
        "logreg_max_iter": _e(300, f"{_PKG}: logistic-regression iterations"),
        "library_size_scale": _e(1e4, "counts-per-10k library-size normalisation before log1p", False),
        "f1_adjusted_resolvable": _e(0.8, f"{_PKG}: chance-corrected per-type F1, (F1 - prevalence) / "
                                          "(1 - prevalence), for 'resolvable' (docs/ROADMAP.md 8h); 0.833 raw F1 "
                                          "with 6 balanced types, 0.9 with 2"),
        "f1_adjusted_weak": _e(0.5, f"{_PKG}: chance-corrected per-type F1 below which a type is "
                                    "'not_resolvable'"),
    },
    "conformal": {
        "alpha": _e(0.1, "miscoverage rate: 90 % target coverage, the usual default", False),
        "min_class_calibration": _e(10, f"{_PKG}: with fewer calibration cells a class uses the marginal "
                                        "threshold (class_conditional=True)"),
    },
    "purity": {
        "detection_threshold": _e(0, f"{_PKG}: a marker is detected when its count is > 0, i.e. one count "
                                     "(docs/ROADMAP.md 8e)"),
        "nmp_max_neg_frac": _e(0.05, "SegTraQ markers_from_reference rarity rule (github.com/LazDaria/SegTraQ): "
                                     "a negative is detected in <= 5 % of the type's reference cells", False),
        "nmp_min_other_frac": _e(0.25, "SegTraQ markers_from_reference rarity rule: ... and in >= 25 % of "
                                       "another type's", False),
        "ncp_coexpr_threshold": _e(0.1, f"{_PKG}: a gene pair co-detected in <= 10 % of the cells detecting "
                                        "either is non-coexpressed"),
        "ncp_max_genes": _e(400, f"{_PKG}: tractability cap on shared genes (top by reference mean)"),
        "ncp_max_pairs": _e(2000, f"{_PKG}: tractability cap on candidate pairs"),
    },
    "spatial": {
        "coherence_n_neighs": _e(15, f"{_PKG}; Shang & Zhou 2022 use 10 neighbours (docs/ROADMAP.md #2)"),
        "pas_threshold": _e(0.2, f"{_PKG}: a cell with < 20 % same-label neighbours is abnormal; Shang & Zhou "
                                 "2022 use <= 40 % (docs/ROADMAP.md #2)"),
        "sanity_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for neighborhood_sanity"),
        "contamination_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for neighbor_contamination"),
    },
    "confidence": {
        "stability_drop_frac": _e(0.2, f"{_PKG} (SpatialScribe QC Layer-5): share of transcripts dropped"),
        "stability_reps": _e(5, f"{_PKG}: thinning repetitions"),
    },
    "signal_qc": {
        "moran_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for Moran's I"),
        "snr_eps": _e(1e-9, f"{_PKG}: floor added to both means; it sets the ratio when controls are near zero "
                            "(docs/ROADMAP.md 8g)"),
        "single_gene_entropy": _e(0.1, f"{_PKG}: normalised detection entropy below which one gene dominates "
                                       "a cell"),
        "control_quantile": _e(0.95, f"{_PKG}: genes whose Moran's I exceeds this quantile of the controls'"),
    },
    "ontology": {
        "min_labels": _e(4, f"{_PKG}: labels with a CL id needed by ontology_consistency (6 label pairs)"),
        "min_cells": _e(10, f"{_PKG}: cells a label needs for its embedding centroid to be used"),
        "n_perm": _e(999, f"{_PKG}: label permutations for the Mantel-style p-value (floor 0.001)"),
    },
    # Recorded, not wired: aqi.py keeps these literals. tests/test_thresholds.py checks that its parameter
    # defaults agree; the body literals (w_coh_range, softmin_floor, default_target_depth,
    # coverage_limited_no_dict_frac) are recorded only.
    "aqi": {
        "min_type_n": _e(50, f"{_PKG} (aqi.py): a predicted type needs this many cells to be scored"),
        "p_softmin": _e(-4.0, f"{_PKG} (aqi.py): power of the soft-min over the purity and fidelity terms"),
        "w_coh_range": _e([0.85, 1.0], f"{_PKG} (aqi.py): range of the coherence multiplier"),
        "softmin_floor": _e(0.01, f"{_PKG} (aqi.py): terms are clipped to at least this before the soft-min"),
        "default_target_depth": _e(50.0, f"{_PKG} (aqi.py): panel depth when median_depth is not given"),
        "coverage_min_frac": _e(0.03, f"{_PKG} (aqi.py): share of confident cells that makes a lineage expected"),
        "coverage_min_cells": _e(50, f"{_PKG} (aqi.py): confident cells that make a lineage expected"),
        "coverage_margin": _e(0.25, f"{_PKG} (aqi.py): top-minus-second marker z-score of a confident call"),
        "reference_coverage_min_frac": _e(0.03, f"{_PKG} (aqi.py): reference share of an abundant lineage"),
        "coverage_limited_no_dict_frac": _e(0.5, f"{_PKG} (aqi.py): no_dict_frac above which coverage_limited "
                                                 "is flagged"),
    },
}


def value(module: str, name: str):
    """``THRESHOLDS[module][name]["value"]``: what the metric modules use as a default."""
    return THRESHOLDS[module][name]["value"]
