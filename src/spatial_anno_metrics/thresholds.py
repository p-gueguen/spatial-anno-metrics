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
        "clisi_n_neighbors": _e(90, "scib-metrics benchmark kNN size for cLISI (perplexity 90 / 3 = 30)", False),
        "isolated_f1_cuts_per_label": _e(2, f"{_PKG}: isolated-label F1 sweeps Ward cuts into 2 .. this x n_labels "
                                            "clusters (scib sweeps Leiden resolutions 0.1-2)"),
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
        "n_pcs": _e(30, "scTypeEval 1.0.0 README recipe (run_pca ndim = 30, run_dissimilarity reduction = TRUE): "
                        "ISCG on the first 30 centred PCs of the pseudobulks, the space the 0.74 guide comes from",
                    False),
        "rcm_marker_n": _e(500, "SingleR classic markers (Aran et al. 2019; scTypeEval's RCM classifier): "
                                "de.n = round(500 * (2/3)^log2(N labels)) per label pair", False),
    },
    "fidelity": {
        "overlap_n_top": _e(25, f"{_PKG}: top DE genes per label in marker_gene_overlap"),
    },
    "external": {
        "ecs_alpha": _e(0.9, "Gates & Ahn 2019 element-centric similarity; clusim's default", False),
        "hierarchical_partial": _e(0.5, f"{_PKG}: credit for a shared root under "
                                        "credit='lineage' (the pre-0.7.0 rule; the default credit='lca' "
                                        "scores by LCA depth)"),
        "spari_alpha": _e(0.8, "Yan, Feng & Luo 2025 (Biometrics, doi 10.1093/biomtc/ujaf127): spARI default",
                          False),
        "spari_max_cells": _e(5000, f"{_PKG}: cells sampled for spARI (the pair sums are O(n^2))"),
    },
    "panel": {
        "target_depth": _e(50.0, f"{_PKG}: counts per cell to thin the reference to"),
        "max_cells_per_type": _e(250, f"{_PKG}: reference cells sampled per type"),
        "cv_folds": _e(3, f"{_PKG}: cross-validation folds"),
        "min_shared_genes": _e(25, f"{_PKG}: fewer genes shared with the panel -> skip"),
        "min_cells_per_type": _e(10, f"{_PKG}: reference types with fewer cells are dropped"),
        "f1_resolvable": _e(0.5, f"{_PKG}: raw-F1 cut for 'resolvable' of the pre-0.7 tiers "
                                 "(null_adjusted=False; docs/ROADMAP.md 8h: lenient)"),
        "f1_weak": _e(0.3, f"{_PKG}: raw F1 below which a type is 'not_resolvable' in the pre-0.7 tiers "
                           "(null_adjusted=False)"),
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
        "ece_n_bins": _e(15, "Guo et al. 2017 (arXiv 1706.04599): 15 equal-width confidence bins", False),
        "clustered_quantiles": _e([0.5, 0.6, 0.7, 0.8, 0.9], "Ding et al. 2023 (arXiv 2306.09335) s. 3: score "
                                  "quantile levels (plus 1 - alpha) that embed a class for clustered conformal",
                                  False),
        "novelty_fdr": _e(0.1, "BH level for conformal novelty p-values (Bates et al. 2023 use 0.1)", False),
        "novelty_frac_train": _e(0.5, f"{_PKG}: reference share that fits the novelty scorer; the rest "
                                      "calibrates the p-values"),
        "novelty_n_components": _e(20, f"{_PKG}: PCA components for the reconstruction-error score"),
    },
    "consensus": {
        "min_share": _e(0.5, f"{_PKG}: vote share a label must exceed to be the consensus (strict majority; "
                             "popV takes the plurality, which ancestor votes would always hand to the root)"),
    },
    "purity": {
        "detection_threshold": _e(0, f"{_PKG}: a marker is detected when its count is > 0, i.e. one count "
                                     "(docs/ROADMAP.md 8e)"),
        "nmp_max_neg_frac": _e(0.05, "SegTraQ markers_from_reference rarity rule (github.com/LazDaria/SegTraQ): "
                                     "a negative is detected in <= 5 % of the type's reference cells", False),
        "nmp_min_other_frac": _e(0.25, f"{_PKG}: ... and in >= 25 % of another type's (SegTraQ's "
                                       "markers_from_reference uses min_pos_frac=0.1 plus a DE step; "
                                       "its 0.25 is MECR's detection floor)"),
        "ncp_coexpr_threshold": _e(0.1, f"{_PKG}: a gene pair co-detected in <= 10 % of the cells detecting "
                                        "either is non-coexpressed"),
        "ncp_max_genes": _e(400, f"{_PKG}: tractability cap on shared genes (top by reference mean)"),
        "ncp_max_pairs": _e(2000, f"{_PKG}: tractability cap on candidate pairs"),
        "null_draws": _e(20, f"{_PKG}: draws of the detection-matched null for crisp_purity / mecr "
                             "(docs/ROADMAP.md 8e)"),
        "nmp_salas_max_share": _e(0.005, "Salas et al. 2025 NMP as txsim negative_marker_purity_reads: "
                                         "(gene, type) is negative when the type holds < 0.5 % of the gene's "
                                         "type-balanced reference expression", False),
        "nmp_salas_min_cells": _e(10, "txsim negative_marker_purity_reads: types with >= 10 cells in both "
                                      "section and reference", False),
        "nmp_n_permutations": _e(20, "Salas et al. 2025 Methods: 0 = mean NMP over 20 random permutations of "
                                     "the cell-type labels", False),
        "coexpr_min_cells": _e(20, "txsim coexpression_similarity: a gene must be detected in >= 20 cells of "
                                   "each modality", False),
        "coexpr_target_sum": _e(1e4, f"{_PKG}: log1p(target_sum * counts / total) before the correlations "
                                     "(txsim reads a lognorm layer; OpenProblems normalises to the median "
                                     "total)"),
    },
    "spatial": {
        "coherence_n_neighs": _e(10, "Shang & Zhou 2022 (10.1038/s41467-022-34879-1) PAS: the 10 nearest "
                                     "neighbours (0.6.x used 15)", False),
        "pas_threshold": _e(0.4, "Shang & Zhou 2022 PAS: abnormal when the label differs from >= 6 of 10 "
                                 "neighbours, i.e. same-label fraction <= 0.4 (0.6.x: < 0.2)", False),
        "sanity_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for neighborhood_sanity"),
        "contamination_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for neighbor_contamination"),
        "contamination_n_perm": _e(199, f"{_PKG}: permutations for neighbor_contamination's p_value "
                                        "(smallest p = 1/200)"),
    },
    "confidence": {
        "stability_drop_frac": _e(0.2, f"{_PKG} (SpatialScribe QC Layer-5): share of transcripts dropped"),
        "stability_reps": _e(5, f"{_PKG}: thinning repetitions"),
    },
    "signal_qc": {
        "moran_n_neighs": _e(6, f"{_PKG}: spatial kNN graph for Moran's I"),
        "snr_eps": _e(1e-9, f"{_PKG}: opt-in (eps=1e-9 restores 0.6; default eps=None): floor added to both "
                            "means, it sets the ratio when controls are near zero (docs/ROADMAP.md 8g)"),
        "single_gene_entropy": _e(0.1, f"{_PKG}: normalised detection entropy below which one gene dominates "
                                       "a cell"),
        "control_quantile": _e(0.95, f"{_PKG}: genes whose Moran's I exceeds this quantile of the controls'"),
        "spatialqm_moran_n_neighs": _e(20, "SpatialQM getMorans (Plummer et al. 2025, commit 36e1a59): "
                                           "20 nearest neighbours weighted 1/d, row-standardised", False),
        "spatialqm_pseudocount": _e(0.1, "SpatialQM getMeanSignalRatio (commit 36e1a59): pseudocount added to "
                                         "the per-feature means before log10", False),
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
