"""Tests for thresholds.py - one table of every hard-coded constant, with its source."""
from __future__ import annotations

import importlib
import inspect

import pytest

anndata = pytest.importorskip("anndata")

# (module, function, parameter) -> (THRESHOLDS module, name). aqi.py keeps its literals (only its PMP call may
# change), so its rows guard the table against drifting from the code.
DEFAULTS = {
    ("internal", "internal_validity", "n_neighbors"): ("internal", "n_neighbors"),
    ("internal", "internal_validity", "subsample"): ("internal", "subsample"),
    ("battery", "annotation_quality", "n_neighbors"): ("battery", "n_neighbors"),
    ("battery", "annotation_quality", "subsample"): ("battery", "subsample"),
    ("isc", "inter_sample_consistency", "min_cells"): ("isc", "min_cells"),
    ("isc", "inter_sample_consistency", "min_samples"): ("isc", "min_samples"),
    ("isc", "isc_regime", "thr_g"): ("isc", "iscg"),
    ("isc", "isc_regime", "thr_l"): ("isc", "iscl"),
    ("fidelity", "marker_gene_overlap", "n_top"): ("fidelity", "overlap_n_top"),
    ("external", "element_centric_similarity", "alpha"): ("external", "ecs_alpha"),
    ("external", "hierarchical_accuracy", "partial"): ("external", "hierarchical_partial"),
    ("panel", "panel_resolvability", "target_depth"): ("panel", "target_depth"),
    ("panel", "panel_resolvability", "max_cells_per_type"): ("panel", "max_cells_per_type"),
    ("panel", "panel_resolvability", "cv"): ("panel", "cv_folds"),
    ("conformal", "conformal_prediction_sets", "alpha"): ("conformal", "alpha"),
    ("purity", "nmp", "max_neg_frac"): ("purity", "nmp_max_neg_frac"),
    ("purity", "nmp", "min_other_frac"): ("purity", "nmp_min_other_frac"),
    ("purity", "ncp", "coexpr_threshold"): ("purity", "ncp_coexpr_threshold"),
    ("purity", "ncp", "max_genes"): ("purity", "ncp_max_genes"),
    ("purity", "ncp", "max_pairs"): ("purity", "ncp_max_pairs"),
    ("spatial", "spatial_coherence", "n_neighs"): ("spatial", "coherence_n_neighs"),
    ("spatial", "spatial_coherence", "pas_threshold"): ("spatial", "pas_threshold"),
    ("spatial", "neighborhood_sanity", "n_neighs"): ("spatial", "sanity_n_neighs"),
    ("spatial", "neighbor_contamination", "n_neighs"): ("spatial", "contamination_n_neighs"),
    ("confidence", "label_stability", "drop_frac"): ("confidence", "stability_drop_frac"),
    ("confidence", "label_stability", "reps"): ("confidence", "stability_reps"),
    ("signal_qc", "moran_signal", "n_neighs"): ("signal_qc", "moran_n_neighs"),
    ("aqi", "annotation_quality_index", "min_type_n"): ("aqi", "min_type_n"),
    ("aqi", "annotation_quality_index", "p_softmin"): ("aqi", "p_softmin"),
    ("aqi", "_lineage_coverage", "min_frac"): ("aqi", "coverage_min_frac"),
    ("aqi", "_lineage_coverage", "min_cells"): ("aqi", "coverage_min_cells"),
    ("aqi", "_lineage_coverage", "margin"): ("aqi", "coverage_margin"),
    ("aqi", "_reference_coverage", "tau"): ("aqi", "reference_coverage_min_frac"),
    ("internal", "clisi", "n_neighbors"): ("internal", "clisi_n_neighbors"),
    ("internal", "isolated_label_scores", "subsample"): ("internal", "subsample"),
    ("isc", "inter_sample_consistency", "n_pcs"): ("isc", "n_pcs"),
    ("purity", "crisp_purity", "n_null"): ("purity", "null_draws"),
    ("purity", "mecr", "n_null"): ("purity", "null_draws"),
    ("purity", "nmp", "n_permutations"): ("purity", "nmp_n_permutations"),
    ("purity", "coexpression_similarity", "min_cells"): ("purity", "coexpr_min_cells"),
    ("spatial", "neighbor_contamination", "n_perm"): ("spatial", "contamination_n_perm"),
    ("external", "spatial_ari", "alpha"): ("external", "spari_alpha"),
    ("external", "spatial_ari", "max_cells"): ("external", "spari_max_cells"),
    ("conformal", "calibration_error", "n_bins"): ("conformal", "ece_n_bins"),
    ("conformal", "conformal_novelty", "fdr"): ("conformal", "novelty_fdr"),
    ("conformal", "conformal_novelty", "n_components"): ("conformal", "novelty_n_components"),
    ("conformal", "conformal_novelty", "frac_train"): ("conformal", "novelty_frac_train"),
    ("consensus", "consensus_entropy", "min_share"): ("consensus", "min_share"),
    ("panel", "panel_resolvability", "f1_resolvable"): ("panel", "f1_adjusted_resolvable"),
    ("panel", "panel_resolvability", "f1_weak"): ("panel", "f1_adjusted_weak"),
    ("ontology", "ontology_consistency", "n_perm"): ("ontology", "n_perm"),
    ("ontology", "ontology_consistency", "min_cells"): ("ontology", "min_cells"),
}


def test_every_entry_has_a_value_a_source_and_a_provisional_flag():
    from spatial_anno_metrics.thresholds import THRESHOLDS

    for module, entries in THRESHOLDS.items():
        for name, e in entries.items():
            assert set(e) == {"value", "source", "provisional"}, (module, name)
            assert isinstance(e["provisional"], bool) and e["source"].strip(), (module, name)
    assert THRESHOLDS["isc"]["iscg"]["value"] == 0.74 and THRESHOLDS["isc"]["iscl"]["value"] == 0.45
    assert THRESHOLDS["isc"]["iscg"]["provisional"] is True


def test_function_defaults_are_the_table_values():
    from spatial_anno_metrics.thresholds import THRESHOLDS

    wrong = []
    for (mod, fn, param), (tm, tn) in DEFAULTS.items():
        f = getattr(importlib.import_module(f"spatial_anno_metrics.{mod}"), fn)
        default = inspect.signature(f).parameters[param].default
        if default != THRESHOLDS[tm][tn]["value"] or type(default) is not type(THRESHOLDS[tm][tn]["value"]):
            wrong.append((mod, fn, param, default, THRESHOLDS[tm][tn]["value"]))
    assert not wrong, wrong


def test_isc_thresholds_echo_the_table():
    from spatial_anno_metrics.isc import ISC_THRESHOLDS
    from spatial_anno_metrics.thresholds import THRESHOLDS

    assert ISC_THRESHOLDS == {k: THRESHOLDS["isc"][k]["value"] for k in ("iscg", "iscl")}
