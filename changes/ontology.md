# changes/ontology.md - work package "ontology" (ROADMAP #25, 8h)

## CHANGELOG

### Added

- `ontology_consistency(adata, label_key, obo_path, label_to_cl, *, embedding="X_pca", n_perm=999, min_cells=10,
  seed=0)` -> `{status, rho, p_value, n_perm, n_labels, most_inconsistent, per_label: {label: {cl_id, loo_delta,
  nearest_embedding, nearest_ontology}}, excluded: {label: reason}}`. It asks whether the label centroids sit in
  the embedding the way the Cell Ontology relates the labels, and needs no ground truth.
  - `rho` is the Spearman correlation, over label pairs, between the distances of the embedding centroids and the
    `is_a` path lengths of the labels' CL terms. The path length is the shortest walk up to a common ancestor, over
    a DAG with several parents.
  - `p_value` comes from a one-sided Mantel test that permutes which label holds which term. With K labels there
    are only K! distinct permutations.
  - `loo_delta` is how much `rho` rises when that label is left out. The label with the largest positive
    `loo_delta` is `most_inconsistent`.
  - The OBO file is read with the standard library: `is_a` lines only, with obsolete terms and `[Typedef]` stanzas
    skipped.
  - It skips (with `excluded`) when fewer than 4 labels have a CL id and at least 10 cells, or when every ontology
    distance between them is the same.
- `panel_resolvability` `per_type` keys:
  - `f1_null`: the type's share of the scored cells. This is the base-rate F1, which is the label-permutation null.
  - `f1_adjusted`: `(f1 - f1_null) / (1 - f1_null)`.
- `panel_resolvability` keyword-only parameters `f1_resolvable`, `f1_weak` and `null_adjusted`.
- `THRESHOLDS["ontology"]` with `min_labels` (4), `min_cells` (10) and `n_perm` (999), all provisional.
- `THRESHOLDS["panel"]` gains `f1_adjusted_resolvable` (0.8) and `f1_adjusted_weak` (0.5).
- `tests/data/cl_tiny.obo`: a 14-term CL slice (4 KB) that the tests use.

### Changed (breaking)

- `panel_resolvability` tiers now use the chance-corrected F1 (ROADMAP 8h).
  - Old: `resolvable` at raw F1 >= 0.5 and `not_resolvable` below 0.3.
  - New: `resolvable` at `f1_adjusted` >= 0.8 and `not_resolvable` below 0.5. The confused-fraction cuts (0.25 / 0.30)
    are unchanged.
  - `f1_resolvable=0.5, f1_weak=0.3, null_adjusted=False` reproduces the old tiers exactly. This holds for all
    288 per-type tiers of the sweep `panel_check` grid.
  - `macro_f1`, `n_types`, the per-type `f1` values and every other key are unchanged.
  - On the sweep grid with no private markers at depth 50 (macro-F1 0.71), `frac_resolvable` falls from 1.00 to
    0.11.
  - `annotation_quality_index` reads `frac_resolvable` only for `components.completeness`, which is a side flag.
    `aqi` and `aqi_core` cannot move. The 12 `tests/test_aqi.py` fixture results are byte-identical, because the
    fixture's 3 types have F1 1.0.

## README rows

Replace the `panel_resolvability` row and add an ontology row to the "Reference-free" table:

| `panel_resolvability(reference, label_key, panel_genes)` | can a target panel resolve a reference's types **at the panel's own depth**? (depth-thinned CV classifier → per-type F1 + confuser; tiers on the chance-corrected F1, `(F1 − prevalence) / (1 − prevalence)` ≥ 0.8 resolvable, < 0.5 not) | — |
| `ontology_consistency(adata, label_key, obo_path, label_to_cl)` | do label centroids in the embedding follow **Cell Ontology** distances? Spearman over label pairs + Mantel permutation p; names the label whose removal raises it most (a "macrophage" inside the epithelium) | Cell Ontology (`is_a`) |

## Catalog

docs/cell-annotation-quality-metrics.md: add a new subsection **§3i** after §3h, add a §5 master-table row, and
update the existing panel-adequacy row.

New subsection after §3h (`<a name="3i-ontology"></a>`):

> ### 3i. Ontology label-relation consistency (needs a label -> CL id map, not ground truth)
>
> When the labels are Cell Ontology terms, their relations are a prior on the embedding: a T cell should sit
> nearer a B cell than a keratinocyte. `ontology.ontology_consistency` takes the labels' embedding centroids and
> the `is_a` path length between their CL terms (the shortest walk up to a common ancestor, read from a
> user-supplied OBO file with the standard library). It then reports the Spearman correlation over label pairs
> (`rho`) with a one-sided Mantel permutation p. Each label also gets a leave-one-out `loo_delta`, the rise in
> `rho` without it. The label with the largest rise is `most_inconsistent`, and `nearest_embedding` vs
> `nearest_ontology` shows why: a "macrophage" whose nearest centroid is epithelial while its nearest ontology
> term is the monocyte.
>
> Read it as a layout score, not a verdict on any label. Raw CL path lengths are coarse: on cl-basic 2026-06-08,
> macrophage-epithelial cell and macrophage-T cell are both 4 steps. Siblings that a good embedding separates
> also lower `rho`. With K labels the p cannot go below 1/K!.
>
> On the Xenium FFPE breast section, with marker-argmax labels for 10 types (12,031 confident cells, PCA 30):
>
> - clean labels: `rho` 0.28 (p 0.08);
> - planting the "macrophage" label on 1,500 luminal cells: `rho` 0.12, and `most_inconsistent` = macrophage;
> - a shuffled label -> CL map: `rho` -0.16.
>
> The order is right, but the effect is modest. ROADMAP #25 cites Genome Biol 2025,
> [10.1186/s13059-025-03781-6](https://doi.org/10.1186/s13059-025-03781-6), for the idea. This package did not
> check that paper's statistic against this one.

Replacement §5 row:

| Panel adequacy (coverage/confusable) | I | M+panel | can the panel resolve the type; tiers on chance-corrected F1 | resolvable | `panel.py` (`panel_resolvability`) |

New §5 row:

| Ontology label-relation consistency | I | L+emb+CL map | do label centroids follow Cell Ontology distances | high rho, low p | `ontology.py` (`ontology_consistency`) |

## Roadmap

- **#25 Ontology label-relation consistency: done.**
  - Code: `src/spatial_anno_metrics/ontology.py` (`ontology_consistency`, `_read_obo`, `_ontology_distances`).
  - Tests: `tests/test_ontology.py` (6 tests) and the fixture `tests/data/cl_tiny.obo`.
    - A planted macrophage centroid inside the epithelium is named `most_inconsistent`. `rho` falls from 0.69 to
      0.19.
    - A random layout gives p 0.89.
    - 4/4 mutants (single `is_a` parent, flipped `loo_delta` sign, lower-tail p, obsolete terms kept) are caught.
  - Real data: the Xenium breast check above. It is modest, and the planted label is named only on confident cells.
    With all 40k marker-argmax cells, the macrophage comes second, behind the plasma cells.
- **8h `panel_resolvability` tiers lenient: done.**
  - Code: `panel.py` uses chance-corrected F1 tiers. The thresholds are `panel.f1_adjusted_resolvable` and
    `panel.f1_adjusted_weak`.
  - Tests: `tests/test_panel.py::test_panel_resolvability_no_private_markers_is_not_resolvable`, `..._old_tiers_reproducible`
    and `..._two_types_at_chance_are_not_resolvable`. All three fail on wp/base.
  - Sweep `panel_check` (3 seeds), mean `frac_resolvable` old -> new:
    - no private markers: 1.00 -> 0.11 at depth 50 and 0.56 -> 0.00 at depth 20;
    - 25 % of private markers kept, depth 20: 0.94 -> 0.72;
    - all markers, depth 5: 1.00 -> 0.33;
    - at depth >= 20 with >= 50 % of markers, and at depth 200, it stays 1.00.
  - Permutation-null check: with permuted labels, the per-type F1 is 0.173 ± 0.026 (n = 90), against a prevalence
    null of 0.167. Every type was `not_resolvable`.

## Integrator notes (files this package may not edit)

- `tests/test_validate.py::_NAN_LABEL_CALLS`: add `"ontology_consistency": lambda sam, a: sam.ontology_consistency(a,
  "cell_type", <path to tests/data/cl_tiny.obo>, {"A": "CL:0000084", "B": "CL:0000236"})`. The NaN check runs
  first, so the map does not matter. `tests/test_ontology.py::test_excluded_labels_and_skips` already checks the
  NaN refusal.
- `tests/test_thresholds.py::DEFAULTS`: add these rows:
  - `("panel", "panel_resolvability", "f1_resolvable"): ("panel", "f1_adjusted_resolvable")`
  - `("panel", "panel_resolvability", "f1_weak"): ("panel", "f1_adjusted_weak")`
  - `("ontology", "ontology_consistency", "n_perm"): ("ontology", "n_perm")`
  - `("ontology", "ontology_consistency", "min_cells"): ("ontology", "min_cells")`
- `THRESHOLDS["panel"]["f1_resolvable"]` / `["f1_weak"]` (0.5 / 0.3) are no longer read by the code. Their source
  text could say "raw-F1 cuts of the pre-0.7 tiers (`null_adjusted=False`)". This package's rule was append-only,
  so it left them as they were.
- `sam.__all__` grows from 33 to 34 names (the contract text says 33).
- `benchmarks/synthetic/results/check_panel_resolvability.csv` gains the columns `frac_not_resolvable`,
  `T_CD8_f1_adjusted` and `T_CD8_tier` at the next full sweep.
