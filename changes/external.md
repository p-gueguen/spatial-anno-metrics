# changes/external.md - work package "external" (ROADMAP #9, #26, #27; ECS + deconvolution goldens)

## CHANGELOG

### Added

- `spatial_ari(*, y_true, y_pred, coords, alpha=0.8, max_cells=5000, seed=0)` -> `{status, spri, spari, n_cells,
  subsampled}`: spatially aware Rand index and adjusted Rand index (Yan, Feng & Luo, *Biometrics* 2025,
  [10.1093/biomtc/ujaf127](https://doi.org/10.1093/biomtc/ujaf127)). A disagreeing pair of cells gets partial
  credit by distance: a wrong merge of nearby cells, or a wrong split of distant cells, costs less. Implemented
  from the paper's definition; the authors' R package is GPL and was read, not copied. It reproduces the values
  printed in the package vignette for its example data to all 7 digits (spRI 0.9260064 / spARI 0.6716089 and
  0.9170200 / 0.6317260). `alpha=0` gives the ARI exactly. Above `max_cells` cells, a seeded random subset is
  scored, because the pair sums are O(n^2).
- `deconvolution_metrics` result keys `aupr`, `ssim` and `ssim_per_type`:
  - `aupr`: presence/absence average precision over all spot x type entries (Spotless, Sang-aram et al.,
    *eLife* 2024, [10.7554/elife.88431](https://doi.org/10.7554/elife.88431)); NaN when every entry is present,
    or none is.
  - `ssim_per_type` and their mean `ssim`: per-type SSIM after max-scaling (Li et al., *Nat Methods* 2022,
    [10.1038/s41592-022-01480-9](https://doi.org/10.1038/s41592-022-01480-9)). Ported from
    QuKunLab/SpatialBenchmarking (BSD-2-Clause, attribution in the code). NaN for an all-zero true column.
- `deconvolution_rank_score(results, metrics=("rmse", "jsd_per_spot", "aupr"))` -> `{status, metrics, mean_rank,
  ranks, order}`: ranks several deconvolution methods on one dataset per metric (1 = best, ties share the mean
  rank) and averages the ranks. The default metrics are the Spotless three. Li et al.'s accuracy score is
  `metrics=("pearson", "ssim", "rmse", "jsd")`.
- `hierarchical_accuracy` result keys:
  - `credit`;
  - `hierarchical_precision` and `hierarchical_recall`, micro-averaged hP and hR;
  - `coarse_call_frac`, the share of predictions that are an internal node of the tree;
  - `coarse_call_accuracy`, the share of those coarse calls that are the truth or one of its ancestors (NaN
    when there are no coarse calls).
- Value tests against the committed goldens:
  - `element_centric_similarity` vs clusim 0.4, 5 pairs x 2 alphas, at 1e-12;
  - `deconvolution_metrics` vs OpenProblems `task_spatial_decomposition`, `r2` on 3 pairs and `jsd` on 2, at
    1e-12.
- `THRESHOLDS["external"]`: `spari_alpha` (0.8, the authors' default) and `spari_max_cells` (5000, provisional).

### Changed (breaking)

- `hierarchical_accuracy(..., partial=0.5)`: the flat partial credit becomes credit by LCA depth. The new
  default is `credit="lca"`: the per-cell hierarchical F-measure (Kiritchenko et al. 2006),
  `hP = (lca+1)/(depth(pred)+1)` and `hR = (lca+1)/(depth(true)+1)`.
  - An honest coarse call (the parent) now scores 2/3, while a wrong sibling subtype still scores 1/2. Before,
    both scored 0.5.
  - On a two-level map, nothing else changes: an exact call scores 1, a wrong sibling 1/2 (= the old 0.5) and a
    wrong lineage 0.
  - `credit="lineage", partial=0.5` reproduces the old scores exactly.
  - `hierarchy` may now be any depth. `lineage_accuracy` means "shares a root"; for a two-level map that is the
    old same-parent fraction.

### Fixed

- `hierarchical_accuracy` scored an honest coarse call ("T" for a CD4 T cell) the same as a confident wrong
  sibling subtype (ROADMAP #9).

## README rows

Replacement rows for the external table (README.md lines 69-74), plus two new rows:

| `element_centric_similarity(a, b)` | exact vectorized ECS for hard partitions (resolution-agnostic); matches clusim to 1e-12 | Gates et al. 2019 |
| `spatial_ari(y_true=, y_pred=, coords=)` | **spRI / spARI**: Rand index / ARI where a disagreeing pair gets partial credit by spatial distance (wrong merge of near cells, wrong split of far cells cost less); `alpha=0` gives ARI | Yan, Feng & Luo, *Biometrics* 2025 |
| `hierarchical_accuracy(y_true=, y_pred=, hierarchy=)` | credit by LCA depth on a label tree of any depth (hierarchical F, Kiritchenko 2006): parent call 2/3, wrong sibling 1/2; micro hP/hR; accuracy of coarse (non-leaf) calls. `credit="lineage"` keeps the old flat partial credit | Zhu et al. 2026; Kiritchenko et al. 2006 |
| `deconvolution_metrics(true_prop=, pred_prop=)` | soft per-spot proportion matrices: **R² / JSD** (OpenProblems), RMSE, per-spot JSD, Pearson, **presence/absence AUPR** (Spotless), **per-type SSIM** (Li et al. 2022) | OpenProblems `task_spatial_decomposition`; Spotless 2024; Li et al. 2022 |
| `deconvolution_rank_score(results)` | rank-aggregated score of several methods on one dataset: mean rank over RMSE / per-spot JSD / AUPR (Spotless) or any chosen metrics | Spotless 2024; Li et al. 2022 |

## Catalog

docs/cell-annotation-quality-metrics.md:

- **§2.3 Partition agreement, row "spARI"**: replace it with
  `| **spARI** (spatially aware ARI) | Rand index / ARI where a pair the two labelings disagree on scores f(d) = α·exp(−d²) (together in prediction, apart in truth) or h(d) = α·(1−exp(−d²)) (apart in prediction, together in truth), coordinates min-max scaled; chance-adjusted with the authors' expectation (Yan, Feng & Luo, Biometrics 2025, doi 10.1093/biomtc/ujaf127; α = 0.8) | ≤1 | A spatially near mislabel costs less than a far one; α = 0 is ARI exactly. O(n²) pairs, so ``spatial_ari`` scores a seeded 5000-cell subset by default. Asymmetric in (truth, prediction). |`
- **§2.6, bullet "Hierarchical annotation accuracy"**: add, after "...instead of a flat miss.":
  "Credit is by distance on the label tree: the hierarchical F-measure of Kiritchenko et al. 2006 from the
  lowest-common-ancestor depth. An honest coarse call (the parent) earns 2/3 and a wrong sibling subtype 1/2, so
  the score prefers 'T cell' over a confident wrong 'CD8 T' for a CD4 T cell. A flat 0.5 for both cannot tell
  them apart. Report **accuracy on coarse (non-leaf) calls** next to it: a method that hedges to the parent is
  right more often at lower resolution."
- **§2.6, bullet "Cell-type composition / proportion accuracy", setting (ii)**: append "Presence/absence
  **AUPR** (Spotless, Sang-aram et al. eLife 2024: predicted proportions as scores for 'type present in spot')
  and per-type **SSIM** after max-scaling (Li et al. Nat Methods 2022) are also returned. To compare several
  deconvolvers on one dataset, ``deconvolution_rank_score`` averages per-metric ranks (Spotless: RMSE, JSD,
  AUPR; Li et al.'s accuracy score: PCC, SSIM, RMSE, JSD)."
- **§5 master table**:
  - The hierarchical row says "credit by LCA depth (parent call 2/3, wrong sibling 1/2) + coarse-call accuracy".
  - The deconvolution row adds "AUPR↑ / SSIM↑".
  - New row: `| **spARI** | E | L+truth+coords | ARI with distance-weighted credit for disagreeing pairs | ≤1, ↑ | external.py |`
- **§6 code block**:
  - `hierarchical_accuracy(*, y_true, y_pred, hierarchy, partial=0.5, credit="lca") -> dict  # credit by LCA depth`
  - `spatial_ari(*, y_true, y_pred, coords, alpha=0.8) -> dict  # spRI / spARI`
  - `deconvolution_rank_score(results, metrics=("rmse", "jsd_per_spot", "aupr")) -> dict`
- **§9 References**:
  - Replace "spARI - spatially-aware ARI" with "Yan, Feng & Luo, 'Spatially aware adjusted Rand index for
    evaluating spatial transcriptomics clustering', *Biometrics* 2025, doi 10.1093/biomtc/ujaf127 (preprint
    10.1101/2025.03.25.645156)".
  - Add: Sang-aram et al., Spotless, *eLife* 2024, doi 10.7554/elife.88431.
  - Add: Li et al., *Nat Methods* 2022, doi 10.1038/s41592-022-01480-9.
  - Add: Kiritchenko et al. 2006, hierarchical precision/recall (Canadian AI 2006, LNCS 4013).

## Roadmap

- **#9 Hierarchical partial credit: done.** `hierarchical_accuracy(credit="lca")` (default) gives credit by LCA
  depth through `_hierarchy.lca_depth`, and returns `coarse_call_frac` and `coarse_call_accuracy`.
  - Evidence: `src/spatial_anno_metrics/external.py` `hierarchical_accuracy`; tests
    `test_hierarchical_coarse_call_beats_a_wrong_sibling` (parent 2/3 > sibling 1/2),
    `test_hierarchical_credit_is_kiritchenko_hf_from_lca_depth` (3-level tree, hand values),
    `test_hierarchical_coarse_call_accuracy_counts_wrong_coarse_calls_and_is_nan_without_any`, and
    `test_hierarchical_lineage_credit_reproduces_the_old_default` (exact equality with the old code).
  - Sweep: 0/8208 values changed. No synthetic scenario emits an internal-node label. With 30% / 100% of T cells
    relabelled to the coarse "T", the old and `credit="lineage"` scores are 0.954167 / 0.846667 (diff 0), and the
    new default gives 0.969444 / 0.897778.
- **#26 Deconvolution: done.** `aupr`, `ssim`, `ssim_per_type` and `deconvolution_rank_score`.
  - Tests: `test_deconvolution_ssim_and_aupr` checks against the paper's SSIM formula and sklearn average
    precision; `test_deconvolution_undefined_aupr_and_ssim_are_nan` covers the NaN cases;
    `test_deconvolution_rank_score_orders_methods` covers ordering, ties, skips and unknown metrics.
  - Sweep check_deconvolution.csv, mean over 3 seeds at noise 0 / 0.1 / 0.25 / 0.5 / 1.0: aupr
    1.000 / 0.997 / 0.975 / 0.919 / 0.257, ssim 1.000 / 0.995 / 0.966 / 0.825 / 0.018.
  - Open: no reference-run golden for SSIM or AUPR. Executing the Li et al. code was not permitted in this
    session, and precrec (Spotless's AUPR) is not installed.
- **#27 spARI: done.** `spatial_ari`.
  - Evidence: `test_spatial_ari_matches_the_pairwise_definition` (an independent pair-by-pair loop, rel 1e-12),
    `test_spatial_ari_is_ari_at_alpha_zero_and_rewards_near_mistakes` (equals sklearn ARI at alpha=0; near
    mistakes score above far ones at equal ARI), and `test_spatial_ari_subsamples_large_inputs_and_validates`.
  - External check, not committed (the data are GPL): the package's `spARI_example_data`, read with base-R
    `load`, gives values identical to the vignette's printed digits for both clusterings.

## Notes for the integrator

- tests/test_thresholds.py `DEFAULTS` needs two rows (the contract asks for them; that file is not mine):
  `("external", "spatial_ari", "alpha"): ("external", "spari_alpha")` and
  `("external", "spatial_ari", "max_cells"): ("external", "spari_max_cells")`.
- `sam.__all__` grows from 33 to 35 names: `spatial_ari` and `deconvolution_rank_score`. `eval_metrics.py`
  does not re-export them; add them there if the shim should stay complete.
- `spatial_ari` takes arrays, like `external_scores`, so it does not belong in `test_validate._NAN_LABEL_CALLS`
  (which calls `(adata, "cell_type")`). It still routes labels through `check_labels`.
- The `hierarchical_partial` source string in thresholds.py still says "ROADMAP #9 proposes credit by LCA
  depth". Under the append-only rule I left it alone. Suggested text: "credit for a shared root under
  `credit='lineage'` (the pre-0.7.0 rule)".
