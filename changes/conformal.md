# changes/conformal.md - work package "conformal" (ROADMAP #20, #22, #24)

## Integrator: files this package could not edit

- `src/spatial_anno_metrics/__init__.py`: `consensus.py` is a new module. Add `consensus` to the
  `from . import (...)` line, add `from .consensus import *  # noqa: F403`, and put `consensus` in the `__all__`
  tuple after `conformal`. With `calibration_error` and `conformal_novelty` (already in `conformal.__all__`),
  `sam.__all__` grows from 33 to 36 names. Until then, `consensus_entropy` is reachable as
  `spatial_anno_metrics.consensus.consensus_entropy`.
- `tests/test_thresholds.py` `DEFAULTS` rows (each checked to pass on this branch):
  ```python
  ("conformal", "calibration_error", "n_bins"): ("conformal", "ece_n_bins"),
  ("conformal", "conformal_novelty", "fdr"): ("conformal", "novelty_fdr"),
  ("conformal", "conformal_novelty", "n_components"): ("conformal", "novelty_n_components"),
  ("conformal", "conformal_novelty", "frac_train"): ("conformal", "novelty_frac_train"),
  ("consensus", "consensus_entropy", "min_share"): ("consensus", "min_share"),
  ```
- `tests/test_validate.py` `_NAN_LABEL_CALLS` row, once `consensus` is exported (raises the expected
  "obs['cell_type'] has 2 missing" on this branch):
  `"consensus_entropy": lambda sam, a: sam.consensus_entropy(a, ["cell_type", "cell_type"]),`.
  `calibration_error` takes label arrays, like `conformal_prediction_sets`, so it is not in that table; its NaN
  refusal is tested in `tests/test_conformal.py`.
- `docs/cell-annotation-quality-metrics.md` section 5 master table: the conformal row still says
  `eval_metrics.py`; the code is in `conformal.py`.

## CHANGELOG

### Added

- `calibration_error(proba, classes, y_true, *, n_bins=15)` -> `{status, ece, classwise_ece, per_class_ece,
  accuracy, mean_confidence, n}`: top-label expected calibration error (Guo et al. 2017,
  [arXiv 1706.04599](https://arxiv.org/abs/1706.04599), 15 equal-width bins) and class-wise ECE (one-vs-rest per
  class, averaged). Check it before reading an annotator's confidence as a probability. Class-wise ECE matches
  values built from sklearn's `calibration_curve` to 1e-12.
- `conformal_prediction_sets(..., class_conditional="clustered", *, n_clusters=None, frac_clustering=None,
  seed=0)`: clustered Mondrian conformal (Ding et al. 2023, [arXiv 2306.09335](https://arxiv.org/abs/2306.09335)),
  written from the paper. It pools classes with similar score distributions, so rare classes get a class-like
  threshold without needing many calibration cells each. The result adds `clusters` (class -> cluster; -1 is the
  null cluster, which uses the marginal threshold). On 40 classes with 30 calibration cells each (half easy, half
  hard), CovGap is 0.024, against 0.101 marginal and 0.059 per-class, and the easy and hard classes never share
  a cluster (test). With 6 classes, Ding's
  heuristic gives 1 cluster when the rarest class is small, which is the marginal threshold; use
  `class_conditional=True` there.
- `conformal_prediction_sets` result key `covgap` (when `y_query` is given): the mean over the query's true
  types of `|coverage_y - (1 - alpha)|`, Ding et al.'s CovGap as a fraction (they report x100).
- `conformal_novelty(X_query, X_ref, *, method="pca", fdr=0.1, n_components=20, frac_train=0.5, seed=0)` ->
  `{status, pvalues, qvalues, novel, score, summary}`: split-conformal novelty p-values in the reference's
  feature space, BH-adjusted (Bates et al. 2023, *Ann Stat*,
  [10.1214/22-AOS2244](https://doi.org/10.1214/22-AOS2244)). The score is the PCA reconstruction error or
  `IsolationForest` (`method="iforest"`). Half of the reference fits the score and the other half calibrates it.
  The p-values and BH are tested against a hand-computed case and scipy's `false_discovery_control`.
  Exchangeability matters: in the synthetic check, an un-thinned 400-count reference flags 85-100 % of known
  150-count cells. Thin the reference to the section's per-cell depths first.
- New module `consensus.py`: `consensus_entropy(adata, label_keys, *, hierarchy=None, min_share=0.5)` ->
  `{status, consensus, consensus_share, entropy, ontology_entropy, summary}`. It gives per-cell agreement
  between annotators with ontology-path voting (popV, Ergen et al. 2024, *Nat Genet*,
  [10.1038/s41588-024-01993-3](https://doi.org/10.1038/s41588-024-01993-3)):
  - every label votes for itself and for its ancestors in the `{child: parent}` map (`_hierarchy.ancestors`);
  - the consensus is the deepest label with more than `min_share` of the votes, or None when no label has that
    share;
  - `entropy` is the flat entropy of the votes;
  - `ontology_entropy` is the entropy of the votes projected to the consensus depth, so a CD4 / CD8 / T split
    scores 0 where the flat entropy scores 1.

  Both entropies are divided by log(number of annotators). The function writes nothing to `adata`.
  - popV propagates only OnClass's votes and takes the plurality. Here every annotator's votes propagate, so the
    plurality would always be the root; that is why the rule is a strict majority.
- `THRESHOLDS["conformal"]`: `ece_n_bins` (15), `clustered_quantiles` ([0.5..0.9], Ding et al.), `novelty_fdr`
  (0.1), `novelty_frac_train` (0.5, provisional), `novelty_n_components` (20, provisional).
  `THRESHOLDS["consensus"]`: `min_share` (0.5, provisional).
- Synthetic benchmark `check_conformal.csv`:
  - new rows: `mondrian="clustered"`, and `held_out="B"`, where B is left out of the reference so the section's
    B cells are a truly novel type;
  - new columns: `held_out`, `covgap`, `pct_empty_known`, `pct_empty_novel`, `ece`, `classwise_ece`, `nov_auc`,
    `nov_flag_known`, `nov_flag_novel`, `nov_fdp`, `nov_flag_known_unmatched`;
  - the 48 old rows and 8 old columns are identical to before.

### Changed

- `conformal_prediction_sets`: `class_conditional` also accepts `"clustered"`, and any other value now raises
  `ValueError`. False and True give identical results; the sweep's 48 old conformal rows match wp/base exactly.

## README rows

Replacement for the `conformal_prediction_sets` row of the "Per-cell confidence" table, plus three new rows:

| `conformal_prediction_sets(proba_query, classes, proba_cal, y_cal, alpha=0.1)` | split-conformal **prediction sets** from any annotator's class probabilities: set size 1 = confident, > 1 = ambiguous, 0 = novel. `class_conditional=True` (Mondrian) keeps coverage on rare types, and `class_conditional="clustered"` pools classes with similar scores when many classes have few calibration cells. An optional `hierarchy` collapses an ambiguous subtype set to a lineage call. With `y_query`, it reports coverage, per-type coverage and CovGap. Coverage holds only when calibration and query cells are exchangeable (a dissociated reference vs a spatial section is not) | split conformal (LAC); clustered conformal, Ding et al. 2023 |
| `calibration_error(proba, classes, y_true)` | **ECE** and **class-wise ECE** of an annotator's probabilities: is the confidence a probability? 0 = calibrated | Guo et al. 2017 |
| `conformal_novelty(X_query, X_ref, method="pca")` | conformal **novelty p-values** in the reference's feature space (PCA reconstruction error or isolation forest), BH-adjusted to a false discovery rate. Valid only when known cells look like the reference: thin the reference to the section's per-cell depths first | Bates et al. 2023 |
| `consensus_entropy(adata, label_keys, hierarchy=None)` | per-cell **consensus** across several annotations with ontology-path voting (a label also votes for its ancestors), plus flat and ontology-aware **entropy**: a subtype split under an agreed parent does not count as disagreement | popV (Ergen et al. 2024) |

## Catalog

For `docs/cell-annotation-quality-metrics.md` section **3d. Per-cell confidence / uncertainty**, replace the
"Conformal label-set size" bullet with these bullets:

- **Conformal label-set size**: distribution-free. Set size 1 = confident, >1 = ambiguous, 0 = novel; it is the
  cleanest formalization of "abstain". Per-class (Mondrian) thresholds keep coverage on rare types. With many
  classes and few calibration cells each, clustered Mondrian (Ding et al. 2023) pools classes whose score
  distributions look alike. CovGap (mean |per-class coverage - target|) is the number to compare them on.
  `conformal_prediction_sets` in `conformal.py`.
- **Calibration error (ECE / class-wise ECE)** (Guo et al. 2017): before reading a softmax or an RCTD weight as
  a probability, check that cells given 0.8 are right about 80 % of the time. `calibration_error` in
  `conformal.py`.
- **Conformal novelty p-values** (Bates et al. 2023): an empty conformal set is a weak novelty signal, because
  depth loss empties sets for known types too. Synthetic check, B left out of the reference, 3 seeds:
  - At 1/4 section depth, 18 % of known cells get an empty set, against 100 % of the novel B cells.
  - The novelty p-values (PCA reconstruction error, BH 10 %, reference thinned to the section's per-cell
    depths) separate B at AUC 0.99. They flag 1.2 % of known cells and 85 % of B cells, at a false-discovery
    proportion of 0.10.
  - Without the depth matching, 87-100 % of known cells are flagged at every depth.

  `conformal_novelty` in `conformal.py`.
- **Consensus entropy with ontology-path voting** (popV, Ergen et al. 2024): when several annotators label the
  same cells, each label also votes for its ancestors, and the consensus is the deepest label with a majority.
  The entropy of the votes, projected to that depth, measures disagreement without penalising subtype splits.
  `consensus_entropy` in `consensus.py`.

Section 5 master table, new rows (and fix the conformal row's code location to `conformal.py`):

| Calibration error (ECE / class-wise ECE) | E | scores+truth | is confidence a probability? | low | `conformal.py` (`calibration_error`) |
| Conformal novelty p-value | I | R+features | novel / out-of-reference cells, FDR-controlled | high p = known | `conformal.py` (`conformal_novelty`) |
| Consensus entropy (ontology voting) | I | several L (+hierarchy) | agreement across annotators | low | `consensus.py` (`consensus_entropy`) |

## Roadmap

- **#20 conformal novelty p-values: done.**
  - `conformal_novelty` in `conformal.py` (PCA reconstruction error or isolation forest; BH).
  - Tests: `test_conformal_pvalues_and_bh_hand_and_scipy` and
    `test_conformal_novelty_flags_novel_and_controls_fdr[pca|iforest]`.
  - Synthetic answer to "novel type or depth shift?" (`check_conformal.csv`, `held_out="B"`, 3 seeds, PCA, BH
    10 %, reference thinned to the section's per-cell depths):

    | section depth | AUC | known cells flagged | B cells flagged | FDP |
    |---|---|---|---|---|
    | 1 | 0.999 | 0.6 % | 99 % | 0.05 |
    | 1/2 | 0.998 | 1.1 % | 99 % | 0.08 |
    | 1/4 | 0.990 | 1.2 % | 85 % | 0.10 |

    For comparison, empty conformal sets (reference calibration, marginal) at the same depths: 0.8 %, 3.2 % and
    17.7 % of known cells, and 100 % of B cells at every depth.
  - Limit: the depth-matching step is the user's. With the 400-count reference as it is, 85-100 % of known
    cells are flagged at every depth.
  - Isolation forest is weaker. On the Gaussian test it needs a larger shift. On the depth-matched sim
    reference, 3 seeds, it flagged no cell at all (AUC 0.84, 0.63 and 0.42 at depth 1, 1/2 and 1/4), so PCA is
    the default.
- **#22 consensus entropy with ontology-path voting: done.** `consensus_entropy` in `consensus.py`.
  - Tests in `tests/test_consensus.py`, with hand-computed values for 6 vote patterns, including Treg -> T_CD4
    propagation and a 3-root no-consensus case.
  - Three mutants each fail a test: no ancestor votes, no depth projection, `>=` majority.
  - Not run on real multi-annotator data.
- **#24 calibration: done.**
  - `calibration_error` (ECE and class-wise ECE): class-wise ECE matches sklearn's `calibration_curve` at
    1e-12, and a calibrated vs sharpened annotator gives ECE < 0.02 vs > 0.1.
  - `covgap` key.
  - `class_conditional="clustered"` (Ding et al. 2023), with
    `test_clustered_mondrian_lowers_covgap_and_finds_the_difficulty_groups` and the null-cluster test.
  - In the 6-type sweep, clustered equals Mondrian when every class has plenty of calibration cells. It falls
    back to about the marginal threshold when the rarest class has about 14 cells: T_CD8 at 2 %, section
    calibration, T_CD8 coverage 0.77 marginal, 0.93 Mondrian, 0.77 clustered.
