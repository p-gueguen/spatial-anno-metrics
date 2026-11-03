# changes/signal.md - work package "signal" (ROADMAP 8g, #8 Moran null, #21; SpatialQM golden repair)

## CHANGELOG

### Added

- `moran_signal(..., method="spatialqm")`: SpatialQM `getMorans`. Per group (genes; negative-control probes),
  it keeps the cells with a count in the group, applies `scater::logNormCounts`, and builds a graph of 20 nearest
  neighbours weighted 1/d, row-standardised. It needs no squidpy and writes nothing to `adata`. On the Xenium
  golden it reproduces all 341 per-feature values to 5.3e-15, with NaN at the same 7 features. The full 167,780-cell
  section takes 3.2 s.
  - `genes=` chooses what is scored; cells and library sizes still come from the whole group, so a feature's I
    does not depend on which others are requested (ERBB2 stays 0.3679 = `getMorans` in any subset).
  - A group with <= 20 cells holding a count gets NaN for all its features, counted in `n_*_excluded`; the other
    group is still scored. spdep's `knearneigh(k = 20)` stops there, so SpatialQM gives no value either.
- `moran_signal` result keys `n_genes_excluded`, `n_controls_excluded` (features with undefined I) and
  `per_feature` (pandas Series of I), plus `method`. `method="spatialqm"` also returns `n_cells_genes` and
  `n_cells_controls`.
- `signal_to_noise(..., method="spatialqm")`: `mean_signal_ratio_log10` and `max_signal_ratio_log10` from
  per-feature means. It matches `getMaxRatio` exactly. `getMeanSignalRatio` differs by 8.2e-5, because we do not
  copy R's recycling of the noise means; SpatialQM's own value moves by 2.9e-4 with the order of the noise rows.
- `signal_to_noise` result keys `pooled_mean_signal_ratio`, `n_cells_with_controls`, `frac_cells_with_controls`;
  keyword-only `eps=`.
- `detection_entropy(..., method="spatialqm")`: `getEntropy`, the Shannon entropy in bits of the count-value
  histogram over every matrix entry (`entropy_bits`).
- `tx_per_area` result key `mean_tx_per_area` (= `getTxPerArea`), next to the median.
- `tx_per_cell(adata, *, layer, key_added="sam_tx_per_cell")`: gene transcripts per cell, mean (= `getTxPerCell`)
  and median.
- `transcript_qc(transcripts, negative_features, *, feature_key="feature_name", cell_key="cell_id")`, from a
  transcripts DataFrame (no spatialdata dependency):
  - `frac_tx_in_cells` (FTC, = `getCellTxFraction`). It accepts Xenium <= 1.2 ids (`-1`) and >= 1.3 ids
    (`"UNASSIGNED"`); SpatialQM returns 1 for the old ids.
  - `fdr` = negative transcripts / all transcripts x genes detected / negatives detected. SpatialQM's
    `getGlobalFDR` is `fdr / 100`. `fdr` is None when no negative feature is detected.
  - Bytes feature names and cell ids (Xenium writes `feature_name` as parquet `binary`) are decoded as UTF-8.
    On row group 0 of the breast `transcripts.parquet` (1e6 rows): fdr 0.006015, 4,324 negative transcripts,
    313 genes detected.
- `moran_signal`, `transcript_qc` and every `method="spatialqm"` call return their result and write nothing
  to `uns['sam']` (documented in the module docstring).
- `decomposition_purity(adata, weights_key, *, key_added="sam_decomposition_purity")`: 1 - the largest
  normalised deconvolution weight per cell, as in Spatial Touchstone. It is ported from SegTraQ
  `sp.rctd_contamination` (LazDaria/SegTraQ `05e90b6`, MIT, attributed). Weights come from `obsm`. Rows with a
  non-finite or negative weight, or a non-positive total, give NaN.
- thresholds `signal_qc.spatialqm_moran_n_neighs` (20) and `signal_qc.spatialqm_pseudocount` (0.1), both from
  SpatialQM (`provisional=False`).

### Changed

- Controls are flagged by feature type, never by name. `var['control']` wins when present. Otherwise every
  feature whose `var['feature_type']` (or `feature_types`) is not `"Gene Expression"` is a control.
  `moran_signal(method="spatialqm")` uses only `"Negative Control Probe"`. Before, without `var['control']`
  there were no controls.
- **`signal_to_noise` (default `method="cell"`) no longer returns the eps floor (breaking, values):**
  - Per-cell ratios are defined only for cells with a control count; the others are NaN (old:
    log2(x / 1e-9)). The medians are taken over the defined cells.
  - It skips (and writes nothing) when no cell has a control count.
  - `frac_cells_signal_below_control` is now taken over the cells with a control count (old: all cells).
  - On the Xenium crop: median 27.956 -> 6.491 (224 of 2000 cells have a control count; pooled 9.511).
  - `eps=1e-9` reproduces the 0.6 values exactly.
- `run_signal_qc` adds a `tx_per_cell` entry.

### Fixed

- `moran_signal` summaries were NaN whenever a feature had no count. On the crop, 139 features (1 gene,
  138 controls) made the medians and `separation` NaN and `frac_genes_above_control_p95` 0.0; now they are
  0.10637 / -0.00058 / 0.10695 / 0.929.
- Docstrings no longer claim SpatialQM ports for different quantities. The default recipes say how they differ;
  `method="spatialqm"` computes SpatialQM's quantities.
- SpatialQM golden:
  - `expected.json` records the sha256 of the three fixtures it was computed from (`meta.input.sha256`), and
    the test pins `expected.json`'s own sha256. All 18 planted corruptions fail.
  - `make_inputs.py outs` gzips the committed `cells.csv` byte for byte, so the regenerated `expected.json` is
    identical from venv310 and venv312.

## README rows

Replacement rows for README.md lines 62-63, plus one new row for the deconvolution table (after line 74):

| `moran_signal(adata, method=)` | **Moran's I gene-vs-control**: separation of gene and control medians, NaN-safe (undefined features counted). `method="squidpy"` (default): raw counts, squidpy kNN graph; `"spatialqm"`: SpatialQM `getMorans` (log-norm, 20-NN 1/d, negative-control probes) | SpatialQM (CSO) |
| `signal_to_noise`, `sparsity`, `detection_entropy`, `tx_per_area`, `tx_per_cell`, `run_signal_qc` | SNR vs controls (per cell where a control count exists + pooled; `method="spatialqm"` = getMeanSignalRatio / getMaxRatio), zero fraction, detection entropy (`method="spatialqm"` = getEntropy), transcripts per area and per cell (mean and median) | SpatialQM (CSO) |
| `transcript_qc(transcripts, negative_features)` | from a transcripts table: fraction of transcripts in cells (FTC) and the negative-codeword false discovery rate (FDR) | SpatialQM (CSO) |
| `decomposition_purity(adata, weights_key)` | 1 - largest normalised deconvolution weight per cell (higher = more mixed): segmentation / reference-gap readout | Spatial Touchstone (Plummer 2025); SegTraQ `rctd_contamination` |

## Catalog

Section **3h. Signal-quality QC (annotation-independent; SpatialQM)** in docs/cell-annotation-quality-metrics.md:
replace the bullet list with:

- **Which features are controls.** They are flagged by feature type, never by name: in the 10x Xenium breast
  panel, 8 of the 28 negative-control probes are named `antisense_<GENE>`, and SpatialQM's name regexes count them
  as genes (`getMaxRatio` 4.333 vs 4.001, FDR 0.593 % vs 0.664 %). Rule: `var['control']` if present, else
  `var['feature_type'] != "Gene Expression"`.
- **Moran's I gene-vs-control** (`getMorans`). Real genes are spatially autocorrelated and controls sit at
  I ~ 0; the separation of the medians is the score. Features with undefined I (no count) are left out and
  counted. With shuffled coordinates the separation collapses to ~0 (test). Two methods:
  - `"squidpy"` (default): raw counts, binary 6-NN graph, all controls.
  - `"spatialqm"`: SpatialQM's recipe, reproduced to 5e-15. Its control median rests on the few cells with a
    probe count (150 of 2000 on the crop), so it sits near E[I] = -1/(n-1).
- **Signal-to-noise** (`getMeanSignalRatio` / `getMaxRatio`):
  - Default: per cell, log2(mean gene / mean control), defined only where the cell has a control count. 81.9 %
    of the full Xenium section's cells and 88.8 % of the crop's have none. The result therefore reports
    `n_cells_with_controls` and the pooled ratio, and never an eps floor (ROADMAP 8g).
  - `method="spatialqm"`: SpatialQM's log10 dataset scalars.
- **Sparsity** (`getSparsity`): zero fraction, matches to 3e-16.
- **Detection entropy:** per-cell normalised entropy by default; `method="spatialqm"` = `getEntropy` (bits of
  the count-value histogram of the whole matrix).
- **Transcripts per area / per cell** (`getTxPerArea`, `getTxPerCell`): mean (SpatialQM's) and median.
- **Transcript QC** (`transcript_qc`, #21): FTC and the negative-codeword FDR from a transcripts table.
  SpatialQM divides its FDR by 100; we report the fraction.

Add to section **2.6** (or next to the deconvolution bullet):

- **Decomposition purity** (`decomposition.decomposition_purity`): 1 - the largest normalised deconvolution
  weight per cell (RCTD full mode). Higher = more of the cell's expression assigned outside its dominant type.
  Spatial Touchstone used this as a segmentation-error readout. It depends on the reference; the code follows
  SegTraQ `rctd_contamination` (MIT).

In §6 (implementation status) and the §5 master table, add `tx_per_cell`, `transcript_qc` and
`decomposition_purity`. In §9, add: SegTraQ, github.com/LazDaria/SegTraQ (MIT), `sp.rctd_contamination`.

Also for the integrator (files this package may not edit):

- `tests/test_thresholds.py` DEFAULTS: no new row needed (the new constants are not parameter defaults).
- `thresholds.py` `signal_qc.snr_eps`: since 0.7.0 this is the opt-in 0.6 value (`eps=None` by default). Its
  source text could say "opt-in (eps=1e-9 restores 0.6)" (append-only rule: not edited here).
- `sam.__all__` grows by 3: `tx_per_cell`, `transcript_qc`, `decomposition_purity`. That is automatic through
  the module `__all__`s; the contract text "33 names" becomes 36 if no other package adds names.
- `benchmarks/synthetic/sweep.py` `_sig` reads `median_mean_signal_ratio`, which now rests on cells with a
  control count. No adapter edit is needed. Its baselines change (see Roadmap 8g).

## Roadmap

- **8g SNR eps floor: done.** No eps floor by default. The per-cell value is NaN without a control count, and the
  result reports `pooled_mean_signal_ratio`, `n_cells_with_controls` and `frac_cells_with_controls`. `eps=1e-9`
  is the old recipe, exact.
  - Evidence: `signal_qc.py` `signal_to_noise`; tests `test_signal_to_noise_has_no_eps_floor`,
    `test_signal_to_noise_skips_when_no_cell_has_a_control_count`, `test_signal_to_noise_old_eps_reproduces_0_6`,
    and golden `test_signal_to_noise_cell_does_not_sit_on_the_eps_floor` (crop: 224 of 2000 cells, median 6.491,
    pooled 9.511; 0.6 gave 27.956).
  - Sweep values (mean of 3 seeds), old -> new:
    - baseline `sig.snr` 30.33 -> 5.11;
    - `ambient` 0.02 / 0.20: 5.54 / 3.26 -> 4.89 / 3.25;
    - `depth_loss` 0.94: 26.71 -> 1.13.
- **#8 Moran shuffled-coordinates null: done (that item only).** Test
  `test_moran_shuffled_coordinates_null[squidpy|spatialqm]`: the separation is > 0.5 on the structured section
  and < 0.05 in absolute value with shuffled coordinates. The other #8 gaps belong to other packages.
- **#21: done.**
  - Decomposition purity: `decomposition.py` `decomposition_purity`; `tests/test_decomposition.py` (hand values
    0 / 0.75 / 0.5, NaN rows, DataFrame weights, skips).
  - TPC / FTC / FDR: `signal_qc.py` `tx_per_cell` and `transcript_qc`. Golden tests `test_tx_per_cell` and
    `test_transcript_qc` reproduce TPC 128.261, FTC 0.999798 and FDR 5.92514e-05 (`fdr / 100`) with SpatialQM's
    name rule. With the feature-type rule: TPC 128.1905, FDR 6.64377e-05. Bytes names and string `-1` /
    `"UNASSIGNED"` ids: `test_transcript_qc_reads_bytes_and_string_ids`.
  - Moran spatialqm edge cases: `test_moran_spatialqm_small_control_group_is_nan_not_a_crash[0|1|2|20]`,
    `test_moran_spatialqm_feature_subset_matches_full_run`, golden
    `test_moran_spatialqm_feature_subset_keeps_getmorans_values`.
- **#7 SpatialQM parity: done (signal QC row).** Every SpatialQM scalar we offer is compared in
  `tests/goldens/signalqc_spatialqm/test_signalqc_spatialqm_metrics.py`. The golden is tied to its inputs by
  sha256, and the regeneration is exact on both Pythons.
- **#13 tutorial: partial (scale check only).** On the full section (167,780 cells):
  - `moran_signal` takes 16.8 s for squidpy (separation 0.152) and 3.2 s for spatialqm (separation 0.115).
  - `signal_to_noise`: 18.1 % of cells have a control count; median 7.15, pooled 9.18 (0.6: 28.51).
  - `tx_per_cell` median 164, the same as the 10x summary parse in SOURCE.md.
  - No tutorial written.
