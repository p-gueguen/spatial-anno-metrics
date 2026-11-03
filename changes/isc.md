# isc (package wp/isc)

## CHANGELOG

### Changed
- `inter_sample_consistency`: ISCG is now computed in scTypeEval's default space, cosine on the first 30 centred PCs of
  the pseudobulks (`reduction="pca"`, `n_pcs=30`, capped at `min(n_profiles, n_genes) - 1`; exact SVD). It matches
  scTypeEval 1.0.0's README ISCG (irlba) within 1.4e-6 on the Yerly 2022 golden, the space the 0.74 guide comes
  from. **Breaking**: ISCG values move, by +0.012 to +0.073 per label on the golden (mean +0.030 / +0.050); on the
  synthetic sweep by -0.012 to +0.025 (scenario means). old -> new: raw-space ISCG -> `reduction=None`
  (reproduces 0.6.x per-label ISCG to 2.3e-14).
- `inter_sample_consistency`: headline `iscg` / `iscl` are means over labels (paper, scTypeEval
  `mean_over_types`), like `integrated`. **Breaking** when labels have unequal sample counts: golden 0.855469 ->
  0.854361 (ISCG, raw space), over-split ISCL 0.569124 -> 0.575964. old -> new: `iscg` / `iscl` (pseudobulk
  means) -> `iscg_profile_mean` / `iscl_profile_mean`.
- `inter_sample_consistency`: the reciprocal-classification match follows scTypeEval's `match_pairs` rule: a mutual
  best hit between two labels counts as a match only when both labels occur in both samples (scTypeEval stores 0.5
  and reads it as 1), and a sample with fewer than 2 scored profiles matches nothing. **Breaking** only when a
  scored label is missing from some samples and is part of a cross-label mutual hit (not triggered on the golden or
  in the synthetic sweep: 0 ISCL values moved).
- `inter_sample_consistency`: pseudobulks are aggregated as `indicator.T @ X`; a sparse matrix is no longer
  densified (values equal to 2.3e-14).
- `isc_regime` / `regime`: marked experimental; the result carries `regime_status = "experimental"`.

### Added
- `inter_sample_consistency(..., reduction="pca" | None, n_pcs=30, rcm_genes=None | "markers" | <gene list>)`.
  `rcm_genes="markers"` runs the reciprocal classifier on SingleR's classic marker union of each training sample
  (de.n = round(500 x (2/3)^log2(N)), positive differences; no fine-tuning or pruning); a gene list (e.g. HVGs)
  restricts it to those genes. ISCG always uses every gene.
- Result keys `iscg_profile_mean`, `iscl_profile_mean`, `regime_status`, `reduction`, `n_pcs`, `rcm_genes`.
- `thresholds.THRESHOLDS["isc"]["n_pcs"]` (30) and `["rcm_marker_n"]` (500). Integrator: add the DEFAULTS row
  `("isc", "inter_sample_consistency", "n_pcs"): ("isc", "n_pcs")` to `tests/test_thresholds.py` (not my file;
  `tests/test_isc.py::test_iscg_reduction_pca_default_and_raw_space` asserts the same default meanwhile).
- ISC golden: tests for the 30-PC ISCG (within 2e-6 of scTypeEval), the label-mean headline (= scTypeEval
  `mean_over_types`) and the markers classifier vs SingleR (0 / 1,128 and 28 / 1,891 differing pairs).

### Fixed
- `benchmarks/parity/isc_yerly2022/run_sctypeeval.R` calls `scTypeEval:::purge_label` instead of copying its regex
  (GPL hygiene); re-run, `expected.json` is byte-identical (md5 eb118de101ca5411f716c8c238d75204).
- `benchmarks/parity/isc_yerly2022/make_golden.py` ran against the removed 0.6 API (`counts_layer=`, 4-tuple
  `_pseudobulk_profiles`); now uses `layer=` and the new `reduction=` (re-run: 30-PC gap 1.4e-6 / 8.5e-7,
  `counts.npz` / `obs.csv` byte-identical).

## README rows

| `inter_sample_consistency(adata, label_key, sample_key, *, layer=None, reduction="pca", n_pcs=30, rcm_genes=None)` | **ISC on two axes** (Garnica et al. 2026), checked against scTypeEval 1.0.0: `iscg` (2-label silhouette on pseudobulk cosine in 30 centred PCs = cross-sample compactness; `reduction=None` for the raw space), `iscl` (silhouette on reciprocal-classification match = cross-sample replicability; `rcm_genes="markers"` or a gene list restricts the classifier), `integrated` = mean over labels of ISCG x ISCL; headline values are means over labels; per-label `regime` via `isc_regime` (consistent / merge / split / remove, **experimental**), `dropped_types` below 5 samples x 10 cells | scTypeEval / Garnica et al. 2026 |

## Catalog

Section **3b. Inter-sample consistency - ISC (scTypeEval)**, replace the "**In this package**" paragraph with:

**In this package**: `inter_sample_consistency(adata, label_key, sample_key, *, min_cells=10, min_samples=5, layer=None, reduction="pca", n_pcs=30, rcm_genes=None)` returns `iscg`, `iscl`, `integrated` (each a mean over labels, as the paper), per-type dicts, the pseudobulk-weighted means (`iscg_profile_mean`, `iscl_profile_mean`), a per-type `regime` from `isc_regime()` flagged `regime_status = "experimental"`, and `dropped_types` (labels below the minimums are reported, not folded in as 1-cell pseudobulks). Parity with scTypeEval 1.0.0 on a Yerly 2022 subsample (`tests/goldens/isc_yerly2022`): ISCG equals scTypeEval's default 30-PC ISCG within 1.4e-6 (irlba vs exact SVD) and, with `reduction=None`, its `reduction = FALSE` ISCG within 1.3e-14; given scTypeEval's RCM distances, ISCL equals its ISCL within 1e-16; the RCM rule (0 for a mutual best hit, 0.5 read as 1 when either label is absent from the other sample, no match for a sample with < 2 profiles) is scTypeEval's. The classifier differs: nearest profile by Spearman over all genes (or `rcm_genes`) in place of SingleR with fine-tuning and pruning. On the authors' labels both make the same mutual hits on 1,128 / 1,128 pseudobulk pairs; on a randomly over-split label the two halves are near-identical and 30 / 1,891 pairs differ (28 with `rcm_genes="markers"`), so per-half ISCL is classifier noise (up to 0.30 apart) while the mean over labels agrees within 0.0002. Do not run ISCL on toy data with flat, low-noise background genes: Spearman ranks there let one foreign cell lift its type's genes above the background, and ISCL collapses (8 types x 64 genes: 0.92 -> 0.14 at 10 % relabelling; 1.00 with a graded per-gene baseline). Regime calls are experimental: on the synthetic sweep only over-splitting reads right (`merge`, 6 / 6); merged labels read `consistent` (Fibro+Epi ISCG 0.984, ISCL 1.000), 20 % relabelling reads `split` for 13 / 18 labels, and a per-sample batch effect `split` for 18 / 18. ISCL on a targeted spatial panel (a few hundred genes) is not validated by the paper; run the over-partitioning positive control (`tests/test_isc.py`, task-2 replica) on your own data before trusting a threshold.

## Roadmap

- **4** -> **partial**. Done: `rcm_genes="markers"` (SingleR classic marker union) and `rcm_genes=<HVG list>`
  (`src/spatial_anno_metrics/isc.py` `_classic_markers`, `_rcm_dissimilarity`; `tests/test_isc.py::
  test_rcm_genes_restricts_the_classifier`, `test_classic_markers_union_of_top_positive_differences`). Real data:
  our Spearman classifier matches SingleR's mutual hits on 1,128 / 1,128 pairs (golden `annotation`), 1,861 /
  1,891 on the over-split (markers: 1,863 / 1,891). The 64-gene collapse is a fixture artefact, not the gene set:
  in `_isc_adata(n_types=8, n_genes=64)` every gene belongs to one type's block over a flat N(0.2, 0.1) background,
  so in a pseudobulk the other 56 genes rank at random and one foreign cell (1 / 38 cells, +0.13 on its 8 genes vs a
  0.016 standard error) lifts its block to ranks 49-56; two pseudobulks that hold one cell of each other's type
  become mutual best hits (e.g. T0/S1 with one T3 cell <-> T3/S0 with one T0 cell; 23 cross-label vs 21 same-label
  cross-sample mutual hits of 80 same-label pairs). ISCL 0.924 / 0.266 / 0.138 / 0.073 at 0 / 5 / 10 / 20 %
  noise (10 seeds), identical with `rcm_genes="markers"` (de.n = 148 >= 64 genes, so the union is every gene);
  1.000 at 0 and 10 % with Pearson or with a fixed lognormal per-gene offset. Integrated formula settled: mean over
  labels of ISCG x ISCL (paper), not the README's geometric mean (golden SOURCE.md). Open: SingleR fine-tuning
  and pruning are not ported.
- **7** (ISC row) -> **done**: 30-PC ISCG (`test_iscg_default_30_pcs_matches_sctypeeval`, 2e-6), label-mean
  headline (`test_headline_is_the_mean_over_labels`), markers vs SingleR (`test_rcm_markers_vs_singler`), plus the
  existing raw-space / ISCL / RCM-gap tests; generator calls `scTypeEval:::purge_label`, `expected.json`
  byte-identical on re-run.
- **8d** -> **partial** (flagged, not fixed). Re-run after this package (sweep `--seeds 3 --only label_noise merge
  over_split batch`, `/srv/GT/analysis/pgueguen/sam_dev/runs/wf3-isc/sweep_after/check_isc_regimes.csv`): calls
  identical to before the change. clean 18 / 18 `consistent`; over-split (Epi halves) 6 / 6 `merge` (right);
  merge 1 (Fibro+Epi) 15 / 15 `consistent`, Fibro+Epi ISCG 0.984 ISCL 1.000 (wrong, expected `split`); merge 3
  9 / 9 `consistent`; label noise 20 % 13 / 18 `split` (ISCG 0.690, ISCL 0.978); batch sd 1.0 18 / 18 `split`
  (ISCG 0.618, ISCL 1.000). ISCL at 40 % noise 0.840. Regime calls now carry `regime_status = "experimental"`
  (`isc.py`, `tests/test_isc.py::test_regime_calls_flagged_experimental`) and the numbers are in the docstring.
- **15** -> **partial**: pseudobulk is `indicator.T @ X` (sparse in, no densify; `isc.py::_pseudobulk_profiles`,
  `tests/test_isc.py::test_pseudobulk_sparse_equals_dense`). The other hotspots in the row (fidelity, `_embedding`,
  `label_stability`, `neighbor_contamination`) belong to other packages.
