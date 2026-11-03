# wp/purity - fragment for the integrator

Package: `purity.py`, `fidelity.py` (+ the `purity` block of `thresholds.py`, the `nmp_txsim` / `coexpr_txsim`
goldens and generators). Roadmap items #1, #8 (ncp test), #8c, #8e, #10, #15.

## Decision: `nmp` default is `score="salas"`

`nmp` now defaults to the cited definition, Salas et al. *Nat Methods* 2025 as txsim / OpenProblems compute it
(`negative_marker_purity_reads`, scaled so 0 = mean of 20 label permutations, 1 = clean). The v0.6.0 score stays
available as `score="raw"`. Reasons:

- The function is named after Salas's metric; with the raw score as default the name promised a number no
  other tool produces (roadmap #1). The Salas score now matches txsim to <= 8.7e-19 on the golden.
- It carries its own null (random labels), which the raw score lacks; on the synthetic sweep it detects spillover
  (1.000 -> 0.968 at 0.4) and label noise (1.000 -> 0.515 at 0.4), stays flat under genuine co-expression
  (>= 0.9996) and under panel growth (1.000 at +960 genes), where the raw score drifts with panel size
  (0.0014 -> 0.0001).
- Kept for `score="raw"`: a per-cell value (`obs['sam_nmp']`), which Salas's per-section score does not have.

Cost: the default call writes nothing to `adata` and returns `nmp` (higher = cleaner), not `mean_nmp`.
`benchmarks/synthetic/sweep.py::_nmp` reads `mean_nmp` and gets NaN for all 96 rows now: change it to
`r = sam.nmp(a, ref, "label", "cell_type"); {"nmp": r.get("nmp")}` with orientation `+1` (or pass `score="raw"`
to keep the old column).

## CHANGELOG

### Added
- `coexpression_similarity(adata, reference, *, label_key=None, ref_label_key="cell_type", min_cells=20,
  layer="counts")`: txsim's co-expression similarity (`1 - |rho_ref - rho_section| / 2` over gene pairs, on
  `log1p(1e4 * counts / total)`), beside `ncp`; with `label_key`, per cell type as txsim's `by_celltype=True`.
  Golden-tested against txsim (max |diff| 6.7e-16).
- `nmp(..., score="salas")`: Salas et al. 2025 negative marker purity as txsim `negative_marker_purity_reads`,
  scaled against `n_permutations=20` label permutations (`seed=0`); returns `nmp`, `nmp_unscaled`, `nmp_null`,
  `per_type`, `per_gene`, `n_negative_pairs`, `types_scored`, `types_absent_from_reference`,
  `types_below_min_cells`. Golden-tested against txsim (overall, per type, per gene, 20 permuted values).
- `crisp_purity` and `mecr` report a detection-matched null (each marker gene's detections reassigned to random
  cells, `n_null=20` draws, `seed=0`): `crisp_purity_null` / `mecr_null` and `vs_null` (observed co-detection /
  null co-detection: 0 = exclusive lineages, 1 = no more exclusive than independent markers). `n_null=0` skips it.
- `pmp`, `marker_program_fidelity`, `marker_gene_overlap` report `labels_without_markers` and `no_dict_frac`
  (labels with no on-panel marker set are not scored); `crisp_purity` and `mecr` report
  `lineages_without_panel_markers`.
- Thresholds (`purity`): `null_draws` 20, `nmp_salas_max_share` 0.005, `nmp_salas_min_cells` 10,
  `nmp_n_permutations` 20, `coexpr_min_cells` 20, `coexpr_target_sum` 1e4.

### Changed
- **Breaking:** `nmp` default `score`: raw per-cell fraction (`mean_nmp`, lower = cleaner, writes
  `obs['sam_nmp']`) -> `score="salas"` (`nmp`, higher = cleaner, writes nothing). Pass `score="raw"` for the old
  score; on the golden it reproduces 0.001030 / 0.004628 exactly, and on the synthetic sweep 78/96 values are
  bit-identical; the other 18 are the NaN fix below (15 over-split / merge values, 3 merge skips).
- **Breaking:** `nmp(score="raw")`: cells of a type absent from the reference, of a type with no negative gene,
  or with no counts on the shared genes are NaN, not 0.0; `mean_nmp` is the mean over the defined cells (new
  `n_defined`); a section with no scorable cell is a skip. `types_without_negatives` -> two lists,
  `types_absent_from_reference` and `types_without_negatives` (in the reference, no gene passes).
  Sweep: over-split / merged labels no longer pull `mean_nmp` towards 0 (`over_split` level 3: 0.000565 ->
  0.001349 for seed 0, the clean value is 0.0014); `merge` level 3 is a skip (was 0.0 on all 3 seeds).
- `ncp` pair counts are one sparse `P.T @ P` instead of a Python loop over pairs; same values (sweep identical,
  a test checks it against the loop).
- `marker_program_fidelity` / `marker_gene_overlap` slice marker columns and compute column moments on sparse
  matrices instead of densifying `X`; CRISP / MECR keep the marker detections sparse. Xenium breast section
  (167,780 cells x 541 genes): fidelity peak memory 1132 -> 57 MB, overlap 2184 -> 650 MB and 12.5 -> 3.3 s,
  same values. The CRISP / MECR null adds about 2 s each there.
- `marker_sets` keys are compared as str (like the labels) in `crisp_purity`, `mecr`, `pmp`,
  `marker_program_fidelity`, `marker_gene_overlap`: integer cluster ids with an int-keyed dict now score instead
  of giving NaN silently; `0` and `"0"` both present raises `ValueError`.

### Fixed
- `nmp(score="raw")` detection fractions are detecting cells / type cells, so the inclusive cuts
  (`<= max_neg_frac`, `>= min_other_frac`) no longer depend on storage: with COL1A1 in exactly 11/220 T reference
  cells a CSR reference gave 0.000 (T unscorable) and a dense one 0.200; both give 0.200 now.
- Attribution of the raw negative rule: only the `<= 5 %` cut is SegTraQ's (`markers_from_reference`); the
  `>= 25 %` other-type floor is a package choice (SegTraQ uses `min_pos_frac=0.1` plus a DE step; its 0.25 is
  MECR's detection floor). `thresholds.py` source text and `provisional` flag corrected.

### Tests
- `test_ncp_...` no longer accepts `status == "skipped"`; new tests pin the `nmp` boundaries (`<=`, `>=`, max over
  >= 3 other types), the strict Salas 0.5 % cut, NaN for unscorable cells, sparse == dense, the CRISP/MECR null
  (positive control: background genes used as markers give `vs_null` 1.0 +- 0.1), int keys, coverage. 18 mutants
  of the new code each fail at least one test.

## README rows

Replace the `crisp_purity` / `mecr` / `pmp`, `nmp` and `ncp` rows, and add a `coexpression_similarity` row:

```
| `crisp_purity(adata, marker_sets)`, `mecr(adata, marker_sets)`, `pmp(adata, marker_sets, label_key)` | marker co-expression **contamination purity**: CRISP per-cell impurity + dataset purity (rule and formula from CRISP, Rose et al. 2026, detected at >= 1 count), MECR (mutually-exclusive co-detection rate over lineage pairs; a marker two sets share is not co-detection, nested sets are skipped), both against a detection-matched null (`vs_null`: 0 = exclusive, 1 = as random as independent markers); PMP (per-cell marker purity, panel-size-invariant). None of them tells spillover from genuine co-expression: use `nmp` | Rose et al. 2026; Hartman & Satija 2024; SpatialScribe QC L3 |
| `nmp(adata, reference, label_key)` | negative marker purity, Salas et al. 2025 as txsim computes it (default `score="salas"`): reads on (gene, type) pairs the reference does not express, scaled so 0 = random labels and 1 = clean; tolerates co-expression the reference has. `score="raw"`: per-cell fraction of counts on genes the reference detects in <= 5 % of the type's cells and >= 25 % of another type's (NaN when a type cannot be scored) | Salas et al. 2025; txsim; raw `<= 5 %` cut from SegTraQ |
| `ncp(adata, reference)` | non-coexpression preservation: fraction of reference-defined non-coexpressed gene pairs that stay non-coexpressed in the section (MECR-style, reference-derived pairs) | Salas et al. 2025 (NCP, as worded in the paper) |
| `coexpression_similarity(adata, reference)` | gene-gene correlation agreement with the reference, `1 - |rho_ref - rho_section| / 2` averaged over pairs; `label_key=` per cell type | txsim |
```

The `marker_program_fidelity` / `marker_gene_overlap` rows: append "; reports `labels_without_markers`". README
line 89 comment: `obs['crisp_impure']` -> `obs['sam_crisp_impure']`. README line 122 (SegTraQ credits): "the
`nmp` negative-marker rule (`markers_from_reference`, without its DE step)" -> "the `<= 5 %` cut of
`nmp(score="raw")` (`markers_from_reference`; its other-type floor is 0.1 plus a DE step, ours 0.25)". Add to the
credits: CRISP (Rose et al. 2026, entry in Catalog below) and txsim (`nmp` Salas score and
`coexpression_similarity`, re-implemented from the definitions: txsim has no licence, run only for the goldens).

Count: `sam.__all__` grows from 33 to 34 names (`coexpression_similarity`).

## Catalog

`docs/cell-annotation-quality-metrics.md`:

**§3c table** (line 204-205), replace the MECR and CRISP rows:
```
| **MECR** (Mutually Exclusive Co-expression Rate) | over lineage marker pairs, `#(both>0)/#(either>0)`, averaged; each pair is judged on the markers the two sets do not share (CD3E in both CD4 T and CD8 T sets is not co-detection), nested sets (T vs CD4 T) are skipped; reported with `mecr_null` and `vs_null = mecr / mecr_null` | lower cleaner | `purity.py` (mecr); SPLIT compare |
| **CRISP purity** | per-cell: impure if it detects markers of ≥2 lineages; purity = 1 − N_impure / N_marker-positive (Rose et al. 2026). Our variant detects at >0 counts (CRISP default: >2) and judges each lineage pair on its unshared markers - not comparable to published CRISP scores. Reported with a detection-matched null and `vs_null` = observed / null impurity | higher cleaner | `purity.py` (`crisp_purity`) |
```

**§3c, new paragraph after the table:**
```
**The null, and what it cannot do.** At the ≥ 1-count rule, CRISP purity is 0.80 on a perfectly clean synthetic section (off-target counts at 0.1 % of the marker rate) and falls with depth alone, so each of `crisp_purity` and `mecr` is reported against a detection-matched null: every marker gene keeps its number of detecting cells, reassigned to random cells (20 draws). `vs_null` = observed co-detection / null co-detection reads 0.22 (CRISP) and 0.12 (MECR) on the clean section, rises with spillover (0.36 / 0.18 at 0.4) and ambient RNA (0.98 / 0.70 at 0.2), and stays stable under label changes. Genuine co-expression raises it as much as spillover (0.35 / 0.17 when Myeloid truly expresses Fibroblast markers at 40 %): no reference-free metric can tell the two apart. NMP with a matched reference can, because the reference carries the co-expression: Salas NMP stays ≥ 0.9996 in that scenario while it drops to 0.968 at spillover 0.4 (`benchmarks/synthetic`, 3 seeds).
```

**§3f NMP bullet** (line 241), replace the `In purity.py ...` sentences with:
```
In `purity.py`, `nmp` computes it as Salas et al. define it and txsim implements it (`negative_marker_purity_reads`, matched to ≤ 1e-18 in `tests/goldens/nmp_txsim`): per type, the type-balanced share of each gene's mean count; a (gene, type) pair is negative when the reference type holds < 0.5 % of the share; NMP = 1 − the section's excess share on those pairs, scaled so 0 = the mean over 20 random label permutations and 1 = no more than the reference (the paper's scaling; its baseline comes from one binning reference, ours from the section itself). Higher cleaner. Because the reference carries genuine co-expression, it is the purity metric to trust where cell types really share markers. `score="raw"` keeps the v0.6.0 per-cell score: the fraction of a cell's counts on genes the reference detects in ≤ 5 % of its type's cells (SegTraQ `markers_from_reference`'s cut) and ≥ 25 % of another type's (a package choice: SegTraQ uses 0.1 plus a DE step, a 10-cell minimum and drops MT-/RPL/RPS genes); cells it cannot score are NaN. The earlier bottom-decile definition picked near-silent genes on large panels: 10% fibroblast spillover into T cells scored 0.004, same as clean, on a synthetic 5,000-gene panel.
```

**§3f, new bullet after NCP:**
```
- **Co-expression similarity** (txsim) — `1 − |ρ_reference − ρ_section| / 2` per gene pair (Pearson on log-normalised counts), averaged; optionally within each cell type. In `purity.py` (`coexpression_similarity`), matched to txsim in `tests/goldens/coexpr_txsim`. On the synthetic sweep it barely moves with spillover (0.979 → 0.978 at 0.4); prefer NCP / NMP for contamination.
```

**§3f NCP bullet**: append "`purity.py` `ncp` follows the paper's wording (fraction of pairs); the authors' code averages a clipped excess instead (`tests/goldens/coexpr_txsim/SOURCE.md`)."

**Line 414**: `Mitchel *Nat Genet* 2025` -> `Mitchel *Nat Genet* 2026`.

**§9 line 439**: replace `Center for Spatial OMICs (CRISP)` with
```
Rose et al. *bioRxiv* 2026, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947) (CRISP); Mitchel et al. *Nat Genet* 2026, [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4) (segmentation admixture, cellAdmix)
```

**§9, new bullets after the SpatialQM entry:**
```
- **CRISP** - Rose, J. R., Rose, E. S., Assumpção, J. A. F., Pathak, H., Peck, H. E., Sasser, L. E., Patel, C. J., Vanover, D. & Santangelo, P. J. "CRISP enables comparisons of image-based spatial transcriptomic segmentation quality across ten organs", *bioRxiv* 2026, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947); code [github.com/santangelo-lab/CRISP](https://github.com/santangelo-lab/CRISP) (R `crispRSeg`, Python `crispySeg`, MIT; Santangelo lab, Emory University and Georgia Tech). Source of the CRISP purity definition, built on MECR (Hartman & Satija 2024): a cell is impure when markers of ≥2 mutually exclusive cell types pass the detection threshold, and `purity = 1 - N_impure / N_marker-positive`. `crisp_purity` keeps that rule and formula with two changes: a marker counts as detected at >0 counts (CRISP: `> threshold`, default 2, compound markers summed first, threshold recommended at about half the mean marker expression), and each lineage pair is judged on the markers the two sets do not share (nested pairs skipped; CRISP expects curated disjoint lists). On disjoint marker sets it equals CRISP run at threshold 0; its values are not comparable to published CRISP scores. `[preprint]`
- **Segmentation admixture** - Mitchel, J., Gao, T., Petukhov, V., Cole, E. & Kharchenko, P. V. "Impact and correction of segmentation errors in spatial transcriptomics", *Nat Genet* 58, 434-444 (2026), [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4); code [github.com/kharchenkolab/cellAdmix](https://github.com/kharchenkolab/cellAdmix). Re-analysing imaging-based data from multiple tissues and platforms, it shows that segmentation errors put other cell types' marker transcripts into cells and that this admixture dominates cell-state analyses (niche differential expression, neighbour influence, ligand-receptor); cellAdmix removes it by factorizing subcellular molecular neighbourhoods, "analogous to doublet filtering". It is the phenomenon MECR / CRISP / PMP read out; transcript-level correction is out of scope here. `[peer-reviewed]`
- **txsim** - [github.com/theislab/txsim](https://github.com/theislab/txsim) (commit `bfdfaa90`; no licence file, so run for the goldens and never copied); cited by OpenProblems as [10.1101/2023.02.13.528102](https://doi.org/10.1101/2023.02.13.528102). Reference implementation of Salas et al.'s NMP (`negative_marker_purity_reads`) and of `coexpression_similarity`; `nmp` (default) and `coexpression_similarity` re-implement both from their definitions and match them on the synthetic goldens. `[tool]`
```

Not added (unverified, per the literature check): the claim that Mitchel et al. show cell-type classification is "comparatively robust".

## Roadmap

- **1** -> done. `nmp(score="salas")` is the default and matches txsim reads (overall, per type, per gene, 20
  permuted values; max |diff| 8.7e-19; `tests/test_purity.py::test_nmp_salas_matches_txsim_reads`); NaN, not 0,
  for unscorable types (`test_nmp_raw_unscorable_cells_are_nan_and_excluded_from_the_mean`); the two type lists
  say why; `coexpression_similarity` added beside `ncp` (`test_coexpression_similarity_matches_txsim`, 6.7e-16).
  Moldia's reads criterion (an older relative-to-mean rule) is recorded in the golden, not implemented.
- **8** (the `ncp` item only) -> done: `test_ncp_scores_clean_as_preserved_and_contamination_as_lost` asserts
  `status == "ok"`, 1.0 clean and 0.0 contaminated; `test_ncp_threshold_is_inclusive` pins `<=`.
- **8c** -> done for `pmp`, `marker_program_fidelity`, `marker_gene_overlap` (`labels_without_markers`,
  `no_dict_frac`) and `crisp_purity` / `mecr` (label-free, so `lineages_without_panel_markers`); tests
  `test_marker_metrics_report_labels_and_lineages_without_markers`, `test_fidelity_reports_labels_without_markers`.
  `neighbor_contamination` (spatial.py) is another package's.
- **8e** -> partial. The detection-matched null and `vs_null` are reported (`test_crisp_mecr_report_a_detection_matched_null`,
  positive control `test_null_positive_control_independent_markers_have_no_excess`). It calibrates depth / detection
  rates, but a clean section still reads `vs_null` 0.22 rather than 0, and it cannot separate genuine co-expression
  from spillover (coexpression 0.4: 0.35 vs spillover 0.4: 0.36); the docs now say NMP is the one to trust there,
  and the sweep backs it (Salas NMP >= 0.9996 under co-expression).
- **10** -> done for the catalog text (CRISP = Rose et al. 2026, Santangelo lab; Mitchel 2026) and the
  `crisp_purity` docstring; the README credits and catalog edits above are for the integrator. The
  spatial-scribe `cell-annotation-qc.md` carries the same two errors (another repository).
- **15** -> partial (this package's part done): fidelity / overlap no longer densify `X`, CRISP / MECR detections
  stay sparse, `ncp` pair loop vectorised (`P.T @ P`). `coexpression_similarity` densifies cells x shared genes
  (marked `ponytail:`). Pseudobulk, `_embedding`, `label_stability`, `neighbor_contamination` are other packages'.

## Integrator notes (files this package may not edit)

- `tests/test_thresholds.py` DEFAULTS rows: `("purity", "crisp_purity", "n_null"): ("purity", "null_draws")`,
  `("purity", "mecr", "n_null"): ("purity", "null_draws")`, `("purity", "nmp", "n_permutations"):
  ("purity", "nmp_n_permutations")`, `("purity", "coexpression_similarity", "min_cells"): ("purity",
  "coexpr_min_cells")` (all checked: equal value and type).
- `tests/test_validate.py::_NAN_LABEL_CALLS`: add `"coexpression_similarity": lambda sam, a:
  sam.coexpression_similarity(a, a.copy(), label_key="cell_type", ref_label_key="truth")` (checked: raises the
  expected ValueError).
- `benchmarks/synthetic/sweep.py::_nmp`: see the decision above.
