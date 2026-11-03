# Roadmap - field scan, 2026-11-03

Where `spatial-anno-metrics` v0.6.0 stands against the field, and what to do next, ranked. Built from
four parallel scans (competing tools, 2024-26 metric literature, validation design, packaging) plus
direct checks against the code at `81aa4ea`, updated for the SegTraQ commit `7008d41` (v0.6.0: shared-marker-safe
MECR / CRISP, new NMP negatives, `neighbor_contamination`). Each item carries its source; items marked
**[checked]** were re-verified by reading the code, re-running a probe, or reading the paper.

**What is ours alone** (no other tool found doing it): depth-thinned `panel_resolvability`,
ISC regime calls in Python, PMP, the lineage-coverage flag, annotator-agnostic `label_stability`.
**Where we only re-implement**: ISC (scTypeEval is the reference), NMP/NCP (txsim / OpenProblems iST),
signal QC (SpatialQM is a superset), conformal (scConform adds hierarchy), `neighbor_contamination` and the NMP
negative rule ([SegTraQ](https://github.com/LazDaria/SegTraQ), which also covers the transcript-level segmentation QC we do not). The package's value is the
consolidation plus the ours-alone items - so the re-implementations must match their sources, and the
ours-alone items must be validated.

---

## P0 - correctness and credibility (do before anyone else uses it)

| # | Item | Why / evidence | Effort | Touches |
|---|---|---|---|---|
| 1 | **NMP has the cited name but a different metric.** v0.6.0 fixed the worst part: negatives were the reference's bottom-decile genes (near-silent on large panels, so 10 % spillover scored like clean cells); they are now SegTraQ's rule, ≤ 5 % of the type's reference cells and ≥ 25 % of another type's, and `types_without_negatives` names the types it cannot score. Still open: those types score 0 (look clean) rather than NaN, and the output is a raw count fraction, lower = cleaner. Salas 2025: negative marker = gene in < 0.5 % of a type's cells, scaled so 0 = random labels, higher = purer; txsim and OpenProblems `task_ist_preprocessing` implement that, so our `nmp` values cannot be compared with theirs. | [checked] `purity.py` `nmp` at `7008d41`; Salas *Nat Methods* 2025 [10.1038/s41592-025-02617-2](https://doi.org/10.1038/s41592-025-02617-2); [txsim](https://github.com/theislab/txsim); [task_ist_preprocessing](https://github.com/openproblems-bio/task_ist_preprocessing) | M | `purity.py`, catalog §3f |
|   | Fix: add the Salas-scaled score beside ours (or rename ours `nmp_raw`), return NaN - not 0 - for a type in `types_without_negatives`, golden-test against the Moldia code ([Xenium_benchmarking](https://github.com/Moldia/Xenium_benchmarking), data Zenodo [11124988](https://doi.org/10.5281/zenodo.11124988)). Add txsim's `coexpression_similarity` beside `ncp`. | | | |
| 2 | **PAS has the cited name but a different threshold.** Ours: abnormal if < 20 % of k = 15 neighbours share the label. Shang & Zhou 2022: abnormal if the label differs from ≥ 6 of 10 neighbours (same-label ≤ 40 %). Our PAS reads systematically lower and cannot be compared with published values. | [checked] `spatial.py:30-54`; [10.1038/s41467-022-34879-1](https://doi.org/10.1038/s41467-022-34879-1) | S | `spatial.py`, catalog §3e |
|   | Fix: default `k=10, pas_threshold=0.4` (strictly: abnormal if differs from ≥ 6), keep ours as an option; add CHAOS (same paper) beside it. | | | |
| 3 | **`internal_validity` cannot run at ≥ 1e5 cells**, and `annotation_quality` always calls it: Ward clustering on every cell, O(N²) memory. Measured: 20 k cells 9 s / 3.6 GB, 40 k cells 67 s / 13 GB. | [checked] `eval_metrics.py:106` (Ward on `E`, before the `sel` subsample at :122) | S | `eval_metrics.py` |
|   | Fix: run Ward on the existing subsample. Same pass: vectorise neighbourhood purity (:88) and avg-similarity (:118), `cdist` for orbital-medoid (:98), sparse `pos.T @ pos` for `ncp` pair loops. | | | |
| 4 | **ISCL collapses at low label noise on the 8-type fixture**: 0.927 → 0.159 at 10 % random relabelling (10 seeds), while ISCG stays 0.997 (expected: uniform relabelling shifts every pseudobulk alike). Likely cause: Spearman nearest-profile over all genes; the paper uses SingleR on HVGs. Synthetic fixture only (64 genes) - real-data behaviour untested. | [checked] re-ran probe, `eval_metrics.py:198-224` | M | `eval_metrics.py` |
|   | Fix: restrict RCM to HVGs / marker genes as the paper does; then parity-test against scTypeEval (below, #7). Also `integrated` = mean of per-type ISCG × ISCL products (:292), while the paper says "product" and the scTypeEval README "geometric mean"; the parity test settles which, then document it. | | | |
| 5 | **AQI validation and tuning**: tracked outside this repo. | - | L | - |
| 6 | **Docstring contradicts code**: `annotation_quality` says "with none of them AQI is a pure coherence index"; the core is {C, M} only, so without `marker_sets` AQI is `None`. | [checked] `eval_metrics.py:902-905` vs :848-858 | S | `eval_metrics.py` |
| 7 | **Parity tests against the original implementations** - none exist today. Commit small JSON goldens, generated offline, so CI stays pure-Python: ECS vs clusim (already matches to 2e-16), deconvolution vs OpenProblems (exact), ISCG/ISCL vs scTypeEval on the Yerly 2022 set (0.44 GB, CC BY, ISC Zenodo [18921437](https://doi.org/10.5281/zenodo.18921437)), signal QC vs SpatialQM on a 2 k-cell CC BY Xenium crop. scTypeEval is GPL-3: commit outputs, never code. | validation scan | M | `tests/goldens/`, `scripts/make_goldens.R` |
| 8 | **Close single-control test gaps** (v0.6.0 added 9 tests for shared markers, NMP spillover and `neighbor_contamination`; none of the gaps below): `neighborhood_sanity` positive control; `label_stability` negative; Moran shuffled-coordinates null; `ncp` test currently accepts `status == "skipped"`; assert `avg_similarity`; catch `Exception` (not just ImportError) around the squidpy import - in `ps_pgueguen` squidpy raises a dask `NotImplementedError` and errors test collection. | validation scan (`tests/test_purity.py:111`, `tests/test_spatial.py:11`) | S | `tests/` |
| 9 | **Hierarchical partial credit**: a wrong sibling subtype gets the same 0.5 as an honest coarse call. GPTCelltype / CASSIA give partial credit only one ontology step away. Fix: credit by distance on a user DAG (LCA depth), plus "accuracy on non-leaf calls". | [checked] `eval_metrics.py:537-553`; GPTCelltype [10.1038/s41592-024-02235-4](https://doi.org/10.1038/s41592-024-02235-4); Genome Biol 2025 [10.1186/s13059-025-03781-6](https://doi.org/10.1186/s13059-025-03781-6) | S | `eval_metrics.py` |
| 10 | **Unverified attributions** left in the catalog: CRISP is credited to the CSO, but a CRISP paper is Rose et al. 2026 (Emory / Georgia Tech, [10.64898/2026.04.16.718947](https://doi.org/10.64898/2026.04.16.718947)) - check whether our `crisp_purity` is that metric, now that v0.6.0 judges each lineage pair on unshared markers only; Mitchel et al. is cited as *Nat Genet* 2025, DOI [10.1038/s41588-025-02497-4](https://doi.org/10.1038/s41588-025-02497-4) (year to confirm). | tools scan; not re-verified | S | catalog §9 |

## P0b - found by the synthetic sweep ([benchmarks/synthetic](../benchmarks/synthetic/README.md), 15 scenarios × 5 seeds, 18 edge cases)

| # | Item | Evidence | Effort |
|---|---|---|---|
| 8a | **Validate labels once at the entry point.** NaN labels crash 8 of 14 metric groups with a `TypeError` deep inside; `neighborhood_sanity` crashes on a single label, and an unused categorical level changes its result. | `results/edge_cases.md` | S |
| 8b | **Warn on non-count input.** With log-normalised `X` and no `counts` layer, ISC, fidelity, PMP, NMP, `neighbor_contamination`, confidence and signal QC silently return different numbers. | edge case `lognorm_X_no_counts` | S |
| 8c | **Disclose labels without a marker set** (the AQI's `no_dict_frac`, for every marker metric). Over-split and merged labels are skipped silently, so PMP / NMP / fidelity / `neighbor_contamination` look clean on both errors. | scenarios `over_split`, `merge` | S |
| 8d | **ISC regime calls fail 3 of 4 cases**: merge reads `consistent` (ISCG rises 0.936 → 0.978), 20 % label noise and a batch effect read `split`; only over-split → `merge` is right. ISCL stays ≥ 0.90 at 40 % noise. Do not ship regime calls before the scTypeEval parity test (#7); see #4. | `results/check_isc_regimes.csv` | M |
| 8e | **CRISP / MECR need a null.** `crisp_purity` is 0.81 on perfectly clean data (off-target 0.1 % of the marker rate, "≥ 1 count" rule), and all reference-free purity metrics flag genuine co-expression (CRISP 0.81 → 0.69 at 5 %). Report against a per-section background (shuffled-marker or leak-rate null), and say that NMP, which carries the reference's co-expression, is the one to trust there. | scenarios `coexpression`, baseline | M |
| 8f | **Spatial coherence / PAS see random errors only**: blind to a compact mislabelled patch, *rise* when types are merged (0.904 → 0.961), fall for correctly labelled infiltrating cells. Document as a randomness detector, not an error detector; PAS at the 20 % threshold stays ≤ 0.012 even for random half-splits (with #2). | scenarios `patch_error`, `merge`, `infiltration` | S |
| 8g | **SNR is set by the eps floor when controls are near zero** (30 at baseline, 5.6 after 0.02 ambient counts per gene): report the control count it rests on, or NaN below a minimum. | scenario `ambient` | S |
| 8h | **`panel_resolvability` tiers are lenient**: `frac_resolvable` = 1.0 at depth 50 with no private markers (macro-F1 0.72). | `results/check_panel_resolvability.csv` | S |

## P1 - release and adoption

| # | Item | Why / evidence | Effort |
|---|---|---|---|
| 11 | **Publish to PyPI** via trusted publishing on a `v*` tag (`pypa/gh-action-pypi-publish`, OIDC); tag v0.6.0. Name is free. Then bioconda (`noarch: python`). | [checked] PyPI 404; no tags; [docs.pypi.org/trusted-publishers](https://docs.pypi.org/trusted-publishers/) | S |
| 12 | **Docs site**: Sphinx + myst-nb on ReadTheDocs (cookiecutter-scverse stack); numpydoc docstrings (only `conformal_prediction_sets` has a Parameters block), ruff `D` with `convention="numpy"`. | [cookiecutter-scverse](https://github.com/scverse/cookiecutter-scverse) | M |
| 13 | **One tutorial on public data** (10x Xenium breast FFPE, Janesick 2023, [GSE243280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243280)): signal QC → reference-free battery → spatial coherence, executed in the docs build - doubles as the 1e5-cell scale test after #3. | packaging scan | M |
| 14 | **API consistency before 1.0** (breaking, so batch into 0.7): `pmp` returns None, `crisp_purity`/`mecr` return floats, everything else dicts; four different "skipped" shapes; `(pred, truth)` vs sklearn's `(y_true, y_pred)`; `k` vs `n_neighs`; layer handling hard-coded to `"counts"` in purity, `layer=` elsewhere, silent fallback to X on a typo; fixed unprefixed `obs` columns (`pmp`, `nmp`, `neighbor_contamination`, `sparsity`…) that can clobber user columns - add `key_added=` and store summaries in `adata.uns["sam"]`. `spatial_coherence` ignores `k` when a graph already exists but reports it. | packaging scan, `confidence.py:33`, `spatial.py:23` | M |
| 15 | **Dense-matrix hotspots**: `_dense(adata.X)` in `marker_program_fidelity` (:310), `marker_gene_overlap` (:348), pseudobulk (:161), `_embedding` dense PCA fallback (:56), `label_stability` copies `adata` per rep; `neighbor_contamination` and the MECR / CRISP helper densify cells × marker genes (fine at 300 markers, 2.4 GB float64 at 1e6 cells). `nmp` went sparse in v0.6.0. 1e6 × 5 k float64 = 40 GB. Slice marker columns sparsely; `indicator.T @ X` for pseudobulk; `TruncatedSVD`. | packaging scan | M |
| 16 | **Release hygiene**: CHANGELOG, Zenodo DOI per release, single-source the version (`importlib.metadata`), `py.typed`, `Repository`/`Issues` URLs, dependency floors (`numpy>=1.24`, `anndata>=0.10`, `scikit-learn>=1.3`, `scipy>=1.10`), one lowest-direct CI run, CONTRIBUTING + issue templates. Consider `requires-python>=3.11` (current numpy/scipy/sklearn need it; anndata 0.13 and squidpy 1.8 need 3.12). | packaging scan; [SPEC 0](https://scientific-python.org/specs/spec-0000/) | S |
| 17 | **SpatialData input**: thin `from_spatialdata(sdata, table_key)`; needed anyway for transcript-level metrics (#21). | [10.1038/s41592-024-02212-x](https://doi.org/10.1038/s41592-024-02212-x) | S-M |
| 18 | **scverse ecosystem listing** needs releases, a PyPI/conda install and API docs (#11-12); the rest is met. | [scverse/ecosystem-packages](https://github.com/scverse/ecosystem-packages) | - |

## P1 - metrics worth adding (ranked by evidence × fit to no-ground-truth spatial)

| # | Metric | Source | Effort |
|---|---|---|---|
| 19 | **Baseline for `neighbor_contamination`** (v0.6.0 ships the per-cell measure, not the null). It counts a neighbour's markers in a cell whether they came from the neighbour or from the cell itself, so it cannot separate spillover from real co-expression: [checked] synthetic two-block grid, B→A 0.037 when every A cell truly expresses a B marker vs 0.040 for real border spillover. Same-type cells *without* that neighbour separate them: 0.037 vs 0.034 (co-expression) and 0.040 vs 0.001 (spillover). Fix: per (source, target) report the excess over non-adjacent same-type cells, plus a label-permutation null; then the SPLIT idea (secondary signal vs the source's share of neighbours) and CellSPA's negative-marker-vs-distance curve. | SegTraQ `contamination_strength`; SPLIT, Bilous *Nat Methods* 2026 [10.1038/s41592-026-03089-8](https://doi.org/10.1038/s41592-026-03089-8); [CellSPA](https://sydneybiox.github.io/CellSPA/articles/CellSPA.html) | S |
| 20 | **Conformal novelty p-values** in reference feature space (isolation forest / PCA reconstruction error, BH-adjusted) - set size 0 from softmax is a weak OOD signal. | López-De-Castro *Bioinformatics* 2025 [10.1093/bioinformatics/btaf521](https://doi.org/10.1093/bioinformatics/btaf521); Bates *Ann Stat* 2023 [10.1214/22-AOS2244](https://doi.org/10.1214/22-AOS2244); scDiagnostics | S-M |
| 21 | **Decomposition purity** = 1 − max RCTD weight per cell, the segmentation-error readout Spatial Touchstone validated; and complete the SpatialQM port (FDR from blank codewords, TPC, FTC). | Plummer *Nat Biotechnol* 2025 [10.1038/s41587-025-02811-9](https://doi.org/10.1038/s41587-025-02811-9) | S |
| 22 | **Multi-annotator consensus entropy with ontology-path voting** (votes count for ancestors) - an agreement signal that does not penalise subtype splits. | popV *Nat Genet* 2024 [10.1038/s41588-024-01993-3](https://doi.org/10.1038/s41588-024-01993-3) | S |
| 23 | **scib `isolated_labels` / cLISI** - `internal_validity` has nothing aimed at rare labels. Wrap scib-metrics or port. | [scib-metrics](https://scib-metrics.readthedocs.io/) | S |
| 24 | **Calibration check (ECE / class-wise ECE)** before reading native confidence as a probability; **CovGap** + clustered-Mondrian conformal for rare classes. | Guo 2017 [arXiv 1706.04599](https://arxiv.org/abs/1706.04599); Ding 2023 [arXiv 2306.09335](https://arxiv.org/abs/2306.09335) | S |
| 25 | **Ontology label-relation consistency**: do per-label centroid distances follow Cell Ontology distances? Catches a "macrophage" next to epithelium with no truth. | Genome Biol 2025 [10.1186/s13059-025-03781-6](https://doi.org/10.1186/s13059-025-03781-6) | M |
| 26 | **Deconvolution**: presence/absence AUPR, per-type SSIM, rank-aggregated score. | Spotless *eLife* 2024 [10.7554/elife.88431](https://doi.org/10.7554/elife.88431); Li *Nat Methods* 2022 [10.1038/s41592-022-01480-9](https://doi.org/10.1038/s41592-022-01480-9) | S |
| 27 | **spARI** (label sets vs each other, spatially tolerant). | *Biometrics* 2025 [10.1093/biomtc/ujaf127](https://doi.org/10.1093/biomtc/ujaf127) | S-M |

## P2 - later

- Validation harness T1-T9: the synthetic half now exists (`benchmarks/synthetic`); still to do on real
  sections (label noise, over-partitioning, merge, type swap, contamination injection,
  segmentation dilation, transcript dropout, coordinate shuffle, nuisance) scored by monotonicity,
  minimum detectable perturbation and constant-fit - the ISC benchmark's design (catalog §7). Feeds #5.
- Threshold provenance: one dict of every constant (ISC 0.74/0.45, panel F1 0.5/0.3, coverage 0.03,
  `min_type_n` 50, …) with source and `provisional: true`, echoed in results; report raw,
  null-anchored and "noise-equivalent" values instead of universal pass/fail.
- Transcript-level admixture (SegTraQ border admixture / z-layer overlap, cellAdmix, ovrlpy signal integrity) as
  per-cell covariates: take them from SegTraQ on the SpatialData (#17) rather than re-implementing.
- Contribute PMP / MECR / CRISP to OpenProblems `task_ist_preprocessing`; `task_label_projection` has
  no calibration metric - an open slot for #20/#24.
- Publication: JOSS needs 6 months of public history (repo created 2026-07-08 → not before ~Jan 2027)
  plus docs, CONTRIBUTING and a comparison to scTypeEval / scib-metrics / SpatialQM / OpenProblems.
  A *Bioinformatics* application note would want #5 as its figure.

## Rules to add to catalog §8

- **Internal metrics on a label-trained embedding are circular** (scANVI, fine-tuned foundation
  models): compute silhouette / kNN purity / ISC on a label-agnostic embedding. ("Metric Mirages",
  [10.1101/2024.04.02.587824](https://doi.org/10.1101/2024.04.02.587824), preprint.)
- **Internal separability mostly measures dataset difficulty**: dataset identity explained ~84 % of
  annotation-κ variance vs ~4 % for the tool ([10.64898/2026.08.28.747622](https://doi.org/10.64898/2026.08.28.747622),
  preprint) - so internal metric levels do not transfer across sections.
- ISC 0.74 / 0.45 still rests on one preprint; no independent 2024-26 replication found.

## Done on branch `docs/quickfixes` (rebased onto `worktree-feat+segtraq-purity` / `7008d41`)

- README: GitHub install (not on PyPI), `conformal_prediction_sets` row.
- Catalog: 22 dead relative links repointed to the public `spatial-scribe` repo; master-table code
  locations for every shipped metric; banner rewritten.
- Attributions: SpatialQM *is* the Plummer et al. 2025 Spatial Touchstone tool (the "name collision"
  warning was wrong); SPLIT is Bilous et al. 2026. Both checked on PubMed.
- `requires-python >= 3.10` (the `[spatial]` extra does not import on 3.9; core passes 50/50 there).
- `CITATION.cff` (validated), version 0.6.0.
- Catalog §5 master table: `neighbor_contamination` row; README and catalog §3e say it cannot separate
  spillover from co-expression (the #19 probe).
