# Roadmap - after 0.7.0

What is left after the 0.7.0 release, checked against the code on `release/0.7.0` (2026-11-03). Item numbers
are those of the 2026-11-03 field scan, so references elsewhere (catalog, CHANGELOG, commit messages) still
resolve. Each open item says why it is open and where the evidence is.

**What is ours alone** (no other tool found doing it): depth-thinned `panel_resolvability`,
ISC regime calls in Python, PMP, the lineage-coverage flag, annotator-agnostic `label_stability`.
**Where we re-implement**: ISC (scTypeEval), NMP / co-expression similarity (txsim / OpenProblems iST), signal
QC (SpatialQM), conformal (scConform), `neighbor_contamination` and decomposition purity
([SegTraQ](https://github.com/LazDaria/SegTraQ)). 0.7.0 gives each of these a value test against its source
(`tests/goldens/`); the ours-alone items still lack real-data validation (P2).

---

## Open

### Release (maintainer actions)

| # | Item | Why it is open | Evidence |
|---|---|---|---|
| 11 | **Publish to PyPI**, then bioconda (`noarch: python`); mint the Zenodo DOI (was #16) | Needs the maintainer's PyPI trusted-publisher setup (pending publisher: project `spatial-anno-metrics`, owner `p-gueguen`, workflow `publish.yml`, environment `pypi`) and a pushed `v0.7.0` tag. The Zenodo DOI needs the GitHub-Zenodo integration switched on before that release. | `.github/workflows/publish.yml` is in place and checks tag = pyproject version; no tag exists (`git tag -l` empty); `.zenodo.json` ready, `CITATION.cff` has no `doi:` yet |
| 18 | **scverse ecosystem listing** (PR to scverse/ecosystem-packages) | Needs the release and the PyPI install (#11); API docs (#12) are done. | - |

### Metrics and code

| # | Item | Why it is open | Evidence |
|---|---|---|---|
| 5 | **AQI validation and tuning** | Tracked outside this repo. | - |
| 8d | **ISC regime calls** are still wrong in 3 of 4 synthetic cases: merge reads `consistent`, 20 % label noise `split` (21/30 labels), a batch effect `split` (30/30); only over-split is right. | Flagged, not fixed: results carry `regime_status="experimental"`. The headline ISCG/ISCL now match scTypeEval (#7), so the calls, not the scores, are the problem; they need real multi-sample data with a known answer before they can be trusted. | `isc.py` (Notes, with the numbers); `tests/test_isc.py::test_regime_calls_flagged_experimental`; `benchmarks/synthetic/results/check_isc_regimes.csv` |

## P2 - later (external PRs & publication)

- **OpenProblems contribution**: PMP / MECR / CRISP to `task_ist_preprocessing`; `task_label_projection` has no
  calibration metric, an open slot for `calibration_error` / `conformal_novelty`.
- **Publication**: JOSS needs 6 months of public history (repo created 2026-07-08, so not before ~Jan 2027),
  plus a comparison to scTypeEval / scib-metrics / SpatialQM / OpenProblems.

---

## Done in 0.7.0

Details, with the `old -> new` call for every breaking change, are in [CHANGELOG.md](../CHANGELOG.md).

- **P0**: #1 NMP (Salas score by default, txsim golden; `score="moldia"` matched to Moldia `Xenium_benchmarking`
  reads golden; `score="cells"` matched to `txsim_cells` / `moldia_seg_cells` golden and its 20 permutations; raw
  score NaN for unscorable types; `coexpression_similarity`; `ncp(score="moldia")` matched to Moldia
  `negative_marker_purity_coexpression` golden), #2 PAS to Shang & Zhou 2022 + CHAOS,
  #3 `internal_validity` at 1e5+ cells (1.64e5 cells: 2,951 s / 211 GB -> 20 s / 0.74 GB), #4 ISCL classifier genes
  (`rcm_genes`) and `integrated` documented and golden-tested, #6 docstring, #7 goldens (clusim, OpenProblems +
  QuKunLab `GenesMetrics.py` SSIM + Spotless `precrec` Davis-Goadrich PR AUC, scTypeEval, SpatialQM, txsim + Moldia),
  #8 single-control test gaps (including the 25 % `weak` confusion cut in `panel_resolvability` and ECE bin-edge side
  in `calibration_error`), #9 LCA-depth hierarchical credit, #10 attributions (CRISP = Rose et al. 2026, Mitchel et
  al. *Nat Genet* 2026).
- **P0b**: 8a label validation, 8b non-count input + NaN/Inf rejection in `get_counts`, 8c dictionary coverage
  (`no_dict_frac`), 8e CRISP / MECR null, 8f spatial-metric docs, 8g SNR without the eps floor, 8h chance-corrected
  panel tiers, `thresholds.thresholds_for(module, **overrides)` helper.
- **P1**: #12 docs site (Sphinx + myst-nb, ReadTheDocs, numpydoc enforced by ruff `D`), #13 Xenium breast
  tutorial (167,780 cells, 86 s), #14 API consistency (including `skip` instead of `KeyError` when spatial metrics
  run without `obsm['spatial']` or `obsp`), #15 sparse paths across all count and feature metrics (including
  `coexpression_similarity`, `conformal_novelty`, `per_cell_confidence` and `label_stability`), #16 release
  hygiene (including `.gitattributes` locking golden and `.obo` line endings; except the Zenodo DOI, moved to #11),
  #17 `from_spatialdata` + `io.segtraq_covariates` passthrough/summary for SegTraQ transcript-level admixture
  columns, #19 contamination baseline + SPLIT (`split_slope`, `split_corr`) and CellSPA (`distance_corr`,
  `distance_curve`) readouts, #20 conformal novelty, #21 decomposition purity + SpatialQM TPC / FTC / FDR,
  #22 consensus entropy (tested on Yerly 2022 BCC annotations), #23 isolated labels + cLISI, #24 ECE / CovGap /
  clustered Mondrian + real-data perturbation harness (`benchmarks/real_perturbation.py` on Yerly 2022 BCC and Xenium
  breast crop) + 5-seed synthetic sweep (`benchmarks/synthetic/sweep.py`), #25 ontology consistency (`distance="path"`
  and `distance="rwr"` Random Walk with Restart after `scGraph-OntoRWR`, tested on Yerly 2022 BCC data),
  #26 deconvolution AUPR (`aupr_method="sklearn" | "precrec"` + `aupr_precrec`) / SSIM / rank score, #27 spARI.
- **Catalog §8**: the label-trained-embedding circularity rule (and `embedding_label_trained=True`), the
  dataset-difficulty variance rule (Wardhana, Zeng & Lu 2026), and the single-preprint caveat on the ISC 0.74 / 0.45
  guide thresholds.
