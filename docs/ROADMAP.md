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
| 8d | **ISC regime calls** are still wrong in 3 of 4 synthetic cases: merge reads `consistent`, 20 % label noise `split` (13/18 labels), a batch effect `split` (18/18); only over-split is right. | Flagged, not fixed: results carry `regime_status="experimental"`. The headline ISCG/ISCL now match scTypeEval (#7), so the calls, not the scores, are the problem; they need real multi-sample data with a known answer before they can be trusted. | `isc.py:231-237` (Notes, with the numbers); `tests/test_isc.py::test_regime_calls_flagged_experimental`; `benchmarks/synthetic/results/check_isc_regimes.csv` |
| 15 | **Dense-matrix hotspots** (remainder) | Three call sites still densify cells x genes: `label_stability` thins a dense float32 copy per rep (making it sparse changes the RNG stream and what `annotate_fn` receives); `coexpression_similarity` densifies cells x shared genes; `conformal_novelty` densifies both inputs. `_embedding` densifies only on scikit-learn < 1.4. | `confidence.py:94`; `purity.py:605-606` (`ponytail:` note); `conformal.py:315-316`; `_utils.py:47` |
| 19 | **`neighbor_contamination` follow-ups** | The baseline, excess and permutation p shipped. Left: the SPLIT readout (secondary signal vs the source's share of neighbours) and CellSPA's negative-marker-vs-distance curve. The p is optimistic on real data (Xenium breast: 75 of 100 source-target pairs at p <= 0.01), so the docs say to read the excess. | `tests/test_spatial.py::test_neighbor_contamination_baseline_separates_spillover_from_coexpression` (spillover excess 0.040, p 0.001; co-expression excess 0.004, p 0.19) |
| §8 | **Two catalog §8 rules not yet written** | (a) Internal separability mostly measures dataset difficulty: dataset identity explained ~84 % of annotation-kappa variance vs ~4 % for the tool ([10.64898/2026.08.28.747622](https://doi.org/10.64898/2026.08.28.747622), preprint; figures to re-check before quoting). (b) ISC 0.74 / 0.45 rests on one preprint with no independent 2024-26 replication: said in the `inter_sample_consistency` docstring, not in §8. | `docs/cell-annotation-quality-metrics.md` §8 (no hit for `747622` or "difficulty"); `isc.py:237` |

### Small gaps found while checking 0.7.0

- **Parity not checked against the original code**: deconvolution SSIM is tested against the paper's formula
  only (running Li et al.'s `GenesMetrics.py` was blocked); AUPR against sklearn average precision, not
  Spotless's `precrec`; NMP against txsim, not Moldia's older reads criterion.
- **Test mutants that survive**: dropping the 25 % `weak` confusion cut in `panel_resolvability` (needs a type
  with adjusted F1 >= 0.8 and confused fraction 0.25-0.30; none stable found in 72 configurations); the
  ECE bin-edge side in `calibration_error` (no value lands on an edge in the tests).
- **Untested on real data**: `consensus_entropy` (no multi-annotator section yet); `ontology_consistency` is
  an own statistic (`is_a` path, Spearman, Mantel), not compared with scGraph-OntoRWR.
- `tests/goldens/signalqc_spatialqm/cells.csv` is pinned by sha256 but has no `.gitattributes` `-text`
  entry, so a CRLF checkout would fail that golden.

## P2 - later (deferred)

- **Real-data perturbation harness** (T1-T9 on real sections: label noise, over-partitioning, merge, type swap,
  contamination injection, segmentation dilation, transcript dropout, coordinate shuffle, nuisance), scored
  by monotonicity, minimum detectable perturbation and constant-fit - the ISC benchmark's design (catalog §7).
  The synthetic half exists (`benchmarks/synthetic`, rerun for 0.7.0).
- **Threshold provenance, remainder**: `THRESHOLDS` (every constant with source and `provisional`) shipped;
  echoing the thresholds used into every result (only ISC does) and reporting noise-equivalent values are not
  done.
- **Transcript-level admixture** (SegTraQ border admixture / z-layer overlap, cellAdmix, ovrlpy) as per-cell
  covariates: take them from SegTraQ on the SpatialData table (`from_spatialdata`) rather than re-implement.
- **OpenProblems contribution**: PMP / MECR / CRISP to `task_ist_preprocessing`; `task_label_projection` has no
  calibration metric, an open slot for `calibration_error` / `conformal_novelty`.
- **Publication**: JOSS needs 6 months of public history (repo created 2026-07-08, so not before ~Jan 2027),
  plus a comparison to scTypeEval / scib-metrics / SpatialQM / OpenProblems.

---

## Done in 0.7.0

Details, with the `old -> new` call for every breaking change, are in [CHANGELOG.md](../CHANGELOG.md).

- **P0**: #1 NMP (Salas score by default, txsim golden; raw score NaN for unscorable types;
  `coexpression_similarity`), #2 PAS to Shang & Zhou 2022 + CHAOS, #3 `internal_validity` at 1e5+ cells
  (1.64e5 cells: 2,951 s / 211 GB -> 20 s / 0.74 GB), #4 ISCL classifier genes (`rcm_genes`) and `integrated`
  documented and golden-tested, #6 docstring, #7 goldens (clusim, OpenProblems, scTypeEval, SpatialQM, txsim),
  #8 single-control test gaps, #9 LCA-depth hierarchical credit, #10 attributions (CRISP = Rose et al. 2026,
  Mitchel et al. *Nat Genet* 2026).
- **P0b**: 8a label validation, 8b non-count input, 8c dictionary coverage (`no_dict_frac`), 8e CRISP / MECR
  null, 8f spatial-metric docs, 8g SNR without the eps floor, 8h chance-corrected panel tiers.
- **P1**: #12 docs site (Sphinx + myst-nb, ReadTheDocs, numpydoc enforced by ruff `D`), #13 Xenium breast
  tutorial (167,780 cells, 86 s), #14 API consistency, #15 sparse paths (except the remainder above), #16 release
  hygiene (except the Zenodo DOI, moved to #11), #17 `from_spatialdata`, #19 contamination baseline, #20
  conformal novelty, #21 decomposition purity + SpatialQM TPC / FTC / FDR, #22 consensus entropy, #23 isolated
  labels + cLISI, #24 ECE / CovGap / clustered Mondrian, #25 ontology consistency, #26 deconvolution AUPR / SSIM
  / rank score, #27 spARI.
- **Catalog §8**: the label-trained-embedding circularity rule (and `embedding_label_trained=True`).
