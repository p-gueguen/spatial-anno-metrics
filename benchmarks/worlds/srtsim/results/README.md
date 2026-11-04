# SRTsim world - the synthetic sweep on realistic spatial counts

The toy world in `benchmarks/synthetic/sim.py` has stripes, Poisson counts and markers 1,000× their off-target
rate. This world keeps the same contract (`world.simulate` / `world.reference`, the same six type names) but its
counts come from **SRTsim** (Zhu, Shang & Zhou, *Genome Biol* 2023,
[10.1186/s13059-023-02879-z](https://doi.org/10.1186/s13059-023-02879-z)) fitted to a real **10x Xenium FFPE human
breast** section (Janesick et al., *Nat Commun* 2023, Rep1, CC BY 4.0). All 15 sweep scenarios run on it unchanged.

```bash
python benchmarks/worlds/srtsim/prepare.py crops --seeds 5        # 1. crops + SRTsim inputs          6 s
Rscript benchmarks/worlds/srtsim/generate.R <run> 5 6             # 2. SRTsim fit + count (R)          95 s on 6 cores
python benchmarks/worlds/srtsim/prepare.py cache --seeds 5        # 3. seed_<k>.h5ad + marker sets     4 s
python benchmarks/worlds/srtsim/run.py --seeds 3 --workers 16     # the sweep on this world           101 s on 16 cores
python benchmarks/worlds/srtsim/realism.py                        # realism.csv / realism.png
python benchmarks/worlds/srtsim/compare.py --toy3 <toy rerun>/verdicts.csv   # verdict_changes.csv
```

The cache (5 seeds, 8.9 MB of h5ad) lives in `SAM_SRTSIM_CACHE` (default
`/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/srtsim/`), not in git.

## How it was fitted

**Data.** The section (167,780 cells, 313 genes, 28 negative-control probes) with the authors' supervised
annotation, `GSM7780153_Xenium_R1_Fig1-5_supervised.csv.gz` from GEO GSE243280 (md5 `197fda5b...`; barcodes match
the section's `cell_id` order exactly). It is used as the truth, mapped onto the six sim types:

| sim type | author clusters |
|---|---|
| T_CD4 | CD4+_T_Cells |
| T_CD8 | CD8+_T_Cells |
| B | B_Cells |
| Myeloid | Macrophages_1, Macrophages_2, IRF7+_DCs, LAMP3+_DCs |
| Fibro | Stromal |
| Epi | Invasive_Tumor, Prolif_Invasive_Tumor, DCIS_1, DCIS_2 |

The other clusters are left out: inside the crops that is 18.9 % of cells (myoepithelial 7.9 %, endothelial
5.0 %, unlabeled 4.4 %, hybrids, perivascular, mast). Their positions stay empty.

**Crops = samples.** A crop is a 400 × 400 µm window whose lower-left corner sits on a 25 µm grid, holds
at least 600 cells of the six types, and has every type at ≥ 4 % and ≥ 12 cells (4,606 windows qualify). Each
world seed draws 6 non-overlapping windows at random (seeded), one per sample. That gives 5 seeds × 6 crops,
27,605 cells in all ([`crops.csv`](crops.csv)). Crops are stacked in y with a 100 µm gap, so a sample is a band,
as in the toy world. The rule selects immune-rich windows: over all crops the mix is Fibro 33.5 %, Epi 17.7 %,
T_CD4 15.4 %, Myeloid 14.8 %, T_CD8 10.9 %, B 7.8 %. `simulate()` keeps the 400 cells nearest each crop centre
(2,400 cells; Fibro 37 %, T_CD4 16 %, Myeloid 15 %, Epi 13 %, T_CD8 10 %, B 8 %, mean of seeds 0-2).

**SRTsim.** Version 0.99.7 from the CRAN archive (removed from CRAN 2025-12-17). Its NAMESPACE imports Morpho,
which pulls in rgl, and rgl needs libGLU, which the compute nodes lack. So `generate.R` sources the package's own
`SRTsim_class.R`, `createSRT.R`, `SRTmodels.R`, `SRTfit.R` and `SRTcount.R`, unmodified, from the md5-checked
tarball. The reference-based mode is run per crop:

- `srtsim_fit(marginal = "auto_choose", sim_scheme = "domain", min_nonzero_num = 2, maxiter = 500)`. Domain = cell
  type. For each feature and type, SRTsim picks Poisson when there is no overdispersion, else NB, ZINB (LRT,
  p < 0.05) or ZIP (AIC). Over 5 seeds × 6 crops × 6 types × 341 features: NB 31,213, Poisson 8,149, ZIP 960,
  ZINB 744. 20,314 slots (33 %) are not fitted because they have fewer than 2 non-zero cells, and SRTsim
  simulates those as 0.
- `srtsim_count(breaktie = "random")`, depth ratio 1. For each feature and type it draws n values from the fitted
  marginal and gives the k-th smallest to the cell with the k-th smallest real count, breaking ties at random.
  Seeded with `set.seed(1000 * seed + crop)`.
- `srtsim_fit(sim_scheme = "tissue")`: parameters only. They feed the uninformative `xg*` genes of `panel_growth`.

**Reference** (`world.reference`). This is `srtsim_count_single` applied to resampled cells. For each type,
`n_per_type` real cells of that type are drawn with replacement from the six crops. Each gene is drawn from the
crop's fitted marginal with the mean × `rr` and placed by rank of the drawn cells' real counts. `rr` = 400 / 150,
the toy's reference-to-section depth ratio. There are no controls and no coordinates. Check on seed 0: per-type
gene means are 0.985-1.011 × (section mean × `rr`), with log-mean correlation 0.986-0.997.

**Markers** ([`marker_sets.json`](marker_sets.json), [`marker_table.csv`](marker_table.csv)). Four per type, as in
the toy world, derived from the simulated counts of all 5 seeds. Counts are depth-normalised; specificity is
log2((mean in type + 0.01) / (max mean in another type + 0.01)), and a gene must be detected in ≥ 20 % of the
type's cells. T_CD4 and T_CD8 share the two genes most specific to T cells as a whole (IL7R, CD3E) and add their
own two:

| type | markers | fold over the next type (real data) |
|---|---|---|
| T_CD4 | IL7R, CD3E, CTLA4, KLRB1 | 0.8-1.6 |
| T_CD8 | IL7R, CD3E, KLRD1, CD8A | 0.8-6.0 |
| B | MS4A1, CD19, BANK1, CD79A | 8.2-9.8 |
| Myeloid | LYZ, FCER1A, HAVCR2, CD163 | 3.4-4.0 |
| Fibro | GJB2, PDGFRB, CRISPLD2, LRRC15 | 10.9-18.1 |
| Epi | ELF3, TACSTD2, MLPH, KRT8 | 27-36 |

The fold for a shared T marker is against the other T subtype. **Controls:** the 28 negative-control probes
(`feature_type == "Negative Control Probe"`, 8 of them named `antisense_*`) are renamed `NegCtrl0..27`, with the
10x id kept in `var['feature_id']`, and flagged in `var['control']`. Codewords and blanks are dropped.

**Scenario knobs** (see `world.py`):
- `weights` is read relative to `sim.WEIGHTS`: T_CD8 cells are dropped at random by the toy's factors, so at
  `rare_type` 0.01, T_CD8 is 0.55-0.75 % of 2,138-2,195 cells.
- `infiltrate` moves T_CD8 cells to uniform positions inside their own crop.
- `extra_bg` adds genes drawn i.i.d. per cell from a random real gene's tissue-level marginal.
- `coexpr` gives Myeloid cells Poisson(s × Fibro's mean × the cell's relative depth) extra counts on Fibro's
  markers, in the section and in the reference. At s = 0.4, Myeloid's mean Fibro-marker count goes from 0.11 to
  0.63.
- `depth` is not honoured. No scenario sets it.

`run.py` adapts two sweep helpers that hard-code the toy:
- edge case `all_zero_marker_gene` zeroes B's first marker instead of `B_m0`;
- `deconv_check` rescales coordinates to the toy density, so a 5 × 5 spot holds about 25 cells.

## Realism: SRTsim vs the real counts of the same cells

All 27,605 cached cells, genes only ([`realism.csv`](realism.csv), [`realism_markers.csv`](realism_markers.csv),
[`realism.png`](realism.png)). The toy column is 3 seeds of `sim.simulate` with its own markers.

| check | real | SRTsim | toy |
|---|---|---|---|
| library size: median / CV | 116 / 0.73 | 116 / 0.72 | 149 / 0.31 |
| library size: 5th-95th percentile | 35-368 | 37-361 | 90-246 |
| zero fraction (cells × genes) | 0.818 | 0.820 | 0.480 |
| mean-variance: slope of log var on log mean | 1.23 | 1.19 | 2.07 |
| Moran's I (6-NN within crop), median over genes / markers | 0.114 / 0.187 | 0.110 / 0.182 | 0.461 / 0.755 |
| marker specificity, median own / max other type | 6.5× | 6.4× | 874× |
| marker-positive cells with markers of ≥ 2 lineages | 82.7 % | 83.1 % | 44.0 % |
| gene-gene \|r\| within Epi, top-50 genes (median) | 0.107 | 0.095 | 0.050 |
| negative-control probe counts per 1,000 cells | 108 | 68 | 487 |

Agreement per gene and per cell, real vs SRTsim:
- per-gene log mean r = 0.991, log variance r = 0.992 (within type: mean r = 0.983);
- zero fraction r = 0.9999, mean |difference| = 0.003;
- per-cell library size Spearman ρ = 0.997, KS statistic 0.007;
- Moran's I per gene r = 0.998 (markers 0.999), median ratio SRTsim / real 0.99;
- gene-gene correlation matrices within Epi r = 0.996.

**What that means.**
- **Close to the real section.** SRTsim's reference mode re-draws every value from a smooth fitted marginal but
  keeps each cell's rank within its type. The section is therefore the real crop with re-sampled noise. Depth,
  sparsity, mean-variance, spatial autocorrelation, marker specificity and much of the gene-gene structure carry
  over, through the shared ranks.
- **Far from the toy world**, on every axis: 4× lower spatial autocorrelation, 136× lower marker specificity, nearly
  twice as many marker-positive cells mixing lineages (83 % vs 44 %; segmentation spillover is in the real counts
  and so in the fit), and a heavy-tailed library size.
- **Two departures.** Control probes lose 37 % of their counts: a probe with < 2 non-zero cells in a type and crop
  is not fitted and becomes 0. And the truth is the authors' annotation, so wherever that annotation is wrong, so
  is "truth".

## Runtime

| step | time |
|---|---|
| `prepare.py crops` | 6 s (rerun: byte-identical inputs) |
| `generate.R` | 95 s wall on 6 cores, 437 CPU-s: domain fit 188 s, count 156 s, tissue fit 93 s; 2.5 GB peak (SLURM MaxRSS) ([`fit_timing.csv`](fit_timing.csv)) |
| `prepare.py cache` | 4 s (rerun: identical marker sets) |
| `run.py` (15 scenarios × 3 seeds, 18 edge cases, 4 checks) | 101 s on 16 cores (the toy sweep at 3 seeds: 84 s) |

## Verdicts vs the toy world

[`matrix.md`](matrix.md), [`verdicts.csv`](verdicts.csv), [`values.csv`](values.csv),
[`verdict_changes.csv`](verdict_changes.csv), [`verdict_summary.md`](verdict_summary.md). The toy is
`benchmarks/synthetic/results` (5 seeds). The toy was also rerun at 3 seeds on this commit, so that changes caused
only by the seed count can be told apart.

| verdict (cells with a claim) | toy, 5 seeds | SRTsim, 3 seeds | toy, 3 seeds |
|---|---|---|---|
| ✅ detects | 134 | 112 | 134 |
| ❌ blind | 61 | 74 | 58 |
| ⚠️ false alarm | 22 | 16 | 21 |
| ⛔ wrong way | 3 | 6 | 3 |
| 〰️ erratic | 35 | 41 | 38 |
| ✔ stable | 280 | 282 | 276 |

133 of 765 (scenario, metric) cells change verdict: 23 move to detects or stable, 43 move away from it, and 67 are
other changes. Only 3 of the 133 also differ between the toy at 5 and at 3 seeds.

**Read this first: the move rule is much stricter here.** A move must exceed 3 × the seed SD at strength 0. In
this world each seed is a different set of tissue crops, so that SD is between-section variability, not
draw-to-draw noise. The SD sets the threshold in 71 % of cells, against 22 % in the toy, and the median 3·SD is
8.9 % of the metric's range, against 0.6 %. Several "lost" detections are real moves just under the bar. The fixed-crop control
(`sweep.py --world srtsim --fixed-crops 0`, [`fixed_crops_0/`](fixed_crops_0/)) keeps one crop set and varies only
the SRTsim draw: 9 of the 30 lost toy detections come back, and 21 (all the sibling-confusion ones but the
isolated-label scores) stay lost. See [`../../README.md`](../../README.md#srtsim-control-fixed-crops).

**Null scenarios.**
- **The toy's false alarms under `coexpression` disappear**, because the real baseline already mixes lineages.
  - On clean data, `crisp_purity` is 0.196 (toy 0.813), `mecr` 0.284 and `pmp` 0.628.
  - At 40 % co-expression: CRISP 0.196 → 0.182 (toy 0.813 → 0.685), MECR 0.284 → 0.298 (toy 0.042 → 0.065),
    PMP 0.628 → 0.617 (toy 0.995 → 0.950). These changes are within the noise.
  - Silhouette, isolated-label ASW and `onto.*` turn ✔ as well.
- **`rare_type`:** silhouette (0.153 → 0.153), `integrated`, Ward-PropMatch, the isolated-label scores and
  `coex.global` turn ✔. `iv.orbital_medoid` becomes a new ⚠️ (0.818 → 0.809).
- **`infiltration`:** coherence and PAS turn ✔ (coherence 0.452 → 0.436; toy 0.920 → 0.734). Real tissue is
  already mixed (clean coherence 0.452, PAS 0.537), so scattering T_CD8 inside its crop adds little. CHAOS stays ⚠️.
- **`panel_growth` raises 8 new false alarms**:
  - ISCG 0.525 → 0.266, ISCL 0.893 → 0.582;
  - Ward-PropMatch 0.748 → 0.532, neighbourhood purity 0.795 → 0.706, Davies-Bouldin 2.23 → 3.16, orbital medoid;
  - cLISI 0.961 → 0.903.

  Here the added genes are truly uninformative and add noise to the PCA. The toy's `xg*` genes carry a mild
  per-type programme (log-sd 0.3), which is why the toy's silhouette-type scores *rose* with panel size.

**Label errors.**
- **`sibling_confusion` goes from ✅ to ❌ for 12 metrics**:
  - internal validity (silhouette 0.153 → 0.143);
  - the isolated-label scores;
  - PMP (0.628 → 0.614), NMP (0.978 → 0.956), `neighbor_contamination`;
  - `ext.hier_acc` (1 → 0.961, threshold 0.05).

  In this section CD4 and CD8 T cells hardly differ on the panel. Their own markers are 1.5-1.6× (CTLA4, KLRB1)
  and 4.6-6.0× (CD8A, KLRD1) over the next type. `panel_resolvability` names T_CD4 as T_CD8's confuser in 48/48
  runs, and T_CD8's F1 is 0.78 at depth 50 with every marker (toy: 1.00).
- **Lost to the stricter bar** (moved the right way; fraction of the threshold in brackets):
  - `over_split`: `iso.f1` 0.744 → 0.596 (0.94), Ward-PropMatch (0.84);
  - `patch_error`: `iso.f1` (0.64), Ward-PropMatch (0.76), PMP 0.628 → 0.571 (0.46);
  - `type_swap`: `mgo` 0.653 → 0.514 (0.96).
- **New detections:**
  - `ext.comp_l1` catches `label_noise` (0 → 0.254) and `type_swap` (0 → 0.481): the composition is unbalanced
    (Fibro 37 %), so even symmetric errors shift proportions.
  - cLISI (0.961 → 0.857) and ISCL (0.893 → 0.753) catch `patch_error`: neither starts at its ceiling of 1 here.
- **Merging, wrong way:** merging still makes labels look better, and now by more than the threshold:
  - neighbourhood purity 0.795 → 0.864 ⛔;
  - CHAOS 0.052 → 0.042 ⛔;
  - `ns.n_suspicious` 1 → 0 ⛔;
  - coherence 0.452 → 0.572 ⛔ (also in the toy).

  ISCG under `merge` is no longer ⛔: it rises 0.525 → 0.578, which is under the threshold, so it reads ❌.
- **`onto.rho`** has a real embedding to read here: 0.598 on clean labels (toy 0.255), as the toy README expected.
  It still misses label noise (0.598 → 0.660, ❌) and goes ⛔ on `patch_error` (0.598 → 0.790).

**Degraded data.**
- `sig.moran_gap` now catches `ambient` (0.092 → 0.053, ✅ 0.2), where the toy was ❌. SNR still does
  (3.83 → 1.63), and `frac_above_ctrl` stays ❌.
- ISCL catches `batch` (0.893 → 0.676); in the toy it sat at 1.000.
- Spillover:
  - `neighbor_contamination` (0.023 → 0.037), NMP (0.978 → 0.767), CRISP (0.196 → 0.092) and MECR
    (0.284 → 0.380) still detect.
  - PMP (0.628 → 0.533, 0.78 of the threshold) and `ncp.moldia` (0.986 → 0.977, 0.99) drop to ❌.
  - `coex.by_type` catches it at 0.4 (0.9673 → 0.9656, threshold 0.0012).

**Dedicated checks** ([`check_*.csv`](.)):
- **Conformal (α = 0.1).** Calibrated on a labelled 30 % of the section, coverage is still 0.90. Calibrated on the
  reference, sets *under*-cover: 0.748 at full depth and 0.581 at a quarter of it (toy: 0.990 and 0.807); 19 % of
  known cells get an empty set. For 2 % T_CD8, its coverage is 0.405 marginal and 0.510 Mondrian. The
  reference-trained classifier's ECE is 0.089 (toy 0.001).
- **Novelty, with B held out.** Conformal novelty ranks B above known cells (AUC 0.83 at full depth, 0.60 at a
  quarter; toy 0.999) but flags none of them at BH 10 %. In the main scenarios `nov.frac_novel` is 0 except under
  ambient noise (0.31 at 0.05, 0.94 at 0.2), and it no longer reacts to batch (0 at log-sd 1.5; toy 0 → 0.998).
- **ISC regimes.** On *correct* labels, 14 of 17 type-seeds are called `split` (ISCG 0.53 on average, ISCL 0.89);
  the toy has 18/18 `consistent`. Crops from different tissue regions are not consistent replicates, so ISC's
  regime call cannot be read as a label-quality verdict across regions.
- **Panel resolvability.** Dropping the 20 private markers barely matters with 313 genes: macro-F1 at depth 50 is
  0.83 with none kept and 0.89 with all (toy 0.71 and 1.00).
- **Deconvolution metrics** stay monotone: R² 1 → 0.41 at 50 % noise, AUPR 1 → 0.87, SSIM 1 → 0.64.

## Limits

- **A semi-synthetic copy.** Rank matching ties every simulated cell to a real one (library size ρ = 0.997). The
  spillover, segmentation errors and annotation errors of the real section are inherited, and "truth" is the
  authors' annotation, not a known generative label.
- **Seeds are different tissue.** That makes the thresholds stricter (above). The toy had 5 seeds; this sweep has 3.
- **Crops were selected to hold every type**, which favours immune-rich regions. Myoepithelial, endothelial,
  perivascular, mast, hybrid and unlabelled cells are left out, which leaves holes in the layout.
- **The simulation departs from the data in places.** Genes are simulated independently within a type, so gene-gene
  structure survives only through shared ranks (within Epi, median |r| 0.095 vs 0.107). Low-count features are set
  to 0, which costs the control probes 37 % of their counts.
- **Knobs on real counts are approximations**: weights drop cells, infiltration moves cells within the crop, and
  co-expression is added on top of the fitted counts.
- **The reference shares the section's spillover**, because it is drawn from the same fits. A real dissociated
  reference would be cleaner.
