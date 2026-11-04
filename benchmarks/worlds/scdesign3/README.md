# scDesign3 world - the synthetic sweep on realistic counts

`benchmarks/synthetic/sim.py` is a toy world: six types in stripes, Poisson counts, markers 1,000× their
off-target rate. This folder builds a world with the same contract from **scDesign3** (Song et al., Nat
Biotechnol 2024) fitted to a real 10x Xenium section, and runs the unchanged sweep (`sim.SCENARIOS`, every
metric, every check) on it.

```bash
python build.py prepare                         # crop + author labels -> cache/input
sbatch generate.sbatch                          # the scDesign3 fit (generate.R), 4-9 min on 32 cores
python build.py pack                            # cache/fit -> cache/params.h5ad
python build.py realism                         # results/realism.*
python run.py --world bands --seeds 3 --workers 16                       # all 15 scenarios + checks
python run.py --world dirichlet --only merge over_split label_noise batch # multi-sample variant
python compare.py --toy <toy sweep dir, 3 seeds>                         # results/verdict_*.csv, isc_regimes.csv
```

The cache (`SAM_SCDESIGN3_CACHE`, default `/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/scdesign3/`)
is not in git. `world.py` loads `params.h5ad` once per process. `run.py` swaps the world into the `sim`
module before the sweep forks its workers; nothing in `benchmarks/synthetic/` is edited.

## How it was fitted

**Data.** The 10x Xenium FFPE human breast cancer section, replicate 1 (Janesick et al. 2023, CC BY; 167,780
cells, 313 genes). The cell types are the **authors' supervised annotation** (GEO GSM7780153,
`Xenium_R1_Fig1-5_supervised.csv.gz`, one row per cell, matched on `cell_id`), not the tutorial's.

**Crop.** x < 1,300 µm, y ≥ 4,000 µm: invasive tumour, stroma and the immune aggregate in the bottom-left
corner. Six author types map onto the toy names, and cells of the other 14 types (DCIS, myoepithelial,
endothelial, hybrids, unlabelled and so on, 37 % of the crop) are left out, so they leave holes in the
tissue:

| toy name | author type | cells |
|---|---|---|
| T_CD4 | CD4+_T_Cells | 1,257 |
| T_CD8 | CD8+_T_Cells | 915 |
| B | B_Cells | 643 |
| Myeloid | Macrophages_1 | 1,033 |
| Fibro | Stromal | 2,822 |
| Epi | Invasive_Tumor | 2,044 |

That is 8,714 cells. The panel is the 304 genes detected in ≥ 0.5 % of crop cells (9 dropped: AHSP, CRHBP,
LEP, MPO, MUC6, PDCD1, SCGB2A1, TCF15, UCP1), plus 28 negative-control probes.

**Model** (`generate.R`, scDesign3 1.7.6, mgcv 1.9.4, R 4.6.0). Each gene gets a negative binomial GAM:

    log mu = cell_type + log(cell_area) + s(x, y, bs = "gp", k = 100),    sigma ~ 1

That is a type effect, a cell-size effect, and a per-gene Gaussian-process smooth over the real centroids.
Gene-gene dependence comes from a Gaussian copula per cell type (`corr_by = cell_type`, distributional
transform). All 304 fits converged. Per gene: theta median 0.99 (10-90 %: 0.12-3.96); the spatial smooth
uses median 30 effective degrees of freedom (max 82 of 99, so k = 100 does not saturate); the log-area
coefficient is 0.71 (0.45-0.97); median deviance explained is 0.28. NB, not ZINB: the zero fraction per gene
already matches (below).

**Controls.** The 28 negative-control probes are not fitted. They are drawn as flat Poisson at each probe's
crop mean (0.006 counts per probe per cell on average, 0.0007-0.088).

**Markers** (`uns['marker_sets']`, derived from the fitted type means by `build.derive_markers`). A gene's
specificity for a type is log2 of its mean in that type over its maximum in any other type. Each type takes
its 4 most specific genes among those with a mean ≥ 0.1. As in the toy world, T_CD4 and T_CD8 share 2
T-lineage genes, scored against the non-T types.

| type | markers | log2 specificity |
|---|---|---|
| T lineage | IL7R, CD3E | 1.6, 1.6 |
| T_CD4 | + KLRB1, CCR7 | 0.67, 0.67 |
| T_CD8 | + KLRD1, KLRC1 | 2.8, 2.6 |
| B | MZB1, CD79A, CD19, MS4A1 | 3.3-2.9 |
| Myeloid | ITGAX, LYZ, HAVCR2, CLEC9A | 2.7-2.0 |
| Fibro | PDGFRB, PDGFRA, LRRC15, POSTN | 3.1-2.4 |
| Epi | MLPH, ELF3, EPCAM, ABCC11 | 3.7-3.4 |

In real Xenium data CD4 T cells have no strongly private gene on this panel: their best two are 1.6× their
level in any other type.

**Sampling.** `world.py` samples the fit the way `scDesign3::simu_new` does: z ~ N(0, R_type), then u = Φ(z),
then the NB quantile at the cell's fitted mean and the gene's dispersion (gamlss `qNBI`). That gives any seed
without R. The R draw and the port agree on every realism statistic (table below), including gene-gene
dependence.

**Two scDesign3 pitfalls, both handled in `generate.R`:**
- **`simu_new` gives each cell type another type's copula** when `corr_group` is a factor. It looks the
  copula up as `copula_list[[x]]` with `x` the factor *code*, while the list is in first-appearance order.
  On this fit Epi cells drew T_CD4's copula, T_CD4 drew Epi's, B drew Fibro's, and only Myeloid matched. The
  draw kept 4 % of the fitted dependence (z-correlation regressed on the fitted correlations: slope 0.04,
  against 0.64 for the port). Passing `corr_group` as character fixes it, and the fixed R draw matches the
  port. `fit_copula` labels its list correctly.
- **The copula fit was not reproducible.** Its distributional transform jitters every count inside forked
  `mcmapply` workers. Under Mersenne-Twister, two identical fits gave copula entries that agreed at r = 0.19
  (T_CD4) to 0.82 (Fibro); the median entry differed by 0.015-0.034, against a typical |correlation| of
  ~0.05. With `RNGkind("L'Ecuyer-CMRG")`, two fits on different nodes agree to 1.9e-12 (copulas) and
  1.7e-10 (means), and the R draw is identical. The copula is still noisy for the small types: that noise is
  estimation noise, now fixed per seed.

**Runtime.** The fit takes 4 min 16 s on 32 cores of fgcz-c-052 and ~9 min on fgcz-c-041. Of the c-052
run, `fit_marginal` is 147 s, `fit_copula` 50 s and `extract_para` 43 s. Peak memory is 91-96 GB summed over
the forked workers, so the job asks for 128 GB. In Python, a 2,400-cell section takes ~0.9 s on first call
(0.25 s to load the cache) and a 1,800-cell reference 0.45 s. The sweep on 16 cores takes 5 min 03 s for
all 15 scenarios × 3 seeds plus the checks and edge cases, against 1 min 22 s for the toy world. The
Dirichlet run (4 scenarios plus the checks) takes 4 min 28 s.

## Contract and kwargs

`simulate()` and `reference()` take the same arguments as `sim.py`, and every kwarg is honoured. Unknown
kwargs raise.

- **Section.** A random subsample of the 8,714-cell pool (`n_cells`, default 2,400) at its real centroids,
  rescaled to one cell per unit², the toy world's density. Cell areas are real. Natural composition: Fibro
  .32, Epi .23, T_CD4 .14, Myeloid .12, T_CD8 .105, B .074.
- **`weights`**: stratified subsampling. `None` is the natural composition, not `sim.WEIGHTS`. `run.py` sets
  `sim.WEIGHTS` to the natural shares, so `rare_type` rescales the other types around them. **`rare_type`
  starts at the crop's own T_CD8 share**: its levels are 0.105, 0.08, 0.04, 0.02, 0.01 instead of 0.15 ...
  0.01.
- **`n_samples`**: horizontal bands, as in the toy world. The tissue is not striped, so the bands differ in
  composition: Epi is 48 % of band s0 and 3 % of s5 (seed 0).
- **`sample_mode="dirichlet"`**: the multi-sample variant (below).
- **`infiltrate`**: that share of T_CD8 cells moves to uniform random positions. Their counts keep the
  location effect of their origin.
- **`extra_bg`**: genes `xg*` with levels drawn from the fitted genes' means, a per-type programme
  lognormal(0, 0.3), and dispersions drawn from the fitted genes. They come from their own stream, so the
  other genes' counts do not change (checked: 0 differing entries).
- **`coexpr`**: Myeloid cells get `coexpr` × Fibro's fitted mean added on Fibro's 4 markers, in both the
  section and the reference. The baseline is not zero: Myeloid already averages 0.70 counts per Fibro-marker
  gene against Fibro's 4.16, from the real data.
- **`depth`**: `None` keeps the real depth (median 128 counts per cell); a number rescales every mean.
- **Reference**: `n_per_type` pool cells per type drawn with replacement, means × 400/150 (as the toy
  world's reference), fresh copula draws, no controls.
- **Scenarios supported: all 15** of `sim.SCENARIOS`. `run.py` also rewires the two edge cases that name toy
  genes (`all_zero_marker_gene`: B's first marker; `extra_marker_type`: two control probes).

## Realism ([`results/realism.csv`](results/realism.csv), [`realism_genes.csv`](results/realism_genes.csv), [`realism.png`](results/realism.png))

The real crop, the R `simu_new` draw and the `world.py` draw (seed 0) use the same 8,714 cells and 304
genes. The toy world is `sim.simulate(0)`: 2,400 cells, its own 82 non-control genes and markers.

| | real crop | scDesign3, R | scDesign3, world.py | toy world |
|---|---|---|---|---|
| library size, median | 129 | 127 | 128 | 149 |
| library size, IQR / median | 1.09 | 1.10 | 1.09 | 0.42 |
| library size, CV | 0.74 | 0.76 | 0.76 | 0.32 |
| library size, KS vs real | - | 0.017 | 0.015 | - |
| zero fraction | 0.798 | 0.794 | 0.794 | 0.495 |
| per-gene zero fraction, r vs real (mean abs diff) | - | 0.999 (0.005) | 0.999 (0.005) | - |
| mean-variance slope (log10) | 1.23 | 1.27 | 1.27 | 1.46 |
| median variance / mean | 1.89 | 1.94 | 1.91 | 1.34 |
| per-gene variance, r of log vs real (median ratio) | - | 0.989 (1.04) | 0.989 (1.02) | - |
| marker mean, own type / max other (median) | 6.5× | 6.9× | 6.3× | 620× |
| marker score AUC, type vs rest (mean; min) | 0.94; 0.90 | 0.91; 0.86 | 0.91; 0.86 | 1.00; 1.00 |
| gene-gene correlation inside Epi, r vs real | - | 0.74 | 0.71 | - |
| same, inside Fibro | - | 0.88 | 0.88 | - |
| Moran's I per gene, median | 0.102 | 0.058 | 0.056 | 0.065 |
| Moran's I per gene, r vs real | - | 0.974 | 0.975 | - |

- Depth, sparsity, overdispersion and marker specificity match the real section. The toy world's markers
  are ~100× too specific (620× vs 6.5×), it has fewer zeros (0.50 vs 0.80), and its library sizes vary
  less than half as much (CV 0.32 vs 0.74).
- **Gene-gene dependence is weaker than real.** Mean |correlation| is 0.044 inside Epi against 0.069 real
  (Fibro 0.051 vs 0.063). Its pattern follows the real one (r 0.71 / 0.88).
- **Spatial autocorrelation is half the real one** at the median (Moran's I 0.056 vs 0.102), although it
  ranks the genes the same way (r 0.975). The smooth captures regional trends; the cell-to-cell part
  (segmentation bleed from neighbours, local niches) is not modelled, so simulated spillover comes only from
  the `spillover` scenario.

## Results: which verdicts change ([`results/verdict_changes.csv`](results/verdict_changes.csv))

The baseline is the toy world at the **same 3 seeds and commit**
([`results/toy_3seeds/`](results/toy_3seeds/)). Seed count alone flips 19 of 765 toy cells (3 vs the
committed 5 seeds), which is the noise floor for this comparison.

| verdict (cells with a claim) | toy, 5 seeds (committed) | toy, 3 seeds | scDesign3, 3 seeds |
|---|---|---|---|
| ✅ detects | 134 | 134 | 121 |
| ❌ blind | 61 | 58 | 49 |
| ⚠️ false alarm | 22 | 21 | 11 |
| ⛔ wrong way | 3 | 3 | 10 |
| 〰️ erratic | 35 | 38 | 53 |
| ✔ stable | 280 | 276 | 274 |

112 of the 535 cells with a claim change verdict (10 of them also flip between 3 and 5 toy seeds); 169 of
765 change overall. Full matrix: [`results/bands/matrix.md`](results/bands/matrix.md).

**Clean baselines shift a lot.** The metrics' "perfect" values are not reached on correct labels in a real
tissue:

| metric | toy | scDesign3 |
|---|---|---|
| silhouette | 0.407 | 0.188 |
| `iv.neighborhood_purity` | 1.000 | 0.792 |
| isolated-label F1 | 1.00 | 0.77 |
| ISCG / ISCL | 0.936 / 1.000 | 0.621 / 0.895 |
| spatial coherence | 0.918 | 0.430 |
| PAS (lower = better) | 0.024 | 0.557 |
| CHAOS | 0.036 | 0.064 |
| `crisp_purity` | 0.800 | 0.109 |
| `mecr` (lower = better) | 0.045 | 0.372 |
| `pmp` / `nmp` | 0.996 / 1.000 | 0.682 / 0.768 |
| `fid.mean_auc` / Cohen's d | 0.995 / 6.25 | 0.900 / 2.01 |
| confidence margin / entropy | 10.6 / 0.014 | 2.83 / 0.56 |
| `onto.rho` | 0.26 | 0.56 |
| `cal.ece` | 0.001 | 0.078 |

More than half the cells of a correctly labelled real section are "abnormal" by PAS, and CRISP calls 89 % of
them impure.

**Verdicts that change:**

1. **Sibling confusion goes undetected without ground truth.** Silhouette, `integrated`, Ward-PropMatch,
   orbital medoid, the isolated-label scores, `pmp`, `nmp`, fidelity and `neighbor_contamination` all go
   from ✅ to ❌ at 75 % of T_CD8 called T_CD4 (e.g. silhouette 0.188 → 0.184, `pmp` 0.682 → 0.681). The two
   subtypes barely differ in the data (CD4 private markers 1.6×). Spatial coherence and PAS move ⛔ the
   wrong way (0.430 → 0.462, 0.557 → 0.499): the mixed T aggregate becomes more uniform in label.
   `mgo.mean_overlap` now detects it (0.82 → 0.72, ❌ → ✅).
2. **Spatial label metrics reward compact errors and merges.** A 50 % patch of Fibro called Epi raises
   coherence 0.430 → 0.499 and lowers PAS 0.557 → 0.493 (toy: flat, ❌; now ⛔). Merging three pairs raises
   coherence 0.430 → 0.592 and lowers PAS 0.557 → 0.366 and CHAOS 0.064 → 0.052 (PAS and CHAOS ❌ → ⛔). So
   does `iv.neighborhood_purity` (0.792 → 0.874, ❌ → ⛔). ISC now catches the patch (ISCG 0.621 → 0.545, ❌ → ✅).
3. **Their toy false alarms disappear where real tissue is already mixed.** Scattering every T_CD8 cell
   leaves coherence at 0.430 → 0.419 and PAS at 0.557 → 0.551 (toy ⚠️ 0.918 → 0.732; now ✔). CHAOS still
   ⚠️ (0.064 → 0.067).
4. **Co-expression is no longer a false alarm for CRISP / MECR** (0.109 → 0.104, 0.372 → 0.375; toy ⚠️).
   Real Myeloid cells already carry Fibro-marker counts, so the added co-expression is small against the
   baseline. `pmp` (0.682 → 0.661) and `neighbor_contamination` (0.0225 → 0.0279) stay ⚠️.
5. **Spillover: `nmp` and `ncp` lose it.** `nmp` falls 0.768 → 0.666 but its threshold is 0.188, because
   `over_split` drives it to 0.18 and widens its range; `ncp` stays at 0.999. CRISP (0.109 → 0.038), MECR,
   `pmp` (0.682 → 0.587) and `neighbor_contamination` (0.0225 → 0.0397, from 5 %) still catch it.
6. **Ambient noise: Moran's gap now detects it** (0.056 → 0.038, from 0.1; toy ❌). `frac_above_ctrl`
   moves 0.782 → 0.637 but stays under its threshold (0.37).
7. **Rare type and panel growth.** ISCG becomes a ⚠️ false alarm under `rare_type` (0.628 → 0.665). Under
   `panel_growth`, silhouette, isolated-label ASW, average similarity and Calinski-Harabasz rise (look better)
   with 960 uninformative genes: silhouette 0.188 → 0.248 (⚠️; the toy world gives ⚠️ at 5 seeds and ↗ at 3, so
   this is the same finding). `cal.ece` improves 0.078 → 0.011 (⚠️). Under `rare_type` the silhouette family is now ✔.
8. **`onto.rho` gets a real lineage structure** (clean 0.56 vs 0.26). It still does not detect label errors:
   it *rises* under label noise (0.56 → 0.65, ⛔) and is blind to sibling confusion and patch errors
   (`onto.rwr_rho` catches the patch, 0.507 → 0.426).
9. **ISCL is blind to batch and merge here too** (0.895 → 0.809 at log-sd 1.5, 0.895 → 0.985 at 3 merges),
   while ISCG catches batch (0.621 → 0.151, from 0.25).
10. **`nov.frac_novel` no longer reacts to batch** (0 → 0; toy 0 → 0.998), but still to ambient noise
    (0 → 0.84).

## Multi-sample variant and ISC regimes ([`results/isc_regimes.csv`](results/isc_regimes.csv))

`sample_mode="dirichlet"`: 6 separate pieces of 400 cells. Each piece takes its composition from
**Dirichlet(20 × natural shares)**, i.e. α = (2.9, 2.1, 1.5, 2.4, 6.5, 4.7) for (T_CD4, T_CD8, B, Myeloid,
Fibro, Epi). Each piece is a stratified subsample of the pool with its own count draw, placed side by side
with a gap so no spatial neighbours cross pieces. Seed 0 gives T_CD8 1.5-21 % and Myeloid 2.3-21 % across
pieces (`uns['world']['sample_composition']`). The `bands` world is a third regime, with strong natural
variation (Epi 3-48 % across bands). The toy world has constant composition.

The paper's fixed-threshold quadrant (`quadrant` in `check_isc_regimes.csv`; the call `regime` returned before
it was calibrated) over 3 seeds. For `merge` and `over_split` only the labels the perturbation made are counted.
The calibrated `regime` is compared across worlds in [`../README.md`](../README.md#isc-regime-calls): it reads
`consistent` for 16/18 (bands) and 14/16 (Dirichlet) clean labels and `merge` for 10 of 11 halves, with false
`remove` and `merge` calls described there.

| case | toy | scDesign3 bands | scDesign3 Dirichlet |
|---|---|---|---|
| clean (6 types) | consistent 18/18 | consistent 6, split 10, remove 2 | consistent 6, split 8, remove 2 (B scored in 1/3 seeds) |
| merge 1: Fibro+Epi | consistent 3/3 (ISCG 0.98) | consistent 3/3 (0.79) | consistent 3/3 (0.93) |
| merge 3: the three merged labels | consistent 9/9 | split 6, consistent 3 | split 5, consistent 4 |
| over_split 1: Epi halves | **merge** 6/6 | **merge** 3, consistent 2 (one half unscored) | **merge** 6/6 |
| label_noise 20 % | split 13, consistent 5 | split 11, remove 7 | split 6, remove 6, consistent 6 |
| batch log-sd 1.0 | split 18/18 | split 13, consistent 3, remove 2 | split 12, consistent 2, remove 2 |

- **Varying composition does not let ISC see a merge.** The merged Fibro+Epi label reads `consistent` in
  every world: ISCG 0.98 / 0.79 / 0.93 and ISCL 1.00 / 1.00 / 1.00. Even when Epi ranges 3-48 % across
  samples, the mixed pseudobulks stay closer to each other than to any other type.
- **The `split` calls on merged immune labels are not detections.** B+Myeloid and T_CD4+T_CD8 read `split`,
  but so do their unmerged parts on clean data: B, Myeloid and T_CD8 read `split` in 3/3 seeds in `bands`
  (ISCG 0.48-0.55 < 0.74, ISCL ≥ 0.91). In the realistic worlds, a correct small immune label already sits
  in the `split` quadrant.
- **Over-splitting still reads right**: `merge` for 6/6 halves under Dirichlet composition, 3/5 in bands.
- So on realistic data the quadrant is less specific than on the toy: clean labels read
  `consistent` for 6/18 (bands) and 6/16 (Dirichlet) against 18/18.

## Dedicated checks (results/bands/check_*.csv)

- **Conformal sets (α = 0.1, full depth, T_CD8 15 %).** Calibrated on the reference, they under-cover: 0.847
  overall and 0.556 on T_CD8 (toy 0.990 / 0.975). Mondrian brings T_CD8 to 0.743. Calibrated on 30 % of the
  section: 0.914, and Mondrian gives T_CD8 0.924. ECE is 0.090 (toy 0.001).
- **Novel type (B left out of the reference).** Conformal novelty p-values rank the novel cells only
  modestly (AUC 0.83 at full depth, 0.62 at a quarter; toy 0.999). At BH 10 % they flag **none** of them
  (toy 99 %).
- **`panel_resolvability`.** At depth 200 macro-F1 is 0.90 and `frac_resolvable` 0.67 with every marker kept
  (toy 1.00 / 1.0). Removing private markers barely matters here (0.88 with none), because real types differ
  across many genes, not just the 4 markers.
- **Deconvolution metrics** stay monotone: R² 1.00 → 0.29 at 50 % noise (toy 0.70), SSIM 1.00 → 0.65.

## Limits

- **One crop of one section and six types.** The other 37 % of cells are removed, so the layout has holes,
  and no DCIS, myoepithelial or endothelial cells are present to confuse with.
- **The world inherits the real data's annotation errors and segmentation bleed** in its fitted means.
  "Truth" is the authors' label: it is not independently validated, and the clean baselines above include
  whatever of it is wrong.
- **Spatial autocorrelation is half the real one, and within-type gene-gene correlation 0.64-0.81× the
  real one**; the copula is noisy for the small types (B: 643 cells).
- **Samples are not biological replicates.** Bands are strips of one tissue; Dirichlet pieces are
  subsamples of one pool. They share every cell-level effect, so the ISC results are about composition, not
  donor variation.
- **Infiltrating cells keep their origin's spatial expression effect**, and the controls are flat Poisson,
  not fitted.
- **3 seeds.** The toy comparison shows 19 of 765 cells flip with the seed count alone.
- `rare_type` starts at 10.5 % (natural T_CD8), not 15 %.
