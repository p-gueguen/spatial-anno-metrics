# Simulation worlds for the synthetic sweep

The toy world in [`../synthetic/sim.py`](../synthetic/sim.py) has six types in stripes, Poisson counts and markers
~1,000× their off-target rate. Its detection thresholds are therefore optimistic. This folder holds three realistic
worlds with the same contract (`simulate` / `reference` / `marker_sets`, the same six type names). The 15 scenarios,
the metrics and the scoring in [`../synthetic/sweep.py`](../synthetic/sweep.py) run on them unchanged.

| world | counts | layout | samples | per-world README |
|---|---|---|---|---|
| `toy` | Poisson, crisp markers | stripes | 6 bands | [`../synthetic/README.md`](../synthetic/README.md) |
| `scdesign3` | scDesign3 NB GAM + Gaussian copula, fitted to one crop of the 10x Xenium breast section (Janesick et al. 2023, Rep1) | the crop's real centroids | 6 bands (`--sample-mode bands`) or 6 Dirichlet-composition pieces (`dirichlet`) | [`scdesign3/README.md`](scdesign3/README.md) |
| `srtsim` | SRTsim reference-based draws (per type and crop, rank-matched to the real counts) on the same section | real centroids of six 400 µm crops | one crop per sample; **each seed is a different crop set** | [`srtsim/results/README.md`](srtsim/results/README.md) |
| `sccube` | scCube VAE trained on Yerly 2022 BCC scRNA-seq, on the 306 Xenium breast-panel genes, Xenium depth | scCube patterns `stripes`, `domains`, `fine`, `mixed` | 6 tiles of 400 cells | [`sccube/README.md`](sccube/README.md) |

## Run

```bash
PYTHONPATH=src python benchmarks/synthetic/sweep.py --world toy                     # the default, as before
PYTHONPATH=src python benchmarks/synthetic/sweep.py --world scdesign3 [--sample-mode dirichlet] [--alpha0 20]
PYTHONPATH=src python benchmarks/synthetic/sweep.py --world srtsim [--fixed-crops 0]
PYTHONPATH=src python benchmarks/synthetic/sweep.py --world sccube --pattern domains
python benchmarks/worlds/compare_worlds.py          # results/cross_world_*.csv|md from the committed results
```

`--seeds`, `--workers`, `--only` and `--out` work as before. Without `--out`, a world writes to its own
`benchmarks/worlds/<world>/results[/<variant>]`, never to the toy results. `worlds.install` (in
[`__init__.py`](__init__.py)) does the swap for every world. It patches the `sim` module the sweep reads, the
`rare_type` levels (scDesign3: they start at the crop's own T_CD8 share), the two edge cases that name toy genes,
and SRTsim's deconvolution check, which rescales the coordinates to the toy's density. The sweep then forks its
workers (`fork` context), so they inherit the patch. Each world's `run.py` is a thin wrapper around it.
BLAS / OpenMP threads are set to 1 per worker before numpy loads (`worlds.limit_threads`). Unset, a 128-core node
started 128 threads in every worker.

The caches are not in git: `SAM_SCDESIGN3_CACHE`, `SAM_SRTSIM_CACHE`, `SAM_SCCUBE_CACHE`, each defaulting to
`/srv/GT/analysis/pgueguen/sam_dev/data/sim_worlds/<world>/`. The per-world READMEs say how to rebuild them.
`--world toy` writes the same files as before byte for byte. Run with 5 seeds on this commit and on the previous
one (305a48f) under Python 3.12, all 10 output files are identical, apart from traceback line numbers in
`errors.csv`. Under Python 3.10 `sig.moran_gap` varies between two runs of the same code by up to 5e-14, so
byte identity cannot hold there. The three `run.py` wrappers also match the scripts on their original branches
byte for byte (1 seed, `label_noise`, all edge cases and checks). The committed `../synthetic/results/` predate
later `src/` changes: a rerun moves values by up to 4e-5 (Python 3.10, no verdict change) and changes 1 of 765
verdicts under Python 3.12 (`batch` × `iso.f1`, ↓ → 〰️, a cell without a claim).

## Verdicts, world by world

Each cell: **toy / scDesign3 (bands) / SRTsim / scCube**. The scCube entry is one symbol when its four patterns
agree, else the counts (`✅3 ❌1`). Only cells with a claim are shown. Symbols are as in
[`../synthetic/README.md`](../synthetic/README.md): ✅ detects · ❌ blind · ⛔ wrong way · 〰️ erratic · ✔ stable
· ⚠️ false alarm · ↗ drifts up (looks better under a null). The toy uses 5 seeds and the realistic worlds 3. Per-cell data,
with each world's clean value: [`results/cross_world_verdicts.csv`](results/cross_world_verdicts.csv).

| metric | label_noise | sibling_confusion | type_swap | patch_error | over_split | merge | spillover | depth_loss | ambient | batch | coord_shuffle | coexpression | rare_type | infiltration | panel_growth |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `cal.ece` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ⚠️3 ↗1 | ✔ / ✔ / ✔ / ✔ | ✔ / ⚠️ / ✔ / ✔ |
| `clisi` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ | ❌ / 〰️ / 〰️ / 〰️ | ❌ / ✅ / ✅ / ✅ | ✅ / 〰️ / 〰️ / ✅ | ❌ / ❌ / ❌ / ❌ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ⚠️ / ✔ |
| `coex.by_type` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / 〰️ / ✅ | ❌ / ✅ / ❌ / ✅ | ❌ / 〰️ / ❌ / 〰️ | ❌ / ❌ / ✅ / ✅ |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ⚠️ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ↗ |
| `coex.global` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ❌ / ❌ / ❌ / ✅ |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ↗ / ✔3 ↗1 |
| `conf.mean_entropy` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |
| `conf.mean_margin` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |
| `crisp_purity` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✅ / ✅ / ✅ / ✅ |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.ari` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.avg_bio` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.comp_l1` | ❌ / ✅ / ✅ / ❌ | ✅ / ✅ / ✅ / ✅ | ❌ / ✅ / ✅ / ❌ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.ecs` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.hier_acc` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.macro_f1` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ext.spari` | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `fid.mean_auc` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ❌ / ✅ / ❌ / 〰️ | ❌ / 〰️ / ❌ / ❌2 〰️2 |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ⚠️3 ✔1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `fid.mean_cohens_d` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ✅2 ❌1 〰️1 | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / 〰️ / ❌2 〰️2 | 〰️ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ↗3 ✔1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `isc.integrated` | ✅ / ✅ / ✅ / ✅ | ❌ / 〰️ / ❌ / ❌3 〰️1 | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ✅ / ❌ / 〰️ | ✅ / ✅ / ✅ / ✅ | ⛔ / 〰️ / ❌ / ❌3 〰️1 |  |  |  | ✅ / ✅ / ✅ / ✅ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ⚠️ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ⚠️ / ⚠️ |
| `isc.iscg` | ✅ / ✅ / ✅ / ✅ | ❌ / 〰️ / ❌ / ❌ | ❌ / 〰️ / ❌ / ❌ | ❌ / ✅ / ❌ / ✅ | ❌ / 〰️ / ❌ / ❌ | ⛔ / 〰️ / ❌ / 〰️ |  |  |  | ✅ / ✅ / ✅ / ✅ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ⚠️ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ⚠️ / ⚠️ |
| `isc.iscl` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / 〰️ / ❌3 〰️1 | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ✅ / ✅ / 〰️ | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ |  |  |  | ❌ / ❌ / ✅ / ❌ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔3 ⚠️1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ⚠️ / ✔ |
| `iso.asw` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ⚠️ / ⚠️ |
| `iso.f1` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / ❌ / 〰️ | ✅ / ❌ / ❌ / ✅3 ❌1 | ✅ / ✅ / ❌ / ✅ | ✅ / ❌ / ❌ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔3 ⚠️1 | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ✔3 ⚠️1 |
| `iv.avg_similarity` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ⚠️ / ⚠️ |
| `iv.calinski_harabasz` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅3 ❌1 | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ⚠️3 ✔1 | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ⚠️ / ⚠️ |
| `iv.davies_bouldin` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ⚠️ / ✔ |
| `iv.integrated` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ↗ / ⚠️ / ⚠️ |
| `iv.neighborhood_purity` | ✅ / ✅ / ✅ / ✅ | ✅ / 〰️ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ❌ / ⛔ / ⛔ / ⛔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ↗ / ↗ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ⚠️ / ✔ |
| `iv.orbital_medoid` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ✅ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ❌ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ⚠️ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ⚠️ / ⚠️3 ✔1 |
| `iv.silhouette` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ⚠️ / ⚠️ |
| `iv.ward_propmatch` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / ❌ / 〰️ | ✅ / ✅ / ❌ / ✅ | ✅ / ✅ / ❌ / ✅ | ✅ / 〰️ / ❌ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️3 ✔1 | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ⚠️ / ⚠️ |
| `mecr` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✅ / ✅ / ✅ / ✅ |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `mgo.mean_overlap` | ❌ / ✅ / ❌ / 〰️3 ⛔1 | ❌ / ✅ / ❌ / ❌ | ✅ / ✅ / ❌ / ✅ | ❌ / ❌ / ❌ / ❌ | ❌ / ❌ / 〰️ / 〰️ | ❌ / ❌ / 〰️ / ❌ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔3 ⚠️1 | ✔ / ↗ / ✔ / ✔3 ⚠️1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ncp` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | 〰️ / ❌ / ❌ / ❌3 ✅1 |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ncp.moldia` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✅ / ❌ / ❌ / ❌ |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ↗ / ↗ |
| `neighbor_contamination` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | 〰️ / 〰️ / 〰️ / 〰️ | ✅ / ❌ / ❌ / ❌3 〰️1 | ❌ / 〰️ / 〰️ / ❌2 〰️2 | ❌ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ✅ / ✅ |  |  |  | ✔ / ✔ / ✔ / ⚠️ | ✔ / ⚠️ / ✔ / ⚠️3 ✔1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ↗ / ✔2 ↗2 |
| `nmp` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ✅ / ❌ | ❌ / ✅ / ❌ / ❌ | ❌ / 〰️ / ❌ / 〰️ | ✅ / ❌ / ✅ / ✅3 ❌1 |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `nov.frac_novel` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `ns.n_suspicious` | ❌ / 〰️ / 〰️ / ❌3 〰️1 | 〰️ / ❌ / ❌ / 〰️2 ✅2 | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ | ❌ / 〰️ / ⛔ / ❌ |  |  |  |  | 〰️ / ✅ / ✅ / 〰️3 ✅1 | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ⚠️ | ✔ / ↗ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ |
| `onto.rho` | 〰️ / ⛔ / ❌ / 〰️3 ❌1 | ❌ / ❌ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ❌ / ⛔ / ⛔ | ❌ / 〰️ / ❌ / 〰️ | 〰️ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `onto.rwr_rho` | 〰️ / ⛔ / ❌ / 〰️3 ❌1 | ❌ / ❌ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ✅ / ⛔ / ⛔ | ❌ / 〰️ / ❌ / 〰️ | 〰️ / 〰️ / 〰️ / 〰️ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `pmp` | ✅ / ✅ / ✅ / ✅ | ✅ / ❌ / ❌ / ✅ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ❌ / ✅ | ❌ / 〰️ / ❌ / ✅ | ❌ / 〰️ / 〰️ / 〰️ | ✅ / ✅ / ❌ / ✅ |  |  |  | ✔ / ✔ / ✔ / ✔ | ⚠️ / ⚠️ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `pmp.no_dict_frac` | ❌ / ❌ / ❌ / ❌ | ❌ / ❌ / ❌ / ❌ | ❌ / ❌ / ❌ / ❌ | ❌ / ❌ / ❌ / ❌ | ✅ / ✅ / ✅ / ✅ | ✅ / ✅ / ✅ / ✅ |  |  |  |  | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |
| `sc.chaos` | ✅ / ✅ / ✅ / ✅ | ❌ / ❌ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️3 ❌1 | ❌ / ❌ / ❌ / ❌3 ⛔1 | ✅ / ✅ / ✅ / ✅ | ❌ / ⛔ / ⛔ / ⛔3 ❌1 |  |  |  |  | ✅ / ✅ / ✅ / ✅ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ✔3 ⚠️1 | ⚠️ / ⚠️ / ⚠️ / ⚠️ | ✔ / ✔ / ✔ / ✔ |
| `sc.mean_coherence` | ✅ / ✅ / ✅ / ✅ | 〰️ / ⛔ / ❌ / ❌ | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ⛔ / ❌ / ❌3 ⛔1 | ✅ / ✅ / ✅ / ✅ | ⛔ / ⛔ / ⛔ / ⛔ |  |  |  |  | ✅ / ✅ / ✅ / ✅ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ✔2 ↗2 | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ |
| `sc.pas` | ✅ / ✅ / ✅ / ✅ | 〰️ / ⛔ / ❌ / ❌3 〰️1 | 〰️ / 〰️ / 〰️ / 〰️ | ❌ / ⛔ / ❌ / ❌3 ⛔1 | ✅ / ✅ / ✅ / ✅ | ❌ / ⛔ / ❌ / ❌2 ⛔2 |  |  |  |  | ✅ / ✅ / ✅ / ✅ | ✔ / ✔ / ✔ / ✔ | ✔ / ↗ / ✔ / ✔2 ↗2 | ⚠️ / ✔ / ✔ / ⚠️ | ✔ / ✔ / ✔ / ✔ |
| `sig.frac_above_ctrl` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  | ❌ / ❌ / ❌ / 〰️3 ✅1 |  | ✅ / ✅ / ✅ / ✅ |  |  |  |  |
| `sig.moran_gap` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  | ❌ / ✅ / ✅ / 〰️3 ✅1 |  | ✅ / ✅ / ✅ / ✅ |  |  |  |  |
| `sig.snr` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  | ✅ / ✅ / ✅ / ✅ |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |
| `sig.tx_per_area` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  | ✅ / ✅ / ✅ / ✅ |  |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |
| `sig.zero_fraction` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  | ✅ / ✅ / ✅ / ✅ |  |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |
| `stab.flip_rate` | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ | ✔ / ✔ / ✔ / ✔ |  |  |  |  | ✔ / ✔ / ✔ / ✔ |  |  |  |  |

Verdict counts over the 535 cells with a claim ([`results/cross_world_counts.csv`](results/cross_world_counts.csv)):

| world | ✅ | ❌ | ⚠️ | ⛔ | 〰️ | ✔ | ↗ |
|---|---|---|---|---|---|---|---|
| toy, 5 seeds | 134 | 61 | 22 | 3 | 35 | 280 | 0 |
| scDesign3 bands | 121 | 49 | 11 | 10 | 53 | 274 | 17 |
| SRTsim (crops vary with the seed) | 112 | 74 | 16 | 6 | 41 | 282 | 4 |
| SRTsim, fixed crops (control, below) | 123 | 38 | 30 | 7 | 65 | 263 | 9 |
| scCube `stripes` | 125 | 48 | 27 | 5 | 55 | 270 | 5 |
| scCube `domains` | 129 | 47 | 28 | 5 | 52 | 270 | 4 |
| scCube `fine` | 130 | 42 | 28 | 6 | 55 | 267 | 7 |
| scCube `mixed` | 129 | 41 | 29 | 9 | 54 | 266 | 7 |

The toy at 3 seeds instead of 5 flips 13-19 of the 535 cells (the noise floor quoted in the scDesign3 and scCube
READMEs).

## Clean baselines

The value of each metric on correctly labelled, undamaged sections (`label_noise` at strength 0, mean over seeds;
[`results/cross_world_baselines.csv`](results/cross_world_baselines.csv) has all 51). The "perfect" values of the
toy are not reached on correct labels in any realistic world:

| metric | toy | scDesign3 | SRTsim | scCube stripes | domains | fine | mixed |
|---|---|---|---|---|---|---|---|
| PAS (lower = better) | 0.022 | 0.557 | 0.538 | 0.062 | 0.217 | 0.508 | 0.508 |
| spatial coherence | 0.920 | 0.430 | 0.452 | 0.796 | 0.654 | 0.472 | 0.416 |
| CHAOS (lower = better) | 0.036 | 0.064 | 0.052 | 0.032 | 0.040 | 0.045 | 0.061 |
| CRISP purity | 0.813 | 0.109 | 0.196 | 0.760 | 0.758 | 0.758 | 0.758 |
| MECR (lower = better) | 0.042 | 0.372 | 0.284 | 0.038 | 0.038 | 0.038 | 0.038 |
| PMP | 0.996 | 0.682 | 0.628 | 0.895 | 0.894 | 0.894 | 0.894 |
| NMP (Salas) | 1.000 | 0.768 | 0.978 | 0.997 | 0.998 | 0.998 | 0.998 |
| `neighbor_contamination` (lower = better) | 0.0003 | 0.022 | 0.023 | 0.009 | 0.008 | 0.012 | 0.012 |
| ISCG / ISCL | 0.938 / 1.000 | 0.621 / 0.895 | 0.525 / 0.893 | 0.892 / 0.970 | 0.898 / 0.982 | 0.898 / 0.982 | 0.898 / 0.982 |
| silhouette | 0.404 | 0.188 | 0.153 | 0.302 | 0.297 | 0.297 | 0.297 |
| kNN neighbourhood purity | 1.000 | 0.792 | 0.795 | 0.883 | 0.883 | 0.883 | 0.883 |
| isolated-label F1 | 1.000 | 0.770 | 0.744 | 0.892 | 0.878 | 0.878 | 0.878 |
| cLISI | 1.000 | 0.968 | 0.961 | 0.996 | 0.997 | 0.997 | 0.997 |
| marker AUC (`fid.mean_auc`) | 0.995 | 0.900 | 0.897 | 0.951 | 0.951 | 0.951 | 0.951 |
| `onto.rho` | 0.255 | 0.560 | 0.598 | 0.618 | 0.609 | 0.609 | 0.609 |
| ECE (reference classifier) | 0.001 | 0.078 | 0.089 | 0.039 | 0.039 | 0.039 | 0.039 |

The spatial label metrics read the tissue's geometry as much as the labels. In scCube, where the same cells are laid
out four ways, PAS on correct labels runs from 0.06 (`stripes`) to 0.51 (`fine`, `mixed`). In the two worlds fitted
to the real section, more than half the cells of a correct annotation read as "abnormal". The purity metrics read the
tissue's own lineage mixing: CRISP purity is 0.11-0.20 on the real-section worlds and 0.76 in scCube, whose VAE
markers are cleaner than real ones.

## Findings

**Holds in every realistic world** (scDesign3, SRTsim, all four scCube patterns):

1. **Merging two labels looks like an improvement to the neighbourhood metrics.** Spatial coherence rises in all
   seven worlds, the toy included (scDesign3 0.43 → 0.59, SRTsim 0.45 → 0.57, scCube `mixed` 0.42 → 0.58). The kNN
   neighbourhood purity, ❌ in the toy, turns ⛔ in all six (0.79 → 0.87 scDesign3, 0.88 → 0.96 scCube). CHAOS turns
   ⛔ in five of six and PAS in three. The toy's merge detections by silhouette, `iv.integrated`, orbital medoid and
   isolated-label ASW all become 〰️.
2. **Conformal novelty never fires.** With B left out of the reference, the novelty p-values still rank B above known
   cells (AUC 0.83 scDesign3, 0.83 SRTsim, 0.90 scCube at full depth; 0.60-0.71 at a quarter of it). At BH 10 % they
   flag **0 %** of the B cells in every realistic world (toy: AUC 1.00, 97 % flagged). `nov.frac_novel` no longer
   reacts to batch either (0 → 0 everywhere; toy 0 → 0.998).
3. **NMP and average similarity miss sibling confusion.** At 75 % of T_CD8 called T_CD4, NMP is ❌ in all six
   (scDesign3 0.768 → 0.800, scCube 0.998 → 0.989; toy ✅ 1.000 → 0.968), and so is `iv.avg_similarity`. The marker
   AUC misses it everywhere, the toy included. Whether the rest of the reference-free metrics miss it depends on the
   world (below).
4. **`ncp(score="moldia")` loses spillover** in all six (SRTsim 0.986 → 0.977, scCube 0.666 → 0.652; toy ✅
   1.000 → 0.962). CRISP, MECR and `neighbor_contamination` catch spillover in all seven worlds.
5. **cLISI catches a compact mislabelled patch**, which the toy could not show because cLISI sat at its ceiling
   there: scDesign3 0.968 → 0.888, SRTsim 0.961 → 0.857, scCube 0.997 → 0.974.
6. **The toy's false alarms of silhouette and isolated-label ASW under `rare_type` and `coexpression`, and of
   `iv.integrated` under `rare_type`, vanish** (toy silhouette 0.404 → 0.310 under `rare_type`; ✔ in all six).
7. **CHAOS false-alarms under infiltration in every world.** ISCG catches batch in every world; ISCL catches it only
   in SRTsim (0.893 → 0.676). `onto.rho` finally has a lineage structure to read (clean 0.56-0.62 against 0.26 in
   the toy) but still does not detect label errors: it rises under `label_noise` everywhere and under `patch_error`
   in SRTsim and every scCube pattern (⛔, 0.61 → 0.72).

**World-specific:**

- **Sibling confusion is invisible without ground truth only on the real Xenium T-cell markers.** In the two worlds
  fitted to the section (scDesign3, SRTsim), CD4 T cells' own markers are 1.5-1.6× their level in the next type.
  Nearly every label-only and marker-purity metric that caught the toy's 75 % confusion goes ❌ in both. That
  covers silhouette (scDesign3 0.188 → 0.184), PMP (0.682 → 0.681), `neighbor_contamination` (0.0225 → 0.0224),
  isolated-label F1 and ASW, and the Ward and integrated scores. The one exception is orbital medoid, in SRTsim
  only (0.818 → 0.773). Hierarchical accuracy goes ❌ too (1 → 0.96, threshold 0.05). What still catches it is the
  ground-truth partition and accuracy metrics, and the reference-based per-type co-expression similarity
  (`coex.by_type`, 0.967 → 0.964 against a threshold of 0.001, ✅ in every world). With scCube's scRNA-seq
  markers (CCL5 / NKG7 against LTB / IL7R) silhouette (0.297 → 0.247), PMP (0.894 → 0.838), isolated-label F1 and
  `neighbor_contamination` still detect it. The SRTsim fixed-crop control (below) confirms this is the data, not
  seed noise: silhouette 0.172 → 0.165 against a threshold of 0.009.
- **NMP's spillover readout.** NMP catches spillover in SRTsim (0.978 → 0.767) and the three patchy scCube patterns
  (0.998 → 0.73-0.89). It misses it in scDesign3 (0.768 → 0.666; its clean value is far from 1, and `over_split`
  stretches its range, so the threshold is 0.19) and in scCube `stripes` (0.997 → 0.952, threshold 0.069). NMP stays
  flat under genuine co-expression in every world.
- **The toy's co-expression and infiltration false alarms vanish only where the tissue is already mixed.** In
  scDesign3 and SRTsim, 40 % co-expression leaves CRISP purity at 0.109 → 0.104 and 0.196 → 0.182, and MECR inside the
  noise. Scattering all T_CD8 cells leaves coherence at 0.43 → 0.42 and 0.45 → 0.44 (toy 0.92 → 0.73). In every
  scCube pattern both false alarms stay (CRISP 0.76 → 0.68; coherence 0.80 → 0.52 in `stripes`), because its markers
  are cleaner than real and its domains are purer.
- **A compact mislabelled patch reads as better** to coherence and PAS in scDesign3 (0.43 → 0.50, 0.557 → 0.493) and
  scCube `mixed` (0.42 → 0.45). `neighbor_contamination`, which caught the patch in the toy, is blind to it in
  scDesign3, SRTsim and three scCube patterns.
- **Ambient RNA**: Moran's signal gap detects it in scDesign3 (0.056 → 0.038), SRTsim (0.092 → 0.053) and scCube
  `domains` (toy ❌).
- **Panel growth**: 960 uninformative genes turn ISCG into a false alarm in SRTsim (0.525 → 0.266) and in every
  scCube pattern (0.898 → 0.677), where the added genes are real low-F genes or i.i.d. draws. scDesign3's added genes
  carry a mild per-type programme, and there ISCG stays put (0.621 → 0.642).
- **`neighbor_contamination` false-alarms when only coordinates are shuffled** in all four scCube patterns
  (0.0085 → 0.018), and under co-expression in scDesign3 and three scCube patterns.
- **ISC regime calls on correct labels**: toy 18/18 `consistent`; scDesign3 bands 6/18 `consistent` (immune types
  read `split`); SRTsim 14 of 17 `split`; scCube only B (naive B plus plasma cells, two real programmes) `split`.

## SRTsim control: fixed crops

In the SRTsim world each seed draws a different set of six tissue crops. The move rule (3 × the seed SD at strength
0, or 5 % of the metric's range) therefore measures between-section variability there. The control
(`sweep.py --world srtsim --fixed-crops 0`, [`srtsim/results/fixed_crops_0/`](srtsim/results/fixed_crops_0/)) keeps
the crops of seed 0 for all 3 seeds and redraws only the SRTsim counts. The redraw is `srtsim_count` (domain
scheme, random tie breaks) ported to Python on the cached fits. Against the R draw of the same crops it agrees on
median library size (119 vs 118), zero fraction (0.820 vs 0.820) and per-gene log mean (r 0.9995). The run took 116 s
on 16 cores.

| | toy | SRTsim, crops vary | SRTsim, fixed crops |
|---|---|---|---|
| claim cells where the seed SD sets the threshold | 19 % | 69 % | 19 % |
| median 3 × seed SD, as % of the metric's range | 0.3 % | 8.9 % | 0.9 % |
| ✅ / ❌ / ⚠️ / 〰️ | 134 / 61 / 22 / 35 | 112 / 74 / 16 / 41 | 123 / 38 / 30 / 65 |

- **Of the 30 toy detections that SRTsim lost, 9 come back with fixed crops**, so between-crop variance is behind
  those: `iso.f1` under `over_split`, `patch_error` and `sibling_confusion`; PMP under `patch_error` (0.641 → 0.578,
  threshold 0.023; varying crops: 0.12) and `spillover` (0.641 → 0.552); Ward-PropMatch and `coex.by_type` under
  `patch_error`; isolated-label ASW under `sibling_confusion`; `mgo` under `type_swap`.
- **21 stay lost, so the metrics or the data are behind them.** These are 10 `sibling_confusion` cells (silhouette,
  PMP, NMP, `neighbor_contamination`, average similarity, Calinski-Harabasz, neighbourhood purity, Ward,
  `iv.integrated`, hierarchical accuracy), 6 `merge` cells (the silhouette family and the isolated-label scores) and 3 `over_split`
  cells. The rest are `ncp.moldia` under spillover and `neighbor_contamination` under `patch_error`.
- **The stricter bar was also hiding false alarms**: they rise from 16 to 30 and erratic cells from 41 to 65.
  Examples are ISCG under `rare_type` (0.583 → 0.612, threshold 0.024) and coherence under infiltration (0.494 →
  0.474, threshold 0.018). That coherence move is a tenth of the toy's (0.92 → 0.73). So the vanishing infiltration
  false alarm in SRTsim (above) is mostly real but not complete.
- Fixed crops are one crop set, so the clean values shift a little (coherence 0.452 → 0.494, PMP 0.628 → 0.641).
  Conformal novelty still flags 0 % of B (AUC 0.83).

## Caveats

- **SRTsim: a seed is a different tissue section.** Its varying-crop verdicts use a stricter bar than the other
  worlds (above). Read its ❌ as "not past between-section variability", and use the fixed-crop run for "not past
  draw noise".
- **scDesign3 and SRTsim are fitted to the same Xenium breast section** and share its author annotation as truth.
  Their agreement is therefore not independent evidence, and both inherit the section's segmentation spillover and
  annotation errors. scDesign3 uses one crop of it. Its bands are strips of one tissue and its Dirichlet pieces are
  subsamples of one pool, so neither is a biological replicate. Its spatial autocorrelation is about half the real
  one (Moran's I median 0.056 vs 0.102).
- **scCube puts a breast panel on skin-cancer cells** (306 of 313 Xenium breast genes on Yerly 2022 BCC
  scRNA-seq). Its markers are cleaner than real (specificity 87 vs 63) and its cells more variable within type.
  Every tile is a 400-cell sample, so its "domains" are narrow. One trained VAE serves all seeds. The three GRF
  patterns share their expression, so their non-spatial verdicts are identical.
- **3 seeds** in every realistic world against the toy's 5. Six types only: no myoepithelial, endothelial or DCIS
  cells to confuse with (left out of the real-section worlds).
- The dedicated checks (`check_*.csv`) and the scDesign3 Dirichlet variant (4 scenarios,
  [`scdesign3/results/dirichlet/`](scdesign3/results/dirichlet/)) are discussed in the per-world READMEs and are not
  in the table above.
