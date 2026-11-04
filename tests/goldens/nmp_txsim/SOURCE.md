# nmp_txsim - negative marker purity (Salas et al. 2025), reference outputs

Parity golden for the Salas et al. negative marker purity (NMP), from the two codebases that implement
it: txsim and the authors' benchmark repository. Generated 2026-11-03. Our `nmp` is a different, raw
score, so this golden records the gap rather than an equality (see "Ours vs reference").

## Files

| File | Content |
|---|---|
| `reference.npz` | matched reference: `counts` (int32, 600 cells x 82 genes), `genes`, `cells`, `labels` (cell type) |
| `section_clean.npz` | clean section: `counts` (int32, 600 x 102), `genes`, `cells`, `labels` |
| `section_spillover.npz` | spillover section: `counts` (int32, 2400 x 102), `genes`, `cells`, `labels` |
| `label_permutations.npz` | `clean` (20 x 600) and `spillover` (20 x 2400), int16: permutation `k` relabels cell `i` with `labels[perm[k, i]]` |
| `expected.json` | reference outputs, `results[section][...]`, below; NaN would be stored as `null` (none occur) |

`expected.json` keys per section: `txsim_reads` and `txsim_cells` (`overall`, `per_celltype`, `per_gene`),
`moldia_seg_reads` and `moldia_seg_cells` (`overall` only), `txsim_reads_permuted` and `txsim_cells_permuted`
(`values` = the 20 permuted-label NMPs in permutation order, `mean`, `scaled_overall`).

## Data

Synthetic, made by `benchmarks/synthetic/sim.py` at commit `f71512f` (this repository, MIT); no external
data, no DOI. Counts are integer-valued in `sim.py` and stored as int32; the arrays are identical when
generated with numpy 2.2.6 and 2.5.4 (a seed-1 section differs, as a check that the comparison can fail).

- `reference` = `sim.reference(seed=0, n_per_type=100)`: 100 cells per type, 6 types, NegCtrl probes dropped.
- `section_clean` = `sim.simulate(seed=0, n_cells=600)`: 22 markers, 60 background genes, 20 NegCtrl probes.
- `section_spillover` = `sim.spillover(0.2, 0)`: each cell keeps a Binomial(0.8) share of its counts and
  gains Poisson(0.2 x the mean counts of its 6 nearest neighbours). `spillover()` simulates the default
  2400 cells, so this is not the clean section perturbed.
- Section labels are `obs['label']`, equal to the truth in both sections; reference labels are `obs['cell_type']`.

## Upstream

| | txsim | Moldia `Xenium_benchmarking` |
|---|---|---|
| URL | https://github.com/theislab/txsim | https://github.com/Moldia/Xenium_benchmarking |
| Commit | `bfdfaa9094c1f465d79b8b4ade855de8f2d2fb78` (main, 2025-06-27, version 0.1.2); same tree as `dev` `0f62644b`, which OpenProblems `task_ist_preprocessing` installs | `50249f7d03d9cfc1453c4b81ef3b1053c4d23052` (main, 2026-04-10) |
| File | `txsim/metrics/_negative_marker_purity.py` | `notebooks/5_segmentation_benchmark/metrics.py` (its README: "the functions to calculate the reported metrics") |
| Licence | none (no LICENSE file or licence metadata; GitHub reports none): run, never copied | none: run, never copied |

Paper: Salas et al., *Nat Methods* 2025, "Optimizing Xenium In Situ data utility by quality assessment and
best-practice analysis workflows", [10.1038/s41592-025-02617-2](https://doi.org/10.1038/s41592-025-02617-2)
(PMC11978515), Methods "Benchmarking segmentation algorithms".

## Exact calls

Both objects carry float64 counts in `X`, `layers['counts']` and `layers['raw']` (the Moldia functions read
`'raw'`), labels in `obs['celltype']`, and are first subset to the reference's genes that are in the section,
in reference order, as OpenProblems `task_ist_preprocessing` (`src/metrics/similarity/script.py` at `31019864`,
MIT) does before calling txsim (txsim indexes the reference by the section's genes, so the NegCtrl probes
would raise):

```python
txsim.metrics.negative_marker_purity_reads(sp, sc, key="celltype", layer="counts", pipeline_output=False)
txsim.metrics.negative_marker_purity_cells(sp, sc, key="celltype", layer="counts", pipeline_output=False)
metrics.negative_marker_purity_reads(sp, sc, key="celltype", pipeline_output=True)   # Moldia
metrics.negative_marker_purity_cells(sp, sc, key="celltype", pipeline_output=True)   # Moldia
```

Permuted baseline: the same two txsim calls with `pipeline_output=True` and
`obs['celltype'] = labels[perm[k]]` for each of the 20 stored permutations
(`np.random.default_rng(k).permutation(n_cells)`, `k = 0..19`; stored so a test does not depend on numpy's stream).

## Definitions, as implemented

`m[g,c]` = mean raw count of gene `g` over the cells of type `c`; `f[g,c]` = fraction of type-`c` cells
with a count > 0; suffixes `sp` (section) and `sc` (reference). Only types present in both with
>= 10 cells in each are used (fewer than 2 such types gives NaN). Checked: a plain numpy/pandas
re-implementation of these definitions (not committed) reproduces every number in `expected.json`, the 80
permuted values included, to <= 1.1e-16 under Python 3.10 and 3.12; the reads and Moldia-reads criteria
give different values (1.0 vs 0.999641 on the clean section), so the check tells the definitions apart.

- **txsim reads**: `x[g,c] = m[g,c] / sum_c' m[g,c']`, a cell-type-balanced share. `(g, c)` is a negative pair
  when `x_sc[g,c] < 0.005`. `X` = mean of `x` over the negative pairs. NMP = `1 - max(0, X_sp - X_sc)`; it is
  exactly 1 when the section is no worse than the reference (the clean section here). `per_celltype` and
  `per_gene` are means over a type's (gene's) negative pairs of `max(0, x_sp - x_sc)`: an excess, 0 = clean,
  clipped per pair, so not `1 -` the overall.
- **txsim cells**: `(g, c)` negative when `f_sc[g,c] < 0.005` (fewer than 0.5 % of the reference type's cells
  detect it). NMP = `1 - max(0, mean f_sp - mean f_sc)` over negative pairs. `per_celltype` / `per_gene` average
  `max(0, f_sc - f_sp)`, the reverse order of the overall, so they are 0 whenever `f_sc = 0`: here all three
  negative pairs (an upstream quirk, recorded as is). With 100 reference cells per type, `f_sc < 0.005` means
  no reference cell detects the gene, so only 3 pairs qualify on these data, by chance.
- **Moldia seg reads**: as txsim reads, but a pair is negative when `m_sc[g,c] / mean_c' m_sc[g,c'] < 0.005`
  (relative to the mean over types, i.e. a share below `0.005 / n_types`), the criterion txsim's code marks
  as old ("lead to too few negative markers"). Its `pipeline_output=False` branch calls `.cip` (a typo) and
  raises, so only the overall value is recorded.
- **Moldia seg cells**: the txsim cells algorithm (same values). Its cell-type filter checks the reference
  counts twice instead of the section's; no effect here (every type has >= 55 section cells).
- **Paper (Methods)**: a negative marker for a set of negative cell types is a gene "expressed in less than
  0.5% of cells in each of the negative cell types" (the cells criterion); NMP is computed on the
  cell-type-balanced read shares `x` (the reads formula); neither function combines exactly these two. Then
  "we scaled the NMP values such that 0 corresponds to the mean NMP of 20 random permutations of cell-type
  assignments for the binning reference, while the maximum is still 1". That scaling is in neither codebase:
  `permut` has 0 hits in `txsim/` and in the Moldia metric files (`notebooks/5_segmentation_benchmark/*.py`,
  `xb/*.py`), while `random` has 14 hits in `txsim/` and `np.random.permutation` is found in 6 Moldia
  notebooks. `scaled_overall` here is `(NMP - b) / (1 - b)` with `b` = the mean over the 20 permutations of
  the same section (our reading of "0 = permuted mean, maximum still 1"); the paper takes `b` from one
  reference segmentation and applies it to every method.

## Ours vs reference

**Since 0.7.0** `nmp(score="salas")`, the default, re-implements txsim reads from the definitions above and
matches it: overall, the 6 `per_celltype` and 22 `per_gene` values and the 20 stored-permutation values, max
|diff| 8.7e-19 on both sections under Python 3.10 and 3.12 (`tests/test_purity.py::test_nmp_salas_matches_txsim_reads`,
tolerance 1e-12; it fails when the 0.5 % cut is halved or multiplied by 20, or the per-pair clip is dropped,
but not at x2: no pair here sits between 0.5 % and 1 %, so `test_nmp_salas_share_cut_is_strict` pins the cut). Its
scaling draws its own 20 permutations (seed 0), so `b` differs from the stored ones: 0.833789 / 0.834392,
scaled 1.0000 / 0.9838, the same to 4 decimals as `scaled_overall`. `nmp(score="moldia")` reproduces the Moldia
relative-to-mean reads criterion (`test_nmp_moldia_matches_xenium_benchmarking_reads`), and `nmp(score="cells")`
reproduces the cell-fraction variant `txsim_cells` / `moldia_seg_cells` and its 20 stored permutations
(`test_nmp_cells_matches_txsim_cells_golden`), both at tolerance 1e-12. The v0.6.0 score below is `score="raw"`;
its values are unchanged.

Ours = `spatial_anno_metrics.nmp` at `f71512f` (v0.6.0), called as a user would: the full 102-gene section
and the reference, labels in `obs['cell_type']`. A gene is a negative for type `t` when it is detected in
<= 5 % of `t`'s reference cells and >= 25 % of another type's; each cell scores the fraction of its counts
(over the genes shared with the reference) on its type's negatives; `mean_nmp` is the mean over cells,
lower = cleaner, unscaled.

| section | txsim reads | txsim cells | Moldia reads | Moldia cells | txsim reads, scaled | txsim cells, scaled | ours `mean_nmp` | ours `types_without_negatives` |
|---|---|---|---|---|---|---|---|---|
| clean | 1.000000 | 0.993346 | 0.999641 | 0.993346 | 1.0000 | 0.9654 | 0.001030 | [] |
| spillover | 0.997312 | 0.957904 | 0.997048 | 0.957904 | 0.9838 | 0.7967 | 0.004628 | [] |

Permuted-label means `b`: reads 0.832627 (clean), 0.833563 (spillover); cells 0.807821, 0.792949. Every
score moves the right way from clean to spillover. The unscaled upstream scores sit within 0.05 of 1; the
scaling is what spreads them.

Per cell type (txsim reads excess vs our per-type mean of `obs['nmp']`, both 0 = clean):

| section | cell type | txsim reads per-type excess | ours mean nmp |
|---|---|---|---|
| clean | B | 0.000245 | 0.000621 |
| clean | Epi | 0.000328 | 0.001089 |
| clean | Fibro | 0.000165 | 0.000754 |
| clean | Myeloid | 0.000220 | 0.000873 |
| clean | T_CD4 | 0.000663 | 0.001503 |
| clean | T_CD8 | 0.000471 | 0.001295 |
| spillover | B | 0.004425 | 0.007219 |
| spillover | Epi | 0.001677 | 0.003306 |
| spillover | Fibro | 0.003695 | 0.005668 |
| spillover | Myeloid | 0.004252 | 0.006481 |
| spillover | T_CD4 | 0.000993 | 0.002282 |
| spillover | T_CD8 | 0.002409 | 0.003800 |

The two rank the six types identically on the spillover section; on the clean section (all values < 0.0016)
B, Myeloid and Fibro change places.

## Regenerate

1. Clone both repositories at the commits above into one directory `<dir>`.
2. Make a venv with txsim's import-time dependencies, resolved as of txsim's commit (txsim itself is imported
   from the clone; cellpose and planktonspace are left out, as txsim imports them only inside segmentation and
   annotation functions):

   ```bash
   uv venv <venv> --python 3.12
   uv pip install --python <venv>/bin/python --exclude-newer 2025-06-28 setuptools wheel squidpy alphashape \
       shapely scikit-image geopandas rasterio anndata scanpy numpy pandas scipy seaborn matplotlib tqdm \
       tifffile networkx scikit-learn
   ```

   Used here: Python 3.12.13, numpy 2.2.6, scipy 1.16.0, pandas 2.3.0, scanpy 1.11.2, anndata 0.11.4,
   scikit-learn 1.7.0, squidpy 1.6.5, matplotlib 3.10.3, seaborn 0.13.2, alphashape 1.3.1.
3. `MPLBACKEND=Agg <venv>/bin/python benchmarks/parity/nmp_txsim/make_golden.py --upstream <dir>` rewrites this
   directory (the npz files are byte-reproducible: fixed zip timestamps) and prints the tables above. The script
   stops if a clone is not at its pinned commit.
