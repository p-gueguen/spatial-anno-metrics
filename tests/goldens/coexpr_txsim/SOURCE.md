# coexpr_txsim - co-expression similarity (txsim) and NCP (Salas et al. 2025), reference outputs

Parity golden for the two reference-based co-expression scores of the spatial-QC literature: txsim's
`coexpression_similarity` (correlation structure, section vs reference) and the negative co-expression
purity (NCP) of Salas et al., from the authors' code. Generated 2026-11-03. We have no
`coexpression_similarity` yet (roadmap item 1 adds it beside `ncp`); our `ncp` is a third definition, so
this golden records the gap.

## Files

| File | Content |
|---|---|
| `reference.npz`, `section_clean.npz`, `section_spillover.npz` | the nmp_txsim inputs, byte-identical: `counts` (int32), `genes`, `cells`, `labels`; see `../nmp_txsim/SOURCE.md` "Data" |
| `expected.json` | reference outputs, `results[section][...]`; NaN would be stored as `null` (none occur) |

`expected.json` keys per section: `txsim_coexpression_similarity` (`overall`, `per_gene`),
`txsim_coexpression_similarity_by_celltype` (`overall`, `per_celltype`, `per_gene`), `moldia_ncp` (`overall`,
`n_negative_pairs`, `per_gene_as_conditioning`, `per_gene_as_target`); top level `lognorm_target_sum`.

## Data

Synthetic, from `benchmarks/synthetic/sim.py` at `f71512f` (this repository, MIT); no external data, no DOI:
`sim.reference(seed=0, n_per_type=100)`, `sim.simulate(seed=0, n_cells=600)` and `sim.spillover(0.2, 0)`
(2400 cells, not the clean section perturbed). Details in `../nmp_txsim/SOURCE.md`.

**Normalisation** (txsim's co-expression reads a log-normalised layer): `lognorm = log1p(1e4 * counts / total)`
in float64, `total` = the cell's counts over all genes of its own file (section: 102 genes incl. NegCtrl;
reference: 82), computed before subsetting to the 82 shared genes. `X` is set to `lognorm` as well, as
OpenProblems does with its `normalized` layer (txsim's gene filter reads `X`; `lognorm > 0` exactly where
counts > 0). OpenProblems normalises with `scanpy.pp.normalize_total(target_sum=None)` (the median total);
the fixed 1e4 keeps the input fully specified.

## Upstream

| | txsim | Moldia `Xenium_benchmarking` |
|---|---|---|
| URL | https://github.com/theislab/txsim | https://github.com/Moldia/Xenium_benchmarking |
| Commit | `bfdfaa9094c1f465d79b8b4ade855de8f2d2fb78` (main, 2025-06-27, version 0.1.2; same tree as `dev` `0f62644b`, which OpenProblems `task_ist_preprocessing` installs) | `50249f7d03d9cfc1453c4b81ef3b1053c4d23052` (main, 2026-04-10) |
| File | `txsim/metrics/_coexpression_similarity.py` | `xb/calculating.py` (`negative_marker_purity_coexpression`, called by `notebooks/3_techniques_comparison/3_4_negative_marker_purity_for_specificity.ipynb`, whose text calls the score NCP, "also refered as NMP through the pipeline") |
| Licence | none (no LICENSE file or licence metadata; GitHub reports none): run, never copied | none: run, never copied |

Papers: txsim is cited by OpenProblems as [10.1101/2023.02.13.528102](https://doi.org/10.1101/2023.02.13.528102);
NCP: Salas et al., *Nat Methods* 2025, [10.1038/s41592-025-02617-2](https://doi.org/10.1038/s41592-025-02617-2)
(PMC11978515), Methods "Benchmarking Xenium against other SRT methods".

## Exact calls

Both objects (labels in `obs['celltype']`, the lognorm above in `X` and `layers['lognorm']`) are first
subset to the reference's genes that are in the section, in reference order, as OpenProblems
`task_ist_preprocessing` (`src/metrics/similarity/script.py` at `31019864`, MIT) does:

```python
txsim.metrics.coexpression_similarity(sp, sc, key="celltype", layer="lognorm", pipeline_output=False)
txsim.metrics.coexpression_similarity(sp, sc, key="celltype", layer="lognorm", by_celltype=True,
                                      pipeline_output=False)
calculating.negative_marker_purity_coexpression(sp, sc, key="celltype", pipeline_output=False, minexp=0.0)
```

Defaults left as they are: `min_cells=20`, `thresh=0`, `correlation_measure="pearson"`.

With `by_celltype=True` txsim loops over the cell types as a Python `set`, whose order follows the string hash
seed; the type-averaged `per_gene` values (and possibly `overall`) then change in the last bit between runs
(seen: 20 of 433 values, up to 2.2e-16). The generator pins `PYTHONHASHSEED=0` and sorts `per_celltype`, so a
re-run gives identical bytes; a comparison needs a tolerance (1e-12 is ample).

## Definitions, as implemented

Checked: a plain numpy/pandas re-implementation of these definitions (not committed), with the lognorm above,
reproduces every number in `expected.json` to <= 3.3e-16 under Python 3.10 and 3.12.

- **txsim `coexpression_similarity`**: genes = shared genes detected (`X > 0`) in >= 20 cells in each modality
  (`scanpy.pp.filter_genes` on `X`; all 82 pass here). `rho` = gene-gene Pearson correlation over cells
  (`np.corrcoef`) of the layer, per modality. `S = 1 - |rho_sc - rho_sp| / 2`, pairs with `|rho_sc| < thresh`
  masked (none at 0), diagonal NaN. `overall` = NaN-mean of the upper triangle; `per_gene` = NaN-mean of each
  row. 1 = the same correlation structure.
  With `by_celltype=True`: per cell type present in both, genes are re-filtered within the type (>= 20 cells in
  each modality) and `S` computed on that type's cells; `per_celltype` = NaN-mean of its upper triangle;
  `overall` = mean of `per_celltype`; `per_gene` = row NaN-mean of the type-averaged `S`.
- **Moldia NCP** (`negative_marker_purity_coexpression`): `K[g, h]` = fraction of cells with `X_g > 0` among
  cells with `X_h > 0` (`minexp = 0`; a column with no positive cell is filled with 1), per modality, on `X`;
  cell types are not used. Negative pairs: `K_sc[g, h] < 0.05` (142 here, all among the 22 markers).
  NCP = `1 - mean over negative pairs of max(0, K_sp - K_sc)`. The per-gene outputs are unclipped means of
  `K_sp - K_sc` over a gene's negative pairs: by column `h` (`per_gene_as_conditioning`, upstream
  `purity_per_gene`) and by row `g` (`per_gene_as_target`; upstream calls it `purity_per_celltype`, but its
  rows are genes). It also prints and draws a clustermap (run with `MPLBACKEND=Agg`).
- **Paper (Methods)** describes NCP as "(1) identifying pairs of non-coexpressed genes in the reference single
  cell dataset, (2) quantifying the presence of these pairs across datasets and (3) computing the percentage of
  these pairs that were not coexpressed in situ": a fraction of pairs. The code computes the mean clipped excess
  above instead.

## Ours vs reference

**Since 0.7.0** `coexpression_similarity` (beside `ncp`) re-implements txsim's from the definitions above: it reads
the counts and log-normalises them exactly as here (1e4, totals over each file's own genes). Overall, `per_gene` and
(with `label_key`) `per_celltype` match, max |diff| 6.7e-16 on both sections under Python 3.10 and 3.12, and the
by-type `per_gene` is NaN for the one gene (`bg28`, clean) that no type scores
(`tests/test_purity.py::test_coexpression_similarity_matches_txsim`, tolerance 1e-12). The `thresh` mask (default 0)
and the Spearman option are not ported. On the synthetic sweep (`benchmarks/synthetic`, 3 seeds) it barely moves
with spillover (0.9788 -> 0.9775 at 0.4) and the by-type score rises slightly (0.9666 -> 0.9671).

Ours = `spatial_anno_metrics.ncp` at `f71512f` (v0.6.0), called as a user would on the full section and the
reference (counts): a shared-gene pair is non-coexpressed when `#(both > 0) / #(either > 0) <= 0.1` (top 400
genes by reference mean, first 2000 such pairs); `ncp` = fraction of the reference's non-coexpressed pairs
still `<= 0.1` in the section. Closest to the paper's wording of NCP; a different number from both upstream
functions, so no equality is expected.

| section | txsim `coexpression_similarity` | txsim, by cell type | Moldia NCP (negative pairs) | ours `ncp` |
|---|---|---|---|---|
| clean | 0.970679 | 0.938020 | 0.998583 (142) | 0.991031 |
| spillover | 0.973491 | 0.951777 | 0.970912 (142) | 0.991031 |

Moldia NCP drops with spillover 0.2; ours keeps the same 221 of 223 pairs in both sections. Ours does move
with stronger contamination (a probe outside this golden, same reference: `sim.spillover(0.4, 0)` 0.968610,
`sim.ambient(0.2, 0)` 0.0), so the flat value is insensitivity at this strength, not a broken call. The two
sections differ in size (600 vs 2400 cells), so the txsim values do not say how co-expression similarity
responds to spillover; they are parity targets only.

## Regenerate

Same clones and venv as nmp_txsim (`../nmp_txsim/SOURCE.md`, "Regenerate"), then

```bash
MPLBACKEND=Agg <venv>/bin/python benchmarks/parity/coexpr_txsim/make_golden.py --upstream <dir>
```

It reuses the nmp_txsim generator's input export (so the three npz files stay byte-identical to nmp_txsim's),
rewrites this directory and prints the table above. The script stops if a clone is not at its pinned commit.
