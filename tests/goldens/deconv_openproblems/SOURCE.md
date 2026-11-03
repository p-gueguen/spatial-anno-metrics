# deconv_openproblems - OpenProblems spatial decomposition metrics, reference outputs

Parity golden for `deconvolution_metrics` against the OpenProblems `task_spatial_decomposition` metric
scripts (`r2`, `jsd`). Generated 2026-11-03. Ours matches `r2` exactly and `jsd` to <= 6e-17, except when a
cell type is absent from every true spot: OpenProblems returns NaN, ours a number (see "Ours vs reference").

## Files

| File | Content |
|---|---|
| `proportions.json` | `{pair: {"true": 40 x 5, "pred": 40 x 5}}`, spots x cell types, rows sum to 1, floats round-trip exactly |
| `expected.json` | `results[pair]` = `{"r2", "jsd"}`; NaN is stored as `null` (`absent_type` `jsd`) |

## Data

Synthetic matrices made by `benchmarks/parity/deconv_openproblems/make_golden.py` (this repository, MIT) with
`np.random.default_rng(0)`; no external data, no DOI.

| pair | true | pred |
|---|---|---|
| `close` | Dirichlet(1) rows | row-normalised `0.9 * true + 0.1 * Dirichlet(1)` |
| `noisy` | the same `true` | an unrelated Dirichlet(1) draw |
| `absent_type` | Dirichlet(1) over types 0-3, type 4 = 0 in every spot | row-normalised `0.95 * true + 0.05 * Dirichlet(1)` (type 4 gets a little mass) |

## Upstream

- `openproblems-bio/task_spatial_decomposition` at `892ab4818b3d28656cd1a05f71b34c4d31af1de2` (main,
  2025-06-11), https://github.com/openproblems-bio/task_spatial_decomposition. Licence MIT (`license: MIT` in
  `_viash.yaml`; the repository has no LICENSE file).
- Files: `src/metrics/r2/script.py`, `src/metrics/jsd/script.py`, executed unchanged except for their
  `## VIASH START ... ## VIASH END` parameter block, which is replaced by our file paths.

## Exact call

For each pair: `true` goes to a solution h5ad as `obsm['proportions_true']`, `pred` to a method-output h5ad as
`obsm['proportions_pred']` (with `uns['dataset_id']`, `uns['method_id']`, which the scripts copy to their
output); each script is run with `par = {'input_method', 'input_solution', 'output'}` and the value read back
from the score h5ad's `uns['metric_values'][0]`. What the scripts compute:

```python
sklearn.metrics.r2_score(prop_true, prop_pred, sample_weight=None, multioutput="uniform_average")    # r2
np.mean(scipy.spatial.distance.jensenshannon(prop_true, prop_pred, axis=0))                          # jsd
```

## Definitions, as implemented

- `r2`: R^2 of each cell type's proportions across spots, averaged over cell types. A constant true column
  (the absent type) scores 1.0 if predicted exactly, else 0.0 (scikit-learn's `force_finite` default).
- `jsd`: Jensen-Shannon distance (natural log) between the true and predicted proportions of each cell type
  across spots, columns renormalised to sum 1 by scipy, averaged over cell types. An all-zero true column
  divides 0 by 0, so that type's distance is NaN, and so is the mean.

## Ours vs reference

Ours = `spatial_anno_metrics.deconvolution_metrics(true, pred)` at `f71512f` (v0.6.0); `r2` and `jsd` are the
same formulas, except that `jsd` averages with `np.nanmean`, so a NaN cell type is dropped instead of turning
the score into NaN. Ours also returns `jsd_per_spot`, `rmse` and `pearson`, which OpenProblems does not have.

| pair | OpenProblems `r2` | ours `r2` | OpenProblems `jsd` | ours `jsd` |
|---|---|---|---|---|
| close | 0.9757759961507725 | 0.9757759961507725 | 0.06708556329896388 | 0.0670855632989639 |
| noisy | -0.9234989842728909 | -0.9234989842728909 | 0.3716543474874577 | 0.3716543474874576 |
| absent_type | 0.7968233676068135 | 0.7968233676068135 | NaN (`null`) | 0.023234464854914334 |

The `jsd` differences on `close` and `noisy` are one floating-point step (summation order of `np.mean` vs
`np.nanmean`).

## Regenerate

```bash
git clone https://github.com/openproblems-bio/task_spatial_decomposition <dir>
git -C <dir> checkout 892ab4818b3d28656cd1a05f71b34c4d31af1de2
uv venv <venv> --python 3.12
uv pip install --python <venv>/bin/python --exclude-newer 2025-06-12 anndata scikit-learn scipy numpy pandas
<venv>/bin/python benchmarks/parity/deconv_openproblems/make_golden.py --upstream <dir>
```

Used here: Python 3.12.13, anndata 0.11.4, scikit-learn 1.7.0, scipy 1.15.3, numpy 2.3.0, pandas 2.3.0 (the
upstream metrics run in the `openproblems/base_python:1` image). The script stops if the clone is not at the
pinned commit, rewrites this directory and prints the table above.

## Test

`tests/test_external.py::test_deconvolution_matches_openproblems_golden` compares `r2` (all three pairs) and `jsd`
(`close`, `noisy`) at rel/abs 1e-12, and for `absent_type` asserts the upstream NaN and that ours equals the
mean JS distance of the four defined cell types. Positive control: adding 1e-9 to one expected value fails it.
The SSIM, AUPR and rank-score keys added in 0.7.0 have no reference outputs here: running the Li et al. 2022
reference code was not permitted in the generating session, so `_ssim_li` is checked against the paper's formula
in `test_deconvolution_ssim_and_aupr` instead.
