# rel_pairwise_txsim - relative pairwise cell-type expression similarity (txsim), reference outputs

Parity golden for `relative_pairwise_ct_expression` against `theislab/txsim` (`relative_pairwise_celltype_expression`), using the same synthetic inputs as `../nmp_txsim/` and `../coexpr_txsim/`.

## Files

| File | Content |
|---|---|
| `reference.npz`, `section_clean.npz`, `section_spillover.npz` | the nmp_txsim inputs, byte-identical (`counts`, `genes`, `cells`, `labels`) |
| `expected.json` | reference `txsim` outputs (`overall`, `per_celltype`, `per_gene`) |

## Upstream

- URL: https://github.com/theislab/txsim (no LICENSE; run only to generate `expected.json`, never copied)
- File: `txsim/metrics/_relative_pairwise_celltype_expression.py`
