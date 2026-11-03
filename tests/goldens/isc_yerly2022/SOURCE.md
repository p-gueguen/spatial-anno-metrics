# ISC parity golden: Yerly 2022 BCC subsample vs scTypeEval

**Status: input only.** `expected.json` (the scTypeEval reference) is not generated yet - run
`benchmarks/parity/isc_yerly2022/run_sctypeeval.R` with scTypeEval installed (see "Regenerate").
Until then this folder holds the input and our own numbers.

| file | content |
|---|---|
| `counts.npz` | raw UMI counts, cells x genes CSR: `data` (uint16), `indices` (uint16), `indptr` (int32), `shape`, `genes` (2000 symbols, original row order). Load: `scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)`. 1.5 MB |
| `obs.csv` | `cell_id` (barcode in the source object), `patient_bcc` (sample), `annotation` (label); row order = matrix rows |
| `expected.json` | not yet generated (schema below) |

## Upstream reference tool

- **scTypeEval 1.1.2**, <https://github.com/carmonalab/scTypeEval>, commit `ca24658e13da5b7d776b4d69a962fac8a7c280ce`
  (master, 2026-07-16; Bioconductor 3.24 devel carries the same version number). Licence **GPL-3**: only outputs may be
  committed here, never its code. Bioconductor 3.23 (R 4.6 release) ships 1.0.0 (tag `v1.0.0`, `05e2c7a`); its ISC
  code path is identical, the only R change being `run_pca()`'s new `aggregation` / `center` / `scale.`
  arguments (same defaults).
- Paper: Garnica J, Andreatta M, Carmona SJ (2026), *Evaluating cell type annotations in single-cell omics in the
  absence of ground truth*, bioRxiv, doi:10.64898/2026.06.10.731285.
- Benchmark code read, not copied (no licence file): <https://github.com/carmonalab/ISC_benchmark_reproducibility>
  @ `8b01161`. Its `data_processing/config/core_datasets.yaml` (id `BCC_Yerly`) sets the label column `annotation`,
  the sample column `patient_bcc`, and excludes `Cancer_cells`.

## Data

- Yerly L. et al., *Integrated multi-omics reveals cellular and molecular interactions governing the invasive niche
  of basal cell carcinoma*, Nat Commun 13, 4897 (2022), doi:10.1038/s41467-022-32670-w, as redistributed by the ISC
  authors: Zenodo record 18921437, doi:10.5281/zenodo.18921437, licence **CC BY 4.0**, file
  `Yerly_2022_35986012.rds` (444,840,639 bytes, md5 `9d0768bc9f0b9f6240ae499028970455`, checked against the
  Zenodo listing). Seurat object, 23,892 genes x 37,018 cells, 7 patients, 10 labels.
- Subsample (`make_input.R`, generated 2026-11-03, R 4.6.0, `set.seed(22)`): `Cancer_cells` dropped (as the
  authors do), then at most 50 cells per (patient, label). **Deliberate edge case**: Mast in BCC3 is cut to 8 cells
  (< `min_cells = 10`), so a scored label is missing from one sample. 2,620 cells:

  | label | BCC1 | BCC2 | BCC3 | BCC4 | BCC5 | BCC6 | BCC7 | scored? |
  |---|---|---|---|---|---|---|---|---|
  | Endothelial, Fibroblast, Melanocyte, MoMacDC, Normal_Kerat, Tcell | 50 | 50 | 50 | 50 | 50 | 50 | 50 | yes, 7 samples |
  | Mast | 50 | 50 | 8 | 50 | 21 | 14 | 38 | yes, 6 samples |
  | Bcell | 32 | 50 | 37 | 42 | 6 | 3 | 4 | no, 4 samples >= 10 cells |
  | PlasmaCell | 50 | 1 | 5 | 50 | 1 | 2 | 6 | no, 2 samples >= 10 cells |

  48 (sample, label) pseudobulks are scored; 2,323 cells pass the filter.
- Genes: the 2000 HVGs scTypeEval's `run_hvg()` defaults would pick, computed with scran directly (scran 1.40.0):
  log1p(count / library size x 1e4) with the library size over all 23,892 genes, on the 2,323 cells that pass the
  filter; README black list removed (`black_list$TCR`, `$Immunoglobulins`, `$Ygenes` from scTypeEval's
  `data/black_list.rda`, 672 genes, 665 present); `modelGeneVar(block = patient_bcc, equiweight = TRUE)`, then
  `getTopHVGs(n = 2000)`. scran 1.40 marks these functions deprecated (scrapper replacements); the golden fixes
  the genes, so a later scran cannot change it.
- **Both tools get exactly these 2000 genes and normalise by their sum** (scTypeEval's default would use all
  genes for the library size; restricting keeps the file < 5 MB and the input identical for both tools).

## Exact reference call (`run_sctypeeval.R`)

The scTypeEval README recipe ("Minimal code to compute the recommended metrics"), HVG step replaced by the saved
genes:

```r
sc <- create_scTypeEval(matrix = counts, metadata = obs)
sc <- run_processing_data(sc, ident = "annotation", sample = "patient_bcc", min_samples = 5, min_cells = 10)
sc <- add_gene_list(sc, gene_list = list(HVG = rownames(counts)))
set.seed(22); sc <- run_pca(sc, gene_list = "HVG", aggregation = "pseudobulk", ndim = 30)
sc <- run_dissimilarity(sc, method = "recip_classif:Match", reduction = FALSE, reciprocal_classifier = "SingleR")
sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = TRUE)   # then again with reduction = FALSE
get_consistency(sc, dissimilarity_slot = ..., consistency_metric = c("silhouette", "2label_silhouette"))
```

- **ISCL** = `silhouette | recip_classif:Match`, **ISCG** = `2label_silhouette | Pseudobulk:Cosine` (README
  `sel.methods`; the paper's Fig. 5 notebook uses the older spellings `RecipClassif:Match`, `2label.silhouette`).
- `expected.json` (planned): `per_type[<obs label>][<metric | space>]` for both metrics in three spaces
  (`recip_classif:Match`, `Pseudobulk:Cosine`, `Pseudobulk:Cosine [reduction=FALSE]`), plus `product` =
  ISCL x ISCG and `geometric_mean` = sqrt(product) (null if either is negative); `mean_over_types` (NA labels
  skipped); `label_map` (obs label -> scTypeEval label); the three dissimilarity matrices with their row labels;
  versions. NaN/NA -> `null`.

Defaults read from the 1.1.2 source:

- **Filter**: a (sample, label) group needs >= `min_cells` = 10 cells; a label needs >= `min_samples` = 5 such
  samples, else it is reported as NA. Labels and sample ids go through `purge_label` (space, `_`, `+`, `-` -> `.`;
  commas dropped), so `Normal_Kerat` is reported as `Normal.Kerat`.
- **Pseudobulk**: sum of raw counts per group, then `Log1p` = log1p(x / column sum x 1e4) over every gene in the
  matrix, then restricted to the gene list; all-zero rows and columns are dropped.
- **ISCG**: `run_dissimilarity()` defaults to `reduction = TRUE`: cosine between the first `ndim` = 30 PCs
  (`irlba::prcomp_irlba`, centred, unscaled, `ndim` capped at min(dim) - 1) of the pseudobulk matrix. Then
  2-label silhouette (a = mean distance to the own label excluding self, b = mean distance to all other labels
  pooled), averaged per label.
- **ISCL**: `recip_classif:Match` (always `reduction = FALSE`): for every ordered sample pair, a sample with itself
  included, SingleR (`de.method = "classic"`: per label pair the top `500 x (2/3)^log2(N)` genes by median
  difference; Spearman on the marker union; quantile 0.8; fine-tuning; pruning, a pruned call becomes `_nores_`)
  is trained on one sample's pseudobulks and predicts the other's. Distance 0 for a mutual best hit, 0.5 for a
  mutual hit where either profile's label is absent from the other sample, 1 otherwise; `get_consistency()` reads
  0.5 as 1. Then `cluster::silhouette`, averaged per label.
- **Combined**: the paper multiplies per label (Fig. 5) and scores a dataset as the mean over labels of that
  product (label-transfer task, R^2 0.98 vs macro-F1); the README says geometric mean.
- The package documentation says a default black list is applied when `black_list = NULL`; the code does not
  (`check_blacklist()` returns the object's own, empty by default). The README passes it explicitly, as done here.

## Ours vs reference

Ours: `spatial_anno_metrics` 0.6.0 (`f71512f`),
`inter_sample_consistency(adata, "annotation", "patient_bcc", counts_layer="counts")` on exactly these files
(`make_golden.py`). Reference columns fill in once `expected.json` exists.

| label | ours ISCG | ours ISCG, 30-PC space | scTypeEval ISCG | ours ISCL | scTypeEval ISCL | ours ISCG x ISCL | scTypeEval product |
|---|---|---|---|---|---|---|---|
| Endothelial | 0.859444 | 0.880916 | pending | 1.000000 | pending | 0.859444 | pending |
| Fibroblast | 0.835680 | 0.861465 | pending | 1.000000 | pending | 0.835680 | pending |
| Mast | 0.801156 | 0.813336 | pending | 1.000000 | pending | 0.801156 | pending |
| Melanocyte | 0.853157 | 0.882166 | pending | 1.000000 | pending | 0.853157 | pending |
| MoMacDC | 0.869293 | 0.901179 | pending | 1.000000 | pending | 0.869293 | pending |
| Normal_Kerat | 0.879024 | 0.928565 | pending | 1.000000 | pending | 0.879024 | pending |
| Tcell | 0.882772 | 0.919472 | pending | 1.000000 | pending | 0.882772 | pending |
| mean over labels | 0.854361 | 0.883871 | pending | 1.000000 | pending | 0.854361 | pending |

Ours drops Bcell (4 samples) and PlasmaCell (2) into `dropped_types`; scTypeEval reports them as NA.
"30-PC space" = our 2-label silhouette on cosine between the first 30 PCs of the centred pseudobulk matrix (exact
SVD), i.e. ours moved into scTypeEval's default ISCG space; it is not the reference.

Gaps expected from the source, to be measured against `expected.json`:

1. **ISCG space.** Ours: cosine on uncentred log-normalised pseudobulks. scTypeEval default: cosine on 30 centred
   PCs. Moving ours into that space raises ISCG by +0.012 to +0.050 per label (mean +0.030). The
   `reduction=FALSE` reference uses our space, so it should equal our ISCG to float precision: that is the
   decisive check of the pseudobulk recipe and the 2-label silhouette. Compare the `reduction=TRUE` values at
   ~1e-4, not 1e-12: irlba (random start, `set.seed(22)`, default tolerance) vs an exact SVD.
2. **ISCL classifier.** Ours: nearest profile by Spearman over every gene in `X`. scTypeEval: SingleR classic
   (pairwise markers, fine-tuning, pruning). **On this input ours gives ISCL = 1 for every label**: all 141
   same-label cross-sample pairs are mutual best hits and 0 cross-label pairs are. ISCL is saturated here, so a
   reference ISCL of 1 would agree without testing the classifier. An over-split label column (the paper's
   over-partitioning task) would make ISCL informative.
3. **Label missing from a sample.** Ours scores any mutual hit 0; scTypeEval scores 1 when either profile's label
   is absent from the other sample (Mast in BCC3 here). Ours found no cross-label mutual hit, so this rule does
   not move our numbers on this input.
4. **Combined score.** Ours `integrated` = mean over labels of ISCG x ISCL = the paper's dataset score (0.854361
   here). The README's geometric mean gives a different scale: mean over labels of sqrt(ISCG x ISCL) = 0.924205.
5. **Headline weighting.** Ours `iscg` / `iscl` average over pseudobulks (a label weighs by its number of samples):
   `iscg` = 0.855469, while the mean over labels is 0.854361 (Mast has 6 pseudobulks, the others 7). scTypeEval
   reports per label; the paper averages over labels.
6. **HVGs and library size.** The same 2000 genes and gene-sum library sizes for both tools by construction; in
   normal use ours scores every gene in `adata` and scTypeEval its HVG list with all-gene library sizes.

## Regenerate

```bash
# 1. input (needs SeuratObject, scran; not scTypeEval)
Rscript benchmarks/parity/isc_yerly2022/make_input.R Yerly_2022_35986012.rds black_list.rda <work>/input
# black_list.rda: https://github.com/carmonalab/scTypeEval/raw/ca24658e13da5b7d776b4d69a962fac8a7c280ce/data/black_list.rda
# 2. golden files + our numbers (checks the npz/csv round trip and the 5 MB cap)
PYTHONPATH=src python benchmarks/parity/isc_yerly2022/make_golden.py <work>/input tests/goldens/isc_yerly2022 <work>/ours.json
# 3. reference (needs scTypeEval 1.1.2 at the commit above, SingleR, jsonlite)
Rscript benchmarks/parity/isc_yerly2022/run_sctypeeval.R <work>/input tests/goldens/isc_yerly2022/expected.json
```

Step 1 is deterministic: two runs on 2026-11-03 wrote byte-identical files.
