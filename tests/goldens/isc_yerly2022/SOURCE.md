# ISC parity golden: Yerly 2022 BCC subsample vs scTypeEval

**Status: complete.** Input, scTypeEval reference for two label columns, and parity tests. Our ISCG equals
scTypeEval's (`reduction = FALSE`) to 1.3e-14; given scTypeEval's own reciprocal-classification distances our ISCL
silhouette equals its ISCL to 1e-16; the one deliberate difference, the classifier inside RCM (ours: Spearman
nearest profile; scTypeEval: SingleR), is pinned pair by pair on an over-split label column.

| file | content |
|---|---|
| `counts.npz` | raw UMI counts, cells x genes CSR: `data` (uint16), `indices` (uint16), `indptr` (int32), `shape`, `genes` (2000 symbols, original row order). Load: `scipy.sparse.csr_matrix((data, indices, indptr), shape=shape)`. 1.5 MB |
| `obs.csv` | `cell_id` (barcode in the source object), `patient_bcc` (sample), `annotation` (label), `annotation_oversplit` (over-split variant, below); row order = matrix rows. 0.15 MB |
| `expected.json` | scTypeEval output, 0.27 MB (schema below) |
| `test_isc_yerly2022.py` | input checks + the parity tests (below) |

`expected.json`: `meta` (versions: scTypeEval, its source, SingleR, scran, irlba, cluster, R, date), `names` (the
metric keys for ISCL, ISCG, ISCG without PCA), and per label column `columns.<column>` with `label_map` (obs label ->
scTypeEval label), `per_type[<obs label>]`, `mean_over_types` (NA labels skipped) and `dissimilarity` (the three
distance matrices, rows named `<sample>_<label>` as scTypeEval names its pseudobulks). `per_type` keys:
`<metric> | <space>` for metric `silhouette` / `2label_silhouette` in space `recip_classif:Match`,
`Pseudobulk:Cosine` (README default, 30 PCs) and `Pseudobulk:Cosine [reduction=FALSE]`; `product` = ISCL x ISCG and
`geometric_mean` = sqrt(product), null when either is negative. NA/NaN -> `null`.

## Upstream reference tool

- **scTypeEval 1.0.0 from Bioconductor 3.23** (`RELEASE_3_23`, git `05e2c7a` = GitHub tag `v1.0.0`) generated
  `expected.json`. GitHub master is **1.1.2** (<https://github.com/carmonalab/scTypeEval>, commit
  `ca24658e13da5b7d776b4d69a962fac8a7c280ce`, 2026-07-16), the source the defaults below were read from. Licence
  **GPL-3**: only outputs are committed here, never its code.
- **1.0.0 vs 1.1.2, checked object by object** (`benchmarks/parity/isc_yerly2022/compare_sctypeeval_src.R`: every
  top-level object of the installed namespace deparsed against the same object sourced from each GitHub tree, plus
  S4 slot types, the `initialize` method and the `data/black_list.rda` md5). The installed package equals `05e2c7a`
  in all 89 objects, 5 classes, `initialize` and the black list. Against 1.1.2 exactly two objects differ,
  `run_pca` and `custom_prcomp`: 1.1.2 adds `aggregation` / `center` / `scale.` arguments whose defaults reproduce
  1.0.0 (PCA on every processed aggregation, centred, unscaled). Every function on the ISC path is identical
  (processing, pseudobulk, Log1p, gene filtering, cosine, the SingleR helper, RCM match, both silhouettes,
  `get_consistency`). Positive control: a copy of the 1.0.0 tree with four one-token edits (`max` -> `min` in
  `compute_2label_silhouette`, `0.5` -> `0.4` in `match_tidy`, one slot type of `dim_red`, one line of
  `initialize`) is reported as exactly those four differences.
- Paper: Garnica J, Andreatta M, Carmona SJ (2026), *Evaluating cell type annotations in single-cell omics in the
  absence of ground truth*, bioRxiv, doi:10.64898/2026.06.10.731285.
- Benchmark code read, not copied (no licence file): <https://github.com/carmonalab/ISC_benchmark_reproducibility>
  @ `8b01161`. Its `data_processing/config/core_datasets.yaml` (id `BCC_Yerly`) sets the label column `annotation`,
  the sample column `patient_bcc`, and excludes `Cancer_cells`. Its over-partitioning task (`SplitCelltype`,
  helper `rand_split_group`) relabels a random half of one cell type's cells in every sample.

## Data

- Yerly L. et al., *Integrated multi-omics reveals cellular and molecular interactions governing the invasive niche
  of basal cell carcinoma*, Nat Commun 13, 4897 (2022), doi:10.1038/s41467-022-32670-w, as redistributed by the ISC
  authors: Zenodo record 18921437, doi:10.5281/zenodo.18921437, licence **CC BY 4.0**, file
  `Yerly_2022_35986012.rds` (444,840,639 bytes, md5 `9d0768bc9f0b9f6240ae499028970455`, checked against the
  Zenodo listing). Seurat object, 23,892 genes x 37,018 cells, 7 patients, 10 labels.
- Subsample (`make_input.R`, R 4.6.0, `set.seed(22)`): `Cancer_cells` dropped (as the authors do), then at most 50
  cells per (patient, label). **Deliberate edge case**: Mast in BCC3 is cut to 8 cells (< `min_cells = 10`), so a
  scored label is missing from one sample. 2,620 cells:

  | label | BCC1 | BCC2 | BCC3 | BCC4 | BCC5 | BCC6 | BCC7 | scored? |
  |---|---|---|---|---|---|---|---|---|
  | Endothelial, Fibroblast, Melanocyte, MoMacDC, Normal_Kerat, Tcell | 50 | 50 | 50 | 50 | 50 | 50 | 50 | yes, 7 samples |
  | Mast | 50 | 50 | 8 | 50 | 21 | 14 | 38 | yes, 6 samples |
  | Bcell | 32 | 50 | 37 | 42 | 6 | 3 | 4 | no, 4 samples >= 10 cells |
  | PlasmaCell | 50 | 1 | 5 | 50 | 1 | 2 | 6 | no, 2 samples >= 10 cells |

  48 (sample, label) pseudobulks are scored; 2,323 cells pass the filter.
- **Over-split variant** `annotation_oversplit` (the paper's over-partitioning task at rate 0.5, re-implemented in
  `make_input.R`, `set.seed(23)`, the benchmark's first split replicate): in every sample a random half of the
  Fibroblast and of the Tcell cells (25 of 50) is relabelled `Fibroblast.Nspl2` / `Tcell.Nspl2`, the other half
  keeps the label. The same cells pass the filter and the HVGs are unchanged; 62 pseudobulks (48 + 14). Each half is
  a random draw from one population, so the right ISCL for a half is low (the paper's "merge" regime). Regenerating
  the input with this column left `counts.mtx`, `genes.txt`, `hvg_rank.txt` and the first three `obs.csv` columns
  byte-identical to the first run (positive control: one changed sample id shows up in the same diff).
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
genes, run once per label column:

```r
sc <- create_scTypeEval(matrix = counts, metadata = obs)
sc <- run_processing_data(sc, ident = <column>, sample = "patient_bcc", min_samples = 5, min_cells = 10)
sc <- add_gene_list(sc, gene_list = list(HVG = rownames(counts)))
set.seed(22); sc <- run_pca(sc, gene_list = "HVG", ndim = 30)
sc <- run_dissimilarity(sc, method = "recip_classif:Match", reduction = FALSE, reciprocal_classifier = "SingleR")
sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = TRUE)   # then again with reduction = FALSE
get_consistency(sc, dissimilarity_slot = ..., consistency_metric = c("silhouette", "2label_silhouette"))
```

- 1.0.0's `run_pca()` has no `aggregation` argument (the call drafted against 1.1.2 passed
  `aggregation = "pseudobulk"` and would stop with "unused argument"). The call above runs on both versions and,
  with 1.1.2's defaults, takes the same code path: single-cell PCA first, then pseudobulk, from one seed. On
  `annotation` irlba warns "You're computing too large a percentage of total singular values" (30 of at most 47
  components for 48 pseudobulks).
- Two runs wrote byte-identical `expected.json` files (18 s for both columns).
- **ISCL** = `silhouette | recip_classif:Match`, **ISCG** = `2label_silhouette | Pseudobulk:Cosine` (README
  `sel.methods`; the paper's Fig. 5 notebook uses the older spellings `RecipClassif:Match`, `2label.silhouette`).

Defaults read from the 1.1.2 source (identical in 1.0.0, see above):

- **Filter**: a (sample, label) group needs >= `min_cells` = 10 cells; a label needs >= `min_samples` = 5 such
  samples, else it is reported as NA. Labels and sample ids go through `purge_label` (space, `_`, `+`, `-` -> `.`;
  commas dropped), so `Normal_Kerat` is reported as `Normal.Kerat`. `run_sctypeeval.R` calls
  `scTypeEval:::purge_label` for its label map instead of copying the rule; re-run on 2026-11-03, it wrote an
  `expected.json` byte-identical to the committed one (md5 `eb118de101ca5411f716c8c238d75204`).
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
`inter_sample_consistency(adata, <column>, "patient_bcc", counts_layer="counts")` on the committed files; the tables
are printed by `make_golden.py`. **Since 0.7.0** the default ISCG space is scTypeEval's (`reduction="pca"`,
`n_pcs=30`: the "ISCG sTE 30 PCs" column, matched to 1.4e-6 / 8.5e-7); the "ISCG ours" column below is now
`reduction=None`, and the headline `iscg` / `iscl` are label means. `make_golden.py` prints both. sTE = scTypeEval 1.0.0. "ISCG sTE 30 PCs" is scTypeEval's default space; "product"
= ISCG x ISCL per label (ours: our ISCG, i.e. our `per_type_integrated`; sTE: its default ISCG); "geo. mean" =
sqrt(sTE product), the README's integrated score.

**`annotation`** (the authors' labels):

| label | ISCG ours | ISCG sTE, reduction=FALSE | ISCG sTE, 30 PCs | ISCL ours | ISCL sTE | product ours | product sTE | geo. mean sTE |
|---|---|---|---|---|---|---|---|---|
| Endothelial | 0.859444 | 0.859444 | 0.880917 | 1.000000 | 1.000000 | 0.859444 | 0.880917 | 0.938572 |
| Fibroblast | 0.835680 | 0.835680 | 0.861465 | 1.000000 | 1.000000 | 0.835680 | 0.861465 | 0.928151 |
| Mast | 0.801156 | 0.801156 | 0.813336 | 1.000000 | 1.000000 | 0.801156 | 0.813336 | 0.901851 |
| Melanocyte | 0.853157 | 0.853157 | 0.882165 | 1.000000 | 1.000000 | 0.853157 | 0.882165 | 0.939236 |
| MoMacDC | 0.869293 | 0.869293 | 0.901179 | 1.000000 | 1.000000 | 0.869293 | 0.901179 | 0.949304 |
| Normal_Kerat | 0.879024 | 0.879024 | 0.928565 | 1.000000 | 1.000000 | 0.879024 | 0.928565 | 0.963621 |
| Tcell | 0.882772 | 0.882772 | 0.919472 | 1.000000 | 1.000000 | 0.882772 | 0.919472 | 0.958891 |
| mean over labels | 0.854361 | 0.854361 | 0.883871 | 1.000000 | 1.000000 | 0.854361 | 0.883871 | 0.939947 |

**`annotation_oversplit`** (Fibroblast and Tcell randomly halved per sample):

| label | ISCG ours | ISCG sTE, reduction=FALSE | ISCG sTE, 30 PCs | ISCL ours | ISCL sTE | product ours | product sTE | geo. mean sTE |
|---|---|---|---|---|---|---|---|---|
| Endothelial | 0.861802 | 0.861802 | 0.919285 | 1.000000 | 1.000000 | 0.861802 | 0.919285 | 0.958793 |
| Fibroblast | 0.765263 | 0.765263 | 0.795751 | 0.293481 | 0.035941 | 0.224590 | 0.028600 | 0.169116 |
| Fibroblast.Nspl2 | 0.760458 | 0.760458 | 0.794789 | -0.054082 | -0.126531 | -0.041127 | -0.100565 | null |
| Mast | 0.805999 | 0.805999 | 0.849795 | 1.000000 | 1.000000 | 0.805999 | 0.849795 | 0.921843 |
| Melanocyte | 0.857120 | 0.857120 | 0.912727 | 1.000000 | 1.000000 | 0.857120 | 0.912727 | 0.955367 |
| MoMacDC | 0.871675 | 0.871675 | 0.929565 | 1.000000 | 1.000000 | 0.871675 | 0.929565 | 0.964140 |
| Normal_Kerat | 0.882922 | 0.882922 | 0.956302 | 1.000000 | 1.000000 | 0.882922 | 0.956302 | 0.977907 |
| Tcell | 0.824684 | 0.824684 | 0.884126 | 0.152438 | 0.455952 | 0.125713 | 0.403119 | 0.634917 |
| Tcell.Nspl2 | 0.797777 | 0.797777 | 0.838404 | -0.208163 | -0.183673 | -0.166068 | -0.153993 | null |
| mean over labels | 0.825300 | 0.825300 | 0.875638 | 0.575964 | 0.575743 | 0.491403 | 0.527204 | 0.797440 (7 labels) |

Both tools report Bcell (4 samples) and PlasmaCell (2) as not scored: ours in `dropped_types`, sTE as NA.
Up to 0.6.x our headline `iscg` / `iscl` were means over pseudobulks, not labels (a label weighs by its number of
samples): `annotation` 0.855469 / 1.000000, `annotation_oversplit` 0.825611 / 0.569124 (kept since 0.7.0 as
`iscg_profile_mean` / `iscl_profile_mean`); `integrated` is the label mean of the products (0.854361 / 0.491403).

**Dataset-level combined score**, three formulas (labels with a negative ISCL have no geometric mean; the mean then
runs over the other 7 of 9):

| formula | `annotation` ours | sTE, 30-PC ISCG | sTE, ISCG reduction=FALSE | `annotation_oversplit` ours | sTE, 30-PC ISCG | sTE, ISCG reduction=FALSE |
|---|---|---|---|---|---|---|
| mean over labels of ISCG x ISCL (paper; ours `integrated`) | 0.854361 | 0.883871 | 0.854361 | 0.491403 | 0.527204 | 0.493365 |
| mean over labels of sqrt(ISCG x ISCL) (README) | 0.924205 | 0.939947 | 0.924205 | 0.779094 | 0.797440 | 0.772034 |
| mean ISCG x mean ISCL (product of label means) | 0.854361 | 0.883871 | 0.854361 | 0.475343 | 0.504143 | 0.475161 |

**What matches, to float precision** (each check run on both columns):

1. ISCG with `reduction = FALSE` = ours: max difference 1.1e-14 (7 labels) and 1.3e-14 (9 labels). This settles
   the pseudobulk recipe (sum, gene-sum library size, log1p x 1e4, all-zero genes dropped) and the 2-label
   silhouette. The cosine distance matrices themselves (48 x 48, 62 x 62) differ by at most 9.8e-15. Positive
   control: one extra UMI in one of the 2,620 cells moves our ISCG by 5.8e-6 against the reference.
2. Our legacy `per_type` (silhouette on pseudobulk cosine) = sTE `silhouette | Pseudobulk:Cosine [reduction=FALSE]`:
   max difference 1.3e-14.
3. ISCL silhouette: fed scTypeEval's RCM distances (0.5 read as 1), our ISCL code returns its ISCL to 8e-17 on all
   labels, including the over-split halves (-0.18 to 0.46). The ISCL gap below is therefore the classifier alone.
   Positive control: with our own RCM distances, which differ on 30 pairs, ISCL moves by up to 0.30.
4. RCM distances: `annotation`: identical on all 1,128 pseudobulk pairs (141 mutual hits, all same-label, no
   cross-label hit). `annotation_oversplit`: identical on 1,861 of 1,891 pairs; the 30 that differ all join two
   halves of one split label (Fibroblast 10, Tcell 20, of the 168 cross-sample pairs among halves of one label);
   15 are mutual hits for us only and 15 for SingleR only, 141 mutual hits in total for each tool.

**Remaining gaps, each explained:**

1. **ISCG space.** Ours: cosine on uncentred log-normalised pseudobulks (= sTE `reduction = FALSE`, check 1). sTE
   default: cosine on 30 centred PCs. The default reads higher by +0.012 to +0.050 per label (mean +0.030) on
   `annotation`, +0.030 to +0.073 (mean +0.050) on `annotation_oversplit`. That is the projection alone: ours moved
   into the 30-PC space with an exact SVD (`make_golden.py`) matches the sTE default to 1.4e-6 and 8.5e-7, the
   residual being irlba's truncated SVD against the exact one. **Closed in 0.7.0**: `reduction="pca"` (default)
   is that exact-SVD space; `test_iscg_default_30_pcs_matches_sctypeeval` holds it within 2e-6.
2. **ISCL classifier.** Ours: nearest profile by Spearman over all 2000 genes. sTE: SingleR classic (pairwise top
   markers, Spearman on their union, 0.8 quantile, fine-tuning, pruning). On the authors' labels both classifiers
   make the same mutual hits on every pair (check 4), so ISCL = 1 everywhere and the classifier is invisible. On
   the over-split column the two halves of a label are near-identical profiles (our smallest Spearman margin
   between a best and second-best hit: 5e-5), so which half a classifier picks is noise: 30 pairs differ and
   per-half ISCL differs by up to 0.30 (Tcell: ours 0.152, sTE 0.456), while the mean over the four halves agrees
   (0.0459 vs 0.0454) and the mean over all labels differs by 0.0002. Through our 0.45 ISCL guide threshold, ours
   reads all four halves as `merge`; sTE's values would read 3 of 4 (Tcell at 0.456 stays `consistent`).
   0.7.0 adds `rcm_genes="markers"` (SingleR's classic marker union per training sample, no fine-tuning or
   pruning): identical to SingleR on `annotation` (0 of 1,128 pairs differ) and 28 of 1,891 differing pairs on
   `annotation_oversplit` (all genes: 30), per-half ISCL 0.143 / -0.095 / 0.088 / -0.229, so the halves stay
   classifier noise; `rcm_genes=<gene list>` restricts the classifier to e.g. HVGs.
3. **Label missing from a sample.** sTE scores a mutual hit 0.5 (read as 1) when either profile's label is absent
   from the other sample; ours scored 0 up to 0.6.x. Not exercised here: neither sTE matrix holds a 0.5 (Mast in
   BCC3 makes no cross-label mutual hit), so this rule moves no number on this input. **Implemented in 0.7.0**
   (`_rcm_dissimilarity` returns 0 / 0.5 / 1 like `recip_classif`, read as 1 before the silhouette) and tested on
   synthetic profiles in `tests/test_isc.py`; also as scTypeEval, a sample with < 2 profiles matches nothing.
4. **Combined score.** The formula dominates: the README's geometric mean sits on a square-root scale, +0.070
   (ours, `annotation`) and +0.288 (ours, `annotation_oversplit`) above the paper's mean of products, and it is
   undefined for a negative ISCL (Fibroblast.Nspl2, Tcell.Nspl2 in both tools), so averaging the defined labels
   drops exactly the worst ones. The product of label means equals the mean of products when ISCL is constant
   (`annotation`) and is 0.016 lower on `annotation_oversplit`. Given the formula, the ISCG space adds +0.030 /
   +0.034 to the mean of products and the classifier moves it by 0.002 (over-split, 0.491403 vs 0.493365). Ours
   follows the paper (mean over labels of the per-label product).
5. **Headline weighting.** Up to 0.6.x ours `iscg` / `iscl` averaged over pseudobulks (above); sTE reports per
   label only and the paper averages over labels. **Closed in 0.7.0**: label means, checked against
   `mean_over_types` in `test_headline_is_the_mean_over_labels`.
6. **HVGs and library size.** The same 2000 genes and gene-sum library sizes for both tools by construction; in
   normal use ours scores every gene in `adata` and scTypeEval its HVG list with all-gene library sizes.

## Tests (`test_isc_yerly2022.py`, offline, ~2 s)

- `test_inputs_parse`: shapes, columns, the scored-label filter, 25 + 25 cells per half per sample.
- `test_iscg_and_legacy_silhouette_match_sctypeeval` (both columns): per-label ISCG and legacy `per_type` within
  1e-12 of sTE `reduction = FALSE`; the same two labels unscored.
- `test_iscl_matches_sctypeeval_on_its_rcm_distances` (both columns): `_rcm_dissimilarity` replaced by sTE's
  matrix; per-label ISCL within 1e-12.
- `test_rcm_classifier_gap_vs_singler`: our RCM matrix equals sTE's on `annotation` and differs on exactly 30 of
  1,891 pairs on `annotation_oversplit`, all between halves of one split label.
- `test_iscg_default_30_pcs_matches_sctypeeval` (both columns): default ISCG within 2e-6 of sTE's 30-PC ISCG per
  label and for the label mean; positive control: every label's 30-PC and `reduction = FALSE` references differ
  by > 0.01.
- `test_headline_is_the_mean_over_labels` (both columns): `iscg` = sTE `mean_over_types`; `iscg_profile_mean` =
  the pre-0.7.0 headline.
- `test_rcm_markers_vs_singler`: `rcm_genes="markers"` differs from SingleR on 0 / 1,128 and 28 / 1,891 pairs.

Each of four one-line mutations of `eval_metrics.py` (now `isc.py`) fails at least one test: library-size scale 1e6, 2-label `b`
x 0.999, one-sided RCM hits counted as matches, RCM rows reversed.

## Regenerate

```bash
# 1. input (needs SeuratObject, scran; not scTypeEval)
Rscript benchmarks/parity/isc_yerly2022/make_input.R Yerly_2022_35986012.rds black_list.rda <work>/input
# black_list.rda: scTypeEval's data/black_list.rda (same md5 in 1.0.0 and 1.1.2), e.g.
# https://github.com/carmonalab/scTypeEval/raw/ca24658e13da5b7d776b4d69a962fac8a7c280ce/data/black_list.rda
# 2. reference (scTypeEval 1.0.0 or 1.1.2, SingleR, jsonlite)
Rscript benchmarks/parity/isc_yerly2022/run_sctypeeval.R <work>/input tests/goldens/isc_yerly2022/expected.json
# 3. golden files + our numbers + the tables above (checks the npz/csv round trip and the 5 MB cap)
PYTHONPATH=src python benchmarks/parity/isc_yerly2022/make_golden.py <work>/input tests/goldens/isc_yerly2022 <work>/ours.json
```

All three steps are deterministic: on 2026-11-03 step 1 reproduced the first run's files byte for byte (plus the
new column), step 2 wrote byte-identical JSON twice, and step 3 rewrote `counts.npz` byte-identical.
