# SpatialQM signal-QC parity golden - Xenium breast crop

**Status (2026-11-03): complete.** Inputs, generators and `expected.json` (SpatialQM's outputs on the crop)
are committed. `test_signalqc_spatialqm_files.py` checks that the files parse, agree with each other, and that
`expected.json` was computed on these inputs (`meta.input`). No metric is compared yet: the owning work
package adds that test, using [Ours vs reference](#ours-vs-reference-on-the-crop) for what should match.

## Upstream

- **SpatialQM**, Center for Spatial OMICs (St. Jude), <https://github.com/Center-for-Spatial-OMICs/SpatialQM>,
  commit `36e1a59d4ca3f8f360a1d826f5e999d9f40c6760` (2026-05-06), version 0.1.0.9000, **MIT**.
  Published with Plummer et al., *Nat Biotechnol* 2025, [10.1038/s41587-025-02811-9](https://doi.org/10.1038/s41587-025-02811-9).
  Installed with pak into the dev R library (R 4.6.0, Bioconductor 3.23). `meta.spatialqm_sha` is the
  installed `RemoteSha`; `make_expected.R` stops unless the source checkout it parses is at that commit.
- **Which code runs, checked on the installed copy.** `R/utils.R` (last edited 2026-05-06, "QC fixes, SPE
  coercion") and `R/utils_update_final.R` (2025-02-20) define the same exported functions. With no `Collate`
  field, `R CMD INSTALL` concatenates `R/*.R` in C-locale order (`DependencyPackages.R`, `coercion-input.R`,
  `utils.R`, `utils_update_final.R`) and the last definition wins. `make_expected.R` compares formals and body
  of every function it calls with both parsed definitions: all 11 equal the `utils_update_final.R` version
  (`meta.function_source`) and differ from the `utils.R` one. The reference is therefore SpatialQM as it
  installs from GitHub; the 2026-05-06 `utils.R` code (Seurat branches, `coerce_to_seurat`) never runs.
- Metric code it calls: `coop::sparsity` (BSD-2), `BioQC::entropy` (GPL-3: run it, never copy it),
  `Voyager::runMoransI` on a `SpatialFeatureExperiment::findSpatialNeighbors` graph (spdep),
  `scater::logNormCounts`. scuttle 1.22.0 warns that the size-factor and normalisation functions behind it are
  deprecated for scrapper, so a later Bioconductor may change `getMorans`.

## Function map (installed code)

Every call except `getMorans` uses the file branch (`seu_obj = NULL`): the installed Seurat branch of
`getMeanSignalRatio`, `getMaxRatio`, `getSparsity`, `getEntropy`, `getTxPerArea`, `getTxPerNuc` and
`getCellTxFraction` is empty. Paths are the crop folder written by `make_inputs.py outs`. Each scalar was
re-derived from the fixtures with independent Python code (<= 4.2e-15 relative) and `getMorans` too (2.1e-14
absolute, NaN at the same 7 features), so the right-hand column is checked, not only read from the code.

| Metric | Exact call | What it computes on these inputs |
|---|---|---|
| Moran's I, gene vs control | `getMorans(seu)`: assay `RNA` = 313 Gene Expression, assay `ControlProbe` = 28 Negative Control Probe, reduction `tissue` = centroids, meta `sample_id`, `platform` | Per group: drop cells with 0 counts in the group (genes: none dropped; probes: 150 of 2000 kept), `scater::logNormCounts` (library-size factors centred at 1, log2(x/sf + 1)), 20 nearest neighbours weighted 1/d and row-standardised (`style = "W"`), `Voyager::runMoransI`. One I per feature, NaN for the 7 features with no count in the kept cells (CRHBP, 6 probes). Seurat renames `_` to `-` (`NegControlProbe-00042`); `expected.json` maps the names back. |
| Mean signal ratio | `getMeanSignalRatio(expMat = <csv>, platform = "CosMx")`, `<csv>` = the same counts, cells x features | The Xenium branch errors ([finding 4](#findings)), so the CosMx branch of the same function reads the counts with `fread`. Rows matching `Neg*\|SystemControl*\|Blank*\|BLANK*\|Unassigned*` = noise (220: 20 NegControlProbe, 41 NegControlCodeword, 159 BLANK), the rest = genes (321, including the 8 `antisense_*` probes). `mean(log10(rowMeans(genes + 0.1)) - log10(rowMeans(noise + 0.1)))`: R recycles the 220 noise means against the 321 gene means (warning suppressed). One log10 scalar. |
| Max ratio | `getMaxRatio(expMat = <csv>, platform = "CosMx")` | `log10(max(gene rowMeans)) - log10(mean(noise rowMeans))`, same split. |
| Sparsity | `getSparsity(expMat = "cell_feature_matrix/matrix.mtx.gz", platform = "Xenium")` | `coop::sparsity(as.matrix(exp))`: fraction of zero entries over all 541 features x 2000 cells. |
| Entropy | `getEntropy(expMat = ..., platform = "Xenium")` | `BioQC::entropy(as.matrix(exp))`: Shannon entropy in bits of the histogram of count values over all 1,082,000 entries (64 distinct values, dominated by the zeros). Not per cell, not over genes. |
| Transcripts per area | `getTxPerArea(cellSegMeta = "cells.csv.gz", platform = "Xenium")` | `mean(total_counts / cell_area)` over cells. Xenium `total_counts` = all 541 features. |
| Transcripts per cell (TPC, #21) | `getTxPerCell(expMat = ..., platform = "Xenium")` | Mean over cells of the column sums over the 321 rows not matching `Unassigned*\|NegControl*\|BLANK*\|Blank*\|SystemControl*\|NegPrb*` (`antisense_*` kept). |
| Fraction of transcripts in cells (FTC, #21) | `getCellTxFraction(tx_file = "transcripts.csv.gz", platform = "Xenium")` | `(n - sum(cell_id == "UNASSIGNED")) / n` over every transcript in the file, no qv filter. |
| False discovery rate (FDR, #21) | `getGlobalFDR(tx_file = "transcripts.csv.gz", platform = "Xenium")` | Negatives = transcripts whose target matches `^(Neg\|Blank\|BLANK)` (case-insensitive; misses `antisense_*`). `expNeg / (expGenes + expNeg) * (numGenes / numNeg) / 100`, numGenes / numNeg = distinct targets with >= 1 transcript; every transcript, no qv filter. It prints the negative transcript count (1223) as "Number of negative control probes". |
| Transcripts per nucleus | `getTxPerNuc(tx_file = ..., platform = "Xenium")` | Transcripts with `cell_id != "UNASSIGNED"`, `overlaps_nucleus == 1` and a target not matching the noise regex, divided by the number of distinct `cell_id` values in the file (2104: 2000 crop cells, 103 cells centred outside the crop, `"UNASSIGNED"`). |
| Cell count | `getNcells(expMat = ..., platform = "Xenium")` | `ncol` of the matrix. |

## Data

- **10x Xenium FFPE Human Breast Cancer, Replicate 1** (Xenium Onboard Analysis `Xenium-1.0.1`, prototype
  instrument, panel PD_260 = 280 predesigned + 33 custom genes), Janesick et al., *Nat Commun* 14:8353
  (2023), [10.1038/s41467-023-43458-x](https://doi.org/10.1038/s41467-023-43458-x); GEO
  [GSE243280](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243280) (GSM7780153 holds the same
  outputs as one 9.9 GB zip).
- **Licence: CC BY 4.0.** The 10x dataset page states "This dataset is licensed under the Creative
  Commons Attribution 4.0 International (CC BY 4.0) license" (as indexed by web search on 2026-11-03;
  the page sits behind a lead-capture form and could not be fetched directly). Attribution: 10x Genomics
  and Janesick et al. 2023.
- Downloaded 2026-11-03 from `https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_<file>`:

  | file | bytes | md5 |
  |---|---|---|
  | `cell_feature_matrix.h5` | 12,148,885 | `3bfd29015fa8bacb821d9c9ee86bedaa` |
  | `cells.csv.gz` | 7,930,261 | `5bbac58bb5edd91bc577603587b6d04e` |
  | `transcripts.parquet` | 669,769,801 | `54fe44356d0284a39cffd6658a167f06` |
  | `metrics_summary.csv` | 726 | `1368d02d9877c56c2d3167b96f48fb8c` |
  | `experiment.xenium` | 1,383 | `dc0c520250cfc8c38f182c824b3dd57f` |

- Section: 167,780 cells; 541 features = 313 Gene Expression, 28 Negative Control Probe (20
  `NegControlProbe_*` + 8 `antisense_*`), 41 Negative Control Codeword, 159 Blank Codeword;
  42,638,083 transcripts (the parquet footer claims 1,000,000 rows; its 43 row groups hold the full count).
  Unassigned transcripts have `cell_id == -1` (Xenium <= 1.2).
- Parsing check against 10x `metrics_summary.csv`: q >= 20 gene transcripts assigned to cells 0.93191
  (10x 0.93207); median gene counts per cell 164 (10x 166). That prototype-run summary does not
  reproduce exactly from the released files; a misread `cell_id` would give 1.0.
- Full-section AnnData for the tutorial and scale test (not in the repo):
  `/srv/GT/analysis/pgueguen/sam_dev/data/xenium_breast/section.h5ad`, 167,780 cells x 541 features,
  `X` = `layers["counts"]`, `obsm["spatial"]` = centroids (um), `var["control"]` = feature_type is not
  Gene Expression, `obs` = the `cells.csv.gz` columns. Written with anndata 0.11.4 and read back with 0.13.4.

## Crop

- **Rule:** centres on a 250 um grid. For each centre, take its 2000 Chebyshev-nearest cells and keep the
  window only if its half-side is <= 400 um and every quadrant holds >= 15 % of the cells. Of those, pick
  the most even epithelial / stromal / immune mix: Shannon entropy of cell classes, where a cell's class is
  the argmax of its summed EPCAM+KRT7+KRT8+CDH1 / LUM+POSTN+PDGFRA+DPT / PTPRC+CD3E+CD68 counts.
- **Window:** centre (252.13251609802245, 2001.3671555936337) um, half-side h = 240.86916015446195 um,
  i.e. |x - cx| < h and |y - cy| < h (x 11.26-493.00, y 1760.50-2242.24 um). This is the ERBB2-high
  tumour edge: a cell in 94 % of its 20 um bins, classes 30 / 32 / 37 % epithelial / stromal / immune.
- **Cells:** the 2000 cells with centroid inside. **Transcripts:** all 321,637 inside the window, any
  qv: 315,916 in crop cells, 5,656 in cells centred outside, 65 unassigned; 80.0 % have qv >= 20 (the
  count matrix holds only those, 256,637 counts).
- `make_inputs.py crop` rerun on 2026-11-03 reproduced all three fixture files byte for byte (md5 below).

## Files

| file | content | md5 |
|---|---|---|
| `counts.npz` | CSR cells x features, 2000 x 541 (h5 counts): `data` int32, `indices` uint16, `indptr` int32, `shape` int64[2], `cell_id` int32, `feature_id` / `feature_name` / `feature_type` str. `np.load(..., allow_pickle=False)`. | `527d72bcfbf6a8a36648382a1eb0fb2c` |
| `cells.csv` | crop rows of `cells.csv.gz`, same order as `counts.npz` rows (um, um^2); `total_counts` = row sum over all 541 features | `dc4f84f763fce622e087646c0d512b7d` |
| `transcripts.npz` | per transcript: `feature_idx` uint16 (into the `counts.npz` features), `cell_id` int32 (-1 = unassigned), `overlaps_nucleus` uint8, `qv` float32; no coordinates (no metric here needs them) | `e8b78016f2eeb727ac4f23898c6d505f` |
| `expected.json` | one key per SpatialQM function in the map above (scalar); `getMorans` = `{feature, type, I}` lists, 341 = 313 `Gene` + 28 `Control`, feature names as in `counts.npz`, NaN as null; `meta` = SpatialQM sha and version, R, `function_source`, `errors` (the Xenium-branch error of the two ratio functions), `calls`, `control_regex_hits_by_feature_type`, package versions, `input` counts, date. Floats at 15 significant digits. | `0b041ab9ac9d2d3de4c7776cafaddb80` |

## Ours vs reference on the crop

Ours = `spatial_anno_metrics.signal_qc` at `f71512f` with default arguments, on an AnnData built from these
files (X = counts of all 541 features, `var["control"]` = the 228 non-gene features, `obs["cell_area"]` and
`obs["total_counts"]` from `cells.csv`; written by `make_inputs.py outs` as `ours.json`). Python 3.10
(squidpy 1.6.5) and 3.12 (squidpy 1.8.4) agree to 4.2e-13 relative over 423 values; the comparator flags a
planted 1e-6 relative change. The Moran gaps are split by walking SpatialQM's recipe to ours one choice at a
time in Python: the walk starts on `expected.json` (2.1e-14) and ends on `ours.json` (3.9e-14).

| Metric | SpatialQM (`expected.json`) | Ours: function, key, value | Gap and why |
|---|---|---|---|
| Moran's I, gene median | 0.06670, NaN-free median of 312 genes (CRHBP NaN) | `moran_signal`, `median_gene_moran`: **NaN** (NaN-free 0.10637) | +0.0397 = +0.0124 from raw counts instead of log-normalised, +0.0272 from 6 binary neighbours instead of 20 weighted 1/d. SpatialQM's cell filter changes nothing here: every crop cell has a gene count. Per-gene Spearman 0.868. The NaN is [finding 1](#findings). |
| Moran's I, control median | -0.00771, 22 of the 28 probes (6 NaN), on the 150 cells with a probe count | `moran_signal`, `median_control_moran`: **NaN** over all 228 non-gene features (NaN-free -0.00058; the same 28 probes -0.00092) | The cell filter. Without autocorrelation E[I] = -1/(n - 1): -0.0067 at n = 150, -0.0005 at n = 2000. Keeping all 2000 cells moves SpatialQM's recipe (raw counts) from -0.00816 to -0.00097. Both sit on their null. Ours also pools 41 codewords and 159 BLANK with the probes. |
| Mean signal ratio | 0.49702 (log10) | `signal_to_noise`, `median_mean_signal_ratio`: 27.956 (log2) | Different quantities. SpatialQM: mean of 321 per-feature log10 ratios, +0.1 pseudocount, recycled noise means. Ours: median over cells of log2((gene mean + 1e-9) / (control mean + 1e-9)); 88.8 % of crop cells have no control count, so the median is the eps floor ([finding 2](#findings)). |
| Max ratio | 4.33268 (log10) | `signal_to_noise`, `median_max_signal_ratio`: 32.705 (log2) | Different quantities: SpatialQM = highest per-feature mean over the mean noise mean; ours = per-cell max count over control mean, eps floor again. |
| Sparsity | 0.902546210720887 | `sparsity`, `zero_fraction`: 0.9025462107208873 | **Equal** (3e-16 relative): same definition, all 541 features. |
| Entropy | 0.675553 bits | `detection_entropy`, `median_detection_entropy`: 0.94072 | Different quantities: SpatialQM = entropy of the count-value histogram of the whole matrix; ours = per-cell entropy over detected features / ln(n detected), median over cells. |
| Transcripts per area | 1.15712454692 (mean) | `tx_per_area`, `median_tx_per_area`: 1.05603 (median) | Median vs mean of the same per-cell `total_counts / cell_area`; the mean of our per-cell values is 1.157124546919905 (4e-15 relative). |
| TPC (#21) | 128.261 | none | Counts over 313 genes + 8 `antisense_*` probes; genes only 128.1905. |
| FTC (#21) | 0.999798 = (321,637 - 65) / 321,637 | none | All qv. 1 on Xenium <= 1.2 ids ([finding 6](#findings)). |
| FDR (#21) | 5.92514e-05 = 1223 / 321,637 x 321 / 206 / 100 | none | Negatives = NegControlProbe + NegControlCodeword + BLANK transcripts, any qv; 321 / 206 = distinct genes / negatives detected. Before the `/ 100` it is 0.005925 (0.59 %). With `antisense_*` as negatives: 6.64377e-05. |
| Transcripts per nucleus | 67.7262 = 142,496 / 2104 | none (not in our API) | Per crop cell: 71.248 ([finding 7](#findings)). |
| Cell count | 2000 | `adata.n_obs`: 2000 | Equal. |

Moran's I per feature along the walk (SpatialQM's recipe, then one change per column; the last column is ours):

| feature | SpatialQM | raw counts | + all 2000 cells | + 6 binary neighbours = ours |
|---|---|---|---|---|
| ERBB2 | 0.36791 | 0.52656 | 0.52656 | 0.60072 |
| KRT7 | 0.59840 | 0.56262 | 0.56262 | 0.63854 |
| EPCAM | 0.59233 | 0.51338 | 0.51338 | 0.59123 |
| LUM | 0.28413 | 0.18736 | 0.18736 | 0.28363 |
| POSTN | 0.21941 | 0.14406 | 0.14406 | 0.22985 |
| PTPRC | 0.33211 | 0.24703 | 0.24703 | 0.27316 |
| CD68 | 0.15666 | 0.10553 | 0.10553 | 0.16736 |
| ACTA2 | 0.20572 | 0.19259 | 0.19259 | 0.26412 |
| NegControlProbe_00042 | -0.00819 | -0.00819 | -0.00053 | -0.00042 |
| antisense_PROKR2 | 0.05009 | 0.03114 | 0.01007 | 0.01630 |
| median, genes | 0.06670 | 0.07913 | 0.07913 | 0.10637 |
| median, 28 probes | -0.00771 | -0.00816 | -0.00097 | -0.00092 |

## Findings

For the owning work package; `src/` is untouched here.

On ours:

1. **`moran_signal` summaries are NaN on real data.** squidpy returns NaN for the 139 features with no count
   in the crop (1 gene, CRHBP, plus 6 / 27 / 105 probes / codewords / BLANK); the zero-count set and the NaN
   set are identical. `np.median` propagates it to `median_gene_moran`, `median_control_moran` and
   `separation`, and `frac_genes_above_control_p95` reads 0.0 because every comparison with a NaN quantile is
   false. NaN-free: genes 0.10637, controls -0.00058, separation 0.10695, fraction 0.929. SpatialQM also
   returns NaN for its 7 zero-count features, but its only summary, the `plotMorans` boxplot, drops them.
2. **SNR sits on the eps floor** (roadmap 8g). 88.8 % of crop cells have no count in any of the 228 control
   features (92.5 % none in the 28 probes). For them the ratio is log2(gene mean) + log2(1 / 1e-9) = -1.782 +
   29.897 at the median (28.115), so the reported median 27.956 is the floor. Over the 224 cells with a
   control count the median is 6.491; pooled over the crop, log2(mean gene count / mean control count) = 9.511.
3. **Three docstrings claim ports of different quantities.** `signal_to_noise` (getMeanSignalRatio /
   getMaxRatio), `detection_entropy` (getEntropy) and `tx_per_area` (median; getTxPerArea is the mean). Only
   `sparsity` reproduces SpatialQM. `moran_signal` is the same statistic with another normalisation, graph and
   control set (table above).

On SpatialQM (installed code at the sha above):

4. **`getMeanSignalRatio` and `getMaxRatio` cannot run on Xenium data.** Their mtx and h5 branches hand a
   `Matrix` object to `base::rowMeans` (the namespace imports only `colSums` and `readMM` from Matrix), which
   stops with `'x' must be an array of at least two dimensions` (`meta.errors` for the crop's mtx; the
   full-section h5 gives the same error); their Seurat branch is empty.
   The shadowed `utils.R` has the same call in its file branch. The reference uses their CosMx branch.
5. **The control regexes miss this panel's 8 `antisense_*` negative-control probes.** All three (ratio, TPC,
   FDR) hit 20 of the 28 probes, 41 of 41 codewords, 159 of 159 BLANK and 0 of 313 genes
   (`meta.control_regex_hits_by_feature_type`). The 8 count as genes in `getMeanSignalRatio`, `getMaxRatio`,
   `getTxPerCell`, `getGlobalFDR` and `getTxPerNuc`, and as controls in `getMorans` (assay by feature type).
   Counted as controls: `getMaxRatio` 4.333 -> 4.001, `getMeanSignalRatio` 0.4970 -> 0.5072, `getTxPerCell`
   128.261 -> 128.1905, `getGlobalFDR` 5.925e-05 -> 6.644e-05. Flag controls by `feature_type`, not by name.
6. `getCellTxFraction` returns 1 on Xenium <= 1.2 transcripts (unassigned = `cell_id == -1`, never
   `"UNASSIGNED"`): 1 instead of 0.99980 on this crop. `make_inputs.py outs` writes the >= 1.3 convention.
7. `getMeanSignalRatio` depends on feature order through the recycling: 0.49696-0.49725 over 1000 random
   orders of the noise rows (0.49702 in file order). `getTxPerNuc` divides by the 2104 distinct `cell_id`
   values, which include `"UNASSIGNED"` and 103 cells centred outside the crop. `getGlobalFDR` returns the
   fraction divided by 100.

## Regenerate

Generated 2026-11-03. From the repo root:

```bash
PY=/srv/GT/analysis/pgueguen/sam_dev/venv310/bin/python   # anndata 0.11: h5ad readable by 0.11 and 0.13
G=benchmarks/parity/signalqc_spatialqm
PYTHONPATH=src $PY $G/make_inputs.py section   # 1. section.h5ad next to the downloaded files
PYTHONPATH=src $PY $G/make_inputs.py crop      # 2. counts.npz, cells.csv, transcripts.npz (this folder)
PYTHONPATH=src $PY $G/make_inputs.py outs      # 3. Xenium-style crop folder + ours.json, from the fixtures only
# 4. SpatialQM into the dev R library (R 4.6; everything else comes from the system library). BioQC 1.40.0
#    came from BiocManager; then, with pak (NMI was archived on CRAN on 2026-05-18):
#    for (r in c("spatialreg", "fastglm", "url::https://cran.r-project.org/src/contrib/Archive/NMI/NMI_2.0.tar.gz",
#                "SpatialFeatureExperiment", "Voyager", "Nanostring-Biostats/InSituType",
#                "Center-for-Spatial-OMICs/SpatialQM@36e1a59"))
#      pak::pkg_install(r, lib = "/srv/GT/analysis/pgueguen/sam_dev/Rlib", dependencies = FALSE, upgrade = FALSE)
#    git clone https://github.com/Center-for-Spatial-OMICs/SpatialQM <src>; git -C <src> checkout 36e1a59d4ca3f8f360a1d826f5e999d9f40c6760
bash -lc "module load Dev/R/4.6.0 && R_LIBS_USER=/srv/GT/analysis/pgueguen/sam_dev/Rlib Rscript $G/make_expected.R \
  /srv/GT/analysis/pgueguen/sam_dev/runs/wf2b-signal/crop_outs tests/goldens/signalqc_spatialqm/expected.json \
  /srv/GT/analysis/pgueguen/sam_dev/data/spatialqm_src/SpatialQM"
```

The re-derivation and the Moran walk are not in the repo:
`/srv/GT/analysis/pgueguen/sam_dev/runs/wf2b-signal/verify_ref.py`,
`/srv/GT/analysis/pgueguen/sam_dev/runs/wf2b-signal/moran_decompose.py`,
`/srv/GT/analysis/pgueguen/sam_dev/runs/wf2b-signal/ours_findings.py`.
