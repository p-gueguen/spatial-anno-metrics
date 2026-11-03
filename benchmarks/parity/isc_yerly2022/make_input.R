# Build the ISC parity input (tests/goldens/isc_yerly2022) from the Yerly et al. 2022 BCC object.
#
# Input : Yerly_2022_35986012.rds, Zenodo 10.5281/zenodo.18921437 (CC BY 4.0),
#         md5 9d0768bc9f0b9f6240ae499028970455
# Output: <out_dir>/counts.mtx (genes x cells, the 2000 HVGs only), genes.txt, obs.csv, hvg_rank.txt
#         -> make_golden.py turns these into tests/goldens/isc_yerly2022/{counts.npz,obs.csv}
#
# Columns and exclusion follow the authors' benchmark config (carmonalab/ISC_benchmark_reproducibility
# 8b01161, data_processing/config/core_datasets.yaml, id BCC_Yerly): ident "annotation",
# sample "patient_bcc", "Cancer_cells" excluded.
# HVGs follow the scTypeEval README recipe (run_hvg defaults): Log1p (counts / library size * 1e4 on all
# genes), README black list (TCR + Immunoglobulins + Ygenes) removed, scran::modelGeneVar blocked by sample
# with equiweight = TRUE, scran::getTopHVGs(n = 2000), on the cells scTypeEval keeps (min_cells = 10 per
# sample-type, min_samples = 5 per type). Calls scran directly: scTypeEval itself is not needed here.
# <black_list.rda> = scTypeEval's data/black_list.rda (data only, loaded, never committed), e.g.
# https://github.com/carmonalab/scTypeEval/raw/ca24658e13da5b7d776b4d69a962fac8a7c280ce/data/black_list.rda
#
# Usage: Rscript make_input.R <Yerly_2022_35986012.rds> <black_list.rda> <out_dir>
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})
a <- commandArgs(TRUE)
stopifnot(length(a) == 3)
out_dir <- a[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(a[1])
counts <- LayerData(obj, assay = "RNA", layer = "counts")
md <- data.frame(patient_bcc = as.character(obj$patient_bcc), annotation = as.character(obj$annotation),
                 row.names = colnames(obj))
rm(obj)
keep <- !is.na(md$annotation) & md$annotation != "Cancer_cells" & !is.na(md$patient_bcc)
md <- md[keep, ]
counts <- counts[, rownames(md)]

# Subsample: at most 50 cells per (patient, label), seed 22 (the authors' processing seed).
# Deliberate edge case: Mast in BCC3 is cut to 8 cells (< min_cells = 10), so a retained type misses one
# sample - exercises the per-group filter and the RCM "type absent from one sample" rule.
set.seed(22)
cap <- 50L
grp <- paste(md$patient_bcc, md$annotation, sep = "|")
picked <- unlist(lapply(split(rownames(md), grp), function(ids) {
  n <- if (identical(md[ids[1], "patient_bcc"], "BCC3") && identical(md[ids[1], "annotation"], "Mast")) 8L else cap
  if (length(ids) > n) sample(ids, n) else ids
}), use.names = FALSE)
picked <- rownames(md)[rownames(md) %in% picked]  # original cell order
md <- md[picked, ]
counts <- counts[, picked]

# Cells scTypeEval keeps for its single-cell slot: (sample, type) groups >= 10 cells, types in >= 5 such samples.
tab <- table(md$annotation, md$patient_bcc)
ok_type <- rownames(tab)[rowSums(tab >= 10) >= 5]
kept <- md$annotation %in% ok_type & tab[cbind(md$annotation, md$patient_bcc)] >= 10

lib <- colSums(counts[, kept])
norm <- counts[, kept]
norm@x <- log1p(norm@x / rep(lib, diff(norm@p)) * 1e4)
e <- new.env(parent = emptyenv())
load(a[2], envir = e)
bl <- unlist(e$black_list[c("TCR", "Immunoglobulins", "Ygenes")], use.names = FALSE)
stopifnot(length(unique(bl)) == 672)
norm <- norm[!rownames(norm) %in% bl, ]
var <- scran::modelGeneVar(norm, block = md$patient_bcc[kept], equiweight = TRUE)
hvg <- scran::getTopHVGs(var, n = 2000)
stopifnot(length(hvg) == 2000)

genes <- rownames(counts)[rownames(counts) %in% hvg]  # original gene order
Matrix::writeMM(counts[genes, ], file.path(out_dir, "counts.mtx"))
writeLines(genes, file.path(out_dir, "genes.txt"))
writeLines(hvg, file.path(out_dir, "hvg_rank.txt"))
write.csv(data.frame(cell_id = rownames(md), patient_bcc = md$patient_bcc, annotation = md$annotation),
          file.path(out_dir, "obs.csv"), row.names = FALSE, quote = FALSE)

cat("cells:", ncol(counts), " kept by scTypeEval filter:", sum(kept), " genes:", length(genes), "\n")
print(table(md$annotation, md$patient_bcc))
cat("scran", as.character(packageVersion("scran")), " SeuratObject", as.character(packageVersion("SeuratObject")),
    " Matrix", as.character(packageVersion("Matrix")), " ", R.version.string, "\n")
