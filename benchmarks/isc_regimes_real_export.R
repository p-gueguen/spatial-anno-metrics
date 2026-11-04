# Export an ISC Zenodo object (10.5281/zenodo.18921437, CC BY 4.0) for benchmarks/isc_regimes_real.py.
#
# Usage: Rscript isc_regimes_real_export.R <rds> <ident_col> <sample_col> <out_dir> [excluded labels, comma-sep]
#                                          [extra obs columns, comma-sep]
# Writes: counts_{i,p,x}.bin (genes x cells CSC of the 2000 HVGs: int32 / int32 / float64), shape.txt,
#         genes.txt, obs.csv (cell_id, sample, label, extra columns, libsize_all), hvg_note.txt.
# HVGs: scTypeEval README recipe (Log1p on the all-gene library size; TCR + Ig + Y black list removed, read from
# scTypeEval's data(black_list); scran::modelGeneVar(block = sample, equiweight = TRUE); getTopHVGs(n = 2000))
# on the cells of scorable groups (>= 10 cells per sample-label, >= 5 samples). Needs SeuratObject, scran and
# scTypeEval (for its black list data only).
# Exports used (md5 of the .rds checked against the Zenodo listing):
#   Yerly_2022_35986012.rds   annotation    patient_bcc  -> yerly_full    excluded Cancer_cells
#   Ganier_2024_38165934.rds  annotation    patient_bcc  -> ganier_full   excluded Cancer_cells  extra scGate_multi
#   Joanito_2022_35773407.rds cell.type     sample.ID    -> joanito_full  excluded iCMS2,iCMS3
#                                                                         extra dataset,sample.origin
# (columns and exclusions from the authors' data_processing/config/core_datasets.yaml)
suppressPackageStartupMessages({library(Matrix); library(SeuratObject)})
a <- commandArgs(TRUE)
rds <- a[1]; ident <- a[2]; samp <- a[3]; out <- a[4]
excl <- if (length(a) >= 5) strsplit(a[5], ",")[[1]] else character(0)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
obj <- readRDS(rds)
cat("class:", class(obj), "\n"); print(dim(obj))
counts <- LayerData(obj, assay = "RNA", layer = "counts")
md <- data.frame(sample = as.character(obj[[samp]][, 1]), label = as.character(obj[[ident]][, 1]),
                 row.names = colnames(obj))
extra <- if (length(a) >= 6) strsplit(a[6], ",")[[1]] else character(0)
for (cc in extra) md[[cc]] <- as.character(obj[[cc]][, 1])
rm(obj); invisible(gc())
keep <- !is.na(md$label) & !is.na(md$sample) & !(md$label %in% excl)
md <- md[keep, ]; counts <- counts[, rownames(md)]
print(table(md$label, md$sample))
tab <- table(md$label, md$sample)
ok_type <- rownames(tab)[rowSums(tab >= 10) >= 5]
kept <- md$label %in% ok_type & tab[cbind(md$label, md$sample)] >= 10
lib <- colSums(counts[, kept])
norm <- counts[, kept]
norm@x <- log1p(norm@x / rep(lib, diff(norm@p)) * 1e4)
e <- new.env(); data("black_list", package = "scTypeEval", envir = e)
bl <- unlist(e$black_list[c("TCR", "Immunoglobulins", "Ygenes")], use.names = FALSE)
norm <- norm[!rownames(norm) %in% bl, ]
var <- scran::modelGeneVar(norm, block = md$sample[kept], equiweight = TRUE)
hvg <- scran::getTopHVGs(var, n = 2000)
# library size over ALL genes is kept as an extra column so the Python side can normalise like scTypeEval
md$libsize_all <- colSums(counts)
X <- as(counts[hvg, ], "CsparseMatrix")
writeBin(as.integer(X@i), file.path(out, "counts_i.bin"), size = 4)
writeBin(as.integer(X@p), file.path(out, "counts_p.bin"), size = 4)
writeBin(as.numeric(X@x), file.path(out, "counts_x.bin"), size = 8)
writeLines(as.character(dim(X)), file.path(out, "shape.txt"))
writeLines(hvg, file.path(out, "genes.txt"))
write.csv(data.frame(cell_id = rownames(md), md), file.path(out, "obs.csv"), row.names = FALSE, quote = TRUE)
writeLines(c(rds, ident, samp, paste(excl, collapse = ","), length(hvg)), file.path(out, "hvg_note.txt"))
cat("done", dim(X), "\n")
