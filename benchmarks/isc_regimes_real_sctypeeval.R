# scTypeEval 1.0.0 ISCG / ISCL on the perturbations of benchmarks/isc_regimes_real.py (Yerly 2022, seed 0).
#
# scTypeEval is GPL-3: this script only calls the package and contains none of its code. Recipe = the README
# "Minimal code to compute the recommended metrics" as in benchmarks/parity/isc_yerly2022/run_sctypeeval.R:
#   ISCL = silhouette | recip_classif:Match (SingleR), ISCG = 2label_silhouette | Pseudobulk:Cosine (30 PCs),
# the HVG step replaced by the 2000 genes of the export. scTypeEval reports scores only (no regime call).
#
# Usage: Rscript isc_regimes_real_sctypeeval.R <dir written by isc_regimes_real.py --sctypeeval-input>
# Writes <dir>/sctypeeval_scores.csv: column, label, iscg, iscl (one row per scored label and perturbation).
suppressPackageStartupMessages({
  library(Matrix)
  library(scTypeEval)
})
d <- commandArgs(TRUE)[1]
read_csc <- function(prefix) {
  shp <- as.integer(readLines(file.path(d, "shape.txt")))
  i <- readBin(file.path(d, paste0(prefix, "_i.bin")), "integer", n = 1e9, size = 4)
  p <- readBin(file.path(d, paste0(prefix, "_p.bin")), "integer", n = 1e9, size = 4)
  x <- readBin(file.path(d, paste0(prefix, "_x.bin")), "double", n = 1e9, size = 8)
  sparseMatrix(i = i + 1L, p = p, x = x, dims = shp)
}
lab <- read.csv(file.path(d, "labels.csv"), colClasses = "character")
genes <- readLines(file.path(d, "genes.txt"))
counts <- read_csc("counts")
dimnames(counts) <- list(genes, lab$cell_id)
rownames(lab) <- lab$cell_id

score <- function(mat, ident) {
  sc <- create_scTypeEval(matrix = mat, metadata = lab)
  sc <- run_processing_data(sc, ident = ident, sample = "sample", min_samples = 5, min_cells = 10)
  sc <- add_gene_list(sc, gene_list = list(HVG = rownames(mat)))
  set.seed(22)
  sc <- run_pca(sc, gene_list = "HVG", ndim = 30)
  sc <- run_dissimilarity(sc, method = "recip_classif:Match", reduction = FALSE, reciprocal_classifier = "SingleR")
  sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = TRUE)
  l <- get_consistency(sc, dissimilarity_slot = "recip_classif:Match", consistency_metric = "silhouette")
  g <- get_consistency(sc, dissimilarity_slot = "Pseudobulk:Cosine", consistency_metric = "2label_silhouette")
  labels <- sort(unique(lab[[ident]]))
  purged <- setNames(labels, scTypeEval:::purge_label(labels))   # scTypeEval's own label clean-up, called
  data.frame(label = purged[g$celltype], iscg = g$measure,
             iscl = l$measure[match(g$celltype, l$celltype)])
}

out <- list()
for (col in grep("^p[0-9]+$", colnames(lab), value = TRUE)) {
  t0 <- Sys.time()
  out[[col]] <- cbind(column = col, score(counts, col))
  cat(col, format(Sys.time() - t0), "\n")
}
for (sd in c("0.5", "1.0")) {
  m <- read_csc(paste0("batch", sd))
  dimnames(m) <- dimnames(counts)
  lab$orig <- lab$p0
  out[[sd]] <- cbind(column = paste0("batch", sd), score(m, "orig"))
}
write.csv(do.call(rbind, out), file.path(d, "sctypeeval_scores.csv"), row.names = FALSE)
