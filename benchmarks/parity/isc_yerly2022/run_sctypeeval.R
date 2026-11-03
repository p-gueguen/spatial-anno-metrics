# Reference ISC on the golden input with scTypeEval -> tests/goldens/isc_yerly2022/expected.json
#
# NOT YET RUN: installing scTypeEval was blocked on 2026-11-03 (see SOURCE.md). Written against the
# scTypeEval 1.1.2 source (github.com/carmonalab/scTypeEval @ ca24658e13da5b7d776b4d69a962fac8a7c280ce, GPL-3);
# this script only calls the package, it contains none of its code.
#
# Recipe = scTypeEval README "Minimal code to compute the recommended metrics", with the HVG step replaced by
# the 2000 genes already in the input (make_input.R picked them with the same scran call):
#   local  ISC (ISCL) = silhouette        | recip_classif:Match (SingleR, always reduction = FALSE)
#   global ISC (ISCG) = 2label_silhouette | Pseudobulk:Cosine   (reduction = TRUE: 30 PCs, the README default)
# plus ISCG with reduction = FALSE (cosine on the log-normalised HVG pseudobulks, our current space).
#
# Usage: Rscript run_sctypeeval.R <make_input_out_dir> <expected.json>
suppressPackageStartupMessages({
  library(Matrix)
  library(scTypeEval)
})
a <- commandArgs(TRUE)
stopifnot(length(a) == 2)

counts <- as(readMM(file.path(a[1], "counts.mtx")), "CsparseMatrix")
obs <- read.csv(file.path(a[1], "obs.csv"), colClasses = "character")
dimnames(counts) <- list(readLines(file.path(a[1], "genes.txt")), obs$cell_id)
rownames(obs) <- obs$cell_id

sc <- create_scTypeEval(matrix = counts, metadata = obs)
sc <- run_processing_data(sc, ident = "annotation", sample = "patient_bcc", min_samples = 5, min_cells = 10)
sc <- add_gene_list(sc, gene_list = list(HVG = rownames(counts)))
set.seed(22)  # irlba start vector
sc <- run_pca(sc, gene_list = "HVG", aggregation = "pseudobulk", ndim = 30)

slots <- list()
consistency <- function(sc, slot) {
  get_consistency(sc, dissimilarity_slot = slot, consistency_metric = c("silhouette", "2label_silhouette"))
}
dmat <- function(sc, slot) {
  m <- as.matrix(sc@dissimilarity[[slot]]@dissimilarity)
  list(labels = rownames(m), matrix = unname(m))
}
sc <- run_dissimilarity(sc, method = "recip_classif:Match", reduction = FALSE, reciprocal_classifier = "SingleR")
sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = TRUE)
cons <- rbind(consistency(sc, "recip_classif:Match"), consistency(sc, "Pseudobulk:Cosine"))
slots[["recip_classif:Match"]] <- dmat(sc, "recip_classif:Match")
slots[["Pseudobulk:Cosine"]] <- dmat(sc, "Pseudobulk:Cosine")
sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = FALSE)
nopca <- consistency(sc, "Pseudobulk:Cosine")
nopca$dissimilarity_method <- "Pseudobulk:Cosine [reduction=FALSE]"
cons <- rbind(cons, nopca)
slots[["Pseudobulk:Cosine [reduction=FALSE]"]] <- dmat(sc, "Pseudobulk:Cosine")

key <- paste(cons$consistency_metric, cons$dissimilarity_method, sep = " | ")
local <- "silhouette | recip_classif:Match"
global <- "2label_silhouette | Pseudobulk:Cosine"
labels <- sort(unique(obs$annotation))
purged <- setNames(gsub(",", "", gsub(" |_|[+]|-", ".", labels)), labels)  # scTypeEval:::purge_label rule
per_type <- lapply(labels, function(l) {
  v <- as.list(setNames(cons$measure[cons$celltype == purged[[l]]], key[cons$celltype == purged[[l]]]))
  p <- v[[local]] * v[[global]]
  v$product <- p
  v$geometric_mean <- if (is.na(p) || v[[local]] < 0 || v[[global]] < 0) NA_real_ else sqrt(p)
  v
})
names(per_type) <- labels
cols <- names(per_type[[1]])
mean_over_types <- lapply(setNames(cols, cols), function(k) mean(vapply(per_type, `[[`, numeric(1), k), na.rm = TRUE))

out <- list(
  tool = "scTypeEval", version = as.character(packageVersion("scTypeEval")),
  commit = packageDescription("scTypeEval")$RemoteSha %||% NA_character_,
  singler = as.character(packageVersion("SingleR")),
  r = R.version.string, generated = format(Sys.Date()),
  names = list(ISCL = local, ISCG = global, ISCG_no_pca = paste0("2label_silhouette | Pseudobulk:Cosine",
                                                                  " [reduction=FALSE]")),
  label_map = as.list(purged),
  per_type = per_type,
  mean_over_types = mean_over_types,
  dissimilarity = slots,
  note = paste("NA = label not scored (< min_samples samples with >= min_cells cells). recip_classif:Match",
               "stores 0.5 for a mutual hit between labels not shared by both samples; get_consistency reads",
               "it as 1. product = ISCL * ISCG (paper Fig. 5, label-transfer task); geometric_mean =",
               "sqrt(product) (README), NA when either is negative.")
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, na = "null", pretty = TRUE), a[2])
print(cons)
