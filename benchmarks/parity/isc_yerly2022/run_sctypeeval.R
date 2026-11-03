# Reference ISC on the golden input with scTypeEval -> tests/goldens/isc_yerly2022/expected.json
#
# Generated with scTypeEval 1.0.0 from Bioconductor 3.23 (git 05e2c7a = GitHub tag v1.0.0; GPL-3: this script
# only calls the package, it contains none of its code). Its R code equals GitHub 1.1.2 (ca24658) except
# run_pca() / custom_prcomp(), which gained aggregation / center / scale. arguments with the 1.0.0 behaviour as
# defaults (SOURCE.md, "Upstream reference tool"); the calls below run unchanged on both versions.
#
# Recipe = scTypeEval README "Minimal code to compute the recommended metrics", with the HVG step replaced by
# the 2000 genes already in the input (make_input.R picked them with the same scran call):
#   local  ISC (ISCL) = silhouette        | recip_classif:Match (SingleR, always reduction = FALSE)
#   global ISC (ISCG) = 2label_silhouette | Pseudobulk:Cosine   (reduction = TRUE: 30 PCs, the README default)
# plus ISCG with reduction = FALSE (cosine on the log-normalised HVG pseudobulks, our space).
# Scored for two label columns: annotation (the authors' labels) and annotation_oversplit (Fibroblast and Tcell
# randomly halved per sample, make_input.R).
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

local <- "silhouette | recip_classif:Match"
global <- "2label_silhouette | Pseudobulk:Cosine"
nopca <- "Pseudobulk:Cosine [reduction=FALSE]"

score_column <- function(ident) {
  sc <- create_scTypeEval(matrix = counts, metadata = obs)
  sc <- run_processing_data(sc, ident = ident, sample = "patient_bcc", min_samples = 5, min_cells = 10)
  sc <- add_gene_list(sc, gene_list = list(HVG = rownames(counts)))
  set.seed(22)  # irlba start vectors (single-cell PCA first, then pseudobulk)
  sc <- run_pca(sc, gene_list = "HVG", ndim = 30)
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
  slots <- list("recip_classif:Match" = dmat(sc, "recip_classif:Match"),
                "Pseudobulk:Cosine" = dmat(sc, "Pseudobulk:Cosine"))
  sc <- run_dissimilarity(sc, method = "Pseudobulk:Cosine", reduction = FALSE)
  np <- consistency(sc, "Pseudobulk:Cosine")
  np$dissimilarity_method <- nopca
  cons <- rbind(cons, np)
  slots[[nopca]] <- dmat(sc, "Pseudobulk:Cosine")

  key <- paste(cons$consistency_metric, cons$dissimilarity_method, sep = " | ")
  labels <- sort(unique(obs[[ident]]))
  purged <- setNames(scTypeEval:::purge_label(labels), labels)  # scTypeEval's own label clean-up, called not copied
  per_type <- lapply(labels, function(l) {
    v <- as.list(setNames(cons$measure[cons$celltype == purged[[l]]], key[cons$celltype == purged[[l]]]))
    p <- v[[local]] * v[[global]]
    v$product <- p
    v$geometric_mean <- if (is.na(p) || v[[local]] < 0 || v[[global]] < 0) NA_real_ else sqrt(p)
    v
  })
  names(per_type) <- labels
  cols <- names(per_type[[1]])
  print(cons)
  list(label_map = as.list(purged), per_type = per_type,
       mean_over_types = lapply(setNames(cols, cols),
                                function(k) mean(vapply(per_type, `[[`, numeric(1), k), na.rm = TRUE)),
       dissimilarity = slots)
}

out <- list(
  meta = list(scTypeEval = as.character(packageVersion("scTypeEval")),
              scTypeEval_source = paste(packageDescription("scTypeEval")$Repository, "git",
                                        packageDescription("scTypeEval")$git_last_commit),
              SingleR = as.character(packageVersion("SingleR")),
              scran = as.character(packageVersion("scran")),
              irlba = as.character(packageVersion("irlba")),
              cluster = as.character(packageVersion("cluster")),
              R = R.version.string, generated = format(Sys.Date())),
  names = list(ISCL = local, ISCG = global, ISCG_no_pca = paste("2label_silhouette |", nopca)),
  columns = list(annotation = score_column("annotation"),
                 annotation_oversplit = score_column("annotation_oversplit")),
  note = paste("NA = label not scored (< min_samples samples with >= min_cells cells). recip_classif:Match",
               "stores 0.5 for a mutual hit between labels not shared by both samples; get_consistency reads",
               "it as 1. product = ISCL * ISCG (paper Fig. 5, label-transfer task); geometric_mean =",
               "sqrt(product) (README), NA when either is negative. scran: version installed alongside, the",
               "HVGs were picked by make_input.R, not here.")
)
writeLines(jsonlite::toJSON(out, auto_unbox = TRUE, digits = NA, na = "null", pretty = TRUE), a[2])
