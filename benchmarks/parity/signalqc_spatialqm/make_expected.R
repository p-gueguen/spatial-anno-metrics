#!/usr/bin/env Rscript
# SpatialQM reference outputs -> tests/goldens/signalqc_spatialqm/expected.json. NOT run in CI.
# Input: the Xenium-style crop folder written by `python make_inputs.py outs`.
# Usage: Rscript make_expected.R <crop_outs_dir> <expected.json> <SpatialQM source checkout>
#
# SpatialQM (MIT), installed with pak from github.com/Center-for-Spatial-OMICs/SpatialQM. The sha comes from
# the installed DESCRIPTION (RemoteSha) and must equal the checkout's HEAD. The checkout is only parsed, to
# record which file's definition of each called function is installed (meta.function_source): R/utils.R
# (2026-05-06) and R/utils_update_final.R (2025-02-20) define the same functions, and the file sourced last
# (utils_update_final.R in C-locale order) wins. What the installed code needs here (details in SOURCE.md):
# - every metric except getMorans uses the file branch (seu_obj = NULL); the installed Seurat branch of
#   getMeanSignalRatio, getMaxRatio, getSparsity, getEntropy, getTxPerArea, getTxPerNuc, getCellTxFraction
#   is empty;
# - getMeanSignalRatio and getMaxRatio error on Xenium input (base::rowMeans gets a Matrix object). The
#   error is kept in meta.errors; the value comes from the same installed functions on their CosMx branch
#   (fread -> base matrix), fed the same counts as a CosMx-style CSV;
# - getMorans takes a Seurat object, built the way readSpatial(seurat = TRUE) / LoadXenium lay it out:
#   RNA = Gene Expression, ControlProbe = Negative Control Probe, reduction "tissue" = cell centroids.
#   Its callees colGraph<-, findSpatialNeighbors, colData<-, rowVars are not imported by SpatialQM, hence
#   the library() calls. Seurat turns "_" into "-" in feature names; expected.json maps them back.
suppressPackageStartupMessages({
  library(SpatialQM); library(Seurat); library(SingleCellExperiment); library(SpatialExperiment)
  library(SpatialFeatureExperiment); library(Voyager); library(MatrixGenerics); library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
outs <- args[1]; out_json <- args[2]; src <- args[3]
sha <- packageDescription("SpatialQM")$RemoteSha
stopifnot(!is.null(sha), identical(system2("git", c("-C", src, "rev-parse", "HEAD"), stdout = TRUE), sha))

mtx <- file.path(outs, "cell_feature_matrix", "matrix.mtx.gz")
tx <- file.path(outs, "transcripts.csv.gz")
cells_csv <- file.path(outs, "cells.csv.gz")
X <- "Xenium"
counts <- as(Matrix::readMM(mtx), "CsparseMatrix")
feat <- data.table::fread(file.path(outs, "cell_feature_matrix", "features.tsv.gz"), header = FALSE)
bc <- data.table::fread(file.path(outs, "cell_feature_matrix", "barcodes.tsv.gz"), header = FALSE)
dimnames(counts) <- list(feat$V2, as.character(bc$V1))
cm <- data.table::fread(cells_csv)
stopifnot(identical(as.character(cm$cell_id), colnames(counts)),
          identical(make.names(feat$V2), feat$V2))  # else data.frame() in the CosMx branch renames features
cosmx_csv <- file.path(tempdir(), "exprMat_cosmx_layout.csv")  # cells x features; the branch drops cell_ID
data.table::fwrite(data.table::data.table(cell_ID = colnames(counts), as.matrix(Matrix::t(counts))), cosmx_csv)

errors <- c(getMeanSignalRatio = tryCatch({getMeanSignalRatio(expMat = mtx, platform = X); NA},
                                          error = conditionMessage),
            getMaxRatio = tryCatch({getMaxRatio(expMat = mtx, platform = X); NA}, error = conditionMessage))
if (anyNA(errors)) stop("the Xenium branch ran: drop the CosMx-branch workaround and update SOURCE.md")
res <- list(
  getNcells = getNcells(expMat = mtx, platform = X),
  getTxPerCell = getTxPerCell(expMat = mtx, platform = X),
  getTxPerArea = getTxPerArea(cellSegMeta = cells_csv, platform = X),
  getTxPerNuc = getTxPerNuc(tx_file = tx, platform = X),
  getCellTxFraction = getCellTxFraction(tx_file = tx, platform = X),
  getGlobalFDR = getGlobalFDR(tx_file = tx, platform = X),
  getMeanSignalRatio = getMeanSignalRatio(expMat = cosmx_csv, platform = "CosMx"),
  getMaxRatio = getMaxRatio(expMat = cosmx_csv, platform = "CosMx"),
  getSparsity = getSparsity(expMat = mtx, platform = X),
  getEntropy = getEntropy(expMat = mtx, platform = X)
)

seu <- CreateSeuratObject(counts = counts[feat$V3 == "Gene Expression", ], assay = "RNA")
seu[["ControlProbe"]] <- CreateAssayObject(counts = counts[feat$V3 == "Negative Control Probe", ])
seu$sample_id <- "Xenium_FFPE_Human_Breast_Cancer_Rep1_crop"
seu$platform <- X
coords <- as.matrix(cm[, c("x_centroid", "y_centroid")])
dimnames(coords) <- list(colnames(seu), c("Tissue_1", "Tissue_2"))
seu[["tissue"]] <- CreateDimReducObject(coords, key = "Tissue_", assay = "RNA")
mor <- getMorans(seu)
dash <- gsub("_", "-", feat$V2)
stopifnot(!anyDuplicated(dash), mor$gene %in% dash)
res$getMorans <- list(feature = feat$V2[match(mor$gene, dash)], type = mor$type, I = mor$value)

defs <- function(f) {  # top-level `name <- function(...)` of one source file; nothing else is evaluated
  e <- new.env()
  for (x in parse(f, keep.source = FALSE))
    if (is.call(x) && (identical(x[[1]], as.name("<-")) || identical(x[[1]], as.name("="))) &&
        is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")))
      assign(as.character(x[[2]]), eval(x[[3]], baseenv()), envir = e)
  e
}
files <- c("utils.R", "utils_update_final.R")
d <- setNames(lapply(file.path(src, "R", files), defs), files)
installed_from <- function(n) {
  f <- get(n, asNamespace("SpatialQM"))
  hit <- vapply(files, function(s) exists(n, d[[s]], inherits = FALSE) &&
    identical(formals(f), formals(d[[s]][[n]])) && identical(body(f), body(d[[s]][[n]])), TRUE)
  paste(files[hit], collapse = " = ")    # both names: the two files define it identically
}
fsrc <- vapply(setdiff(names(res), "meta"), installed_from, "")
stopifnot(nzchar(fsrc))                  # each called function matches one file's definition exactly

regex <- c(getMeanSignalRatio_getMaxRatio = "Neg*|SystemControl*|Blank*|BLANK*|Unassigned*",  # as installed
           getTxPerCell = "Unassigned*|NegControl*|BLANK*|Blank*|SystemControl*|NegPrb*",
           getGlobalFDR = "(?i)^(Neg|Blank|BLANK)")
txd <- data.table::fread(tx)
pv <- function(p) as.character(packageVersion(p))
res$meta <- list(
  spatialqm_sha = sha, spatialqm_version = pv("SpatialQM"), R = R.version.string,
  function_source = as.list(fsrc), errors = as.list(errors),
  calls = list(getMeanSignalRatio = "expMat = <cells x features CSV of the same counts>, platform = \"CosMx\"",
               getMaxRatio = "expMat = <cells x features CSV of the same counts>, platform = \"CosMx\"",
               getMorans = "Seurat object as above", other = "file branch, platform = \"Xenium\""),
  control_regex_hits_by_feature_type = lapply(regex, function(r)
    as.list(tapply(grepl(r, feat$V2, perl = TRUE), feat$V3, sum))),
  packages = lapply(setNames(nm = c("Seurat", "SeuratObject", "Voyager", "SpatialFeatureExperiment",
                                    "SpatialExperiment", "SingleCellExperiment", "spdep", "scater", "scuttle",
                                    "BiocParallel", "BioQC", "coop", "data.table", "Matrix", "fs")), pv),
  input = list(n_cells = ncol(counts), n_features = nrow(counts), n_transcripts = nrow(txd),
               n_transcripts_unassigned = sum(txd$cell_id == "UNASSIGNED"),
               n_features_by_type = as.list(table(feat$V3))),
  date = format(Sys.Date())
)
write_json(res, out_json, auto_unbox = TRUE, digits = NA, na = "null", pretty = TRUE)
cat("wrote", out_json, "\n")
