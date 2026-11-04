# SRTsim fit + count for the SRTsim world (step 2 of 3; see prepare.py and results/README.md).
#
#   Rscript generate.R RUN N_SEEDS N_CORES
#
# RUN holds prepare.py's inputs: genes.csv, seed_<k>_cells.csv (cell, crop, label, x_local, y_local, ...) and
# seed_<k>_counts.mtx (features x cells, the real Xenium counts). For every seed and crop (one crop = one sample):
#   1. createSRT(count_in = crop counts, loc_in = (x_local, y_local, label))
#   2. srtsim_fit(marginal = "auto_choose", sim_scheme = "domain"): a Poisson / NB / ZIP / ZINB marginal per
#      feature and per cell type (domain = the authors' cell type), chosen by SRTsim's LRT / AIC rules;
#      features with < 2 non-zero cells in a type are not fitted (SRTsim simulates them as 0)
#   3. srtsim_count(breaktie = "random"), rrr = NULL (same depth): new counts drawn from each marginal and
#      placed by the rank of the real counts within the type (reference-based, tissue layout kept)
#   4. srtsim_fit(sim_scheme = "tissue") for the parameters only (world.py draws uninformative extra genes)
# Writes seed_<k>_sim_counts.mtx (+ _sim_cells.csv / _sim_genes.csv), seed_<k>_params_{domain,tissue}.csv and
# seed_<k>_timing.csv to RUN. Seeded per crop: set.seed(1000 * seed + crop).
#
# SRTsim 0.99.7 (CRAN archive; removed from CRAN 2025-12-17). Its NAMESPACE imports Morpho -> rgl, which needs
# libGLU, absent on the compute nodes; so unless library(SRTsim) loads, the five files that hold the fit and
# count code are sourced from the md5-checked archive tarball, unmodified.

args <- commandArgs(trailingOnly = TRUE)
run <- args[1]
n_seeds <- as.integer(args[2])
n_cores <- as.integer(ifelse(length(args) >= 3, args[3], 1))

suppressPackageStartupMessages({
  library(methods); library(Matrix); library(S4Vectors); library(magrittr); library(dplyr); library(parallel)
})

tarball_url <- "https://cran.r-project.org/src/contrib/Archive/SRTsim/SRTsim_0.99.7.tar.gz"
tarball_md5 <- "bfd7f8b16cec6788a91708851e1cf5e6"
srtsim_source <- "package"
if (!suppressWarnings(requireNamespace("SRTsim", quietly = TRUE))) {
  tb <- file.path(run, "SRTsim_0.99.7.tar.gz")
  if (!file.exists(tb)) download.file(tarball_url, tb, quiet = TRUE)
  stopifnot(unname(tools::md5sum(tb)) == tarball_md5)
  src <- file.path(run, "SRTsim_src")
  untar(tb, exdir = src)
  for (f in c("SRTsim_class.R", "createSRT.R", "SRTmodels.R", "SRTfit.R", "SRTcount.R")) {
    source(file.path(src, "SRTsim", "R", f))
  }
  srtsim_source <- paste0("sourced from ", basename(tb), " (md5 ", tarball_md5, ")")
} else {
  suppressPackageStartupMessages(library(SRTsim))
}
cat("SRTsim:", srtsim_source, "\n")

params_long <- function(fit, crop) {
  do.call(rbind, lapply(names(fit), function(lab) {
    f <- fit[[lab]]
    p <- as.data.frame(f$marginal_param1)
    data.frame(crop = crop, label = lab, gene = names(f$gene_sel1), pi0 = p$pi0, theta = p$theta, mu = p$mu,
               model = p$model_selected)
  }))
}

genes <- read.csv(file.path(run, "genes.csv"))$gene
for (seed in seq_len(n_seeds) - 1) {
  cells <- read.csv(file.path(run, sprintf("seed_%d_cells.csv", seed)))
  X <- as(Matrix::readMM(file.path(run, sprintf("seed_%d_counts.mtx", seed))), "CsparseMatrix")
  dimnames(X) <- list(genes, cells$cell)
  res <- mclapply(sort(unique(cells$crop)), function(crop) {
    set.seed(1000 * seed + crop)
    cc <- cells[cells$crop == crop, ]
    loc <- data.frame(x = cc$x_local, y = cc$y_local, label = cc$label, row.names = cc$cell)
    obj <- createSRT(count_in = X[, cc$cell], loc_in = loc)
    t0 <- proc.time()[["elapsed"]]
    fd <- srtsim_fit(obj, marginal = "auto_choose", sim_scheme = "domain")
    t1 <- proc.time()[["elapsed"]]
    fd <- srtsim_count(fd, breaktie = "random")
    t2 <- proc.time()[["elapsed"]]
    ft <- srtsim_fit(obj, marginal = "auto_choose", sim_scheme = "tissue")
    t3 <- proc.time()[["elapsed"]]
    list(counts = fd@simCounts[genes, cc$cell], dom = params_long(fd@EstParam, crop),
         tis = params_long(ft@EstParam, crop),
         timing = data.frame(seed = seed, crop = crop, n_cells = nrow(cc), fit_domain = t1 - t0,
                             count = t2 - t1, fit_tissue = t3 - t2, seconds = t3 - t0))
  }, mc.cores = n_cores, mc.set.seed = FALSE)
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) stop("seed ", seed, ": ", paste(sapply(res[bad], as.character), collapse = "\n"))
  S <- do.call(cbind, lapply(res, `[[`, "counts"))
  Matrix::writeMM(as(S, "CsparseMatrix"), file.path(run, sprintf("seed_%d_sim_counts.mtx", seed)))
  write.csv(data.frame(cell = colnames(S)), file.path(run, sprintf("seed_%d_sim_cells.csv", seed)), row.names = FALSE)
  write.csv(data.frame(gene = rownames(S)), file.path(run, sprintf("seed_%d_sim_genes.csv", seed)), row.names = FALSE)
  write.csv(do.call(rbind, lapply(res, `[[`, "dom")), file.path(run, sprintf("seed_%d_params_domain.csv", seed)),
            row.names = FALSE)
  write.csv(do.call(rbind, lapply(res, `[[`, "tis")), file.path(run, sprintf("seed_%d_params_tissue.csv", seed)),
            row.names = FALSE)
  tm <- do.call(rbind, lapply(res, `[[`, "timing"))
  write.csv(tm, file.path(run, sprintf("seed_%d_timing.csv", seed)), row.names = FALSE)
  cat(sprintf("seed %d: %d cells, %.0f s fit+count (sum over crops)\n", seed, ncol(S), sum(tm$seconds)))
}
writeLines(c(srtsim_source, capture.output(sessionInfo())), file.path(run, "generate_session.txt"))
