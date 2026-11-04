# scDesign3 fit of the Xenium breast crop -> parameter cache for world.py.
#
#   Rscript generate.R <cache> <n_cores> [k_gp=60] [n_genes=all]
#
# Reads <cache>/input/{counts.mtx, genes.txt, obs.csv} (written by `python build.py prepare`) and writes
# <cache>/fit/:
#   mean.csv.gz        cells x genes fitted NB mean (extract_para$mean_mat)
#   sigma.csv          per-gene NB dispersion in gamlss NBI form, sigma = 1 / theta (Var = mu + sigma mu^2)
#   copula_<type>.csv.gz  Gaussian-copula correlation matrix per cell type (fit_copula, corr_by = cell type)
#   rdraw.mtx          one scDesign3::simu_new draw (genes x cells), to validate world.py's Python sampler
#   gene_fit.csv       per gene: theta, edf of the spatial smooth, log_area coefficient, deviance explained
#   fit_info.json      parameters, versions, timings, failed genes
#
# Model (scDesign3 1.7.x, Song et al. Nat Biotechnol 2024), one GAM per gene with mgcv:
#   counts ~ NB(mu, theta),  log(mu) = cell_type + log_area + s(spatial1, spatial2, bs = "gp", k = k_gp)
# i.e. a cell-type effect, a cell-size effect, and a per-gene Gaussian-process smooth over the real cell
# centroids (the same smooth for every type: regional expression shifts). Gene-gene dependence: a Gaussian
# copula per cell type, fitted on the distributional-transform PIT of the marginals (DT = TRUE).
#
# Reproducibility: the copula's distributional transform jitters every count at random inside forked
# mcmapply workers. Under the default Mersenne-Twister those draws are not reproducible (two fits gave copula
# entries agreeing at r = 0.19-0.82 across types); L'Ecuyer-CMRG gives each worker a stream derived from the
# seed, so the fit reproduces exactly for the same n_cores.
suppressPackageStartupMessages({
  library(scDesign3); library(SingleCellExperiment); library(Matrix)
})
RNGkind("L'Ecuyer-CMRG")
args <- commandArgs(trailingOnly = TRUE)
cache <- args[1]
n_cores <- as.integer(args[2])
k_gp <- if (length(args) >= 3) as.integer(args[3]) else 60L
n_genes <- if (length(args) >= 4) as.integer(args[4]) else NA_integer_
out <- file.path(cache, if (is.na(n_genes)) "fit" else sprintf("fit_test_%d", n_genes))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
t0 <- Sys.time()
tm <- list()
stamp <- function(what) { tm[[what]] <<- as.numeric(difftime(Sys.time(), t0, units = "secs")); message(what, ": ", round(tm[[what]]), " s") }

counts <- as(readMM(file.path(cache, "input", "counts.mtx")), "CsparseMatrix")
genes <- readLines(file.path(cache, "input", "genes.txt"))
obs <- read.csv(file.path(cache, "input", "obs.csv"), stringsAsFactors = FALSE)
rownames(counts) <- genes
colnames(counts) <- as.character(obs$cell)
if (!is.na(n_genes)) counts <- counts[seq_len(n_genes), ]
sce <- SingleCellExperiment(assays = list(counts = counts),
                            colData = DataFrame(cell_type = factor(obs$cell_type), spatial1 = obs$x,
                                                spatial2 = obs$y, log_area = obs$log_area,
                                                row.names = as.character(obs$cell)))
set.seed(2026)
data <- construct_data(sce = sce, assay_use = "counts", celltype = "cell_type", pseudotime = NULL,
                       spatial = c("spatial1", "spatial2"), other_covariates = "log_area",
                       corr_by = "cell_type")
stamp("construct_data")
mu_formula <- sprintf("cell_type + log_area + s(spatial1, spatial2, bs = 'gp', k = %d)", k_gp)
marginal <- fit_marginal(data = data, predictor = "gene", mu_formula = mu_formula, sigma_formula = "1",
                         family_use = "nb", n_cores = n_cores, usebam = FALSE, parallelization = "mcmapply")
stamp("fit_marginal")
failed <- names(marginal)[vapply(marginal, function(m) !inherits(m$fit, "gam"), logical(1))]
gene_fit <- do.call(rbind, lapply(names(marginal), function(g) {
  f <- marginal[[g]]$fit
  if (!inherits(f, "gam")) return(data.frame(gene = g, theta = NA, edf_spatial = NA, log_area = NA, dev_expl = NA))
  s <- summary(f)
  data.frame(gene = g, theta = f$family$getTheta(TRUE), edf_spatial = s$s.table[1, "edf"],
             log_area = unname(coef(f)["log_area"]), dev_expl = s$dev.expl)
}))
write.csv(gene_fit, file.path(out, "gene_fit.csv"), row.names = FALSE)

set.seed(2027)
copula <- fit_copula(sce = sce, assay_use = "counts", marginal_list = marginal, family_use = "nb",
                     copula = "gaussian", n_cores = n_cores, input_data = data$dat, DT = TRUE)
stamp("fit_copula")
para <- extract_para(sce = sce, marginal_list = marginal, n_cores = n_cores, family_use = "nb",
                     new_covariate = data$newCovariate, data = data$dat)
stamp("extract_para")
# scDesign3 1.7.6 simu_new() looks a group's copula up with copula_list[[x]], x from unique(corr_group): with
# a factor corr_group that is the factor CODE (alphabetical), while copula_list is in first-appearance order,
# so cell types silently get another type's copula (measured on this fit: Epi cells drew T_CD4's copula, the
# draw's gene-gene dependence fell to a slope of 0.04 on the fitted correlations). A character corr_group
# makes [[x]] look up by name. fit_copula() names its list correctly; only the sampler is affected.
dat_chr <- data$dat
dat_chr$corr_group <- as.character(dat_chr$corr_group)
new_chr <- data$newCovariate
if ("corr_group" %in% colnames(new_chr)) new_chr$corr_group <- as.character(new_chr$corr_group)
stopifnot(setequal(names(copula$copula_list), unique(dat_chr$corr_group)))
set.seed(2028)
rdraw <- simu_new(sce = sce, mean_mat = para$mean_mat, sigma_mat = para$sigma_mat, zero_mat = para$zero_mat,
                  quantile_mat = NULL, copula_list = copula$copula_list, n_cores = n_cores, family_use = "nb",
                  input_data = dat_chr, new_covariate = new_chr, filtered_gene = data$filtered_gene)
stamp("simu_new")

mm <- para$mean_mat
gz <- gzfile(file.path(out, "mean.csv.gz"), "w")
write.csv(mm, gz)
close(gz)
write.csv(data.frame(gene = colnames(para$sigma_mat), sigma = para$sigma_mat[1, ]),
          file.path(out, "sigma.csv"), row.names = FALSE)
for (t in names(copula$copula_list)) {
  C <- as.matrix(copula$copula_list[[t]])
  gz <- gzfile(file.path(out, sprintf("copula_%s.csv.gz", t)), "w")
  write.csv(C, gz)
  close(gz)
}
writeMM(as(as.matrix(rdraw), "CsparseMatrix"), file.path(out, "rdraw.mtx"))
writeLines(rownames(rdraw), file.path(out, "rdraw_genes.txt"))
stamp("write")
info <- list(mu_formula = mu_formula, sigma_formula = "1", family = "nb", copula = "gaussian",
             corr_by = "cell_type", DT = TRUE, k_gp = k_gp, n_cells = ncol(sce), n_genes = nrow(sce),
             failed_genes = failed, filtered_gene = as.character(unlist(data$filtered_gene)),
             copula_aic = copula$model_aic, n_cores = n_cores, timings_s = tm,
             scDesign3 = as.character(packageVersion("scDesign3")), mgcv = as.character(packageVersion("mgcv")),
             R = R.version.string, host = Sys.info()[["nodename"]], seeds = c(construct = 2026, copula = 2027, simu = 2028))
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, digits = NA), file.path(out, "fit_info.json"))
message("done")
