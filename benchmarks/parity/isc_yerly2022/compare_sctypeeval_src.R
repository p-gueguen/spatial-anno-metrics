# Compare the INSTALLED scTypeEval namespace with the R/ sources of one or more source trees, object by object.
# Functions are compared by deparsed formals + body (srcref and whitespace ignored), other top-level objects by
# deparse, S4 classes by their slot types, the S4 initialize method by its definition, plus the black-list md5.
# Trees = unpacked GitHub archives, e.g. github.com/carmonalab/scTypeEval/archive/<sha>.tar.gz (GPL-3: read,
# never copied here). Result and positive control: SOURCE.md, "Upstream reference tool".
# Usage: Rscript compare_sctypeeval_src.R <tree> [<tree> ...]
suppressPackageStartupMessages(library(scTypeEval))
ns <- asNamespace("scTypeEval")
canon <- function(x) if (is.function(x)) c(deparse(formals(x)), deparse(body(x))) else deparse(x)

compare_tree <- function(tree) {
  env <- new.env(parent = globalenv())
  classes <- list(); init <- NULL
  files <- sort(list.files(file.path(tree, "R"), "[.]R$", full.names = TRUE))
  for (f in files) for (e in parse(f, keep.source = FALSE)) {
    head <- if (is.call(e)) paste(deparse(e[[1]]), collapse = "") else ""
    if (head %in% c("<-", "=") && is.name(e[[2]])) {
      assign(as.character(e[[2]]), eval(e[[3]], env), envir = env)
    } else if (head %in% c("setClass", "methods::setClass")) {
      e[[1]] <- quote(methods::setClass); e$where <- env
      cl <- e[[2]]
      suppressWarnings(suppressMessages(eval(e, env)))
      classes[[cl]] <- methods::getClassDef(cl, where = env)@slots
    } else if (head %in% c("setMethod", "methods::setMethod") && identical(e[[2]], "initialize")) {
      init <- eval(e[[4]], env)
    }
  }
  src <- ls(env, all.names = TRUE)
  inst <- setdiff(ls(ns, all.names = TRUE), c(".__NAMESPACE__.", ".__S3MethodsTable__.", ".packageName",
                                              grep("^[.]__[CT]__", ls(ns, all.names = TRUE), value = TRUE)))
  common <- intersect(src, inst)
  differ <- common[!vapply(common, function(n) identical(canon(get(n, env)), canon(get(n, ns))), logical(1))]
  cls_differ <- names(classes)[!vapply(names(classes), function(cl)
    identical(classes[[cl]], methods::getClassDef(cl, where = ns)@slots), logical(1))]
  init_inst <- methods::getMethod("initialize", "scTypeEval", where = ns)@.Data
  cat(sprintf("\n== %s\n%d R files, %d top-level objects in source, %d in installed namespace, %d compared\n",
              basename(tree), length(files), length(src), length(inst), length(common)))
  cat("only in source   :", if (length(setdiff(src, inst))) setdiff(src, inst) else "-", "\n")
  cat("only in installed:", if (length(setdiff(inst, src))) setdiff(inst, src) else "-", "\n")
  cat("differ           :", if (length(differ)) differ else "-", "\n")
  cat(sprintf("S4 classes       : %d compared, differ: %s\n", length(classes),
              if (length(cls_differ)) paste(cls_differ, collapse = " ") else "-"))
  cat("initialize method:", if (identical(canon(init), canon(init_inst))) "identical" else "DIFFERS", "\n")
  data_src <- tools::md5sum(file.path(tree, "data", "black_list.rda"))
  data_inst <- tools::md5sum(file.path(find.package("scTypeEval"), "data", "black_list.rda"))
  cat("data/black_list.rda md5:", if (identical(unname(data_src), unname(data_inst))) "identical" else "DIFFERS", "\n")
  invisible(differ)
}
cat("installed scTypeEval", as.character(packageVersion("scTypeEval")), "from",
    packageDescription("scTypeEval")$Repository, "git", packageDescription("scTypeEval")$git_last_commit, "\n")
for (t in commandArgs(TRUE)) compare_tree(t)
