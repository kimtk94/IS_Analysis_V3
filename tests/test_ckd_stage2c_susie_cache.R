#!/usr/bin/env Rscript
# Base-R-only regressions: no coloc installation or scientific analysis runs.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 1L)
root <- normalizePath(args[[1]], mustWork=TRUE)
scope <- new.env(parent=globalenv())
source(file.path(root, "scripts", "run_ckd_stage2c_susie.R"), local=scope)
fixtures <- file.path(root, "tests", "fixtures")
summary <- read.delim(file.path(fixtures, "ckd_susie_cache_summary.tsv"),
                      stringsAsFactors=FALSE)
LD <- as.matrix(read.delim(file.path(fixtures, "ckd_susie_cache_ld.tsv"),
                           row.names=1, check.names=FALSE))
d <- list(beta=summary$beta, varbeta=summary$varbeta, snp=summary$snp,
          position=summary$pos37, MAF=summary$MAF, N=100L, sdY=1, LD=LD, type="quant")
context <- list(summary=summary,
                ordered_ld_metadata=data.frame(snp=summary$snp, sign=c(1, -1, 1)))
config <- scope$susie_model_config(length(d$snp))
versions <- list(R="synthetic-R", coloc="synthetic-coloc", susieR="synthetic-susieR")
baseline <- list(d=d, context=context, config=config, versions=versions)

runs <- 0L
fake_fitter <- function(d, suffix, config) {
  runs <<- runs + 1L
  list(ok=TRUE,
       fit=structure(list(converged=TRUE, observed_beta=d$beta, fit_call=runs,
                          observed_config=config), class="susie"),
       error="")
}
invoke <- function(path, state=baseline, fitter=fake_fitter) {
  scope$fit_susie_cached(state$d, "synthetic", path, state$context,
                         config=state$config, versions=state$versions, fitter=fitter)
}
assert <- function(value, message) {
  if (!isTRUE(value)) stop(message, call.=FALSE)
}

scratch <- tempfile("ckd_susie_cache_tests_")
dir.create(scratch)
# Run inside a function so all temporary model caches are cleaned on failure too.
run_tests <- function() {
  on.exit(unlink(scratch, recursive=TRUE), add=TRUE)
  first_path <- file.path(scratch, "first.rds")
  first <- invoke(first_path)
  repeat_fit <- invoke(first_path)
  assert(runs == 1L && first$cache_status == "miss" && repeat_fit$cache_status == "hit",
         "Identical inputs should reuse exactly one successful fit")
  assert(inherits(readRDS(first_path), "susie"), "Cache must remain a SuSiE fit RDS")

  changes <- list(
    summary_effect=function(x) { x$d$beta[1] <- x$d$beta[1] + 0.05; x },
    effect_variance=function(x) { x$d$varbeta[1] <- x$d$varbeta[1] * 2; x },
    scalar_N=function(x) { x$d$N <- x$d$N + 1L; x },
    per_SNP_N=function(x) { x$context$summary$N[1] <- 99L; x },
    LD=function(x) { x$d$LD[1, 2] <- x$d$LD[2, 1] <- 0.3; x },
    ordered_SNPs=function(x) {
      order <- c(3, 2, 1)
      for (field in c("beta", "varbeta", "snp", "position", "MAF"))
        x$d[[field]] <- x$d[[field]][order]
      x$d$LD <- x$d$LD[order, order]
      x$context$summary <- x$context$summary[order, ]
      x$context$ordered_ld_metadata <- x$context$ordered_ld_metadata[order, ]
      x
    },
    allele_context=function(x) {
      x$context$summary$effect_allele[1] <- "A"
      x$context$summary$other_allele[1] <- "G"
      x
    },
    model_config=function(x) { x$config$maxit <- 500L; x },
    phenotype_scale=function(x) { x$d$sdY <- 2; x },
    R_version=function(x) { x$versions$R <- "different-R"; x },
    coloc_version=function(x) { x$versions$coloc <- "different-coloc"; x },
    susieR_version=function(x) { x$versions$susieR <- "different-susieR"; x }
  )
  for (name in names(changes)) {
    path <- file.path(scratch, paste0(name, ".rds"))
    before <- runs
    invoke(path)
    changed <- changes[[name]](baseline)
    actual <- invoke(path, changed)
    assert(runs == before + 2L && actual$cache_status == "miss",
           paste(name, "must invalidate the previous fit"))
    again <- invoke(path, changed)
    assert(runs == before + 2L && again$cache_status == "hit",
           paste(name, "must reuse its newly matching fit"))
    assert(identical(actual$fit$observed_beta, changed$d$beta),
           paste(name, "returned a stale effect estimate"))
    assert(identical(actual$fit$observed_config, changed$config),
           paste(name, "did not pass fingerprinted settings to the fitter"))
  }

  bad_path <- file.path(scratch, "old_or_invalid.rds")
  saveRDS(structure(list(converged=TRUE), class="susie"), bad_path)
  before <- runs
  assert(invoke(bad_path)$cache_status == "miss" && runs == before + 1L,
         "Legacy fit without fingerprint must be recomputed")
  bad_fit <- readRDS(bad_path)
  bad_fit$converged <- FALSE
  saveRDS(bad_fit, bad_path)
  before <- runs
  assert(invoke(bad_path)$cache_status == "miss" && runs == before + 1L,
         "Nonconverged fit must never be reused")
  writeLines("not an RDS file", bad_path)
  before <- runs
  assert(invoke(bad_path)$cache_status == "miss" && runs == before + 1L,
         "Unreadable cache must be recomputed")

  failed_path <- file.path(scratch, "failed.rds")
  failed <- invoke(failed_path, fitter=function(d, suffix, config)
    list(ok=FALSE, fit=NULL, error="synthetic failure"))
  assert(!failed$ok && !file.exists(failed_path), "Failed fit must not create a reusable cache")
  unconv <- invoke(failed_path, fitter=function(d, suffix, config)
    list(ok=TRUE, fit=list(converged=FALSE), error=""))
  assert(!unconv$ok && !file.exists(failed_path), "Nonconverged new fit must not be cached")
  cat("CKD_SUSIE_CACHE_TEST_PASS changes=", length(changes), "\n", sep="")
}
run_tests()
