#!/usr/bin/env Rscript
# Execute ORIGINAL coloc pinned check_dataset and process.dataset functions.
# Requires only base R; no package installation, no simulation of original R.
# Each test exercises a genuine author R function from pin-verified Git blobs.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2L)
source(args[[1L]], local=.GlobalEnv, echo=FALSE)
source(args[[2L]], local=.GlobalEnv, echo=FALSE)
stopifnot(exists("check_dataset", mode="function"),
          exists("process.dataset", mode="function"),
          exists("approx.bf.estimates", mode="function"))

test <- function(name, passed) {
    if (!isTRUE(passed)) stop(paste0("ORIGINAL COLOC R CONTRACT FAILED: ", name))
    cat(paste0("ORIGINAL_COLOC_TEST_", name, "=PASS\n"))
}

error_message <- function(expression) {
    tryCatch({force(expression); ""},error=function(e) conditionMessage(e))
}

cc <- list(
    type="cc",
    snp=c("snp1","snp2","snp3","snp4"),
    beta=c(0.22,-0.16,0.04,0.12),
    varbeta=c(0.009,0.012,0.008,0.01)
)
ok <- suppressWarnings(check_dataset(cc))
res <- suppressWarnings(process.dataset(cc,"df1"))
test("cc_beta_varbeta_without_s",
     is.null(ok) && nrow(res)==4L &&
     all(is.finite(res[["lABF.df1"]])) &&
     identical(res[["snp"]],cc$snp))

cc_small <- cc
cc_small$s <- 0.07
cc_large <- cc
cc_large$s <- 0.57
suppressWarnings(check_dataset(cc_small))
suppressWarnings(check_dataset(cc_large))
small <- suppressWarnings(process.dataset(cc_small,"df1"))
large <- suppressWarnings(process.dataset(cc_large,"df1"))
test("cc_beta_varbeta_with_s_same_abf",
     isTRUE(all.equal(res[["lABF.df1"]],small[["lABF.df1"]],tolerance=0)) &&
     isTRUE(all.equal(res[["lABF.df1"]],large[["lABF.df1"]],tolerance=0)))

cc_p <- list(
    type="cc",
    snp=cc$snp,
    pvalues=c(0.0001,0.002,0.08,0.01),
    MAF=c(0.12,0.22,0.32,0.42),
    N=1200
)
err <- error_message(suppressWarnings(check_dataset(cc_p)))
test("cc_pvalue_only_missing_s_rejected",
     nzchar(err) && grepl("s, proportion",err,fixed=TRUE))

cc_p$s <- 0.15
ok <- suppressWarnings(check_dataset(cc_p))
test("cc_pvalue_only_supplied_s_allowed",is.null(ok))

qt <- list(
    type="quant",
    snp=cc$snp,
    beta=cc$beta,
    varbeta=cc$varbeta,
    sdY=1.7
)
ok <- suppressWarnings(check_dataset(qt))
qt_result <- suppressWarnings(process.dataset(qt,"df2"))
test("quant_beta_varbeta_sdY_explicit_allowed",
     is.null(ok) && nrow(qt_result)==4L &&
     all(is.finite(qt_result[["lABF.df2"]])))

qt_missing <- qt
qt_missing$sdY <- NULL
err <- error_message(suppressWarnings(check_dataset(qt_missing)))
test("quant_beta_varbeta_missing_sdY_N_MAF_rejected",
     nzchar(err) && grepl("sdY for type quant",err,fixed=TRUE))
