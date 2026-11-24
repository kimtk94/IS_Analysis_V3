#!/usr/bin/env Rscript
# UMOD isolated in-memory SuSiE sensitivity check. No source/output files changed.
suppressPackageStartupMessages(library(susieR))
root <- "/srv/is-analysis"
inp <- file.path(root,"results/ckd/stage2c_susie/susie_input/UMOD.tsv.gz")
ldroot <- file.path(root,"data/ckd/stage2c_ld/ld")
dat <- read.delim(gzfile(inp),check.names=FALSE)
meta <- read.delim(file.path(ldroot,"UMOD.ld.meta.tsv"),check.names=FALSE)
ord <- match(meta$snp,dat$snp)
stopifnot(!anyNA(ord),all(meta$ref_id==readLines(file.path(ldroot,"UMOD.ld.vars"))))
dat <- dat[ord,,drop=FALSE]
n <- nrow(meta)
con <- file(file.path(ldroot,"UMOD.ld.bin"),"rb")
ld <- matrix(readBin(con,"numeric",n=n*n,size=4,endian="little"),n,n,byrow=TRUE)
close(con)
stopifnot(all(is.finite(ld)))
sgn <- as.numeric(meta$effect_vs_ld_major_sign)
ld <- ld * tcrossprod(sgn)
diag(ld) <- 1
stopifnot(max(abs(ld-t(ld)))<1e-5)
target <- match("rs11860566",dat$snp)
stopifnot(!is.na(target))
cat("UMOD_SENSITIVITY_INPUT n=",n," target_row=",target,"\n")
# Predeclare a fixed signal-enriched subset shared by both arms.
absz <- abs(dat$beta_outcome_aligned/dat$se_outcome)
selected <- sort(unique(c(target,head(order(absz,decreasing=TRUE),180L))))
cat("UMOD_SENSITIVITY_SUBSET n=",length(selected)," selection=top_180_abs_egfr_z_plus_outlier; exploratory_only=TRUE\n")
# Fit eGFR summary statistics in memory with identical fit parameters both arms.
for (drop in c(FALSE,TRUE)) {
 keep <- if(drop) setdiff(selected,target) else selected
 z <- dat$beta_outcome_aligned[keep]/dat$se_outcome[keep]
 N <- round(median(dat$n_outcome[keep]))
 fit <- tryCatch(susie_rss(z=z,R=ld[keep,keep,drop=FALSE],n=N,L=3,max_iter=100,
  estimate_residual_variance=FALSE,refine=FALSE,verbose=FALSE),error=function(e)e)
 if (inherits(fit,"error")) {
  cat("UMOD_SENSITIVITY_RESULT arm=",if(drop)"EXCLUDED" else "INCLUDED",
      " status=ERROR message=",conditionMessage(fit),"\n",sep="")
 } else {
  cat("UMOD_SENSITIVITY_RESULT arm=",if(drop)"EXCLUDED" else "INCLUDED",
      " converged=",isTRUE(fit$converged),
      " iterations=",fit$niter,
      " sets=",if(is.null(fit$sets$cs))0 else length(fit$sets$cs),
      " elbo_tail=",tail(fit$elbo,1),"\n",sep="")
 }
 rm(fit);gc(verbose=FALSE)
}
