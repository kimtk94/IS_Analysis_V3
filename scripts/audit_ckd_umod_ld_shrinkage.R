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
# Compare ridge-to-identity shrinkage on a fixed included SNP subset.
keep <- selected
z <- dat$beta_outcome_aligned[keep]/dat$se_outcome[keep]
N <- round(median(dat$n_outcome[keep]))
baseR <- ld[keep,keep,drop=FALSE]
for (lambda in c(0,0.001,0.01,0.05,0.1)) {
 R <- (1-lambda)*baseR
 diag(R) <- 1
 fit <- tryCatch(susie_rss(z=z,R=R,n=N,L=3,max_iter=100,
     estimate_residual_variance=FALSE,refine=FALSE,verbose=FALSE),error=function(e)e)
 if(inherits(fit,"error")) {
   cat("UMOD_LD_SHRINKAGE lambda=",lambda," status=ERROR msg=",conditionMessage(fit),"\n",sep="")
 } else {
   cat("UMOD_LD_SHRINKAGE lambda=",lambda," converged=",isTRUE(fit$converged),
       " iterations=",fit$niter," sets=",if(is.null(fit$sets$cs))0 else length(fit$sets$cs),
       " elbo_tail=",tail(fit$elbo,1),"\n",sep="")
 }
 rm(fit);gc(verbose=FALSE)
}
