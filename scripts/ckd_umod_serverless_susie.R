#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(susieR))
d <- Sys.getenv("CKD_UMOD_INPUT_DIR","input/ckd/umod")
dat <- read.delim(gzfile(file.path(d,"UMOD.tsv.gz")),check.names=FALSE)
meta <- read.delim(file.path(d,"UMOD.ld.meta.tsv"),check.names=FALSE)
vars <- readLines(file.path(d,"UMOD.ld.vars"))
stopifnot(length(vars)==nrow(meta),all(meta$ref_id==vars),!anyDuplicated(meta$snp))
ord <- match(meta$snp,dat$snp)
stopifnot(!anyNA(ord))
dat <- dat[ord,,drop=FALSE]
n <- nrow(meta)
con <- file(file.path(d,"UMOD.ld.bin"),"rb")
ld <- matrix(readBin(con,"numeric",n=n*n,size=4,endian="little"),n,n,byrow=TRUE)
close(con)
stopifnot(all(is.finite(ld)),max(abs(ld-t(ld)))<1e-5)
sgn <- as.numeric(meta$effect_vs_ld_major_sign)
stopifnot(all(sgn %in% c(-1,1)))
ld <- ld*tcrossprod(sgn)
diag(ld)<-1
target<-match("rs11860566",dat$snp)
stopifnot(!is.na(target))
absz<-abs(dat$beta_outcome_aligned/dat$se_outcome)
sel<-sort(unique(c(target,head(order(absz,decreasing=TRUE),180L))))
z<-dat$beta_outcome_aligned[sel]/dat$se_outcome[sel]
N<-round(median(dat$n_outcome[sel]))
R0<-ld[sel,sel,drop=FALSE]
for(lambda in c(0,0.001,0.01,0.05,0.1)){
 R<-(1-lambda)*R0
 diag(R)<-1
 fit<-tryCatch(susie_rss(z=z,R=R,n=N,L=3,max_iter=100,estimate_residual_variance=FALSE,refine=FALSE,verbose=FALSE),error=function(e)e)
 if(inherits(fit,"error")){
  cat("UMOD_LD_SHRINKAGE lambda=",lambda," status=ERROR message=",conditionMessage(fit),"\n",sep="")
 }else{
  cat("UMOD_LD_SHRINKAGE lambda=",lambda," converged=",isTRUE(fit$converged)," iterations=",fit$niter,
  " sets=",if(is.null(fit$sets$cs))0 else length(fit$sets$cs),"\n",sep="")
 }
}
