#!/usr/bin/env Rscript
# Read-only post-fit stability audit. Never promotes exploratory loci.
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<1) stop("Usage: Rscript audit_alcohol_finite_ld_stability.R <sandbox_dir>")
root <- normalizePath(args[[1]],mustWork=TRUE)
rows <- list()
for(locus in c("ADH1B","ALDH2")){
  ids <- c("FINITE_REF_504","FINITE_REF_504_EB_MISMATCH")
  fits <- lapply(ids,function(id) readRDS(file.path(root,"models",locus,paste0(id,".rds"))))
  cs <- lapply(fits,function(f) if(is.null(f$sets$cs)) list() else f$sets$cs)
  sets <- lapply(cs,function(z) lapply(z,as.integer))
  jaccard <- function(a,b){
    if(!length(a) && !length(b)) return(1)
    if(!length(a)||!length(b)) return(0)
    best <- vapply(a,function(x) max(vapply(b,function(y){
      union_n <- length(union(x,y))
      if(!union_n) return(1)
      length(intersect(x,y))/union_n
    },numeric(1))),numeric(1))
    mean(best)
  }
  cs_detail <- lapply(fits,function(f) {
    z <- if(is.null(f$sets$cs)) list() else f$sets$cs
    lapply(seq_along(z),function(j) {
      ix <- as.integer(z[[j]])
      ix <- ix[!is.na(ix) & ix>=1 & ix<=length(f$pip)]
      purity <- if(!is.null(f$sets$purity) && nrow(as.matrix(f$sets$purity))>=j)
        as.numeric(f$sets$purity[j,]) else numeric()
      list(cs_name=names(z)[j],n_snps=length(ix),
           lead_index=if(length(ix)) ix[which.max(f$pip[ix])] else NULL,
           lead_pip=if(length(ix)) max(f$pip[ix]) else NULL,
           purity=purity)
    })
  })
  d <- lapply(fits,function(f) f$R_finite_diagnostics)
  rows[[locus]] <- list(
    locus=locus,model_ids=ids,
    credible_set_counts=vapply(sets,length,integer(1)),
    credible_set_details=cs_detail,
    ser_alpha_dimensions=lapply(fits,function(f) dim(f$alpha)),
    cs_jaccard_A_to_B=jaccard(sets[[1]],sets[[2]]),
    cs_jaccard_B_to_A=jaccard(sets[[2]],sets[[1]]),
    max_abs_pip_delta=max(abs(fits[[1]]$pip-fits[[2]]$pip)),
    mean_abs_pip_delta=mean(abs(fits[[1]]$pip-fits[[2]]$pip)),
    pip_delta_invariant_pass=(mean(abs(fits[[1]]$pip-fits[[2]]$pip)) <= max(abs(fits[[1]]$pip-fits[[2]]$pip)) + 1e-15),
    Q_art=lapply(d,function(x) if(is.null(x$Q_art)) NULL else x$Q_art),
    R_reliability_flag=lapply(d,function(x) x$R_reliability_flag),
    R_sensitivity_flag=lapply(d,function(x) x$R_sensitivity_flag),
    status=if(locus=="ALDH2") "BLOCKED" else "EXPLORATORY",
    causal_inference_validated=FALSE)
}
out <- file.path(root,"ALCOHOL_FINITE_LD_STABILITY_READONLY.json")
write_json(list(status="DIAGNOSTIC_ONLY",results=rows),out,pretty=TRUE,auto_unbox=TRUE,digits=NA,na="null")
cat("STABILITY_DIAGNOSTICS_WRITTEN",out,"\n")
