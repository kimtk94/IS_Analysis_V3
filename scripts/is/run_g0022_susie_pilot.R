#!/usr/bin/env Rscript
# G0022 signed-LD SuSiE-RSS pilot. Causal inference NOT established.
# Assumptions: 1KG EAS 504-person LD; 250kb truncated windows;
# case/control effective N unavailable. Compare 20k and reported total 256274
# as explicit hypothetical N sensitivity, NOT valid effective sample-size estimates.
suppressPackageStartupMessages(library(susieR))
options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>0) args[[1]] else "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
files <- list.files(root,pattern="^input_qc\\.json$",full.names=TRUE,recursive=TRUE)
if(length(files)!=5L) stop(sprintf("Expected exactly 5 real pilot input records; got %d",length(files)))
rows <- list()
for(f in files) {
  folder <- dirname(f)
  signal <- basename(folder)
  zfile <- file.path(folder,"z.tsv")
  rfile <- file.path(folder,"signed_ld.tsv.gz")
  vfile <- file.path(folder,"variants.tsv")
  z <- read.table(zfile,header=TRUE)$z
  R <- as.matrix(read.table(gzfile(rfile),sep="\t",header=FALSE))
  v <- read.table(vfile,sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
  if(length(z)!=nrow(R) || nrow(R)!=ncol(R) || nrow(v)!=length(z)) stop("Mismatched z/R/variants")
  if(any(!is.finite(z)) || any(!is.finite(R))) stop("nonfinite z/R")
  if(max(abs(diag(R)-1))>1e-5 || max(abs(R-t(R)))>1e-5) stop("LD diagonal/symmetry invalid")
  for(n in c(20000,256274)) {
    label <- paste0("N",n)
    status <- "RUN_FAILED"
    top_variant <- ""
    top_pip <- NA_real_
    cs_count <- 0L
    converged <- FALSE
    elapsed <- 0
    errors <- ""
    cs_purity_min <- NA_real_
    fit <- NULL
    t0 <- proc.time()[["elapsed"]]
    fit <- tryCatch(susie_rss(z=z,R=R,n=n,L=5,
            estimate_residual_variance=FALSE,max_iter=150,
            coverage=0.95,check_prior=TRUE,verbose=FALSE),
        error=function(e){errors <<- conditionMessage(e);NULL})
    elapsed <- round(proc.time()[["elapsed"]]-t0,2)
    if(!is.null(fit)) {
      converged <- isTRUE(fit$converged)
      status <- if(converged) "EXPLORATORY_CONVERGED" else "EXPLORATORY_NOT_CONVERGED"
      pips <- susie_get_pip(fit)
      if(length(pips)!=nrow(v)) stop("PIP dimension mismatch")
      top <- which.max(pips)
      top_variant <- v$variant_id[top]
      top_pip <- pips[top]
      cs <- susie_get_cs(fit,Xcorr=R,coverage=0.95,min_abs_corr=0.5)
      cs_count <- if(is.null(cs$cs)) 0L else length(cs$cs)
      if(!is.null(cs$purity) && nrow(cs$purity)>0)
        cs_purity_min <- min(cs$purity$min.abs.corr,na.rm=TRUE)
      out <- v
      out$pip <- pips
      write.table(out,file.path(folder,paste0("susie_",label,"_PIP.tsv")),
        sep="\t",quote=FALSE,row.names=FALSE)
      saveRDS(fit,file.path(folder,paste0("susie_",label,"_pilot.rds")))
    }
    rows[[length(rows)+1L]] <- data.frame(
      signal=signal, hypothetical_n=n,n_snps=length(z),status=status,
      converged=converged,cs_count=cs_count,
      min_abs_corr_purity=cs_purity_min,top_pip_variant=top_variant,
      top_pip=top_pip,elapsed_sec=elapsed,error_message=errors,
      scientific_verdict="PILOT_ONLY_NOT_CAUSAL_OR_FINAL_FINE_MAPPING",
      stringsAsFactors=FALSE)
    cat(signal,label,status,"CS",cs_count,"Top PIP",top_pip,"elapsed",elapsed,"sec\n")
    flush.console()
  }
}
result <- do.call(rbind,rows)
write.table(result,file.path(root,"G0022_SUSIE_RSS_PILOT_SUMMARY.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("SuSiE scenarios:",nrow(result),"; converged:",sum(result$converged),"; with CS:",sum(result$cs_count>0),"\n")
