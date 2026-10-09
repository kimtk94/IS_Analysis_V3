#!/usr/bin/env Rscript
# Diagnose EAS summary-statistics / 1000G EAS LD mismatch in 5 G0022 pilot windows.
# Larger lambda indicates a poorer fit of GWAS z scores to reference LD.
suppressPackageStartupMessages(library(susieR))
argv <- commandArgs(trailingOnly=TRUE)
root <- if(length(argv)>0) argv[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
files <- list.files(root,pattern="^input_qc\\.json$",recursive=TRUE,full.names=TRUE)
if(length(files)!=5L) stop("Require 5 signed-LD input windows")
out <- list()
for(f in files) {
  folder <- dirname(f)
  name <- basename(folder)
  z <- as.numeric(read.table(file.path(folder,"z.tsv"),header=TRUE)[,1])
  R <- as.matrix(read.table(gzfile(file.path(folder,"signed_ld.tsv.gz")),
                            header=FALSE,sep="\t"))
  if(any(!is.finite(R)) || nrow(R)!=length(z) || ncol(R)!=length(z))
    stop("Invalid z/R dimensions")
  for(n in c(20000,256274)) {
    est <- NA_real_
    error <- ""
    t0 <- proc.time()[["elapsed"]]
    est <- tryCatch(estimate_s_rss(z,R,n=n,method="null-mle"),
                    error=function(e){error <<- conditionMessage(e);NA_real_})
    sec <- round(proc.time()[["elapsed"]]-t0,2)
    status <- if(!is.finite(est)) "DIAGNOSTIC_FAILED" else
      if(est>0.2) "VERY_HIGH_LD_MISMATCH_REVIEW" else
      if(est>0.05) "LD_MISMATCH_REVIEW" else
      "LOW_ESTIMATED_MISMATCH_NOT_VALIDATION"
    out[[length(out)+1L]] <- data.frame(
      signal=name,hypothetical_n=n,ld_mismatch_s=est,
      mismatch_status=status,elapsed_sec=sec,error_message=error,
      final_fine_mapping_ready=FALSE,
      stringsAsFactors=FALSE)
    cat(name,n,est,status,"\n");flush.console()
  }
}
report <- do.call(rbind,out)
write.table(report,file.path(root,"G0022_GWAS_LD_CONSISTENCY.tsv"),
            row.names=FALSE,quote=FALSE,sep="\t")
cat("LD consistency diagnostics complete",nrow(report),"\n")
