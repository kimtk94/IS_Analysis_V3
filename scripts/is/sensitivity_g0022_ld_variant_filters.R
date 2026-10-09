#!/usr/bin/env Rscript
# Conditional LD mismatch sensitivity to minor-allele frequency / EAF QC.
# Changing the set of SNPs changes the statistical model; improvements do NOT
# prove allele errors and are not evidence that a subset fixes all discrepancies.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
fns<-list.files(root,pattern="^input_qc\\.json$",full.names=TRUE,recursive=TRUE)
n_table<-read.table(file.path(root,"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv"),
  header=TRUE,sep="\t",stringsAsFactors=FALSE)
if(length(fns)!=5L||nrow(n_table)!=2L)stop("Pilot inputs or sample provenance absent")
filters<-list(
  list(name="BASELINE_MAF_001_EAF_DELTA_015",minmaf=.01,maxdelta=.15),
  list(name="STRICT_MAF_005_EAF_DELTA_005",minmaf=.05,maxdelta=.05),
  list(name="VERY_STRICT_MAF_010_EAF_DELTA_003",minmaf=.1,maxdelta=.03)
)
result<-list()
for(p in fns){
  folder<-dirname(p);name<-basename(folder)
  trait<-if(grepl("_AIS_",name,fixed=TRUE))"AIS" else "AS"
  n<-n_table$working_n_rounded[n_table$trait==trait]
  variants<-read.table(file.path(folder,"variants.tsv"),header=TRUE,sep="\t")
  z<-read.table(file.path(folder,"z.tsv"),header=TRUE,sep="\t")$z
  R<-as.matrix(read.table(gzfile(file.path(folder,"signed_ld.tsv.gz")),sep="\t"))
  if(nrow(variants)!=length(z)||length(z)!=nrow(R))stop("Dimension mismatch")
  rf<-as.numeric(variants$reference_alt_af)
  gf<-as.numeric(variants$gwas_alt_eaf)
  for(conf in filters){
    keep<-is.finite(rf)&is.finite(gf)&pmin(rf,1-rf)>=conf$minmaf &
      abs(rf-gf)<=conf$maxdelta
    k<-sum(keep)
    if(k<30L){
      result[[length(result)+1L]]<-data.frame(signal=name,trait=trait,
        filter=conf$name,retained_snps=k,
        relative_variants_kept=k/length(z),estimated_s=NA_real_,
        status="TOO_FEW_VARIANTS",stringsAsFactors=FALSE)
      next
    }
    zsub<-z[keep]; rsub<-R[keep,keep,drop=FALSE]
    s<-tryCatch(estimate_s_rss(zsub,rsub,n=n,method="null-mle"),
       error=function(e)NA_real_)
    state<-if(!is.finite(s))"ESTIMATION_FAILED" else
      if(s>.2)"VERY_HIGH_MISMATCH" else if(s>.05)"MISMATCH_REVIEW" else
      "LOW_ESTIMATED_MISMATCH_NOT_VALIDATION"
    result[[length(result)+1L]]<-data.frame(signal=name,trait=trait,
      filter=conf$name,retained_snps=k,
      relative_variants_kept=k/length(z),estimated_s=s,
      status=state,stringsAsFactors=FALSE)
    cat(name,conf$name,k,s,state,"\n")
  }
}
output<-do.call(rbind,result)
write.table(output,file.path(root,"G0022_LD_MISMATCH_FILTER_SENSITIVITY.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE)
cat("total filters",nrow(output),"\n")
