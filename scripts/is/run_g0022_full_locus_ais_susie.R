#!/usr/bin/env Rscript
# Entire chr12 G0022 AIS locus, 4,649 full-genotype QC SNPs.
# Exploratory SuSiE-RSS: EAS 1000G panel Nref=504, approximated N_eff=70474,
# missing per-variant effective N and LD-reference uncertainty correction.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1"
q<-jsonlite::fromJSON(file.path(root,"FULL_LOCUS_INPUT_QC.json"))
if(q$status!="FULL_G0022_SIGNED_LD_INPUT_QC_PASS"||q$genotype_n!=504)
  stop("Source/genotype reference QC missing")
v<-read.table(file.path(root,"variants.tsv"),header=TRUE,sep="\t",
  stringsAsFactors=FALSE,check.names=FALSE)
p<-nrow(v)
if(p!=q$matrix_dimension || p>6500L)stop("Unbounded full-locus dimension")
n_effective<-70474L
if(file.info(file.path(root,"ld.rowmajor.f64"))$size!=8*p*p)
  stop("Wrong binary matrix byte size")
avail<-scan(text=sub("^MemAvailable:","",grep("^MemAvailable:",
 readLines("/proc/meminfo"),value=TRUE)),what=double(),nmax=1,quiet=TRUE)
if(length(avail)!=1||avail<2.5*1024*1024)stop("Under 2.5GiB MemAvailable; abort")
values<-readBin(file.path(root,"ld.rowmajor.f64"),what="double",
 n=p*p,size=8,endian="little")
R<-matrix(values,nrow=p,ncol=p,byrow=TRUE)
rm(values);gc(verbose=FALSE)
z<-as.numeric(v$z)
if(!all(is.finite(z))||!all(is.finite(R))||max(abs(diag(R)-1))>1e-6 ||
   max(abs(R-t(R)))>1e-6)stop("Full locus signed LD QC failed")
if(!all(q$lead_ids_present %in% v$variant_id))stop("Missing clump signals")
cat("START full locus",p,"SNPs",q$genotype_n,"EAS samples, N_eff approx",n_effective,"\n")
flush.console()
start<-proc.time()[["elapsed"]]
status<-"NOT_FINISHED";fit<-NULL;err<-""
fit<-tryCatch(susie_rss(z=z,R=R,n=n_effective,L=8,
     estimate_residual_variance=FALSE,max_iter=120,coverage=.95,
     check_prior=TRUE,verbose=FALSE),
   error=function(e){err<<-conditionMessage(e);NULL})
elapsed<-proc.time()[["elapsed"]]-start
result<-list(
  status="FULL_LOCUS_EXPLORATORY_ONLY",lead_ids=q$lead_ids_present,
  n_snps=p,ref_n=504,approximate_n_eff=n_effective,
  n_eff_is_per_snp=FALSE,ld_reference_correction="NOT_AVAILABLE_SUSIER_0_14_2",
  multi_ancestry_replication="NOT_TESTED",
  converged=FALSE,n_credible_sets=0L,credible_set_sizes=integer(),
  credible_set_variant_ids=list(),
  max_pip=NA_real_,max_pip_snp="",
  max_pairwise_cs_abs_cor=NA_real_,
  runtime_sec=round(elapsed,2),error=err,
  final_causal_signal_validated=FALSE,
  final_causal_gene_validated=FALSE)
if(!is.null(fit)){
  result$converged<-isTRUE(fit$converged)
  cs<-susie_get_cs(fit,Xcorr=R,coverage=.95,min_abs_corr=.5)
  if(!is.null(cs$cs)){
    result$n_credible_sets<-length(cs$cs)
    result$credible_set_sizes<-vapply(cs$cs,length,integer(1))
    result$credible_set_variant_ids<-lapply(cs$cs,function(ids)v$variant_id[ids])
  }
  pip<-susie_get_pip(fit)
  if(length(pip)!=p)stop("PIP dimension mismatch")
  imax<-which.max(pip)
  result$max_pip<-pip[imax]
  result$max_pip_snp<-v$variant_id[imax]
  write.table(data.frame(variant_id=v$variant_id,pos=v$pos,pip=pip,z=z),
    file.path(root,"G0022_FULL_AIS_PIP.tsv"),sep="\t",
    row.names=FALSE,quote=FALSE)
  saveRDS(fit,file.path(root,"G0022_FULL_AIS_EXPLORATORY_SUSIE.rds"))
}
jsonlite::write_json(result,file.path(root,"G0022_FULL_AIS_SUSIE_SUMMARY.json"),
  auto_unbox=TRUE,pretty=TRUE,digits=10)
cat("END full locus",result$converged,"CS",result$n_credible_sets,
 "maxPIP",result$max_pip,"duration",elapsed,"sec\n")
