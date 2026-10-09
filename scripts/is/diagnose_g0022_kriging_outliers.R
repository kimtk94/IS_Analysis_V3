#!/usr/bin/env Rscript
# SNP-level conditional-z / allele-flip diagnostics, no automatic allele changes.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0) args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
cases<-list.files(root,pattern="^input_qc\\.json$",recursive=TRUE,full.names=TRUE)
if(length(cases)!=5L)stop("Expected 5 real signed-LD pilot windows")
rows<-list();detail<-list()
for(f in cases){
  directory<-dirname(f);signal<-basename(directory)
  z<-as.numeric(read.table(file.path(directory,"z.tsv"),header=TRUE)$z)
  m<-as.matrix(read.table(gzfile(file.path(directory,"signed_ld.tsv.gz")),sep="\t",header=FALSE))
  v<-read.table(file.path(directory,"variants.tsv"),sep="\t",header=TRUE,stringsAsFactors=FALSE)
  if(nrow(v)!=length(z)||length(z)!=nrow(m)||length(z)!=ncol(m))stop("Mismatched variants/Z/LD")
  kr<-tryCatch(kriging_rss(z,m,n=20000),
      error=function(e)e)
  if(inherits(kr,"error")){
    rows[[length(rows)+1L]]<-data.frame(signal=signal,
     status="KRIGING_FAILED",n_variants=length(z),
     flagged_loglr2_absz2=NA_integer_,top_loglr=NA_real_,
     top_outlier_variant="",note=conditionMessage(kr))
    next
  }
  if(is.null(kr$conditional_dist))stop("kriging output structure changed")
  tab<-as.data.frame(kr$conditional_dist)
  if(nrow(tab)!=length(z)||!("logLR"%in%colnames(tab)))
    stop("Invalid kriging output dimensions")
  loglr<-as.numeric(tab$logLR)
  if(any(!is.finite(loglr)))stop("Nonfinite logLR")
  flag<-loglr>2 & abs(z)>2
  top<-which.max(loglr)
  rows[[length(rows)+1L]]<-data.frame(
    signal=signal,status="DIAGNOSTIC_COMPLETE",n_variants=length(z),
    flagged_loglr2_absz2=sum(flag),
    top_loglr=max(loglr),
    top_outlier_variant=v$variant_id[top],
    note="OUTLIER_IS_NOT_PROVEN_ALLELE_FLIP")
  temp<-data.frame(
    signal=signal,variant_id=v$variant_id,
    pos=v$pos,observed_z=z,
    reference_alt_af=v$reference_alt_af,
    gwas_alt_eaf=v$gwas_alt_eaf,
    logLR=loglr,
    suspected_isolated_flip=as.integer(flag),
    annotation="DO_NOT_AUTO_FLIP_OR_REMOVE")
  detail[[length(detail)+1L]]<-temp
  cat(signal,"flagged",sum(flag),"top",v$variant_id[top],"logLR",max(loglr),"\n")
  flush.console()
}
summary<-do.call(rbind,rows)
write.table(summary,file.path(root,"G0022_KRIGING_OUTLIER_SUMMARY.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE)
if(length(detail)){
  all<-do.call(rbind,detail)
  all<-all[order(all$signal,-all$logLR),]
  write.table(all,file.path(root,"G0022_KRIGING_SNP_DIAGNOSTICS.tsv"),
    sep="\t",quote=FALSE,row.names=FALSE)
}
cat("TOTAL_FLAGGED",sum(summary$flagged_loglr2_absz2,na.rm=TRUE),"\n")
