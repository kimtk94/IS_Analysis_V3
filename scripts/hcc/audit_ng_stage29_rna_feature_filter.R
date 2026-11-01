#!/usr/bin/env Rscript
# Stage29: inspect why exactly 588 RNA feature rows of HCCfull are absent in
# HCCfullz. Aggregate counts and quantiles ONLY; no sample or gene identifiers.
argv<-commandArgs(trailingOnly=TRUE)
if(length(argv)!=1L)stop("Usage: Rscript audit_ng_stage29_rna_feature_filter.R OUT.json")
out<-argv[[1]]
COMMIT<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
BASE<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",COMMIT,"/HCCinputs/")
HASH<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
get<-function(name){
 f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
 suppressWarnings(download.file(paste0(BASE,name),f,method="libcurl",mode="wb",quiet=TRUE))
 sz<-file.info(f)$size
 if(is.na(sz)||sz<100L||sz>16000000L)stop("file size violation")
 dig<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
 if(length(dig)!=1L||sub(" .*","",dig)!=HASH[[name]])stop("checksum drift")
 readRDS(f)
}
a<-get("HCCfull.rds");b<-get("HCCfullz.rds")
if(!is.list(a)||!is.list(b)||is.null(a$T)||is.null(b$T))stop("source T missing")
x<-as.matrix(a$T);y<-as.matrix(b$T)
if(!is.numeric(x)||!is.numeric(y)||nrow(x)!=20356L||nrow(y)!=19768L||
   ncol(x)!=50L||ncol(y)!=50L||anyDuplicated(rownames(x))||
   anyDuplicated(rownames(y))||anyDuplicated(colnames(x))||anyDuplicated(colnames(y)))
 stop("source geometry drift")
present<-rownames(x)%in%rownames(y)
if(sum(present)!=19768L||sum(!present)!=588L)stop("unexpected source RNA filtering")
# Source row-level SD is invariant to sample reordering. Do not impute or
# assume source fullz is an original RNA count dataset.
sx<-apply(x,1,sd)
bad<-!is.finite(sx)
if(any(bad))stop("source row SD NA/Inf")
s0<-sx[!present];s1<-sx[present]
vals<-function(s)list(
 n=length(s),zero_n=sum(s==0),near_zero_lt_1e_12_n=sum(s<1e-12),
 below_1e_6_n=sum(s<1e-6),
 min=min(s),p10=unname(quantile(s,.1)),median=median(s),
 p90=unname(quantile(s,.9)),max=max(s))
excluded<-vals(s0);kept<-vals(s1)
# Test exact variance filter hypothesis, independently without exposing rows.
retained_if_nonzero<-sum(sx>0)
retained_if_gt_eps<-sum(sx>1e-12)
all_excluded_constant<-all(s0==0)
const_genes_kept<-sum(s1==0)
tojson<-function(o){
 num<-function(x)if(!is.finite(x))"null" else sprintf("%.9f",x)
 paste0('{"n":',o$n,',"zero_sd_n":',o$zero_n,
  ',"sd_lt_1e_minus12_n":',o$near_zero_lt_1e_12_n,
  ',"sd_lt_1e_minus6_n":',o$below_1e_6_n,
  ',"sd_min":',num(o$min),',"sd_p10":',num(o$p10),
  ',"sd_median":',num(o$median),',"sd_p90":',num(o$p90),
  ',"sd_max":',num(o$max),'}')
}
payload<-paste0('{"schema":"hcc-stage29-rna-fullz-exclusion-qc-v1",',
 '"source_commit":"',COMMIT,'",',
 '"full_sha256":"',HASH[["HCCfull.rds"]],'",',
 '"fullz_sha256":"',HASH[["HCCfullz.rds"]],'",',
 '"rna_full_feature_n":20356,',
 '"rna_fullz_feature_n":19768,',
 '"strict_matching_feature_n":19768,',
 '"rna_full_only_feature_n":588,',
 '"feature_sample_n":50,',
 '"excluded_rna_feature_sd":',tojson(excluded),',',
 '"retained_rna_feature_sd":',tojson(kept),',',
 '"rna_full_nonzero_sd_feature_n":',retained_if_nonzero,',',
 '"rna_full_gt_1e_minus12_sd_feature_n":',retained_if_gt_eps,',',
 '"all_588_excluded_features_exactly_constant":',
 tolower(as.character(all_excluded_constant)),',',
 '"retained_constant_feature_n":',const_genes_kept,',',
 '"feature_selection_author_code_verified":false,',
 '"per_gene_or_patient_values_exported":false,',
 '"source_rds_files_overwritten":false}')
dir.create(dirname(out),showWarnings=FALSE,recursive=TRUE)
writeLines(payload,out,useBytes=TRUE);cat(payload,"\n")
