#!/usr/bin/env Rscript
# Stage27 per-feature rank order preservation audit: public Suter processed
# HCCfull/HCCfullz matrices, source pinned, 50 biopsy measurements.
# Only aggregate numerical summaries: never log row/sample identifiers or rhos.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage Rscript audit_ng_stage27_fullz_rank_lineage.R OUT.json")
out<-args[[1]]
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
base<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
expected<-c(
  HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
  HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007"
)
get<-function(n){
  p<-tempfile(fileext=".rds")
  on.exit(unlink(p),add=TRUE)
  suppressWarnings(download.file(paste0(base,n),p,mode="wb",method="libcurl",quiet=TRUE))
  sz<-file.info(p)$size
  if(is.na(sz)||sz<100||sz>16000000L)stop("RDS byte guard failed")
  hash<-suppressWarnings(system2("sha256sum",p,stdout=TRUE,stderr=FALSE))
  if(length(hash)!=1L||sub(" .*","",hash)!=expected[[n]])stop("RDS source drift")
  readRDS(p)
}
full<-get("HCCfull.rds")
fullz<-get("HCCfullz.rds")
layers<-c("M","CN","T","P","PP")
if(!is.list(full)||!is.list(fullz)||!all(layers%in%names(full))||!all(layers%in%names(fullz)))
 stop("missing public multimodal list")
fmt<-function(v){
  if(length(v)!=1L||!is.finite(v))return("null")
  sprintf("%.9f",v)
}
jsbool<-function(x)if(isTRUE(x))"true" else "false"
summary_one<-function(k){
  a<-as.matrix(full[[k]])
  b<-as.matrix(fullz[[k]])
  if(!is.numeric(a)||!is.numeric(b)||ncol(a)!=50L||ncol(b)!=50L)
    stop("unexpected modality numeric/sample type")
  if(anyDuplicated(rownames(a))||anyDuplicated(rownames(b))||
     anyDuplicated(colnames(a))||anyDuplicated(colnames(b)))
    stop("nonunique source axis")
  rows<-intersect(rownames(a),rownames(b))
  cols<-intersect(colnames(a),colnames(b))
  if(length(rows)<20L||length(cols)!=50L)stop("unmatched source 50 biopsy IDs")
  a<-a[rows,cols,drop=FALSE]
  b<-b[rows,cols,drop=FALSE]
  if(any(!is.finite(a))||any(!is.finite(b)))
    stop("nonfinite source observations")
  # Do not export per-feature rhos. Both Pearson and Spearman quantify
  # relationship among the SAME SOURCE 50 biopsies for each matching row.
  sp<-vapply(seq_len(nrow(a)),function(i){
    if(length(unique(a[i,]))<2L||length(unique(b[i,]))<2L)return(NA_real_)
    suppressWarnings(cor(a[i,],b[i,],method="spearman"))
  },numeric(1))
  pe<-vapply(seq_len(nrow(a)),function(i){
    if(sd(a[i,])==0||sd(b[i,])==0)return(NA_real_)
    suppressWarnings(cor(a[i,],b[i,],method="pearson"))
  },numeric(1))
  valid<-is.finite(sp)&is.finite(pe)
  if(sum(valid)<20L)stop("too few source matched feature correlations")
  sr<-sp[valid];pr<-pe[valid]
  # Correlation near 1 implies rank preservation, not exact scaling
  # and a high median alone does not establish full source lineage.
  obj<-paste0('{"modality":"',k,'",',
    '"source_common_features_n":',nrow(a),',',
    '"biopsies_n":50,',
    '"rho_finite_features_n":',sum(valid),',',
    '"undefined_rho_feature_n":',nrow(a)-sum(valid),',',
    '"spearman_rank_correlation_median":',fmt(median(sr)),',',
    '"spearman_rank_correlation_p10":',fmt(unname(quantile(sr,.1))),',',
    '"spearman_rank_correlation_p90":',fmt(unname(quantile(sr,.9))),',',
    '"spearman_fraction_ge_0_99":',fmt(mean(sr>=.99)),',',
    '"spearman_fraction_lt_0_9":',fmt(mean(sr<.9)),',',
    '"spearman_fraction_negative":',fmt(mean(sr<0)),',',
    '"pearson_median":',fmt(median(pr)),',',
    '"pearson_fraction_ge_0_99":',fmt(mean(pr>=.99)),',',
    '"feature_or_sample_identifiers_exported":false}')
  obj
}
result<-paste0('{"schema":"hcc-stage27-full-fullz-paired-rank-v1",',
 '"source_repository":"cbg-ethz/HCC",',
 '"source_commit":"',ref,'",',
 '"source_full_sha256":"',expected[["HCCfull.rds"]],'",',
 '"source_fullz_sha256":"',expected[["HCCfullz.rds"]],'",',
 '"layer_metrics":[',paste(vapply(layers,summary_one,character(1)),collapse=","),'],',
 '"exact_source_matrix_generation_reproduced":false,',
 '"original_Figure2b_gene_universe_reproduced":false,',
 '"original_Ng_2022_tumor_normal_relative_quantitation_reproduced":false,',
 '"gene_level_FDR_or_biological_significance_tested":false,',
 '"per_feature_rhos_or_patient_identifiers_exported":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(result,out,useBytes=TRUE)
cat(result,"\n")
