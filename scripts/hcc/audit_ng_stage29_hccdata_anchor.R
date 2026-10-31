#!/usr/bin/env Rscript
# Stage29 independent prebuilt Suter HCCdata.rds (small BN network data)
# cross-check against HCCfull/HCCfullz source 50-biopsy labels and numerical
# profiles. Source is same author repository, NOT independent clinical truth.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript audit_ng_stage29_hccdata_anchor.R OUT.json")
out<-args[[1]]
revision<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
base<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",revision,"/HCCinputs/")
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
 HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
smallBlob<-"45538b8452a1d2e3844b596f6bd50f1086f0bf21"
fetch<-function(n){
  f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
  suppressWarnings(download.file(paste0(base,n),f,mode="wb",method="libcurl",quiet=TRUE))
  sz<-file.info(f)$size
  if(is.na(sz)||sz<100L||sz>16000000L)stop("source file size invalid")
  if(n=="HCCdata.rds"){
    got<-suppressWarnings(system2("git",c("hash-object",shQuote(f)),stdout=TRUE,stderr=FALSE))
    if(length(got)!=1L||got!=smallBlob||sz!=185393L)stop("HCCdata git blob fingerprint drift")
  }else{
    got<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
    if(length(got)!=1L||sub(" .*","",got)!=sha[[n]])stop("source HCCfull hash drift")
  }
  readRDS(f)
}
anchor<-fetch("HCCdata.rds");full<-fetch("HCCfull.rds");fullz<-fetch("HCCfullz.rds")
if(!all(c("T","P","PP")%in%names(anchor))||!all(c("T","P","PP")%in%names(full))||
 !all(c("T","P","PP")%in%names(fullz)))stop("required source modalities missing")
num<-function(x)if(length(x)!=1L||!is.finite(x))"null" else sprintf("%.9f",x)
bool<-function(x)if(isTRUE(x))"true" else "false"
same_gene<-function(a,b){
  # Same feature ID across source products and exact named biopsy columns.
  common<-intersect(rownames(a),rownames(b))
  commonCols<-intersect(colnames(a),colnames(b))
  if(length(common)<10||length(commonCols)!=50L){
    return(list(shared=length(common),axis=length(commonCols),
                finite=0L,rho=NA_real_,positive=NA_real_))
  }
  x<-a[common,commonCols,drop=FALSE];y<-b[common,commonCols,drop=FALSE]
  if(any(!is.finite(x))||any(!is.finite(y)))stop("source nonfinite")
  rho<-vapply(seq_len(nrow(x)),function(i){
    if(sd(x[i,])<1e-12||sd(y[i,])<1e-12)return(NA_real_)
    suppressWarnings(cor(x[i,],y[i,],method="spearman"))
  },numeric(1))
  good<-rho[is.finite(rho)]
  list(shared=length(common),axis=length(commonCols),
       finite=length(good),rho=if(length(good))median(good) else NA_real_,
       positive=if(length(good))mean(good>.95) else NA_real_)
}
rows<-character()
for(k in c("T","P","PP")){
 a<-as.matrix(anchor[[k]]);b<-as.matrix(full[[k]]);c<-as.matrix(fullz[[k]])
 # HCCdata is the authors' small network-object input, NOT a guaranteed
 # quantitative 50-biopsy matrix. Preserve its real structure and report
 # whether same-name exact-gene Spearman is actually computable.
 if(length(dim(a))!=2L||is.null(nrow(a))||is.null(ncol(a))||
    nrow(a)<1L||ncol(a)<1L||
    anyDuplicated(rownames(a))||anyDuplicated(colnames(a)))
     stop("HCCdata structural source matrix invalid")
 has_numeric<-is.numeric(a) && all(is.finite(a))
 comparable<-has_numeric && ncol(a)==50L && !is.null(colnames(a)) &&
    length(intersect(colnames(a),colnames(b)))==50L &&
    length(intersect(colnames(a),colnames(c)))==50L
 if(comparable){
    r<-same_gene(a,b);s<-same_gene(a,c)
 }else{
    # A small network input without 50 quantitative biopsy columns
    # cannot validate per-biopsy labels. This is a scientific HOLD, not
    # a source-download or code-execution failure.
    r<-list(shared=length(intersect(rownames(a),rownames(b))),
       axis=length(intersect(colnames(a),colnames(b))),
       finite=0L,rho=NA_real_,positive=NA_real_)
    s<-list(shared=length(intersect(rownames(a),rownames(c))),
       axis=length(intersect(colnames(a),colnames(c))),
       finite=0L,rho=NA_real_,positive=NA_real_)
 }
 rows<-c(rows,paste0('{"modality":"',k,'",',
 '"small_data_feature_n":',nrow(a),',',
 '"small_data_sample_n":',ncol(a),',',
 '"small_data_source_numeric_matrix":',bool(has_numeric),',',
 '"small_data_exact_50_biopsy_comparison_eligible":',bool(comparable),',',
 '"small_full_common_feature_n":',r$shared,',',
 '"small_full_common_sample_labels_n":',r$axis,',',
 '"small_full_rho_finite_feature_n":',r$finite,',',
 '"small_full_same_label_per_feature_rho_median":',num(r$rho),',',
 '"small_full_fraction_gene_rho_gt_095":',num(r$positive),',',
 '"small_fullz_common_feature_n":',s$shared,',',
 '"small_fullz_common_sample_labels_n":',s$axis,',',
 '"small_fullz_rho_finite_feature_n":',s$finite,',',
 '"small_fullz_same_label_per_feature_rho_median":',num(s$rho),',',
 '"small_fullz_fraction_gene_rho_gt_095":',num(s$positive),',',
 '"biopsy_or_gene_identifiers_exported":false}'))
}
payload<-paste0('{"schema":"hcc-stage29-small-HCCdata-rds-anchor-v1",',
 '"source_commit":"',revision,'",',
 '"source_HCCdata_git_blob_sha":"',smallBlob,'",',
 '"source_full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"source_fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"modalities":[',paste(rows,collapse=","),'],',
 '"HCCdata_is_independent_patient_cohort":false,',
 '"HCCdata_authoritative_biopsy_identity_proven":false,',
 '"source_sample_labels_corrected":false,',
 '"clinical_or_gene_effect_inference_performed":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(payload,out,useBytes=TRUE)
cat(payload,"\n")
