#!/usr/bin/env Rscript
# Stage28: independently infer full/fullz biopsy-column permutation on TRAIN
# features, validate held-out feature rank recovery, check T/P/PP permutation
# agreement. NO source matrices/subject or feature IDs/permutation vector in JSON.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript audit_ng_stage28_crossmodal_permutation.R OUT.json")
out<-args[[1]]
commit<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
base<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",commit,"/HCCinputs/")
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
       HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
readsource<-function(name) {
  f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
  suppressWarnings(download.file(paste0(base,name),f,quiet=TRUE,mode="wb",method="libcurl"))
  size<-file.info(f)$size
  if(is.na(size)||size<100L||size>16000000L)stop("source size out of range")
  digest<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
  if(length(digest)!=1L||sub(" .*","",digest)!=sha[[name]])stop("source SHA-256 drift")
  readRDS(f)
}
fmt<-function(x)if(length(x)!=1L||!is.finite(x))"null" else sprintf("%.9f",x)
bool<-function(x)if(isTRUE(x))"true" else "false"
similarity<-function(x,y) {
  if(!identical(dim(x),dim(y)))stop("source aligned matrices differ in shape")
  sx<-apply(x,1,sd);sy<-apply(y,1,sd)
  good<-is.finite(sx)&is.finite(sy)&sx>1e-12&sy>1e-12
  if(sum(good)<30L)stop("too few valid training features")
  x<-x[good,,drop=FALSE];y<-y[good,,drop=FALSE]
  sx<-sx[good];sy<-sy[good]
  zx<-sweep(sweep(x,1,rowMeans(x),"-"),1,sx,"/")
  zy<-sweep(sweep(y,1,rowMeans(y),"-"),1,sy,"/")
  # Scale across feature mean; score is NOT bounded by 1.
  crossprod(zx,zy)/nrow(zx)
}
derive<-function(x,y)max.col(similarity(x,y),ties.method="first")
rhos<-function(x,y,perm) {
  # Source feature row identity and common 50 columns already aligned.
  # Spearman each gene across 50 biopsies, summarize only; use for HOLDOUT.
  n<-nrow(x);a<-rep(NA_real_,n);b<-rep(NA_real_,n)
  for(i in seq_len(n)){
    if(sd(x[i,])<1e-12||sd(y[i,])<1e-12)next
    a[i]<-suppressWarnings(cor(x[i,],y[i,],method="spearman"))
    b[i]<-suppressWarnings(cor(x[i,],y[i,perm],method="spearman"))
  }
  valid<-is.finite(a)&is.finite(b)
  if(sum(valid)<30L)stop("too few evaluation features")
  v<-a[valid];w<-b[valid]
  list(n=sum(valid),before=median(v),after=median(w),
       q10=unname(quantile(w,.1)),q90=unname(quantile(w,.9)),
       improvement=median(w-v),gain=median(w)-median(v),
       after_gt_0_95=mean(w>.95),before_gt_0_95=mean(v>.95),
       after_negative=mean(w<0))
}
# Deterministic synthetic regression: sample permutation across two layers,
# independently split training/held out by feature, should recover 100% mapping.
set.seed(20261031)
toy<-matrix(rnorm(6000),nrow=600,ncol=10)
perm_seed<-c(3L,1L,8L,6L,2L,10L,5L,7L,9L,4L)
toy_z<-toy[,perm_seed,drop=FALSE]
fitted<-derive(toy[seq(1,599,2),,drop=FALSE],toy_z[seq(1,599,2),,drop=FALSE])
expect<-match(seq_len(10),perm_seed)
if(!identical(fitted,expect))stop("synthetic crossed-column permutation regression FAILED")
eval<-rhos(toy[seq(2,600,2),,drop=FALSE],toy_z[seq(2,600,2),,drop=FALSE],fitted)
if(eval$after<.9999||eval$before>.5)stop("synthetic heldout recovery regression FAILED")
full<-readsource("HCCfull.rds")
fullz<-readsource("HCCfullz.rds")
layers<-c("M","CN","T","P","PP")
permutations<-list()
results<-list()
source_summary<-function(k) {
  a<-as.matrix(full[[k]]);b<-as.matrix(fullz[[k]])
  if(!is.numeric(a)||!is.numeric(b)||ncol(a)!=50L||ncol(b)!=50L)
    stop("unexpected sample cardinality")
  if(is.null(rownames(a))||is.null(rownames(b))||is.null(colnames(a))||
     is.null(colnames(b))||anyDuplicated(rownames(a))||
     anyDuplicated(rownames(b))||anyDuplicated(colnames(a))||
     anyDuplicated(colnames(b)))stop("unsafe duplicated source label")
  idr<-intersect(rownames(a),rownames(b))
  idc<-intersect(colnames(a),colnames(b))
  if(length(idr)<100L||length(idc)!=50L)stop("insufficient source axis overlap")
  a<-a[idr,idc,drop=FALSE];b<-b[idr,idc,drop=FALSE]
  if(any(!is.finite(a))||any(!is.finite(b)))stop("unbounded source data")
  # Deterministic exact independent feature split, no data-dependent choice.
  o<-order(rownames(a))
  a<-a[o,,drop=FALSE];b<-b[o,,drop=FALSE]
  train<-seq(1,nrow(a),by=2L)
  test<-seq(2,nrow(a),by=2L)
  if(length(train)<30L||length(test)<30L)stop("insufficient split")
  p<-derive(a[train,,drop=FALSE],b[train,,drop=FALSE])
  p2<-derive(a[test,,drop=FALSE],b[test,,drop=FALSE])
  isperm<-identical(sort(p),seq_len(50L))
  isperm2<-identical(sort(p2),seq_len(50L))
  perm_agree<-sum(p==p2)
  holdout<-rhos(a[test,,drop=FALSE],b[test,,drop=FALSE],p)
  reverse<-rhos(a[train,,drop=FALSE],b[train,,drop=FALSE],p2)
  # Source sample IDs and inferred mappings remain only in runner process memory.
  permutations[[k]]<<-p
  results[[k]]<<-list(
    feature_n=nrow(a),train_n=length(train),heldout_n=length(test),
    mapping_is_bijection=isperm,
    reverse_mapping_is_bijection=isperm2,
    train_and_reverse_holdout_matches_n=perm_agree,
    unchanged_named_positions_n=sum(p==seq_len(50L)),
    holdout_rho=holdout,reverse_heldout_rho=reverse)
}
for(k in layers)source_summary(k)
encodeRho<-function(x,prefix)paste0(
 '"',prefix,'_finite_feature_n":',x$n,',',
 '"',prefix,'_rho_before_median":',fmt(x$before),',',
 '"',prefix,'_rho_after_median":',fmt(x$after),',',
 '"',prefix,'_rho_after_q10":',fmt(x$q10),',',
 '"',prefix,'_rho_after_q90":',fmt(x$q90),',',
 '"',prefix,'_median_within_feature_rho_delta":',fmt(x$improvement),',',
 '"',prefix,'_frac_rho_gt_095_before":',fmt(x$before_gt_0_95),',',
 '"',prefix,'_frac_rho_gt_095_after":',fmt(x$after_gt_0_95),',',
 '"',prefix,'_frac_rho_negative_after":',fmt(x$after_negative))
encodeLayer<-function(k){
  v<-results[[k]]
  paste0('{"modality":"',k,'",',
  '"features_common_n":',v$feature_n,',',
  '"train_feature_n":',v$train_n,',',
  '"heldout_feature_n":',v$heldout_n,',',
  '"sample_n":50,',
  '"train_map_is_one_to_one":',bool(v$mapping_is_bijection),',',
  '"reverse_heldout_map_is_one_to_one":',bool(v$reverse_mapping_is_bijection),',',
  '"train_vs_reverse_map_agreement_n":',v$train_and_reverse_holdout_matches_n,',',
  '"unchanged_named_positions_n":',v$unchanged_named_positions_n,',',
  encodeRho(v$holdout_rho,"test"),',',
  encodeRho(v$reverse_heldout_rho,"reverse_test"),',',
  '"reordered_data_persisted":false,',
  '"sample_permutation_or_gene_identifiers_exported":false}')
}
pairs<-list(c("T","P"),c("T","PP"),c("P","PP"),c("M","CN"))
same_map<-function(pair){
  a<-pair[[1]];b<-pair[[2]]
  paste0('{"modality_a":"',a,'","modality_b":"',b,'",',
  '"matching_sample_map_positions_n":',sum(permutations[[a]]==permutations[[b]]),',',
  '"both_maps_bijective":',bool(identical(sort(permutations[[a]]),seq_len(50L))&&
    identical(sort(permutations[[b]]),seq_len(50L))),'}')
}
outjson<-paste0(
 '{"schema":"hcc-stage28-crossmodal-trainheldout-permutation-v1",',
 '"source_commit":"',commit,'",',
 '"full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"synthetic_permutation_and_heldout_recovery_passed":true,',
 '"each_modality_reported_independently":true,',
 '"modalities":[',paste(vapply(layers,encodeLayer,character(1)),collapse=","),'],',
 '"crossmodal_comparisons":[',paste(vapply(pairs,same_map,character(1)),collapse=","),'],',
 '"publisher_original_biopsy_identity_reconciled":false,',
 '"original_suter_preprocessing_generation_reproduced":false,',
 '"sample_labels_persistently_corrected":false,',
 '"technical_alignment_is_not_a_biological_replication":true,',
 '"per_sample_mapping_and_per_gene_values_exported":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(outjson,out,useBytes=TRUE)
cat(outjson,"\n")
