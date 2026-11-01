#!/usr/bin/env Rscript
# Stage28b: independently infer HCCfull→HCCfullz biopsy permutation from
# training genes only, then test on held-out genes whether fullz equals
# per-feature R scale(full), not just rank-preserving. No source identifier
# or permutation is written to artifacts.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage Rscript audit_ng_stage28b_numeric_transform.R OUTPUT.json")
out<-args[[1]]
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
src<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
 HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
get<-function(name) {
 f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
 suppressWarnings(download.file(paste0(src,name),f,mode="wb",method="libcurl",quiet=TRUE))
 sz<-file.info(f)$size
 if(is.na(sz)||sz<100L||sz>16000000L)stop("source bytes unexpected")
 digest<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
 if(length(digest)!=1L||sub(" .*","",digest)!=sha[[name]])stop("source sha drift")
 readRDS(f)
}
rowstandard<-function(x){
 mu<-rowMeans(x);sdv<-apply(x,1,sd)
 if(any(!is.finite(sdv))||any(sdv<1e-12))stop("zero/invalid variance")
 sweep(sweep(x,1,mu,"-"),1,sdv,"/")
}
map<-function(x,y) {
 a<-rowstandard(x);b<-rowstandard(y)
 scores<-crossprod(a,b)/nrow(a)
 max.col(scores,ties.method="first")
}
residual<-function(x,y){
 if(!identical(dim(x),dim(y)))stop("matrix shape")
 z<-rowstandard(x)
 r<-abs(z-y)
 list(median=median(r),p95=unname(quantile(r,.95)),
      max=max(r),within1e6=mean(r<1e-6),
      within1e4=mean(r<1e-4),
      all_gene_row_max_lt_1e5=mean(apply(r,1,max)<1e-5))
}
num<-function(v)if(length(v)!=1L||!is.finite(v))"null" else sprintf("%.9f",v)
bval<-function(x)if(isTRUE(x))"true" else "false"
# Synthetic shuffled source with exact row standardized values must pass.
set.seed(281)
toy<-matrix(rnorm(6000),ncol=10)
perm_seed<-c(4L,1L,10L,3L,9L,6L,2L,8L,5L,7L)
toy_z<-rowstandard(toy)[,perm_seed,drop=FALSE]
idx<-seq(1,600,2)
p<-map(toy[idx,,drop=FALSE],toy_z[idx,,drop=FALSE])
if(!identical(p,match(seq_len(10),perm_seed)))stop("source synthetic permutation not recovered")
if(residual(toy[-idx,,drop=FALSE],toy_z[-idx,p,drop=FALSE])$max>1e-10)
 stop("synthetic numeric fullz restoration failed")
full<-get("HCCfull.rds")
fullz<-get("HCCfullz.rds")
one<-function(k) {
 a<-as.matrix(full[[k]]);b<-as.matrix(fullz[[k]])
 if(ncol(a)!=50L||ncol(b)!=50L||!is.numeric(a)||!is.numeric(b))
    stop("source cardinality or type")
 if(anyDuplicated(rownames(a))||anyDuplicated(rownames(b))||
    anyDuplicated(colnames(a))||anyDuplicated(colnames(b)))stop("source axes duplicated")
 rows<-intersect(rownames(a),rownames(b))
 cols<-intersect(colnames(a),colnames(b))
 if(length(rows)<1000L||length(cols)!=50L)stop("no source identity overlap")
 a<-a[rows,cols,drop=FALSE];b<-b[rows,cols,drop=FALSE]
 if(any(!is.finite(a))||any(!is.finite(b)))stop("nonfinite sources")
 ix<-order(rownames(a))
 a<-a[ix,,drop=FALSE];b<-b[ix,,drop=FALSE]
 train<-seq(1,nrow(a),by=2L);test<-seq(2,nrow(a),by=2L)
 p<-map(a[train,,drop=FALSE],b[train,,drop=FALSE])
 if(!identical(sort(p),seq_len(50L)))stop("nonbijective source permutation")
 before<-residual(a[test,,drop=FALSE],b[test,,drop=FALSE])
 after<-residual(a[test,,drop=FALSE],b[test,p,drop=FALSE])
 paste0('{"modality":"',k,'",',
 '"common_features_n":',nrow(a),',',
 '"training_features_n":',length(train),',',
 '"test_features_n":',length(test),',',
 '"sample_n":50,',
 '"mapping_bijective":true,',
 '"unchanged_named_positions_n":',sum(p==seq_len(50L)),',',
 '"heldout_median_abs_z_residual_before":',num(before$median),',',
 '"heldout_median_abs_z_residual_after":',num(after$median),',',
 '"heldout_p95_abs_z_residual_after":',num(after$p95),',',
 '"heldout_max_abs_z_residual_after":',num(after$max),',',
 '"heldout_cell_fraction_abs_lt_1e6_after":',num(after$within1e6),',',
 '"heldout_cell_fraction_abs_lt_1e4_after":',num(after$within1e4),',',
 '"heldout_feature_fraction_row_max_abs_lt_1e5_after":',num(after$all_gene_row_max_lt_1e5),',',
 '"sample_feature_identifiers_exported":false}')
}
layers<-c("T","P","PP")
payload<-paste0('{"schema":"hcc-stage28b-fullz-numeric-reconstruction-v1",',
 '"source_commit":"',ref,'",',
 '"full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"synthetic_numeric_reconstruction_passed":true,',
 '"layers":[',paste(vapply(layers,one,character(1)),collapse=","),'],',
 '"source_original_Ng_tumor_normal_2022_values_reproduced":false,',
 '"metadata_verified_canonical_sample_ids":false,',
 '"source_files_overwritten":false,',
 '"public_artifact_contains_sample_gene_or_permutation_vector":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(payload,out,useBytes=TRUE)
cat(payload,"\n")
