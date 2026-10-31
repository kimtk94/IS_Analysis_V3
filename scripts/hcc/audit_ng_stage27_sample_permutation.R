#!/usr/bin/env Rscript
# Stage27b: source-matched sample-axis permutation check across paired RNA,
# protein, phosphoprotein, CN and mutation feature matrices. Aggregate only.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript audit_ng_stage27_sample_permutation.R OUTPUT.json")
out<-args[[1]]
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
base<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
 HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
get<-function(n){
 f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
 suppressWarnings(download.file(paste0(base,n),f,mode="wb",quiet=TRUE,method="libcurl"))
 if(is.na(file.info(f)$size)||file.info(f)$size<100||file.info(f)$size>16000000L)
  stop("source size changed")
 got<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
 if(length(got)!=1L||sub(" .*","",got)!=sha[[n]])stop("unpinned source RDS")
 readRDS(f)
}
pair_similarity<-function(a,b){
 if(ncol(a)!=ncol(b)||ncol(a)<3||nrow(a)!=nrow(b))stop("alignment incorrect")
 # Rows are same exact feature IDs. Compare within-feature-normalized sample
 # vectors. Matrix multiplication aggregates over thousands of independent
 # features without exporting any sample-pair similarity matrix.
 za<-rowMeans(a);zb<-rowMeans(b)
 va<-apply(a,1,sd);vb<-apply(b,1,sd)
 ok<-is.finite(va)&is.finite(vb)&va>0&vb>0
 if(sum(ok)<20)stop("too few variable matched features")
 ax<-sweep(sweep(a[ok,,drop=FALSE],1,za[ok],"-"),1,va[ok],"/")
 bx<-sweep(sweep(b[ok,,drop=FALSE],1,zb[ok],"-"),1,vb[ok],"/")
 s<-crossprod(ax,bx)/sum(ok)
 best<-max.col(s,ties.method="first")
 diagonal<-diag(s)
 list(
  variable_features_n=sum(ok),
  diagonal_best_sample_n=sum(best==seq_len(ncol(s))),
  diagonal_best_sample_fraction=mean(best==seq_len(ncol(s))),
  best_sample_indices_form_permutation=identical(sort(best),seq_len(ncol(s))),
  diagonal_similarity_median=median(diagonal),
  diagonal_similarity_min=min(diagonal),
  diagonal_similarity_max=max(diagonal),
  max_crosssample_similarity_median=median(apply(s,1,max)),
  max_crosssample_similarity_gt_0_7_fraction=mean(apply(s,1,max)>0.7),
  source_sample_axis_n=ncol(s)
 )
}
# Internal synthetic test: sample order rotation can be detected.
set.seed(27)
toy<-matrix(rnorm(2000),nrow=200,ncol=10)
colnames(toy)<-paste0("sample",1:10)
rownames(toy)<-paste0("feature",1:200)
same<-pair_similarity(toy,toy)
if(same$diagonal_best_sample_n!=10L||same$diagonal_similarity_median<.89)
 stop("identity sample matching fails synthetic test")
rotated<-toy[,c(2:10,1),drop=FALSE]
mis<-pair_similarity(toy,rotated)
if(mis$diagonal_best_sample_n>1L||mis$max_crosssample_similarity_median<.89)
 stop("synthetic sample permutation detection failed")
full<-get("HCCfull.rds");fullz<-get("HCCfullz.rds")
fmt<-function(x)if(!is.finite(x))"null" else sprintf("%.8f",x)
bval<-function(v)if(isTRUE(v))"true" else "false"
make<-function(k){
 x<-as.matrix(full[[k]]);y<-as.matrix(fullz[[k]])
 if(!is.numeric(x)||!is.numeric(y)||ncol(x)!=50L||ncol(y)!=50L)
  stop("source axes unexpected")
 if(anyDuplicated(rownames(x))||anyDuplicated(rownames(y))||
    anyDuplicated(colnames(x))||anyDuplicated(colnames(y)))
  stop("source identities not unique")
 rows<-intersect(rownames(x),rownames(y))
 cols<-intersect(colnames(x),colnames(y))
 if(length(rows)<100L||length(cols)!=50L)stop("same biopsy IDs not recoverable")
 x<-x[rows,cols,drop=FALSE];y<-y[rows,cols,drop=FALSE]
 if(any(!is.finite(x))||any(!is.finite(y)))stop("nonfinite source values")
 o<-pair_similarity(x,y)
 paste0('{"modality":"',k,'",',
 '"features_in_similarity_n":',o$variable_features_n,',',
 '"samples_n":50,',
 '"diagonal_wins_n":',o$diagonal_best_sample_n,',',
 '"diagonal_wins_fraction":',fmt(o$diagonal_best_sample_fraction),',',
 '"best_matches_form_one_to_one_permutation":',bval(o$best_sample_indices_form_permutation),',',
 '"same_label_similarity_median":',fmt(o$diagonal_similarity_median),',',
 '"same_label_similarity_min":',fmt(o$diagonal_similarity_min),',',
 '"same_label_similarity_max":',fmt(o$diagonal_similarity_max),',',
 '"best_offdiagonal_allowed_similarity_median":',fmt(o$max_crosssample_similarity_median),',',
 '"best_similarity_gt_0_7_fraction":',fmt(o$max_crosssample_similarity_gt_0_7_fraction),',',
 '"no_sample_or_feature_identifiers_exported":true}')
}
layers<-c("M","CN","T","P","PP")
outjson<-paste0('{"schema":"hcc-stage27b-paired-source-sample-match-v1",',
 '"source_commit":"',ref,'",',
 '"source_full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"source_fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"synthetic_identity_permutation_controls_passed":true,',
 '"layer_metrics":[',paste(vapply(layers,make,character(1)),collapse=","),'],',
 '"reordering_feature_or_sample_identifiers_proven_as_cause":false,',
 '"exact_original_preprocessing_code_reproduced":false,',
 '"biological_or_clinical_significance_proven":false,',
 '"all_source_patient_and_feature_identifiers_exported":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(outjson,out,useBytes=TRUE)
cat(outjson,"\n")
