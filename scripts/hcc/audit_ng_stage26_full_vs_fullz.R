#!/usr/bin/env Rscript
# Stage26: empirical, source-pinned HCCfull vs HCCfullz matrix lineage.
# No sample, patient, gene identifiers, per-feature or expression matrices in artifact.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("usage: Rscript audit_ng_stage26_full_vs_fullz.R output.json")
out<-args[[1]]
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
base<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
files<-c("HCCfull.rds","HCCfullz.rds")
sha_expected_full<-"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
getsource<-function(filename){
 p<-tempfile(fileext=".rds")
 on.exit(unlink(p),add=TRUE)
 suppressWarnings(download.file(paste0(base,filename),p,quiet=TRUE,mode="wb",method="libcurl"))
 if(is.na(file.info(p)$size)||file.info(p)$size<100L||file.info(p)$size>16000000L)
  stop("unbounded or empty source RDS")
 sha<-suppressWarnings(system2("sha256sum",p,stdout=TRUE,stderr=FALSE))
 if(length(sha)!=1L||!grepl("^[a-f0-9]{64} ",sha))stop("SHA256 unavailable")
 list(object=readRDS(p),sha256=sub(" .*","",sha),bytes=as.integer(file.info(p)$size))
}
a<-getsource(files[[1]]);b<-getsource(files[[2]])
if(!identical(a$sha256,sha_expected_full))stop("full source drift")
keys<-c("M","CN","T","P","PP")
if(!is.list(a$object)||!is.list(b$object)||
   !all(keys%in%names(a$object))||!all(keys%in%names(b$object)))
 stop("source list schema changed")
fmt<-function(x)if(length(x)!=1L||!is.finite(x))"null" else sprintf("%.9f",x)
bval<-function(x)if(isTRUE(x))"true" else "false"
layer<-function(x,y,k){
 if(!(is.matrix(x)||is.data.frame(x))||!(is.matrix(y)||is.data.frame(y)))
   stop("source matrix type changed")
 xa<-as.matrix(x);ya<-as.matrix(y)
 if(!is.numeric(xa)||!is.numeric(ya))stop("non-numeric layer")
 ids_row<-identical(rownames(xa),rownames(ya))
 ids_col<-identical(colnames(xa),colnames(ya))
 if(is.null(rownames(xa))||is.null(rownames(ya))||
    is.null(colnames(xa))||is.null(colnames(ya)))
  stop("missing dimensional identity in source matrix")
 if(anyDuplicated(rownames(xa))||anyDuplicated(rownames(ya))||
    anyDuplicated(colnames(xa))||anyDuplicated(colnames(ya)))
  stop("duplicate source feature or sample identities")
 if(ncol(xa)!=50L||ncol(ya)!=50L)stop("unexpected source sample cardinality")
 common_rows<-intersect(rownames(xa),rownames(ya))
 common_cols<-intersect(colnames(xa),colnames(ya))
 if(length(common_rows)<20L||length(common_cols)!=50L)
  stop("unverifiable common-source feature or patient alignment")
 # Align by exactly matching original feature/sample identities. Row sets may
 # differ because HCCfullz can have its own source feature selection.
 full_features_n<-nrow(xa);fullz_features_n<-nrow(ya)
 full_samples_n<-ncol(xa);fullz_samples_n<-ncol(ya)
 sample_set_identical<-setequal(colnames(xa),colnames(ya))
 feature_set_identical<-setequal(rownames(xa),rownames(ya))
 xa<-xa[common_rows,common_cols,drop=FALSE]
 ya<-ya[common_rows,common_cols,drop=FALSE]
 both<-is.finite(xa)&is.finite(ya)
 # Exact transformation comparison requires same finite positions
 both_finite<-sum(both)
 if(both_finite<1000)stop("insufficient finite pairs")
 xmean<-rowMeans(xa,na.rm=TRUE)
 xsd<-apply(xa,1,sd,na.rm=TRUE)
 yz<-sweep(sweep(xa,1,xmean,"-"),1,xsd,"/")
 mask<-both&is.finite(yz)
 res<-abs(yz[mask]-ya[mask])
 exact<-both & (xa==ya)
 # Distinguish row-z feature scaling vs column-z sample scaling.
 cmean<-colMeans(xa,na.rm=TRUE)
 csd<-apply(xa,2,sd,na.rm=TRUE)
 ycol<-sweep(sweep(xa,2,cmean,"-"),2,csd,"/")
 mc<-both&is.finite(ycol)
 cres<-abs(ycol[mc]-ya[mc])
 list(
  modality=k,features=as.integer(nrow(xa)),samples=as.integer(ncol(xa)),
  source_full_features_n=as.integer(full_features_n),
  source_fullz_features_n=as.integer(fullz_features_n),
  source_full_samples_n=as.integer(full_samples_n),
  source_fullz_samples_n=as.integer(fullz_samples_n),
  source_common_feature_n=as.integer(length(common_rows)),
  feature_set_identical=feature_set_identical,
  sample_set_identical=sample_set_identical,
  row_and_column_names_identical=ids_row&&ids_col,
  both_finite_n=as.integer(both_finite),
  full_nonfinite_n=as.integer(sum(!is.finite(xa))),
  z_nonfinite_n=as.integer(sum(!is.finite(ya))),
  identical_numeric_cells_n=as.integer(sum(exact,na.rm=TRUE)),
  full_feature_sd_zero_n=as.integer(sum(xsd==0,na.rm=TRUE)),
  full_feature_sd_positive_n=as.integer(sum(xsd>0,na.rm=TRUE)),
  fullz_gene_row_mean_median=as.numeric(median(rowMeans(ya,na.rm=TRUE),na.rm=TRUE)),
  fullz_gene_row_sd_median=as.numeric(median(apply(ya,1,sd,na.rm=TRUE),na.rm=TRUE)),
  row_z_valid_cells_n=as.integer(sum(mask)),
  row_z_median_abs_residual=if(length(res))as.numeric(median(res)) else NA_real_,
  row_z_max_abs_residual=if(length(res))as.numeric(max(res)) else NA_real_,
  row_z_frac_abs_lt_1e_5=if(length(res))as.numeric(mean(res<1e-5)) else NA_real_,
  col_z_valid_cells_n=as.integer(sum(mc)),
  col_z_median_abs_residual=if(length(cres))as.numeric(median(cres)) else NA_real_,
  col_z_frac_abs_lt_1e_5=if(length(cres))as.numeric(mean(cres<1e-5)) else NA_real_
 )
}
stat<-lapply(keys,function(k)layer(a$object[[k]],b$object[[k]],k))
field<-function(x)paste0('{',
 '"modality":"',x$modality,'","features":',x$features,',"samples":',x$samples,
 ',"source_full_features_n":',x$source_full_features_n,
 ',"source_fullz_features_n":',x$source_fullz_features_n,
 ',"source_full_samples_n":',x$source_full_samples_n,
 ',"source_fullz_samples_n":',x$source_fullz_samples_n,
 ',"source_common_feature_n":',x$source_common_feature_n,
 ',"feature_set_identical":',bval(x$feature_set_identical),
 ',"sample_set_identical":',bval(x$sample_set_identical),
 ',"row_and_column_names_identical":',bval(x$row_and_column_names_identical),
 ',"both_finite_n":',x$both_finite_n,
 ',"full_nonfinite_n":',x$full_nonfinite_n,
 ',"z_nonfinite_n":',x$z_nonfinite_n,
 ',"identical_numeric_cells_n":',x$identical_numeric_cells_n,
 ',"full_feature_sd_zero_n":',x$full_feature_sd_zero_n,
 ',"full_feature_sd_positive_n":',x$full_feature_sd_positive_n,
 ',"fullz_gene_row_mean_median":',fmt(x$fullz_gene_row_mean_median),
 ',"fullz_gene_row_sd_median":',fmt(x$fullz_gene_row_sd_median),
 ',"row_z_valid_cells_n":',x$row_z_valid_cells_n,
 ',"row_z_median_abs_residual":',fmt(x$row_z_median_abs_residual),
 ',"row_z_max_abs_residual":',fmt(x$row_z_max_abs_residual),
 ',"row_z_frac_abs_lt_1e_5":',fmt(x$row_z_frac_abs_lt_1e_5),
 ',"col_z_valid_cells_n":',x$col_z_valid_cells_n,
 ',"col_z_median_abs_residual":',fmt(x$col_z_median_abs_residual),
 ',"col_z_frac_abs_lt_1e_5":',fmt(x$col_z_frac_abs_lt_1e_5),'}')
encoded<-paste(vapply(stat,field,character(1)),collapse=",")
result<-paste0('{"schema":"hcc-stage26-full-vs-fullz-lineage-v1",',
 '"source_commit":"',ref,'",',
 '"full_sha256":"',a$sha256,'","fullz_sha256":"',b$sha256,'",',
 '"full_bytes":',a$bytes,',"fullz_bytes":',b$bytes,',',
 '"all_source_layer_names_identical":true,',
 '"layer_statistics":[',encoded,'],',
 '"zscore_recognition_is_empirical_not_original_source_code_confirmation":true,',
 '"original_ng_RSEM_TMM_edgeR_reproduced":false,',
 '"original_ng_5481_gene_universe_reproduced":false,',
 '"source_patient_and_gene_identifiers_exported":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(result,out,useBytes=TRUE)
cat(result,"\n")
