#!/usr/bin/env Rscript
# Stage29 source fixed: 588 full RNA rows missing from HCCfullz.
# Check low-variance, row-zero, mean/abundance filtering without disclosing IDs.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("usage: Rscript audit_ng_stage29_filtered_588.R OUT.json")
out<-args[[1]]
commit<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
 HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
sourcebase<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",commit,"/HCCinputs/")
get<-function(n){
 f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
 suppressWarnings(download.file(paste0(sourcebase,n),f,method="libcurl",mode="wb",quiet=TRUE))
 if(is.na(file.info(f)$size)||file.info(f)$size>16000000||file.info(f)$size<100)stop("source file size")
 h<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
 if(length(h)!=1L||sub(" .*","",h)!=sha[[n]])stop("source hash mismatch")
 readRDS(f)
}
a<-as.matrix(get("HCCfull.rds")$T);b<-as.matrix(get("HCCfullz.rds")$T)
if(ncol(a)!=50L||ncol(b)!=50L||nrow(a)!=20356L||nrow(b)!=19768L||
   anyDuplicated(rownames(a))||anyDuplicated(rownames(b))||
   !all(rownames(b)%in%rownames(a)))stop("RNA source feature geometry invalid")
miss<-setdiff(rownames(a),rownames(b))
if(length(miss)!=588L)stop("exact 588 omitted features changed")
sdv<-apply(a,1,sd)
mu<-rowMeans(a)
if(any(!is.finite(sdv))||any(!is.finite(mu)))stop("source numeric nonfinite")
drop<-rownames(a)%in%miss
selected<-!drop
num<-function(x)if(length(x)==1L&&is.finite(x))sprintf("%.9f",x)else "null"
frac<-function(x)num(mean(x))
counts_lowest<-sum(drop[order(sdv)[seq_len(588)]])
counts_highest<-sum(drop[order(mu)[seq_len(588)]])
# Perfect variance-threshold segregation if all retained stdev > all excluded.
strict_sd_separation<-max(sdv[drop])<min(sdv[selected])
strict_mean_separation<-max(mu[drop])<min(mu[selected])
qstats<-function(v)list(min=min(v),med=median(v),max=max(v),q10=unname(quantile(v,.1)),q90=unname(quantile(v,.9)))
put<-function(name,x)paste0('"',name,'":{',
 '"min":',num(x$min),',"median":',num(x$med),',"max":',num(x$max),',
 '"q10":',num(x$q10),',"q90":',num(x$q90),'}')
str<-paste0('{"schema":"hcc-stage29-source-filter-588-transcripts-v1",',
 '"source_commit":"',commit,'",',
 '"full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"full_RNA_features_n":20356,',
 '"fullz_RNA_features_n":19768,',
 '"fullz_RNA_missing_from_full_n":0,',
 '"omitted_RNA_features_n":588,',
 '"retained_RNA_features_n":19768,',
 '"dropped_zero_standard_deviation_n":',sum(sdv[drop]==0),',',
 '"retained_zero_standard_deviation_n":',sum(sdv[selected]==0),',',
 '"dropped_in_bottom_588_by_SD_n":',counts_lowest,',',
 '"dropped_in_bottom_588_by_rowMean_n":',counts_highest,',',
 '"SD_strict_threshold_partition_possible":',tolower(as.character(strict_sd_separation)),',',
 '"mean_strict_threshold_partition_possible":',tolower(as.character(strict_mean_separation)),',',
 put("dropped_SD_summary",qstats(sdv[drop])),',',
 put("retained_SD_summary",qstats(sdv[selected])),',',
 put("dropped_mean_summary",qstats(mu[drop])),',',
 put("retained_mean_summary",qstats(mu[selected])),',',
 '"original_author_filter_rule_source_verified":false,',
 '"filtered_gene_identifiers_or_source_expression_exported":false,',
 '"clinical_significance_evaluated":false}')
dir.create(dirname(out),showWarnings=FALSE,recursive=TRUE)
writeLines(str,out,useBytes=TRUE)
cat(str,"\n")
