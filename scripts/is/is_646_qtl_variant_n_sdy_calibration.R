#!/usr/bin/env Rscript
# Source-pinned coloc 5.2.3 trait-SDY variance calibration using SNP-specific observed QTL n.
# No per-SNP ABF N correction; beta/SE (varbeta) are retained unchanged.
suppressPackageStartupMessages({library(data.table);library(coloc)})
if(as.character(packageVersion("coloc"))!="5.2.3")stop("PINNED_COLOC_5_2_3_REQUIRED")
a<-commandArgs(trailingOnly=TRUE)
if(length(a)!=4L)stop("USAGE: SCRIPT INDEX MASTER INPUT_DIR NEW_OUTPUT_DIR")
index<-fread(a[1]); master<-fread(a[2]); folder<-normalizePath(a[3],mustWork=TRUE); out<-a[4]
if(dir.exists(out))stop("NO_OVERWRITE_EXISTING_OUTPUT")
key<-c("locus","dataset_key","gene_base")
if(nrow(index)!=646L||nrow(master)!=646L||anyDuplicated(index,by=key)||anyDuplicated(master,by=key))stop("ORIGINAL_646_INVALID")
m<-merge(index,master,by=key,sort=FALSE,suffixes=c("_idx","_orig"))
if(nrow(m)!=646L||any(m$status_idx!="READY")||any(m$status_orig!="PASS"))stop("INDEX_MASTER_NOT_ALIGNED")
setorder(m,locus,dataset_key,gene_base)
sdy_est<-getFromNamespace("sdY.est","coloc")
outrows<-vector("list",646L)
for(i in seq_len(646L)){
  z<-m[i]
  name<-paste0(z$locus,"__",z$dataset_key,"__",z$gene_base,".tsv")
  if(basename(z$file)!=name)stop("UNEXPECTED_INDEX_FILE_NAMING")
  f<-file.path(folder,name)
  row<-data.table(locus=z$locus,dataset_key=z$dataset_key,gene_base=z$gene_base,
      source_filename=name,status="MISSING_INPUT",reason="LOCAL_INPUT_NOT_PRESENT",
      n_snps=NA_integer_,qtl_n_first=NA_real_,qtl_n_min=NA_real_,
      qtl_n_fraction_below_first=NA_real_,qtl_an_equals_2n_all=NA,
      sdY_scalar_first=NA_real_,sdY_variant_n=NA_real_,
      sdY_relative_change=NA_real_,
      H3_original=NA_real_,H4_original=NA_real_,
      H3_variantN_sdY=NA_real_,H4_variantN_sdY=NA_real_,
      H4_delta=NA_real_,H4_cond_original=NA_real_,H4_cond_variantN=NA_real_,
      H4_ge_0_5_stable=NA,H4_ge_0_75_stable=NA)
  if(file.exists(f)){
    tryCatch({
      x<-fread(f)
      columns<-c("locus","dataset_key","gene_base","match_key","eqtl_n","eqtl_an",
        "eqtl_beta","eqtl_se","eqtl_maf","gwas_beta","gwas_se","gwas_maf","qtl_n_scalar")
      if(!all(columns%in%names(x))||nrow(x)!=z$nsnps_idx||uniqueN(x$match_key)!=nrow(x))stop("SOURCE_SNP_SCHEMA")
      if(any(x$locus!=z$locus|x$dataset_key!=z$dataset_key|x$gene_base!=z$gene_base))stop("SOURCE_KEY_DISAGREE")
      n<-as.numeric(x$eqtl_n);an<-as.numeric(x$eqtl_an)
      if(anyNA(n)||anyNA(an)||any(n<=0)||any(an<=0))stop("INVALID_QTL_N")
      if(abs(n[1]-z$qtl_n_first)>1e-8||
         abs(n[1]-z$qtl_n_max)>1e-8||
         abs(min(n)-z$qtl_n_min)>1e-8||
         uniqueN(x$qtl_n_scalar)!=1L)stop("ORIGINAL_FIRST_MAX_N_NOT_MATCHING")
      v<-as.numeric(x$eqtl_se)^2
      maf<-as.numeric(x$eqtl_maf)
      if(any(!is.finite(v)|v<=0)||any(!is.finite(maf)|maf<=0|maf>=.5+1e-9))stop("EQTL_SE_MAF_INVALID")
      y0<-suppressWarnings(sdy_est(v,maf,n[1]))
      yv<-suppressWarnings(sdy_est(v,maf,n))
      if(any(!is.finite(c(y0,yv)))||any(c(y0,yv)<=0))stop("ESTIMATED_SDY_INVALID")
      g<-list(beta=as.numeric(x$gwas_beta),varbeta=as.numeric(x$gwas_se)^2,
           snp=as.character(x$match_key),MAF=as.numeric(x$gwas_maf),N=174686,
           s=22664/174686,type="cc")
      e<-list(beta=as.numeric(x$eqtl_beta),varbeta=v,
           snp=as.character(x$match_key),type="quant",sdY=y0)
      pp0<-suppressMessages(suppressWarnings(coloc.abf(g,e,p1=1e-4,p2=1e-4,p12=1e-5)))$summary
      e$sdY<-yv
      ppv<-suppressMessages(suppressWarnings(coloc.abf(g,e,p1=1e-4,p2=1e-4,p12=1e-5)))$summary
      p0<-as.numeric(pp0[c("PP.H3.abf","PP.H4.abf")])
      pv<-as.numeric(ppv[c("PP.H3.abf","PP.H4.abf")])
      if(any(!is.finite(c(p0,pv)))||abs(p0[1]-z[["PP.H3"]])>1e-8||
         abs(p0[2]-z[["PP.H4"]])>1e-8)stop("BASELINE_SDY_EXPLICIT_REPRODUCTION_FAILURE")
      row[,':='(status="PASS",reason="",n_snps=nrow(x),
        qtl_n_first=n[1],qtl_n_min=min(n),
        qtl_n_fraction_below_first=mean(n<n[1]),
        qtl_an_equals_2n_all=all(abs(an-2*n)<1e-8),
        sdY_scalar_first=y0,sdY_variant_n=yv,sdY_relative_change=yv/y0-1,
        H3_original=p0[1],H4_original=p0[2],
        H3_variantN_sdY=pv[1],H4_variantN_sdY=pv[2],
        H4_delta=pv[2]-p0[2],
        H4_cond_original=p0[2]/sum(p0),
        H4_cond_variantN=pv[2]/sum(pv),
        H4_ge_0_5_stable=(p0[2]>=0.5)==(pv[2]>=0.5),
        H4_ge_0_75_stable=(p0[2]>=0.75)==(pv[2]>=0.75))]
    },error=function(err){
       row[,':='(status="FAIL",reason=gsub("[\r\n\t]"," ",conditionMessage(err)))]
    })
  }
  outrows[[i]]<-row
  if(i%%100==0L||i==646L)
    cat("QTL_SDY_PROGRESS",i,"/646 PASS",sum(vapply(outrows[seq_len(i)],function(v)v$status=="PASS",logical(1))),"\n")
}
res<-rbindlist(outrows,fill=TRUE)
pass<-res[status=="PASS"]
if(nrow(res)!=646L||anyDuplicated(res,by=key)||sum(res$status%in%c("PASS","FAIL","MISSING_INPUT"))!=646L)stop("FINAL_ACCOUNTING_FAILED")
summary<-list(
schema="IS_PHASE10_11_VARIANT_N_SDY_CALIBRATION_V1",
original_assays=646L,pass=nrow(pass),fail=sum(res$status=="FAIL"),
missing=sum(res$status=="MISSING_INPUT"),
total_gene_snp_records=sum(pass$n_snps),
qtl_AN_equals_twice_N_assays=sum(pass$qtl_an_equals_2n_all),
median_absolute_delta_H4=if(nrow(pass))median(abs(pass$H4_delta))else NA_real_,
max_absolute_delta_H4=if(nrow(pass))max(abs(pass$H4_delta))else NA_real_,
max_absolute_relative_change_sdY=if(nrow(pass))max(abs(pass$sdY_relative_change))else NA_real_,
original_H4_ge_0_5=sum(pass$H4_original>=.5),
variantN_sdY_H4_ge_0_5=sum(pass$H4_variantN_sdY>=.5),
original_H4_ge_0_75=sum(pass$H4_original>=.75),
variantN_sdY_H4_ge_0_75=sum(pass$H4_variantN_sdY>=.75),
threshold_0_5_crossings=sum(!pass$H4_ge_0_5_stable),
threshold_0_75_crossings=sum(!pass$H4_ge_0_75_stable),
interpretation="Variant-specific n only in estimated common trait sdY; same beta,SE and ABF effect prior; not a per-SNP n-adjusted SE likelihood",
causal_evidence="NOT_ESTABLISHED")
dir.create(out,recursive=TRUE)
fwrite(res,file.path(out,"IS_VARIANT_N_SDY_CALIBRATION_STATUS_646.tsv"),sep="\t",na="NA")
jsonlite::write_json(summary,file.path(out,"IS_VARIANT_N_SDY_CALIBRATION_SUMMARY.json"),auto_unbox=TRUE,pretty=TRUE,digits=16)
cat("QTL_VARIANT_N_SDY_STATUS PASS",summary$pass,"FAIL",summary$fail,"MISSING",summary$missing,"\n")
print(summary)
