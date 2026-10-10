#!/usr/bin/env Rscript
# Sensitivity of the archived original coloc v5.2.3 first-SNP QTL N assumption.
# Scalar-N stress test ONLY, not correct per-variant sample size model.
suppressPackageStartupMessages({library(data.table);library(coloc)})
stopifnot(as.character(packageVersion("coloc"))=="5.2.3")
argv<-commandArgs(trailingOnly=TRUE)
if(length(argv)!=4L)stop("usage: SCRIPT original_index original_master existing_snp_dir new_output_dir")
ip<-normalizePath(argv[1],mustWork=TRUE)
mp<-normalizePath(argv[2],mustWork=TRUE)
dir_in<-normalizePath(argv[3],mustWork=TRUE)
out<-argv[4]
if(dir.exists(out))stop("FAIL_CLOSED_EXISTING_OUTPUT")
idx<-fread(ip);orig<-fread(mp)
key<-c("locus","dataset_key","gene_base")
if(nrow(idx)!=646L||nrow(orig)!=646L||
   anyDuplicated(idx,by=key)||anyDuplicated(orig,by=key))stop("ORIGINAL_INDEX_SCHEMA_OR_KEYS_INVALID")
joined<-merge(idx,orig,by=key,sort=FALSE,suffixes=c("_idx","_orig"))
if(nrow(joined)!=646L||any(joined$status_idx!="READY")||
   any(joined$status_orig!="PASS")||
   any(joined$nsnps_idx!=joined$nsnps_orig))stop("INDEX_MASTER_ALIGNMENT_INVALID")
setorder(joined,locus,dataset_key,gene_base)
n_gwas<-174686;s_gwas<-22664/n_gwas
priors<-list(p1=1e-4,p2=1e-4,p12=1e-5)
rows<-vector("list",646)
for(i in seq_len(nrow(joined))){
  s<-joined[i]
  name<-paste0(s$locus,"__",s$dataset_key,"__",s$gene_base,".tsv")
  if(!identical(basename(s$file),name))stop("INDEX_INPUT_FILE_NAMING_INVALID")
  file<-file.path(dir_in,name)
  row<-data.table(locus=s$locus,dataset_key=s$dataset_key,gene_base=s$gene_base,
                  filename=name,status="MISSING_INPUT",failure_reason="SOURCE_TSV_NOT_LOCALLY_CACHED",
                  n_snps=NA_integer_,n_qtl_min=NA_real_,n_qtl_median=NA_real_,
                  n_qtl_first=NA_real_,n_qtl_max=NA_real_,n_qtl_iqr=NA_real_,
                  fraction_snp_N_below_first=NA_real_,
                  H4_FIRST=NA_real_,H4_MEDIAN=NA_real_,H4_MIN=NA_real_,
                  H3_FIRST=NA_real_,H3_MEDIAN=NA_real_,H3_MIN=NA_real_,
                  H4_cond_FIRST=NA_real_,H4_cond_MEDIAN=NA_real_,H4_cond_MIN=NA_real_,
                  abs_delta_H4_min=NA_real_,abs_delta_H4_median=NA_real_,
                  robust_H4_ge_0_5=NA,robust_H4_ge_0_75=NA)
  if(file.exists(file)){
    tryCatch({
      x<-fread(file)
      req<-c(key,"match_key","eqtl_n","qtl_n_scalar","eqtl_beta","eqtl_se","eqtl_maf","gwas_beta","gwas_se","gwas_maf")
      if(!all(req%in%names(x))||nrow(x)!=s$nsnps_idx||
         uniqueN(x$match_key)!=nrow(x))stop("INVALID_SNP_INPUT")
      if(any(x$locus!=s$locus)||any(x$dataset_key!=s$dataset_key)||
         any(x$gene_base!=s$gene_base))stop("INPUT_COMPOSITE_KEY_MISMATCH")
      n<-as.numeric(x$eqtl_n)
      if(any(!is.finite(n))||any(n<=0))stop("QTL_N_INVALID")
      if(abs(n[1]-s$qtl_n_first)>1e-8||
         abs(max(n)-s$qtl_n_max)>1e-8||
         abs(min(n)-s$qtl_n_min)>1e-8||
         abs(s$qtl_n_first-s$qtl_n_max)>1e-8||
         uniqueN(x$qtl_n_scalar)!=1L||
         abs(x$qtl_n_scalar[1]-n[1])>1e-8)
         stop("ARCHIVED_QTL_N_RANGE_OR_FIRST_MAX_NOT_REPRODUCED")
      a<-c(First=n[1],Median=median(n),Min=min(n))
      d1<-list(beta=as.numeric(x$gwas_beta),varbeta=as.numeric(x$gwas_se)^2,
               snp=as.character(x$match_key),MAF=as.numeric(x$gwas_maf),
               N=n_gwas,s=s_gwas,type="cc")
      h4<-setNames(numeric(3),names(a))
      h3<-h4
      for(j in seq_along(a)){
        d2<-list(beta=as.numeric(x$eqtl_beta),varbeta=as.numeric(x$eqtl_se)^2,
                 snp=as.character(x$match_key),MAF=as.numeric(x$eqtl_maf),
                 N=as.numeric(a[j]),type="quant")
        fit<-suppressMessages(suppressWarnings(coloc.abf(d1,d2,p1=priors$p1,p2=priors$p2,p12=priors$p12)))
        h4[j]<-as.numeric(fit$summary["PP.H4.abf"])
        h3[j]<-as.numeric(fit$summary["PP.H3.abf"])
      }
      if(any(!is.finite(h4))||any(!is.finite(h3)))stop("NONFINITE_RESULT")
      if(abs(h4["First"]-s[["PP.H4"]])>1e-8||
         abs(h3["First"]-s[["PP.H3"]])>1e-8)stop("BASELINE_NOT_REPRODUCED")
      conditional<-h4/(h4+h3)
      if(any(!is.finite(conditional)))stop("H4_CONDITIONAL_NOT_DEFINED")
      row[,':='(n_snps=nrow(x),n_qtl_min=min(n),n_qtl_median=median(n),
                n_qtl_first=n[1],n_qtl_max=max(n),n_qtl_iqr=IQR(n),
                fraction_snp_N_below_first=mean(n<n[1]),
                H4_FIRST=h4["First"],H4_MEDIAN=h4["Median"],H4_MIN=h4["Min"],
                H3_FIRST=h3["First"],H3_MEDIAN=h3["Median"],H3_MIN=h3["Min"],
                H4_cond_FIRST=conditional["First"],H4_cond_MEDIAN=conditional["Median"],H4_cond_MIN=conditional["Min"],
                abs_delta_H4_min=abs(h4["Min"]-h4["First"]),
                abs_delta_H4_median=abs(h4["Median"]-h4["First"]),
                robust_H4_ge_0_5=all(h4>=0.5),robust_H4_ge_0_75=all(h4>=0.75),
                status="PASS",failure_reason="")]
    },error=function(e){
      row[,':='(status="FAIL",failure_reason=gsub("[\\r\\n\\t]"," ",conditionMessage(e)))]
    })
  }
  rows[[i]]<-row
  if(i%%100==0||i==646){
    cat("N_SCALAR_PROGRESS",i,"/646 PASS",sum(vapply(rows[seq_len(i)],function(z)z$status=="PASS",logical(1))),"\n")
  }
}
res<-rbindlist(rows,fill=TRUE)
if(nrow(res)!=646||anyDuplicated(res,by=key)||
   sum(res$status%in%c("PASS","FAIL","MISSING_INPUT"))!=646)stop("STATUS_PARTITION_FAILED")
pass<-res[status=="PASS"]
summary<-data.table(
   metric=c("total_original_test_count","numeric_sensitivity_PASS","numeric_sensitivity_FAIL",
            "inputs_not_cached","n_first_equals_max","n_median_equals_first",
            "n_baseline_H4_ge_0_5","n_min_H4_ge_0_5",
            "n_baseline_H4_ge_0_75","n_min_H4_ge_0_75",
            "n_threshold_0_5_changed","n_threshold_0_75_changed",
            "max_abs_delta_H4_min","median_abs_delta_H4_min"),
   value=as.character(c(nrow(res),nrow(pass),sum(res$status=="FAIL"),
            sum(res$status=="MISSING_INPUT"),
            sum(pass$n_qtl_first==pass$n_qtl_max),
            sum(pass$n_qtl_median==pass$n_qtl_first),
            sum(pass$H4_FIRST>=.5),sum(pass$H4_MIN>=.5),
            sum(pass$H4_FIRST>=.75),sum(pass$H4_MIN>=.75),
            sum((pass$H4_FIRST>=.5)!=(pass$H4_MIN>=.5)),
            sum((pass$H4_FIRST>=.75)!=(pass$H4_MIN>=.75)),
            if(nrow(pass))max(pass$abs_delta_H4_min)else NA_real_,
            if(nrow(pass))median(pass$abs_delta_H4_min)else NA_real_))
)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(res,file.path(out,"IS_646_QTL_N_SCALAR_SENSITIVITY.tsv"),sep="\t",na="NA")
fwrite(summary,file.path(out,"IS_646_QTL_N_SCALAR_SUMMARY.tsv"),sep="\t")
writeLines(c(
 "# IS original-646 QTL scalar-N stress test",
 "",
 "- Uses unmodified archived SNP beta/SE/MAF and coloc v5.2.3 priors, and compares first/max N versus median N versus minimum SNP N.",
 "- These are **three scalar N approximations**. They do not implement per-SNP effective N correction, estimate phenotype variance from donor genotypes, or provide biological replication.",
 "- All original 646 index rows are reported with PASS/FAIL/MISSING_INPUT and source file key; only actual cached input TSVs can be PASS.",
 "- Scalar N sensitivity can be weak because coloc.abf quantitative-trait sdY approximation uses both per-SNP SE and MAF; no change does not establish correct n_eff.",
 "- Variant effect-allele alignment and GTEx cohort-matched LD remain separate unresolved issues.",
 paste0("- PASS=",nrow(pass),", FAIL=",sum(res$status=="FAIL"),
        ", MISSING_INPUT=",sum(res$status=="MISSING_INPUT")," (local snapshot only).")
),file.path(out,"IS_646_QTL_N_SCALAR_INTERPRETATION.md"))
cat("IS_QTL_N_SCALAR_SENSITIVITY",if(sum(res$status=="FAIL")==0)"NO_COMPUTATION_FAILURE"else"FAILURE_PRESENT",
 "PASS",nrow(pass),"FAIL",sum(res$status=="FAIL"),"MISSING",sum(res$status=="MISSING_INPUT"),"\n")
