#!/usr/bin/env Rscript
# Source-level filtration stress test of Japanese alcohol exposure GWAS vs EAS LD.
# No SuSiE unless primary diagnostic passes; selection-dependent s cannot
# repair 504-reference sample size or evidence of overlap/pleiotropy.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
scenarios<-list(
 list(tag="SOURCE_QC_DEFAULT",min_maf=.05,max_ef_diff=.10),
 list(tag="REF_MAF_10_EAFDELTA_05",min_maf=.10,max_ef_diff=.05),
 list(tag="REF_MAF_10_EAFDELTA_03",min_maf=.10,max_ef_diff=.03),
 list(tag="REF_MAF_05_EAFDELTA_03",min_maf=.05,max_ef_diff=.03))
rows<-list()
for(locus in c("ADH1B","ALDH2")){
 f<-file.path(root,locus)
 qc<-fromJSON(file.path(f,"input_qc.json"))
 v<-read.delim(file.path(f,"variants.tsv"),check.names=FALSE)
 p<-nrow(v)
 if(p!=qc$variants||p>1600||qc$reference_EAS_n!=504)
   stop("Source QC gate failed")
 block<-readBin(file.path(f,"ld.f64.rowmajor"),what="double",
    n=p*p,size=8,endian="little")
 R<-matrix(block,nrow=p,ncol=p,byrow=TRUE);rm(block)
 z<-as.numeric(v$source_z)
 original_s<-fromJSON(file.path(f,"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json"))$LD_mismatch_s
 for(sc in scenarios){
   maf<-pmin(v$reference_EAS_ALT_EAF,1-v$reference_EAS_ALT_EAF)
   freq_diff<-abs(v$reference_EAS_ALT_EAF-v$Japanese_original_ALT_EAF)
   use<-which(maf>=sc$min_maf & freq_diff<=sc$max_ef_diff)
   if(length(use)<50){sm<-NA_real_;status<-"BLOCKED_LT_50_SNPS"}
   else{
     sm<-tryCatch(susieR::estimate_s_rss(z[use],R[use,use],
       n=as.numeric(stats::median(v$source_n[use]))),
       error=function(e)NA_real_)
     status<-if(!is.finite(sm))"DIAGNOSTIC_FAIL" else if(sm>.05)
       "BLOCKED_GWAS_LD_MISMATCH" else "EXPLORATORY_ONLY_NOT_CAUSAL"
   }
   rows[[length(rows)+1]]<-data.frame(locus=locus,scenario=sc$tag,
     variant_n=length(use),MAF_min=sc$min_maf,absolute_EAF_diff_max=sc$max_ef_diff,
     LD_mismatch_s=sm,original_full_input_s=original_s,
     diagnostic_gate=status,
     causal_SNP_validated=FALSE,causal_gene_validated=FALSE,
     result_is_proof_of_cohort_matched_LD=FALSE)
   cat("SENSITIVITY",locus,sc$tag,"n",length(use),"s",sm,"status",status,"\n")
   flush.console()
 }
}
out<-do.call(rbind,rows)
write.table(out,file.path(root,"ALCOHOL_ADH1B_ALDH2_GWAS_LD_MISMATCH_FILTER_SENSITIVITY.tsv"),
  sep="\t",quote=FALSE,row.names=FALSE)
write_json(list(status="CONDITIONAL_FREQUENCY_SENSITIVITY_EXPLORATORY_ONLY",
 tested_scenarios=nrow(out),original_fit_may_be_BLOCKED=TRUE,
 no_causal_claim=TRUE),
 file.path(root,"ALCOHOL_GWAS_LD_MISMATCH_FILTER_SENSITIVITY_SUMMARY.json"),
 pretty=TRUE,auto_unbox=TRUE)
