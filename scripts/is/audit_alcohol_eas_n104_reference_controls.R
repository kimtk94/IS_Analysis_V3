#!/usr/bin/env Rscript
# Fixed-seed EAS104 subsampling diagnostic for Japanese GWAS vs 1KG LD.
# Observational panel-size-only stress test, NOT a bootstrap sampling
# distribution of causality. 12 correlated draws, no MR or SuSiE fit.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
out<-list()
for(locus in c("ADH1B","ALDH2")){
 folder<-file.path(root,locus)
 v<-read.delim(file.path(folder,"variants.tsv"),check.names=FALSE)
 q<-fromJSON(file.path(folder,"EAS104_random_reference_controls","random_reference_manifest.json"))
 if(q$locus!=locus||q$source_1000G_EAS_n!=504||
     q$draw_n!=104||q$draws_requested!=12||
     q$reproducibility_seed!=20261010||nrow(q$draws)!=12)
   stop("Unsafe or incomplete fixed control panel manifest")
 p<-nrow(v);z<-v$source_z
 n<-as.numeric(median(v$source_n))
 readR<-function(path){
   if(file.info(path)$size!=8*p*p)stop("Wrong fixed-size LD bytes")
   x<-readBin(path,"double",n=p*p,size=8,endian="little")
   R<-matrix(x,nrow=p,ncol=p,byrow=TRUE)
   if(max(abs(diag(R)-1))>1e-6||max(abs(R-t(R)))>1e-6)stop("Invalid LD")
   R
 }
 records<-list()
 for(i in seq_len(nrow(q$draws))){
   x<-q$draws[i,]
   if(x$status!="COMPLETE"||x$n_SNPs_monomorphic>0)
     stop("Monomorphic control draw cannot be silently treated as normal")
   R<-readR(file.path(folder,"EAS104_random_reference_controls",x$file))
   s<-tryCatch(as.numeric(susieR::estimate_s_rss(z,R,n=n,method="null-mle")),
        error=function(e)NA_real_)
   records[[length(records)+1]]<-data.frame(
      locus=locus,reference="RANDOM_EAS104",
      replicate_index=i,seed=20261010,
      GWAS_effectively_constant_median_n=n,ref_n=104,SNP_n=p,
      diagnostic_s=s,
      population_is_JPT_only=FALSE,
      all_subjects_come_from_EAS=TRUE,
      sampled_replicates_not_independent=TRUE,
      MR_or_causal_finemapping_performed=FALSE)
   cat("EAS_104_CONTROL",locus,i,"s",format(s,digits=8),"\n")
   flush.console()
 }
 out[[locus]]<-do.call(rbind,records)
}
result<-do.call(rbind,out)
write.table(result,file.path(root,"ALCOHOL_EAS104_RANDOM_DOWNSAMPLING_LD_MISMATCH.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE)
other<-read.delim(file.path(root,"ALCOHOL_ALDH2_ADH1B_EAS504_VS_JPT104_LD_MISMATCH.tsv"))
brief<-lapply(c("ADH1B","ALDH2"),function(locus){
 x<-result[result$locus==locus,]
 ja<-other[other$locus==locus &
           other$panel=="JPT104" & other$max_abs_af_difference==.1,]
 e<-other[other$locus==locus &
          other$panel=="EAS504" & other$max_abs_af_difference==.1,]
 if(nrow(ja)!=1||nrow(e)!=1||any(!is.finite(x$diagnostic_s)))
     stop("JPT and EAS104 diagnostic sources incomplete")
 list(locus=locus,
   EAS504_s=e$LD_mismatch_s,JPT104_s=ja$LD_mismatch_s,
   EAS104_repeat_n=nrow(x),
   EAS104_s_median=median(x$diagnostic_s),
   EAS104_s_min=min(x$diagnostic_s),EAS104_s_max=max(x$diagnostic_s),
   JPT104_s_outside_EAS104_random_range=
      ja$LD_mismatch_s<min(x$diagnostic_s)||
      ja$LD_mismatch_s>max(x$diagnostic_s))
})
names(brief)<-c("ADH1B","ALDH2")
report<-list(
 status="REF_N_SENSITIVITY_COMPLETE_EXPLORATORY_ONLY",
 reference_504_total_n=504,downsample_n=104,repeats=12,
 fixed_seed=20261010,estimates=brief,
 independent_Japanese_population_LD_external_found=FALSE,
 large_Japanese_cohort_LD_available=FALSE,
 random_reference_replicates_are_correlated=TRUE,
 does_not_estimate_sampling_distribution_of_true_GWAS_LD=TRUE,
 LD_mismatch_mechanism_is_not_uniquely_identified=TRUE,
 causal_alcohol_pleiotropy_stroke_not_inferred=TRUE)
write_json(report,file.path(root,
 "ALCOHOL_JPT104_EAS104_PANEL_SIZE_CONTROL_SUMMARY.json"),
 pretty=TRUE,auto_unbox=TRUE,digits=9)
cat("PANEL_SIZE_DIAGNOSTIC_PASS",toJSON(brief,auto_unbox=TRUE), "\n")
