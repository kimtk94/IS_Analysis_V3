#!/usr/bin/env Rscript
# Real source-aligned Japanese alcohol GWAS 2-locus (ALDH2 and ADH1B)
# SNP-level z and EAS504 signed LD. Diagnostic before exploratory SuSiE.
# Note the per-variant N differs; reference is external and population differs.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>=1)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
traits<-if(length(args)>=2&&args[[2]]!="all")args[[2]] else c("ADH1B","ALDH2")
for(trait in traits){
  sub<-file.path(root,trait)
  qc<-fromJSON(file.path(sub,"input_qc.json"))
  if(qc$status!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS"||
     qc$reference_EAS_n!=504||!isTRUE(qc$source_beta_SE_per_SNP_verified))
     stop(paste("Source provenance not complete:",trait))
  v<-read.delim(file.path(sub,"variants.tsv"),check.names=FALSE)
  p<-nrow(v)
  if(p!=qc$variants||p>1600||p<50)stop("Invalid input dimension")
  if(file.info(file.path(sub,"ld.f64.rowmajor"))$size!=8*p*p)
      stop("LD double binary size mismatch")
  memory<-as.numeric(strsplit(trimws(gsub("MemAvailable:|kB","",
    grep("^MemAvailable:",readLines("/proc/meminfo"),value=TRUE)))," +")[[1]][1])
  if(!is.finite(memory)||memory < 2*1024*1024)stop("Insufficient MemAvailable; no dense SuSiE")
  ld<-readBin(file.path(sub,"ld.f64.rowmajor"),"double",n=p*p,size=8,
       endian="little")
  R<-matrix(ld,nrow=p,ncol=p,byrow=TRUE);rm(ld)
  if(any(!is.finite(R))||max(abs(R-t(R)))>1e-7||
     max(abs(diag(R)-1))>1e-6)stop("Signed dosage LD invalid")
  z<-as.numeric(v$source_z)
  n<-as.numeric(stats::median(v$source_n))
  if(any(!is.finite(z))||n<100000||max(v$source_n)-min(v$source_n)>60000)
    stop("Real-source z/per-variant N quality gate invalid")
  cat("DIAG_START",trait,"p",p,"n_median",n,"max|z|",max(abs(z)),"\n")
  flush.console()
  s<-tryCatch(susieR::estimate_s_rss(z,R,n=n,method="null-mle"),
    error=function(e)structure(NA_real_,error=conditionMessage(e)))
  mismatch<-as.numeric(s)
  status<-if(!is.finite(mismatch))"BLOCKED_DIAGNOSTIC_ERROR" else
    if(mismatch > .10)"BLOCKED_SEVERE_GWAS_LD_MISMATCH" else
    if(mismatch > .05)"BLOCKED_MODERATE_GWAS_LD_MISMATCH" else
    "EXPLORATORY_MODEL_RUN_ALLOWED_NOT_CAUSAL"
  result<-list(locus=trait,phenotype="Japanese_alcohol_log2_grams_day_plus1",
    reference="1000G_Phase3_EAS504_GRCh37_external",
    GWAS_population="Japanese_Koyanagi_2024",reference_n=504,
    n_ref_lt_GWAS=TRUE,
    snps=p,median_sample_n=n,min_sample_n=min(v$source_n),
    max_sample_n=max(v$source_n),
    maximum_absolute_GWAS_z=max(abs(z)),
    maximum_EAF_mismatch=qc$max_reference_vs_Japanese_EAF_gap,
    LD_mismatch_s=mismatch,LD_gate=status,
    GWAS_LD_representativeness_proven=FALSE,
    per_variant_n_used_by_susie=FALSE,
    n_input_is_median_not_true_joint_GWAS_n=TRUE,
    converged=FALSE,n_credible_sets=0L,
    top_pip=NA_real_,top_snp="",
    cofounding_pleiotropy_and_sample_overlap_excluded=FALSE,
    final_finemapping_validated=FALSE,
    independent_causal_signals_identified=FALSE,
    causal_alcohol_to_BP_AIS_mediation_estimated=FALSE)
  if(status=="EXPLORATORY_MODEL_RUN_ALLOWED_NOT_CAUSAL"){
    fit<-tryCatch(susieR::susie_rss(z=z,R=R,n=n,L=6,max_iter=150,
       check_prior=TRUE,estimate_residual_variance=FALSE,coverage=.95,
       verbose=FALSE),error=function(e){result$susie_error<<-conditionMessage(e);NULL})
    if(!is.null(fit)){
       cs<-susieR::susie_get_cs(fit,Xcorr=R,coverage=.95,min_abs_corr=.5)
       pip<-susieR::susie_get_pip(fit)
       which_top<-which.max(pip)
       result$converged<-isTRUE(fit$converged)
       result$n_credible_sets<-if(is.null(cs$cs))0L else length(cs$cs)
       result$top_pip<-unname(pip[which_top])
       result$top_snp<-as.character(v$ID[which_top])
       result$credible_set_sizes<-if(is.null(cs$cs))list() else
           as.list(vapply(cs$cs,length,integer(1)))
       result$credible_set_SNPs<-if(is.null(cs$cs))list() else
           lapply(cs$cs,function(i)as.character(v$ID[i]))
       write.table(data.frame(variant_id=v$ID,pos=v$pos,z=z,pip=pip),
           file.path(sub,"EXPLORATORY_SUSIE_PIP.tsv"),quote=FALSE,
           row.names=FALSE,sep="\t")
       saveRDS(fit,file.path(sub,"EXPLORATORY_SUSIE.rds"))
    }
  }
  write_json(result,file.path(sub,"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json"),
    pretty=TRUE,auto_unbox=TRUE,digits=11,null="null")
  cat("DIAG_END",trait,"mismatch_s",mismatch,"status",status,
      "CS",result$n_credible_sets,"converged",result$converged,"\n")
  flush.console()
}
