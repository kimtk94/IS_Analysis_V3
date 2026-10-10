#!/usr/bin/env Rscript
# Compare Japanese Koyanagi original alcohol GWAS summary stats to two
# *different* external 1000G LD panels, same variant universe and GWAS z.
# EAS504 and JPT104 s diagnostics. Do NOT infer LD fix from a smaller s
# or causal gene/phenotypic mediation. No SuSiE refit on small JPT.
suppressPackageStartupMessages({library(susieR); library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0) args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
out<-list()
for(locus in c("ADH1B","ALDH2")){
 f<-file.path(root,locus)
 qc<-fromJSON(file.path(f,"input_qc.json"))
 jqc<-fromJSON(file.path(f,"JPT104_sensitivity","input_qc.json"))
 if(qc$status!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS" ||
   jqc$status!="JPT104_LD_SOURCE_RECONSTRUCTED_EAS_BASELINE_CROSSCHECK_PASS" ||
   jqc$max_reconstructed_existing_EAS_signed_LD_absolute_difference>1e-9||
   jqc$JPT_n!=104||jqc$source_original_EAS_n!=504)
   stop("Full genotype reference provenance not verified")
 e<-read.delim(file.path(f,"variants.tsv"),check.names=FALSE)
 j<-read.delim(file.path(f,"JPT104_sensitivity","variants.tsv"),check.names=FALSE)
 if(nrow(j)!=nrow(e)||!identical(as.character(e$ID),as.character(j$ID)))
   stop("JPT/EAS variants and positions do not match")
 p<-nrow(e)
 if(p<50||p>1600)stop("Unsafe LD matrix dimension")
 getR<-function(path){
   if(file.info(path)$size!=8*p*p)stop("LD binary size mismatch")
   x<-readBin(path,what="double",n=p*p,size=8,endian="little")
   R<-matrix(x,nrow=p,byrow=TRUE)
   if(any(!is.finite(R))||max(abs(diag(R)-1))>1e-7||
       max(abs(R-t(R)))>1e-7)stop("Source LD not signed symmetric correlation")
   R
 }
 E<-getR(file.path(f,"ld.f64.rowmajor"))
 J<-getR(file.path(f,"JPT104_sensitivity","ld.f64.rowmajor"))
 z<-e$source_z
 n<-as.numeric(median(e$source_n))
 if(any(!is.finite(z))||n<1e5)stop("GWAS source beta/se N invalid")
 if(any(abs(e$source_z*e$source_se-e$source_beta_ALT)>1e-6))
    stop("GWAS source beta/se sign QC failed")
 afE<-e$reference_EAS_ALT_EAF
 afJ<-j$JPT_ALT_freq
 mafE<-pmin(afE,1-afE)
 mafJ<-pmin(afJ,1-afJ)
 mafG<-as.numeric(e$Japanese_original_ALT_EAF)
 results<-list()
 for(panel in c("EAS504","JPT104")){
   R<-if(panel=="EAS504") E else J
   af<-if(panel=="EAS504") afE else afJ
   maf<-if(panel=="EAS504") mafE else mafJ
   for(threshold in c(.10,.05,.03)){
     # Frequency threshold on original Japanese GWAS EAF vs this panel's
     # *own* frequency, so sample sizes may change, recorded explicitly.
     use<-which(maf>=.05 & abs(af-mafG)<=threshold)
     if(length(use)<40){s<-NA_real_;gate<-"INSUFFICIENT_SNPS"}
     else {
       s<-tryCatch(as.numeric(susieR::estimate_s_rss(z[use],
         R[use,use,drop=FALSE],n=n,method="null-mle")),
         error=function(err) NA_real_)
       gate<-if(!is.finite(s))"ESTIMATE_S_FAILED" else
         if(s>.05)"INTERNAL_QC_MISMATCH_HIGH" else
           "DIAGNOSTIC_ONLY_NOT_AN_INFERENCE_PASS"
     }
     results[[length(results)+1]]<-data.frame(
       locus=locus,panel=panel,reference_n=if(panel=="EAS504")504 else 104,
       Japanese_GWAS_N_median=n,MAF_filter=.05,
       max_abs_af_difference=threshold,SNPs_included=length(use),
       LD_mismatch_s=s,gate=gate,
       smaller_s_does_not_validate_conditional_finemapping=TRUE,
       causal_ALDH2_or_ADH1B_effect_validated=FALSE,
       clinical_alcohol_BP_AIS_mediation_estimated=FALSE)
     cat("JPT_LD_COMPARE",locus,panel,"AF_threshold",threshold,
       "SNPs",length(use),"s",format(s,digits=8),"gate",gate,"\n")
     flush.console()
   }
 }
 out[[locus]]<-do.call(rbind,results)
}
all<-do.call(rbind,out)
write.table(all,file.path(root,"ALCOHOL_ALDH2_ADH1B_EAS504_VS_JPT104_LD_MISMATCH.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)
summary<-list(
 status="JPT104_EAS504_DIAGNOSTIC_COMPARISON_ONLY",
 original_EAS_reference_n=504,Japanese_subset_JPT_n=104,
 GWAS_is_original_Koyanagi_2024_Japanese=TRUE,
 scenarios=nrow(all),all_panels_allele_harmonized=TRUE,
 population_match_improves_LD_automatically=FALSE,
 smaller_s_implies_original_genomewide_model_correct=FALSE,
 LD_reference_104_too_small_for_confident_finemapping=TRUE,
 unmeasured_individual_cohort_overlap_resolved=FALSE,
 ALDH2_causal_gene_resolved=FALSE,ADH1B_causal_variant_resolved=FALSE,
 AIS_mediated_causal_effect_calculated=FALSE)
write_json(summary,file.path(root,
 "ALCOHOL_ALDH2_ADH1B_EAS_JPT_LD_MISMATCH_SUMMARY.json"),
 pretty=TRUE,auto_unbox=TRUE)
cat("JPT_EAS_DIAGNOSTIC_COMPLETE",nrow(all),"\n")
