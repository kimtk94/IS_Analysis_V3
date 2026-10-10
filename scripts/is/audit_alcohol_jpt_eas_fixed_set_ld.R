#!/usr/bin/env Rscript
# Exact same SNP set on original 1KG EAS504 vs true JPT104 LD.
# EAS104 random controls already used all original 637/574 SNPs.
# No per-panel AF filtering: fixes important confounding of prior
# panel-specific source-AF-filter comparison.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
readR<-function(path,p){
 if(file.info(path)$size!=8*p*p)stop("Mismatch binary signed LD size")
 R<-matrix(readBin(path,"double",n=p*p,size=8,endian="little"),
    nrow=p,byrow=TRUE)
 if(any(!is.finite(R))||max(abs(diag(R)-1))>1e-6||
    max(abs(R-t(R)))>1e-6)stop("Broken genotype signed LD matrix")
 R
}
out<-list()
for(locus in c("ADH1B","ALDH2")){
 d<-file.path(root,locus)
 q<-fromJSON(file.path(d,"input_qc.json"))
 jq<-fromJSON(file.path(d,"JPT104_sensitivity","input_qc.json"))
 x<-read.delim(file.path(d,"variants.tsv"),check.names=FALSE)
 y<-read.delim(file.path(d,"JPT104_sensitivity","variants.tsv"),check.names=FALSE)
 if(q$reference_EAS_n!=504||jq$JPT_n!=104||nrow(x)!=nrow(y)||
    !identical(as.character(x$ID),as.character(y$ID))||
    jq$max_reconstructed_existing_EAS_signed_LD_absolute_difference>1e-9)
     stop("Source allele-fixed same-SNP set verification failed")
 n<-median(x$source_n);z<-x$source_z;p<-nrow(x)
 r_e<-readR(file.path(d,"ld.f64.rowmajor"),p)
 r_j<-readR(file.path(d,"JPT104_sensitivity","ld.f64.rowmajor"),p)
 for(pop in c("EAS504","JPT104")){
   refR<-if(pop=="EAS504")r_e else r_j
   s<-as.numeric(susieR::estimate_s_rss(z,refR,n=n,method="null-mle"))
   if(!is.finite(s)||s<0||s>1)stop("Invalid s returned")
   out[[length(out)+1]]<-data.frame(
     locus=locus,reference=pop,
     source_snp_count=p,used_snp_count=p,
     source_genomic_variants_identical_across_references=TRUE,
     allele_frequency_screen_after_SNP_set_fixed=FALSE,
     original_GWAS_median_n=n,reference_n=if(pop=="EAS504")504 else 104,
     LD_mismatch_s=s,
     no_causal_finemap_or_MR_claim=TRUE)
   cat("FIXED_SET",locus,pop,"SNPs",p,"s",format(s,digits=9),"\n")
 }
}
all<-do.call(rbind,out)
write.table(all,file.path(root,"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH.tsv"),
   sep="\t",row.names=FALSE,quote=FALSE)
write_json(list(status="FIXED_SET_IDENTICAL_SNP_LD_REFERENCE_COMPARISON_COMPLETE",
 original_EAS_ref_n=504,JPT_ref_n=104,
 SNP_selection_identical_between_panels=TRUE,
 models=nrow(all),
 cohort_matched_reference_LD_established=FALSE,
 independent_AIS_mediation_or_CSK_causality_established=FALSE),
 file.path(root,"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH_SUMMARY.json"),
 pretty=TRUE,auto_unbox=TRUE)
