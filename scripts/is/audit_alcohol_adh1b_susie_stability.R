#!/usr/bin/env Rscript
# ADH1B source-based exploratory SuSiE stability: L, AF/MAF filter choice.
# No causal target inference; ALDH2 full locus remains separately blocked.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>=1)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
folder<-file.path(root,"ADH1B")
qc<-fromJSON(file.path(folder,"input_qc.json"))
v<-read.delim(file.path(folder,"variants.tsv"),check.names=FALSE)
p<-nrow(v)
if(qc$reference_EAS_n!=504||p!=637||p!=qc$variants)stop("Unexpected source ADH1B SNP set")
vec<-readBin(file.path(folder,"ld.f64.rowmajor"),what="double",n=p*p,size=8,endian="little")
R<-matrix(vec,nrow=p,byrow=TRUE);rm(vec)
n<-as.numeric(median(v$source_n))
if(n!=154559)stop("Unverified GWAS N")
scenarios<-list(
 list(name="BASELINE_MAF05_DELTA10",ma=.05,delta=.10),
 list(name="RESTRICTED_MAF10_DELTA03",ma=.10,delta=.03))
rows<-list()
for(sc in scenarios){
  maf<-pmin(v$reference_EAS_ALT_EAF,1-v$reference_EAS_ALT_EAF)
  freq_delta<-abs(v$Japanese_original_ALT_EAF-v$reference_EAS_ALT_EAF)
  use<-which(maf>=sc$ma & freq_delta<=sc$delta)
  if(length(use)<50)stop("SNP selection unexpectedly too small")
  x<-v[use,];LD<-R[use,use,drop=FALSE];zz<-x$source_z
  estimated_s<-susieR::estimate_s_rss(zz,LD,n=n)
  if(estimated_s>.05)stop("Filtered ADH1B source LD mismatch too high")
  for(L in c(3L,6L,10L)){
    fit<-tryCatch(susieR::susie_rss(z=zz,R=LD,n=n,L=L,max_iter=120,
      coverage=.95,estimate_residual_variance=FALSE,verbose=FALSE),
      error=function(e)NULL)
    if(is.null(fit)){conv<-FALSE;num<-0L;top<-"";pip<-NA_real_;css=""}
    else{
      cs<-susieR::susie_get_cs(fit,Xcorr=LD,coverage=.95,min_abs_corr=.5)
      conv<-isTRUE(fit$converged);num<-if(is.null(cs$cs))0L else length(cs$cs)
      score<-susieR::susie_get_pip(fit)
      i<-which.max(score);top<-as.character(x$ID[i]);pip<-as.numeric(score[i])
      css<-if(is.null(cs$cs))"" else
        paste(vapply(cs$cs,function(ids)length(ids),integer(1)),collapse=";")
    }
    row<-data.frame(locus="ADH1B",scenario=sc$name,
      SNP_count=length(use),L_maxcausal_assumed=L,
      MAF_cutoff=sc$ma,max_allele_freq_discordance=sc$delta,
      LD_mismatch_s=estimated_s,converged=conv,
      exploratory_95pct_credible_sets=num,
      credible_set_sizes=css,top_pip_snp=top,top_pip=pip,
      cross_model_independent_causal_signals_proven=FALSE,
      biological_target_or_mediation_validated=FALSE)
    rows[[length(rows)+1]]<-row
    cat("ADH1B_STABILITY",sc$name,"L",L,"N",length(use),"s",estimated_s,
        "converged",conv,"CS",num,"top",top,"PIP",pip,"\n")
    flush.console()
  }
}
out<-do.call(rbind,rows)
write.table(out,file.path(root,"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.tsv"),
  row.names=FALSE,quote=FALSE,sep="\t")
final<-list(status="ADH1B_EXPLORATORY_STABILITY_ONLY",
  models=nrow(out),nonconverged=sum(!out$converged),
  across_L_and_filters_CS_range=c(min(out$exploratory_95pct_credible_sets),
    max(out$exploratory_95pct_credible_sets)),
  causally_independent_signals_established=0L,
  rs671_not_present_in_ADH1B_region=TRUE,
  AIS_outcome_not_used_for_exposure_finemap=TRUE)
write_json(final,file.path(root,"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY_SUMMARY.json"),
  pretty=TRUE,auto_unbox=TRUE)
print(final)
