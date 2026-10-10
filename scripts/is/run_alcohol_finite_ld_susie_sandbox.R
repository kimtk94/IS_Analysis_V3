#!/usr/bin/env Rscript
# Isolated susieR >= 0.16.6 finite-reference / empirical-Bayes LD mismatch
# sensitivity on original Japanese alcohol GWAS ALDH2 and ADH1B loci.
# This deliberately does NOT update any canonical research results or assert
# causal validity, particularly for earlier blocked ALDH2.
suppressPackageStartupMessages({library(susieR);library(jsonlite)})
vers<-packageVersion("susieR")
if(vers < "0.16.6" || !("R_finite" %in% names(formals(susieR::susie_rss))) ||
   !("R_mismatch" %in% names(formals(susieR::susie_rss))))
  stop("Isolated susieR 0.16.6+ unavailable: do not silently use old method")
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
out<-if(length(args)>1)args[[2]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1"
loc<-if(length(args)>2)args[[3]] else "ADH1B"
if(!loc%in%c("ADH1B","ALDH2"))stop("Unapproved locus")
source_dir<-file.path(root,loc)
q<-fromJSON(file.path(source_dir,"input_qc.json"))
old<-fromJSON(file.path(source_dir,"ALCOHOL_SUSIE_DIAGNOSTIC_AND_GATE.json"))
v<-read.delim(file.path(source_dir,"variants.tsv"),check.names=FALSE)
p<-nrow(v)
if(q$status!="ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS" ||
   !isTRUE(q$source_beta_SE_per_SNP_verified) ||
   !isTRUE(q$LD_input_order_exact_variant_table) ||
   q$reference_EAS_n!=504 || p!=q$variants || p<50 || p>1600 ||
   old$snps!=p)
  stop("Original data source provenance failed")
f<-file.path(source_dir,"ld.f64.rowmajor")
if(file.info(f)$size!=p*p*8)stop("Signed LD matrix binary size mismatches")
R<-matrix(readBin(f,"double",n=p*p,size=8,endian="little"),p,p,byrow=TRUE)
if(any(!is.finite(R))||max(abs(R-t(R)))>1e-8||
   max(abs(diag(R)-1))>1e-6)stop("Signed reference LD not valid")
z<-as.numeric(v$source_z)
N<-as.numeric(median(v$source_n))
if(N<100000 || any(!is.finite(z)) ||
   any(abs(v$source_z*v$source_se-v$source_beta_ALT)>1e-6))
  stop("SNP-level effect direction or GWAS summary N invalid")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(out,"models",loc),recursive=TRUE,showWarnings=FALSE)
params<-list(
 list(id="FINITE_REF_504",ref=504,mismatch="none"),
 list(id="FINITE_REF_504_EB_MISMATCH",ref=504,mismatch="eb"))
cst<-function(fit)if(is.null(fit$sets)||is.null(fit$sets$cs)) 0L else length(fit$sets$cs)
diag_summary<-function(x){
  if(is.null(x))return(list(diagnostics_provided=FALSE))
  out<-list(diagnostics_provided=TRUE,names=names(x))
  for(k in c("B","effective_rank","r_over_B","lambda_bias","B_corrected",
     "Q_art","R_sensitivity_flag","R_reliability_flag")){
    val<-x[[k]]
    if(!is.null(val)&&is.atomic(val)&&length(val)<=30)
       out[[k]]<-val
  }
  if(!is.null(x$per_variable_penalty) && is.numeric(x$per_variable_penalty)){
     a<-as.numeric(x$per_variable_penalty)
     out$penalty_median<-median(a,na.rm=TRUE)
     out$penalty_max<-max(a,na.rm=TRUE)
  }
  out
}
for(param in params){
 id<-param$id
 cat("SANDBOX_START",loc,id,"p",p,"study_n_median",N,
    "susieR",as.character(vers),"\n")
 flush.console()
 caught<-character()
 fit<-tryCatch(withCallingHandlers(
   susieR::susie_rss(z=z,R=R,n=N,L=6,max_iter=80,
     estimate_residual_variance=FALSE,coverage=.95,
     R_finite=param$ref,R_mismatch=param$mismatch,verbose=FALSE),
   warning=function(w){caught<<-c(caught,conditionMessage(w));
     invokeRestart("muffleWarning")}),
   error=function(e){caught<<-c(caught,paste("ERROR:",conditionMessage(e)));NULL})
 result<-list(
  status="EXPLORATORY_SANDBOX_NO_CAUSAL_INFERENCE",
  locus=loc,susieR_version=as.character(vers),reference_panel="1000G_EAS_504",
  source="Koyanagi_2024_Japanese_alcohol_beta_SE_original",
  original_no_correction_s=old$LD_mismatch_s,
  original_s_still_blocks_ALDH2=loc=="ALDH2",
  original_snp_count=p,source_median_n=N,
  reference_actual_B=param$ref,
  R_finite_enabled=TRUE,R_mismatch=param$mismatch,
  source_alleles_signed_and_order_QC=TRUE,
  finite_LD_reference_correction_sandbox_only=TRUE,
  output_is_calibrated_causal_gene_or_IV=FALSE,
  source_cohort_ld_matched=FALSE,
  AIS_MR_mediation_validated=FALSE,
  finished_without_error=!is.null(fit),
  converged=if(is.null(fit))FALSE else isTRUE(fit$converged),
  credible_sets_exploratory=if(is.null(fit))NA_integer_ else cst(fit),
  max_PIP=if(is.null(fit))NA_real_ else max(fit$pip),
  top_variant=if(is.null(fit))"" else v$ID[which.max(fit$pip)],
  warnings=caught,
  diagnostics=if(is.null(fit))list(diagnostics_provided=FALSE)
              else diag_summary(fit$R_finite_diagnostics))
 if(!is.null(fit)){
   path<-file.path(out,"models",loc,paste0(id,".rds"))
   saveRDS(fit,path)
   result$stored_fit_RDS<-path
   write.table(data.frame(variant_id=v$ID,pip=fit$pip),
     file.path(out,"models",loc,paste0(id,".pip.tsv")),
     row.names=FALSE,col.names=TRUE,sep="\t",quote=FALSE)
 }
 outfile<-file.path(out,"models",loc,paste0(id,".qc.json"))
 write_json(result,outfile,pretty=TRUE,auto_unbox=TRUE,digits=10,na="null")
 cat("SANDBOX_END",loc,id,"converged",result$converged,
    "CS",result$credible_sets_exploratory,"maxPIP",result$max_PIP,
    "warnings",length(caught),"diag_flags",
    paste(result$diagnostics$names,collapse=","),"\n")
 flush.console()
}
