#!/usr/bin/env Rscript
# Approximate case-control N sensitivity. Not calibrated per-SNP effective N.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
manifest<-read.table(file.path(root,"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv"),
  sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(manifest)!=2||!all(c("AS","AIS")%in%manifest$trait))stop("Invalid provenance manifest")
if(any(manifest$validation!="SAMPLE_TOTAL_MATCHED_OFFICIAL_YAML"))stop("Invalid sample provenance")
files<-list.files(root,pattern="^input_qc\\.json$",recursive=TRUE,full.names=TRUE)
if(length(files)!=5L)stop("Need 5 input-QC windows for G0022")
outdir<-file.path(root,"effective_n_sensitivity")
dir.create(outdir,showWarnings=FALSE,recursive=TRUE)
results<-list()
for(f in files){
  folder<-dirname(f);name<-basename(folder)
  trait<-if(grepl("_AIS_",name,fixed=TRUE))"AIS" else
    if(grepl("_AS_",name,fixed=TRUE))"AS" else stop("Bad source label")
  neff<-manifest$working_n_rounded[manifest$trait==trait]
  v<-read.table(file.path(folder,"variants.tsv"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
  z<-read.table(file.path(folder,"z.tsv"),header=TRUE,sep="\t")$z
  R<-as.matrix(read.table(gzfile(file.path(folder,"signed_ld.tsv.gz")),sep="\t"))
  if(nrow(R)!=ncol(R)||nrow(R)!=length(z)||nrow(v)!=length(z))stop("Input mismatch")
  if(any(!is.finite(R))||max(abs(R-t(R)))>1e-5)stop("Invalid LD")
  fit<-NULL;err<-""
  fit<-tryCatch(susie_rss(z=z,R=R,n=neff,L=5,
    estimate_residual_variance=FALSE,max_iter=150,
    coverage=.95,check_prior=TRUE,verbose=FALSE),
    error=function(e){err<<-conditionMessage(e);NULL})
  converged<-FALSE;top_variant<-"";top_pip<-NA_real_
  cs_count<-NA_integer_;purity<-NA_real_;cs_sizes<-""
  if(!is.null(fit)){
    converged<-isTRUE(fit$converged)
    pip<-susie_get_pip(fit)
    if(length(pip)!=nrow(v))stop("PIP length mismatch")
    top<-which.max(pip);top_variant<-v$variant_id[top];top_pip<-pip[top]
    cs<-susie_get_cs(fit,Xcorr=R,coverage=.95,min_abs_corr=.5)
    cs_count<-if(is.null(cs$cs))0L else length(cs$cs)
    cs_sizes<-if(cs_count>0)paste(vapply(cs$cs,length,integer(1)),collapse=";") else ""
    if(!is.null(cs$purity)&&nrow(cs$purity)>0){
      purity<-min(cs$purity$min.abs.corr,na.rm=TRUE)
    }
    saveRDS(fit,file.path(outdir,paste0(name,"_neff_",neff,".rds")))
    v$pip<-pip
    write.table(v,file.path(outdir,paste0(name,"_neff_",neff,"_PIP.tsv")),
      row.names=FALSE,quote=FALSE,sep="\t")
  }
  result<-data.frame(signal=name,trait=trait,approximate_n_eff=neff,
    total_study_n=manifest$official_meta_total_n[manifest$trait==trait],
    pilot_snps=length(z),converged=converged,
    purity_filtered_cs_count=cs_count,cs_sizes=cs_sizes,
    min_abs_corr_purity=purity,top_variant=top_variant,
    top_pip=top_pip,error_message=err,
    scientific_status=if(converged)
      "APPROXIMATE_N_SENSITIVITY_PILOT_NOT_VALIDATED"
      else "NONCONVERGED_BLOCKED",stringsAsFactors=FALSE)
  results[[length(results)+1L]]<-result
  cat(name,"n_eff",neff,"converged",converged,"CS",cs_count,
    "PIP",top_pip,"\n");flush.console()
}
summary<-do.call(rbind,results)
write.table(summary,file.path(outdir,"G0022_SUSIE_EFFECTIVE_N_SENSITIVITY.tsv"),
  sep="\t",row.names=FALSE,quote=FALSE)
cat("Completed",nrow(summary),"pilot windows; converged",sum(summary$converged),"\n")
