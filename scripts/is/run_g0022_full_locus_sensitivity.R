#!/usr/bin/env Rscript
# Full 4-Mb G0022 model sensitivity. Experimental only, not validated fine mapping.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1"
q<-jsonlite::fromJSON(file.path(root,"FULL_LOCUS_INPUT_QC.json"))
v<-read.table(file.path(root,"variants.tsv"),header=TRUE,sep="\t",stringsAsFactors=FALSE)
p<-nrow(v)
if(p!=q$genotype_snp_count||p>6500L||q$genotype_n!=504)stop("Input QC failed")
ldpath<-file.path(root,"ld.rowmajor.f64")
if(file.info(ldpath)$size!=p*p*8)stop("Wrong LD binary size")
R<-matrix(readBin(ldpath,what="double",n=p*p,size=8,endian="little"),
          nrow=p,byrow=TRUE)
if(max(abs(diag(R)-1))>1e-5 || max(abs(R-t(R))>1e-5))stop("Invalid LD")
z<-as.numeric(v$z)
if(length(z)!=p||any(!is.finite(z)))stop("Invalid z")
configs<-list(
 list(name="LOWER_N20000_L8",n=20000,L=8,shrink=0),
 list(name="APPROX_NEFF_L3",n=70474,L=3,shrink=0),
 list(name="APPROX_NEFF_L8_LD_SHRINK1PCT",n=70474,L=8,shrink=0.01)
)
rows<-list()
for(cfg in configs){
  workR<-R
  if(cfg$shrink>0){
    workR<-R*(1-cfg$shrink)
    diag(workR)<-1
  }
  fit<-NULL;error<-""
  start<-proc.time()[["elapsed"]]
  fit<-tryCatch(susie_rss(z=z,R=workR,n=cfg$n,L=cfg$L,
     max_iter=120,estimate_residual_variance=FALSE,
     coverage=.95,check_prior=TRUE,verbose=FALSE),
     error=function(e){error<<-conditionMessage(e);NULL})
  elapsed<-proc.time()[["elapsed"]]-start
  converged<-FALSE;cs_n<-NA_integer_;cs_sizes<-"";maxpip<-NA_real_;lead<-""
  if(!is.null(fit)){
    converged<-isTRUE(fit$converged)
    cs<-susie_get_cs(fit,Xcorr=workR,coverage=.95,min_abs_corr=.5)
    cs_n<-if(is.null(cs$cs))0L else length(cs$cs)
    cs_sizes<-if(cs_n)paste(vapply(cs$cs,length,integer(1)),collapse=";") else ""
    pi<-susie_get_pip(fit)
    m<-which.max(pi)
    maxpip<-pi[m];lead<-v$variant_id[m]
    saveRDS(fit,file.path(root,paste0("G0022_FULL_AIS_",cfg$name,"_PILOT.rds")))
  }
  rows[[length(rows)+1L]]<-data.frame(
    config=cfg$name,n_assumed=cfg$n,L=cfg$L,
    ad_hoc_ld_diagonal_shrink=cfg$shrink,
    converged=converged,
    purity_filtered_cs=cs_n,cs_sizes=cs_sizes,
    max_pip=maxpip,max_pip_variant=lead,
    seconds=round(elapsed,2),error=error,
    validated_causal_signal=FALSE,
    caveat="N_EFF_APPROX_REFERENCE_N504_AD_HOC_LD_SHRINK_NOT_FORMAL_FINITE_LD_CORRECTION",
    stringsAsFactors=FALSE)
  cat(cfg$name,"converged",converged,"cs",cs_n,"maxPIP",maxpip,"seconds",elapsed,"\n")
  flush.console()
}
out<-do.call(rbind,rows)
write.table(out,file.path(root,"G0022_FULL_AIS_SENSITIVITY.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE)
cat("Completed full-locus sensitivity",nrow(out),"configurations\n")
