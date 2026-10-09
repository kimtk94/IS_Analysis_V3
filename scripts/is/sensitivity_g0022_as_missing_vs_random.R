#!/usr/bin/env Rscript
# Is improvement due to choosing shared SNPs vs merely reducing variant count?
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
shared<-file.path(root,"matched_as_ais_ld_sensitivity")
ninfo<-read.table(file.path(root,"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv"),
 header=TRUE,sep="\t")
neff<-ninfo$working_n_rounded[ninfo$trait=="AS"]
dirs<-list.dirs(shared,recursive=FALSE)
if(length(dirs)!=2L)stop("Require both AS windows")
set.seed(20261009L)
out<-list()
for(dir in dirs){
  name<-basename(dir)
  allfolder<-file.path(root,name)
  originals<-read.table(file.path(allfolder,"variants.tsv"),header=TRUE,
                        sep="\t",stringsAsFactors=FALSE)
  z<-read.table(file.path(allfolder,"z.tsv"),header=TRUE)$z
  R<-as.matrix(read.table(gzfile(file.path(allfolder,"signed_ld.tsv.gz")),sep="\t"))
  sharedids<-read.table(file.path(dir,"variants.tsv"),header=TRUE,sep="\t")$variant_id
  keep<-which(originals$variant_id %in% sharedids)
  if(length(keep)!=length(sharedids)||length(keep)>nrow(R))stop("Invalid shared ID mapping")
  stats<-function(sel)estimate_s_rss(z[sel],R[sel,sel,drop=FALSE],
                                      n=neff,method="null-mle")
  full<-stats(seq_along(z))
  chosen<-stats(keep)
  trials<-numeric(20)
  for(i in seq_along(trials)){
    random<-sort(sample(seq_along(z),length(keep),replace=FALSE))
    trials[i]<-stats(random)
  }
  out[[length(out)+1]]<-data.frame(
    window=name,full_n=length(z),shared_n=length(keep),
    full_s=full,matched_study_shared_s=chosen,
    random_dropout_replicates=length(trials),
    random_dropout_s_min=min(trials),random_dropout_s_median=median(trials),
    random_dropout_s_max=max(trials),
    random_subset_s_le_matched=sum(trials<=chosen),
    caveat="EXPLORATORY_SELECTION_SENSITIVITY_NOT_CAUSAL",
    stringsAsFactors=FALSE)
  cat(name,"full_s",full,"shared_s",chosen,
      "random_min",min(trials),"random_median",median(trials),"\n")
}
res<-do.call(rbind,out)
write.table(res,file.path(shared,"G0022_AS_MATCHED_SNP_VS_RANDOM_DROPOUT.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE)
