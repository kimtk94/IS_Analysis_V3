#!/usr/bin/env Rscript
# Study z vs same 1000G EAS LD, exactly identical SNP set / same allele orientation.
suppressPackageStartupMessages(library(susieR))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/susie_g0022_pilot"
subdir<-file.path(root,"matched_as_ais_ld_sensitivity")
mf<-read.table(file.path(root,"G0022_GIGASTROKE_EAS_EFFECTIVE_N_APPROX.tsv"),
              header=TRUE,sep="\t",stringsAsFactors=FALSE)
dirs<-list.dirs(subdir,recursive=FALSE,full.names=TRUE)
if(length(dirs)!=2L)stop("Require exactly 2 AS windows")
output<-list()
for(folder in dirs) {
  variants<-readLines(file.path(folder,"variants.tsv"))[-1]
  R<-as.matrix(read.table(gzfile(file.path(folder,"ld.tsv.gz")),sep="\t"))
  if(nrow(R)!=length(variants)||nrow(R)!=ncol(R))stop("LD dimensions invalid")
  for(trait in c("AS","AIS")){
    z<-as.numeric(readLines(file.path(folder,paste0(tolower(trait),"_z.tsv"))))
    if(length(z)!=nrow(R)||any(!is.finite(z)))stop("Z dimensions invalid")
    n<-mf$working_n_rounded[mf$trait==trait]
    s<-estimate_s_rss(z,R,n=n,method="null-mle")
    output[[length(output)+1]]<-data.frame(
      window=basename(folder),trait=trait,
      matched_variants=length(variants),approximate_n_eff=n,
      estimated_s=s,
      threshold_flag=if(s>.2)"VERY_HIGH_LD_MISMATCH" else
        if(s>.05)"LD_MISMATCH" else "LOW_MISMATCH_NOT_PROOF",
      biological_interpretation="NOT_INDEPENDENT_VALIDATION",
      stringsAsFactors=FALSE)
    cat(basename(folder),trait,"s",s,"n_snps",length(variants),"\n")
  }
}
res<-do.call(rbind,output)
write.table(res,file.path(subdir,"G0022_SHARED_SNP_AS_AIS_LD_DIAGNOSTIC.tsv"),
  row.names=FALSE,quote=FALSE,sep="\t")
