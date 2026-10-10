#!/usr/bin/env Rscript
# Direct per-SNP coloc.abf re-estimation for all 646 archived BBJ IS x GTEx v8
# gene-tissue pairs under a PREDECLARED p12 prior grid.
# Keep all 646 assays; never cherry-pick an H4-maximizing prior.
# No cohort-matched GTEx-QTL LD / multi-signal validation is claimed.
suppressPackageStartupMessages(library(coloc))
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(k) {
  ix <- which(args==k)
  if(length(ix)!=1 || ix>=length(args)) stop(paste("Missing",k))
  args[[ix+1]]
}
indir <- getarg("--source-dir")
master_file <- getarg("--master")
replay_file <- getarg("--verified-replay")
out_dir <- getarg("--out-dir")
priors <- c(1e-6,3e-6,1e-5,3e-5,1e-4)
if(!dir.exists(indir) || !file.exists(master_file) || !file.exists(replay_file))
 stop("Missing existing SNP source or baseline master")
if(dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE))>0)
 stop("Refuse to overwrite output")

master <- read.delim(master_file,check.names=FALSE,stringsAsFactors=FALSE)
baseline <- read.delim(replay_file,check.names=FALSE,stringsAsFactors=FALSE)
key <- function(df) paste(df$locus,df$dataset_key,df$gene_base,sep="__")
if(nrow(master)!=646 || nrow(baseline)!=646 || anyDuplicated(key(master)) ||
   anyDuplicated(key(baseline)) || !setequal(key(master),key(baseline)))
 stop("Original 646 assay universe mismatch")
if(any(baseline$status!="PASS") || any(as.numeric(baseline$delta_max)>1e-8))
 stop("Complete direct SNP baseline replay verification is required")
replay_by_key <- baseline[match(key(master),key(baseline)),]
required <- c("locus","dataset_key","gene_base","match_key",
 "gwas_beta","gwas_se","gwas_maf","eqtl_beta","eqtl_se","eqtl_maf",
 "qtl_n_scalar","harmonization")
records <- vector("list",646)
for(i in seq_len(646)) {
  r <- master[i,,drop=FALSE]
  file <- file.path(indir,paste0(key(r),".tsv"))
  if(!file.exists(file))stop(paste("Source missing",file))
  d <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE)
  if(nrow(d)!=as.integer(r$nsnps) || !all(required%in%names(d)))
   stop(paste("Source schema or SNP-N mismatch",basename(file)))
  if(anyDuplicated(d$match_key) || anyNA(d$match_key) ||
     any(d$harmonization!="EXACT_REF_ALT_GRCH38") ||
     !all(d$locus==r$locus & d$dataset_key==r$dataset_key & d$gene_base==r$gene_base))
   stop(paste("Variant or SNP-alignment source QC failed",basename(file)))
  for(col in c("gwas_beta","gwas_se","gwas_maf","eqtl_beta","eqtl_se","eqtl_maf","qtl_n_scalar"))
   d[[col]] <- suppressWarnings(as.numeric(d[[col]]))
  if(any(!vapply(d[c("gwas_beta","gwas_se","gwas_maf","eqtl_beta",
                    "eqtl_se","eqtl_maf","qtl_n_scalar")],
          function(x) all(is.finite(x)),logical(1))))
   stop(paste("Nonfinite source value",basename(file)))
  if(any(d$gwas_se<=0 | d$eqtl_se<=0 | d$gwas_maf<=0 | d$eqtl_maf<=0 |
         d$gwas_maf>0.5 | d$eqtl_maf>0.5 | d$qtl_n_scalar<5))
   stop(paste("Bad SE/MAF/N",basename(file)))
  ds1 <- list(snp=d$match_key,beta=d$gwas_beta,varbeta=d$gwas_se^2,
              MAF=d$gwas_maf,N=174686,s=22664/174686,type="cc")
  ds2 <- list(snp=d$match_key,beta=d$eqtl_beta,varbeta=d$eqtl_se^2,
              MAF=d$eqtl_maf,N=d$qtl_n_scalar[[1]],type="quant")
  one <- vector("list",length(priors))
  for(j in seq_along(priors)) {
   null <- capture.output(
     fit <- suppressWarnings(suppressMessages(coloc.abf(
       dataset1=ds1,dataset2=ds2,p1=1e-4,p2=1e-4,p12=priors[[j]])))
   )
   posterior <- as.numeric(fit$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                        "PP.H3.abf","PP.H4.abf")])
   if(length(posterior)!=5 || any(!is.finite(posterior)) ||
      abs(sum(posterior)-1)>1e-8)
    stop("Invalid posterior probabilities")
   one[[j]] <- data.frame(
     locus=r$locus,dataset_key=r$dataset_key,gene_base=r$gene_base,
     n_snps=nrow(d),p1=1e-4,p2=1e-4,p12=priors[[j]],
     PP_H0=posterior[[1]],PP_H1=posterior[[2]],PP_H2=posterior[[3]],
     PP_H3=posterior[[4]],PP_H4=posterior[[5]],
     H4_gt_H3=posterior[[5]]>posterior[[4]],
     H4_ge_0p5=posterior[[5]]>=0.5,
     H4_ge_0p8=posterior[[5]]>=0.8,
     baseline_H4=as.numeric(r[["PP.H4"]]),
     n_QTL_first=d$qtl_n_scalar[[1]],
     interpretation="PRIOR_SENSITIVITY_ONLY_NOT_CAUSAL_COLOCALIZATION",
     stringsAsFactors=FALSE,check.names=FALSE)
   if(priors[[j]]==1e-5 &&
      max(abs(posterior-as.numeric(r[c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")])))>1e-8)
    stop("Original baseline not reproduced at p12=1e-5")
  }
  records[[i]] <- do.call(rbind,one)
  if(i%%50==0 || i==646)
    cat("P12_PROGRESS",i,"of 646\n",file=stderr())
}
if(!dir.exists(out_dir)) dir.create(out_dir,recursive=TRUE)
grid <- do.call(rbind,records)
out <- file.path(out_dir,"IS_646_ABF_DIRECT_P12_PRIOR_GRID.tsv")
write.table(grid,out,sep="\t",quote=FALSE,row.names=FALSE,na="")
by_prior <- lapply(priors,function(x) {
 sub <- grid[grid$p12==x,,drop=FALSE]
 list(p12=x,n=nrow(sub),H4_ge_0p5=sum(sub$H4_ge_0p5),
      H4_ge_0p8=sum(sub$H4_ge_0p8),H4_gt_H3=sum(sub$H4_gt_H3))
})
base <- grid[grid$p12==1e-5,,drop=FALSE]
by_assay <- split(grid,interaction(grid$locus,grid$dataset_key,grid$gene_base,drop=TRUE))
robust <- lapply(by_assay,function(d) {
 list(assay=paste(d$locus[[1]],d$dataset_key[[1]],d$gene_base[[1]],sep="__"),
      min_H4=min(d$PP_H4),max_H4=max(d$PP_H4),
      robust_H4_gt_H3=all(d$H4_gt_H3),
      robust_H4_ge_0p5=all(d$H4_ge_0p5),
      robust_H4_ge_0p8=all(d$H4_ge_0p8))
})
summary <- list(
 status="DIRECT_SNP_COLOC_ABF_P12_PRIOR_SENSITIVITY_NOT_VALIDATED_CAUSALITY",
 original_assays=646,grid_p12=priors,models_per_assay=5,total_models=nrow(grid),
 coloc_version=as.character(packageVersion("coloc")),
 p1=1e-4,p2=1e-4,
 by_prior=by_prior,
 any_assays_with_H4_ge_0p8_at_any_prior=sum(vapply(robust,function(r)r$max_H4>=0.8,logical(1))),
 assays_robust_H4_ge_0p8_all_priors=sum(vapply(robust,function(r)r$robust_H4_ge_0p8,logical(1))),
 assays_robust_H4_gt_H3_all_priors=sum(vapply(robust,function(r)r$robust_H4_gt_H3,logical(1))),
 baseline_p12_1e5_atleast_H4_0p8=sum(base$PP_H4>=0.8),
 matched_GTex_molecular_QTL_LD=FALSE,
 single_causal_signal_assumption_unverified=TRUE,
 source_dataset_BBJ_not_independent_from_GIGASTROKE_EAS=TRUE,
 causal_genes_established=0,
 warning="Do not choose a p12 prior post hoc to inflate PP.H4, do not treat correlated assay/tissue results as independent associations."
)
if(requireNamespace("jsonlite",quietly=TRUE))
 jsonlite::write_json(summary,file.path(out_dir,"IS_646_DIRECT_P12_SENSITIVITY_MANIFEST.json"),
                     pretty=TRUE,auto_unbox=TRUE)
cat("IS_646_DIRECT_P12_GRID_COMPLETE",nrow(grid),"MODELS",
    "MAX_H4_BASELINE",max(base$PP_H4),"\n")
