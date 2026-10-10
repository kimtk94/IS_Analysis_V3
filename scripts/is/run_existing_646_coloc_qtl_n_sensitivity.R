#!/usr/bin/env Rscript
# Post-recovery GTEx-v8 BBJ-IS coloc ABF sample-N sensitivity only.
# Historical model uses FIRST SNP QTL scalar N despite within-assay variable N.
# No per-SNP N estimator is claimed; test several defensible SCALAR N choices.
# These share the same summary statistics, not independent studies.
suppressPackageStartupMessages(library(coloc))
args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key) {
 idx <- which(args==key)
 if(length(idx)!=1 || length(args)<idx+1)stop(paste("Expected",key))
 args[[idx+1]]
}
source_dir <- get_arg("--source-dir")
master_file <- get_arg("--master")
replay_file <- get_arg("--verified-replay")
out_dir <- get_arg("--out-dir")
if(!dir.exists(source_dir) || !file.exists(master_file) || !file.exists(replay_file))
 stop("Expected complete copied source and archived baseline/replay files")
if(dir.exists(out_dir) && length(list.files(out_dir,all.files=TRUE,no..=TRUE))>0)
 stop("Refusing to overwrite output")
master <- read.delim(master_file,check.names=FALSE,stringsAsFactors=FALSE)
replay <- read.delim(replay_file,check.names=FALSE,stringsAsFactors=FALSE)
key <- function(df)paste(df$locus,df$dataset_key,df$gene_base,sep="__")
if(nrow(master)!=646 || nrow(replay)!=646 || anyDuplicated(key(master)) ||
   anyDuplicated(key(replay)) || !setequal(key(master),key(replay)))
 stop("Replay and original coloc assay keys drifted")
if(any(replay$status!="PASS") || any(as.numeric(replay$delta_max)>1e-8))
 stop("Baseline 646 SNP-level ABF reproduction is not verified")
map <- match(key(master),key(replay))
required <- c("eqtl_beta","eqtl_se","eqtl_maf","eqtl_n",
 "gwas_beta","gwas_se","gwas_maf","match_key","qtl_n_scalar","harmonization")
records <- vector("list",646)
for(i in seq_len(646)) {
 r <- master[i,,drop=FALSE]
 filename <- paste0(key(r),".tsv")
 x <- read.delim(file.path(source_dir,filename),
                 check.names=FALSE,stringsAsFactors=FALSE)
 if(!all(required%in%names(x)) || nrow(x)!=as.integer(r$nsnps))
  stop(paste("Input schema/row count drift:",filename))
 for(c in c("eqtl_beta","eqtl_se","eqtl_maf","eqtl_n",
            "gwas_beta","gwas_se","gwas_maf","qtl_n_scalar"))
  x[[c]] <- as.numeric(x[[c]])
 if(any(!vapply(x[c("eqtl_beta","eqtl_se","eqtl_maf","eqtl_n",
        "gwas_beta","gwas_se","gwas_maf","qtl_n_scalar")],
        function(v)all(is.finite(v)),logical(1))))
  stop(paste("Nonfinite input:",filename))
 if(any(x$eqtl_se<=0 | x$gwas_se<=0 | x$eqtl_n<5 |
        x$eqtl_maf<=0 | x$eqtl_maf>0.5 | x$gwas_maf<=0 | x$gwas_maf>0.5))
  stop(paste("Input variance/frequency failed:",filename))
 if(anyDuplicated(x$match_key) || any(x$harmonization!="EXACT_REF_ALT_GRCH38"))
  stop(paste("SNP allele alignment failed:",filename))
 n_first <- x$qtl_n_scalar[[1]]
 if(abs(n_first-x$eqtl_n[[1]])>1e-8)
  stop(paste("Original sample N no longer matches first-snp eqtl_n:",filename))
 ns <- c(FIRST=n_first,MIN=min(x$eqtl_n),MEDIAN=as.numeric(round(median(x$eqtl_n))),
         MAX=max(x$eqtl_n))
 d1 <- list(snp=x$match_key,beta=x$gwas_beta,varbeta=x$gwas_se^2,
            MAF=x$gwas_maf,N=174686,s=22664/174686,type="cc")
 output <- c()
 for (mode in names(ns)) {
  d2 <- list(snp=x$match_key,beta=x$eqtl_beta,varbeta=x$eqtl_se^2,
             MAF=x$eqtl_maf,N=unname(ns[[mode]]),type="quant")
  null <- capture.output(
   fit <- suppressWarnings(suppressMessages(
     coloc.abf(dataset1=d1,dataset2=d2,p1=1e-4,p2=1e-4,p12=1e-5)))
  )
  h4 <- as.numeric(fit$summary[["PP.H4.abf"]])
  h3 <- as.numeric(fit$summary[["PP.H3.abf"]])
  if(!is.finite(h4) || !is.finite(h3))stop(paste("Nonfinite coloc posterior",filename,mode))
  output[[paste0("h4_",tolower(mode))]] <- h4
  output[[paste0("h3_",tolower(mode))]] <- h3
 }
 original_h4 <- as.numeric(r[["PP.H4"]])
 if(abs(original_h4-output[["h4_first"]])>1e-8)
  stop(paste("Baseline first-N replay disagrees with reference posterior",filename))
 variants <- c(output[["h4_min"]],output[["h4_median"]],output[["h4_max"]])
 delta <- max(abs(variants-original_h4))
 records[[i]] <- data.frame(
  locus=r$locus,dataset_key=r$dataset_key,gene_base=r$gene_base,
  n_variants=nrow(x),n_distinct_per_snp_N=length(unique(x$eqtl_n)),
  qtl_N_first=ns[["FIRST"]],qtl_N_min=ns[["MIN"]],
  qtl_N_median_rounded=ns[["MEDIAN"]],qtl_N_max=ns[["MAX"]],
  archived_H4=original_h4,
  H4_first=output[["h4_first"]],
  H4_min=output[["h4_min"]],
  H4_median=output[["h4_median"]],
  H4_max=output[["h4_max"]],
  H3_first=output[["h3_first"]],
  H3_min=output[["h3_min"]],
  H3_median=output[["h3_median"]],
  H3_max=output[["h3_max"]],
  max_abs_H4_shift_vs_first=delta,
  any_H4_crosses_0p5=any((variants>=0.5)!=(original_h4>=0.5)),
  any_H4_crosses_0p8=any((variants>=0.8)!=(original_h4>=0.8)),
  interpretation="SCALAR_QTL_N_SENSITIVITY_NOT_CAUSAL_OR_PER_SNP_N_MODEL",
  stringsAsFactors=FALSE,check.names=FALSE)
 if(i%%50==0 || i==646)
  cat("N_SENSITIVITY_PROGRESS",i,"of",646,
      "max_H4_delta",max(vapply(records[seq_len(i)],
          function(a)a$max_abs_H4_shift_vs_first,numeric(1))),"\n",file=stderr())
}
if(!dir.exists(out_dir))dir.create(out_dir,recursive=TRUE)
joined <- do.call(rbind,records)
write.table(joined,file.path(out_dir,"IS_646_COLOC_QTL_N_SCALAR_SENSITIVITY.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE,na="")
summary <- list(
 status="SCALAR_QTL_N_SENSITIVITY_ONLY",
 input_assays=646,models_per_assay=4,
 qtl_N_modes=c("FIRST_SOURCE_AS_ORIGINAL","MIN","ROUNDED_MEDIAN","MAX"),
 all_source_genes_reference_unchanged=TRUE,
 assays_with_variable_per_SNP_QTL_N=sum(joined$n_distinct_per_snp_N>1),
 max_abs_H4_shift=max(joined$max_abs_H4_shift_vs_first),
 median_abs_H4_shift=median(joined$max_abs_H4_shift_vs_first),
 assays_with_abs_H4_shift_gt_0p05=sum(joined$max_abs_H4_shift_vs_first>0.05),
 assays_with_abs_H4_shift_gt_0p20=sum(joined$max_abs_H4_shift_vs_first>0.2),
 assays_crossing_H4_0p5=sum(joined$any_H4_crosses_0p5),
 assays_crossing_H4_0p8=sum(joined$any_H4_crosses_0p8),
 p1=1e-4,p2=1e-4,p12=1e-5,
 GTEx_cohort_matched_LD_verified=FALSE,
 no_causal_inference=TRUE,
 note="N mode choice is not an estimate of individual-SNP N sensitivity; all sample-N models are scalar approximations."
)
if(requireNamespace("jsonlite",quietly=TRUE))
 jsonlite::write_json(summary,file.path(out_dir,"IS_646_COLOC_QTL_N_SENSITIVITY_MANIFEST.json"),
                      pretty=TRUE,auto_unbox=TRUE)
cat("IS_646_QTL_N_SENSITIVITY_COMPLETE",summary$input_assays,
    "MAX_DELTA",summary$max_abs_H4_shift,
    "N_SHIFT_GT_0P05",summary$assays_with_abs_H4_shift_gt_0p05,"\n")
