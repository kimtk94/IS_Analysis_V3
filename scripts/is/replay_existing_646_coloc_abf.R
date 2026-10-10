#!/usr/bin/env Rscript
# Read-only independent per-SNP coloc.abf reproduction from previously archived
# BBJ IS + GTEx v8 harmonized summary-statistics inputs.
#
# This is *not* conditional fine-mapping, independent QTL LD confirmation,
# causal gene identification, nor EAS independent stroke replication.
#
# Outputs are written only to a new, empty derivative output directory.
suppressPackageStartupMessages(library(coloc))

args <- commandArgs(trailingOnly=TRUE)
getopt <- function(key, default=NULL) {
  i <- which(args==key)
  if(length(i)==0) return(default)
  if(length(i)!=1 || length(args)<i+1) stop(paste("Invalid arg",key))
  args[[i+1]]
}
source_dir <- getopt("--source-dir")
master_file <- getopt("--master")
outdir <- getopt("--out-dir")
limit_arg <- getopt("--limit", "0")
if(any(vapply(list(source_dir,master_file,outdir),is.null,logical(1)))) {
  stop("Usage: --source-dir DIR --master COLOC_ABF_MASTER.tsv --out-dir EMPTY_DIR [--limit 0|N]")
}
limit <- as.integer(limit_arg)
if(is.na(limit) || limit<0) stop("--limit must be nonnegative")
if(!dir.exists(source_dir)) stop("Source input directory does not exist")
if(!file.exists(master_file)) stop("Archived coloc ABF master missing")
if(dir.exists(outdir) && length(list.files(outdir, all.files=TRUE, no..=TRUE))>0) {
  stop("Output already exists and is nonempty; refuse overwrite")
}
master <- read.delim(master_file,check.names=FALSE,stringsAsFactors=FALSE)
required_master <- c("locus","dataset_key","gene_base","nsnps",
                     "PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")
if(!all(required_master%in%names(master)) || nrow(master)!=646) stop("Historic coloc master schema/count mismatch")
key <- paste(master$locus,master$dataset_key,master$gene_base,sep="::")
if(anyDuplicated(key)) stop("Duplicate original archived locus/dataset/gene")
if(limit>0) master <- head(master,limit)
required_input <- c("locus","dataset_key","gene_base","match_key",
                    "gwas_beta","gwas_se","gwas_maf",
                    "eqtl_beta","eqtl_se","eqtl_maf",
                    "qtl_n_scalar","eqtl_n","harmonization")
results <- vector("list",nrow(master))
for(i in seq_len(nrow(master))) {
  entry <- master[i,,drop=FALSE]
  file <- file.path(source_dir,paste(entry$locus,entry$dataset_key,entry$gene_base,sep="__"))
  file <- paste0(file,".tsv")
  result <- list(locus=entry$locus,dataset_key=entry$dataset_key,
                 gene_base=entry$gene_base,source_file=basename(file),
                 expected_nsnps=entry$nsnps,status="NOT_RUN",
                 delta_max="",maf_diff_gt_0p1="",maf_diff_fraction="",n_qtl_sample_N_distinct="",qtl_N_per_SNP_min="",qtl_N_per_SNP_max="",new_PP_H0="",new_PP_H1="",new_PP_H2="",
                 new_PP_H3="",new_PP_H4="",detail="")
  if(!file.exists(file)) {
    result$status <- "MISSING_INPUT"
    result$detail <- "FILE_NOT_PRESENT_IN_SELECTED_SOURCE_DIR"
  } else {
    tryCatch({
      d <- read.delim(file,check.names=FALSE,stringsAsFactors=FALSE)
      if(!all(required_input%in%names(d))) stop("MISSING_REQUIRED_INPUT_COLUMN")
      if(nrow(d)!=as.integer(entry$nsnps)) stop("VARIANT_COUNT_NOT_MATCHED")
      if(anyDuplicated(d$match_key) || anyNA(d$match_key)) stop("DUPLICATED_OR_MISSING_SNP")
      if(!all(d$locus==entry$locus & d$dataset_key==entry$dataset_key & d$gene_base==entry$gene_base))
        stop("LOCUS_DATASET_GENE_DISCORDANCE")
      if(!all(d$harmonization=="EXACT_REF_ALT_GRCH38")) stop("ALLELE_ALIGNMENT_UNEXPECTED")
      numeric_cols <- c("gwas_beta","gwas_se","gwas_maf",
                        "eqtl_beta","eqtl_se","eqtl_maf","qtl_n_scalar","eqtl_n")
      for(n in numeric_cols) d[[n]] <- as.numeric(d[[n]])
      if(any(!vapply(d[numeric_cols],function(x)all(is.finite(x)),logical(1)))) stop("NONFINITE_SUMMARY_STATS")
      if(any(d$gwas_se<=0 | d$eqtl_se<=0 | d$qtl_n_scalar<5 | d$eqtl_n<5)) stop("INVALID_SE_OR_N")
      if(any(d$gwas_maf<=0 | d$gwas_maf>0.5 |
             d$eqtl_maf<=0 | d$eqtl_maf>0.5)) stop("INVALID_MAF")
      result$maf_diff_gt_0p1 <- sum(abs(d$gwas_maf-d$eqtl_maf)>0.1)
      result$maf_diff_fraction <- result$maf_diff_gt_0p1/nrow(d)
      result$n_qtl_sample_N_distinct <- length(unique(d$eqtl_n))
      result$qtl_N_per_SNP_min <- min(d$eqtl_n)
      result$qtl_N_per_SNP_max <- max(d$eqtl_n)
      if(abs(d$qtl_n_scalar[[1]]-d$eqtl_n[[1]])>1e-8) stop("ORIGINAL_QTL_FIRST_N_DISAGREES_WITH_PER_SNP_SOURCE")
      sample_n <- d$qtl_n_scalar[[1]]
      ds1 <- list(snp=d$match_key, beta=d$gwas_beta,
                  varbeta=d$gwas_se^2,MAF=d$gwas_maf,
                  N=174686, s=22664/174686, type="cc")
      ds2 <- list(snp=d$match_key,beta=d$eqtl_beta,
                  varbeta=d$eqtl_se^2,MAF=d$eqtl_maf,
                  N=sample_n,type="quant")
      fit <- suppressWarnings(suppressMessages(coloc.abf(dataset1=ds1,dataset2=ds2,
                                        p1=1e-4,p2=1e-4,p12=1e-5)))
      values <- as.numeric(fit$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf",
                                         "PP.H3.abf","PP.H4.abf")])
      expected <- as.numeric(entry[c("PP.H0","PP.H1","PP.H2","PP.H3","PP.H4")])
      if(any(!is.finite(values)) || any(!is.finite(expected))) stop("NONFINITE_POSTERIOR")
      delta <- max(abs(values-expected))
      for(j in 0:4) result[[paste0("new_PP_H",j)]] <- values[[j+1]]
      result$delta_max <- delta
      result$status <- if(delta<=1e-8) "PASS" else "POSTERIOR_DISCORDANCE"
      result$detail <- paste0("coloc=",as.character(packageVersion("coloc")),
                              ";qtl_N_first=",sample_n,
                              ";gwas_N=174686;cases=22664")
    },error=function(e){
      result$status <<- "SCIENTIFIC_QC_REVIEW"
      result$detail <<- substr(conditionMessage(e),1,250)
    })
  }
  results[[i]] <- result
  if(i %% 25==0 || i==nrow(master)) {
    tbl <- table(vapply(results[seq_len(i)],function(x)x$status,character(1)))
    cat("PROGRESS",i,"of",nrow(master),paste(paste(names(tbl),tbl,sep=":"),collapse=","),"\n",file=stderr())
  }
}
if(!dir.exists(outdir)) dir.create(outdir,recursive=TRUE)
out <- as.data.frame(do.call(rbind,lapply(results,as.data.frame)),
                     stringsAsFactors=FALSE,check.names=FALSE)
write.table(out,file.path(outdir,"IS_646_GWAS_QTL_SNP_ABF_REPLAY.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE,na="")
qc <- as.list(table(out$status))
summary <- list(status="PER_SNP_ABF_RERUN_DIAGNOSTIC_NO_CAUSAL_CLAIM",
                total=nrow(out),original_expected_tests=646,
                status_counts=qc,
                p1=1e-4,p2=1e-4,p12=1e-5,
                coloc_version=as.character(packageVersion("coloc")),
                GWAS="BBJ_IS_case_control_174686_N_22664_cases",
                molecular_QTL="GTEx_V8_bulk_eqtl_original_first_scalar_N",
                cohort_matched_molecular_QTL_LD_verified=FALSE,
                no_causal_gene_claim=TRUE,
                maf_discordance_is_not_snp_harmonization_proof=TRUE,
                per_snp_QTL_sample_N_range_recorded=TRUE,
                out_file="IS_646_GWAS_QTL_SNP_ABF_REPLAY.tsv")
if(requireNamespace("jsonlite",quietly=TRUE))
 jsonlite::write_json(summary,file.path(outdir,"IS_646_SNP_ABF_REPLAY_MANIFEST.json"),pretty=TRUE,auto_unbox=TRUE)
cat("IS_COLOC_ABF_REPLAY_COMPLETE",nrow(out),
    paste(paste(names(qc),unlist(qc),sep=":"),collapse=","),"\n")
