#!/usr/bin/env Rscript
# Ischemic stroke Phase10/11: fail-closed, one-prior direct SNP replay of
# ALL 646 ORIGINAL legacy GTEx gene-tissue tests (not 80 expanded windows).
# For missing inputs, emit MISSING_INPUT; for failures, emit FAIL, never PASS.
# This runner does not alter source files, inputs, notebooks or canonical data.
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})
if (as.character(packageVersion("coloc")) != "5.2.3")
  stop("COLOC_VERSION_PINNED_5.2.3")
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)!=4L)
  stop("Usage: Rscript script.R ORIGINAL_INDEX.tsv ORIGINAL_ABF_MASTER.tsv INPUT_DIR NEW_OUTPUT_DIR")
idx_path <- normalizePath(argv[[1]],mustWork=TRUE)
master_path <- normalizePath(argv[[2]],mustWork=TRUE)
input_dir <- normalizePath(argv[[3]],mustWork=TRUE)
out_dir <- argv[[4]]
if (dir.exists(out_dir))
  stop("REFUSE_ANY_EXISTING_OUTPUT_DIRECTORY")
idx <- fread(idx_path)
original <- fread(master_path)
if (nrow(idx)!=646L || nrow(original)!=646L)
  stop("LEGACY_646_DENOMINATOR_MISMATCH")
key <- c("locus","dataset_key","gene_base")
needed_idx <- c(key,"file","nsnps","qtl_n_first","qtl_n_min","qtl_n_max","status")
needed_orig <- c(key,"nsnps","qtl_n","status",paste0("PP.H",0:4))
if (!all(needed_idx %in% names(idx)) || !all(needed_orig %in% names(original)))
  stop("SOURCE_INDEX_OR_MASTER_SCHEMA_INVALID")
if (anyDuplicated(idx,by=key) || anyDuplicated(original,by=key))
  stop("SOURCE_COMPOSITE_KEY_NOT_UNIQUE")
joined <- merge(idx,original,by=key,suffixes=c("_index","_old"),sort=FALSE)
if (nrow(joined)!=646L || any(joined$status_index!="READY") ||
    any(joined$status_old!="PASS") ||
    any(joined$nsnps_index!=joined$nsnps_old) ||
    any(abs(joined$qtl_n_first-joined$qtl_n)>1e-8))
  stop("SOURCE_INDEX_MASTER_ALIGNMENT_INVALID")
setorder(joined,locus,dataset_key,gene_base)
N_total <- 174686
N_cases <- 22664
p1 <- 1e-4
p2 <- 1e-4
p12 <- 1e-5
source_fields <- c(key,"match_key","gwas_beta","gwas_se","gwas_maf",
                   "eqtl_beta","eqtl_se","eqtl_maf","eqtl_n","qtl_n_scalar",
                   "harmonization")
rows <- vector("list",nrow(joined))
for (i in seq_len(nrow(joined))) {
  s <- joined[i]
  fname <- basename(s$file)
  # Reject path traversal and unknown index inputs, and require name follows composite key.
  expected <- paste0(s$locus,"__",s$dataset_key,"__",s$gene_base,".tsv")
  if (!identical(fname,expected))
    stop(paste("FILENAME_INDEX_MISMATCH",i,fname,expected))
  path <- file.path(input_dir,fname)
  stat <- "MISSING_INPUT"
  why <- "SOURCE_FILE_NOT_PRESENT"
  n_snps <- NA_integer_
  seen_n_first <- NA_real_
  n_unique <- NA_integer_
  delta_n <- NA_real_
  maf_big <- NA_integer_
  harm_labels <- ""
  pph <- rep(NA_real_,5)
  max_diff <- NA_real_
  has_file <- file.exists(path)
  if (has_file) {
    tryCatch({
      x <- fread(path)
      if (!all(source_fields %in% names(x)))
        stop("SOURCE_COLUMNS_MISSING")
      if (nrow(x)!=s$nsnps_index || uniqueN(x$match_key)!=nrow(x))
        stop("SNP_DENOMINATOR_OR_DUPLICATE_KEY")
      if (uniqueN(x$locus)!=1L || uniqueN(x$dataset_key)!=1L ||
          uniqueN(x$gene_base)!=1L ||
          x$locus[1]!=s$locus || x$dataset_key[1]!=s$dataset_key ||
          x$gene_base[1]!=s$gene_base)
        stop("SOURCE_COMPOSITE_KEY_MISMATCH")
      vals <- as.matrix(x[,.(gwas_beta,gwas_se,gwas_maf,eqtl_beta,eqtl_se,
                            eqtl_maf,eqtl_n,qtl_n_scalar)])
      if (any(!is.finite(vals)) || any(x$gwas_se<=0) ||
          any(x$eqtl_se<=0) || any(x$gwas_maf<=0) ||
          any(x$eqtl_maf<=0) || any(x$gwas_maf>0.5) ||
          any(x$eqtl_maf>0.5))
        stop("INVALID_BETA_SE_MAF_OR_N")
      if (uniqueN(x$qtl_n_scalar)!=1L ||
          abs(as.numeric(x$qtl_n_scalar[1])-s$qtl_n_first)>1e-8 ||
          abs(as.numeric(x$eqtl_n[1])-s$qtl_n_first)>1e-8)
        stop("FIRST_QTL_N_NOT_REPRODUCIBLE")
      if (abs(min(x$eqtl_n)-s$qtl_n_min)>1e-8 ||
          abs(max(x$eqtl_n)-s$qtl_n_max)>1e-8)
        stop("SNP_SPECIFIC_N_RANGE_MISMATCH")
      # Provenance on allele orientation, NOT external cross-cohort validation.
      harm_labels <- paste(sort(unique(as.character(x$harmonization))),collapse=";")
      d1 <- list(beta=as.numeric(x$gwas_beta),
                 varbeta=as.numeric(x$gwas_se)^2,
                 snp=as.character(x$match_key),
                 MAF=as.numeric(x$gwas_maf),
                 N=N_total,s=N_cases/N_total,type="cc")
      d2 <- list(beta=as.numeric(x$eqtl_beta),
                 varbeta=as.numeric(x$eqtl_se)^2,
                 snp=as.character(x$match_key),
                 MAF=as.numeric(x$eqtl_maf),
                 N=as.numeric(s$qtl_n_first),type="quant")
      fit <- suppressMessages(suppressWarnings(
        coloc.abf(dataset1=d1,dataset2=d2,p1=p1,p2=p2,p12=p12)
      ))
      new <- as.numeric(fit$summary[paste0("PP.H",0:4,".abf")])
      old <- as.numeric(unlist(s[,paste0("PP.H",0:4),with=FALSE]))
      if (length(new)!=5L || any(!is.finite(new)) ||
          any(new<0) || abs(sum(new)-1)>1e-8)
        stop("COL0C_H0_H4_INVALID")
      max_diff <- max(abs(new-old))
      pph <- new
      n_snps <- nrow(x)
      n_unique <- uniqueN(x$match_key)
      seen_n_first <- as.numeric(x$qtl_n_scalar[1])
      delta_n <- max(x$eqtl_n)-min(x$eqtl_n)
      maf_big <- sum(abs(x$gwas_maf-x$eqtl_maf)>0.1)
      if (max_diff>1e-8) stop("HISTORICAL_H0_H4_NOT_REPRODUCED")
      stat <- "PASS"
      why <- "H0_H4_REPLAY_MATCH_1e-8"
    },error=function(e) {
      stat <<- "FAIL"
      why <<- gsub("[\r\n\t]"," ",conditionMessage(e))
    })
  }
  rows[[i]] <- data.table(
    locus=s$locus,dataset_key=s$dataset_key,gene_base=s$gene_base,
    filename=fname,status=stat,detail=why,
    expected_nsnps=s$nsnps_index,actual_nsnps=n_snps,unique_snps=n_unique,
    qtl_n_first_expected=s$qtl_n_first,qtl_n_first_actual=seen_n_first,
    original_qtl_n_range=s$qtl_n_max-s$qtl_n_min,
    observed_qtl_n_range=delta_n,maf_diff_gt_0_1=maf_big,
    input_harmonization_labels=harm_labels,
    PP_H0=pph[1],PP_H1=pph[2],PP_H2=pph[3],PP_H3=pph[4],
    PP_H4=pph[5],max_abs_hypothesis_delta=max_diff,
    model="DIRECT_SNP_ABF_ORIGINAL_FIRST_N",
    stringsAsFactors=FALSE
  )
  if (i%%50L==0L || i==nrow(joined)) {
    cat("BATCH_PROGRESS",i,"/",nrow(joined),
        "PASS",sum(vapply(rows[seq_len(i)],function(z) z$status=="PASS",logical(1))),
        "FAIL",sum(vapply(rows[seq_len(i)],function(z) z$status=="FAIL",logical(1))),
        "\n")
    flush.console()
  }
}
all <- rbindlist(rows)
if (nrow(all)!=646L || anyDuplicated(all,by=key))
  stop("FINAL_DENOMINATOR_MISMATCH")
counts <- all[, .N, by=status]
pass <- sum(all$status=="PASS")
failed <- sum(all$status=="FAIL")
missing <- sum(all$status=="MISSING_INPUT")
if (pass+failed+missing!=646L) stop("STATUS_PARTITION_INVALID")
# Use one new isolated output folder with all statuses; do not mislabel PARTIAL.
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
fwrite(all,file.path(out_dir,"IS_LEGACY_646_SNP_REPLAY_STATUS.tsv"),
       sep="\t",na="NA",quote=FALSE)
fwrite(all[status=="PASS"],file.path(out_dir,"IS_LEGACY_SNP_REPLAY_PASS_ONLY.tsv"),
       sep="\t",na="NA",quote=FALSE)
writeLines(c(
  "# Original BBJ–GTEx legacy coloc 646-input direct-SNP replay",
  paste0("- Status: ",if (pass==646L) "FULL_646_REPLAY_PASS" else "PARTIAL_OR_FAILED_REPLAY"),
  paste0("- PASS: ",pass,"; FAIL: ",failed,"; MISSING_INPUT: ",missing,
         "; total: 646"),
  "- PASS is assigned only after an original same-SNP-input H0–H4 coloc.abf 5.2.3 rerun matches archival hypotheses within 1e-8.",
  "- Missing input is never PASS. Errors are reported per gene/tissue; original source and output preserved.",
  "- Original p1=p2=1e-4, p12=1e-5; BBJ case-control N=174686, cases=22664; GTEx QTL N from original FIRST scalar.",
  "- All 646 results are from original FOUR BBJ loci; no expanded 80-window analyses have been performed here.",
  "- No false causal/eQTL-to-stroke gene claim: cohort ancestry differences, multi-signal LD, selected-tissue multiplicity remain.",
  "- QTL sample N varies per SNP; original FIRST N is reproduced as a historical computation, not a validated effective-N model."
),file.path(out_dir,"IS_LEGACY_646_REPLAY_README.md"))
cat("IS_LEGACY_646_BATCH_REPLAY=COMPLETE_INPUT_ACCOUNTING\n")
cat("BATCH_STATUS",if (pass==646L) "FULL_PASS" else "PARTIAL",
    "PASS",pass,"FAIL",failed,"MISSING",missing,"\n")
