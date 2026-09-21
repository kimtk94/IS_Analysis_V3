#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("usage: run_ckd_stage2c_susie.R <susie_input_dir> <ld_dir> <stage2b_default.tsv> <output_dir>")
}
input_dir <- args[[1]]
ld_dir <- args[[2]]
abf_file <- args[[3]]
output_dir <- args[[4]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(susieR))

read_ld <- function(gene, dat) {
  meta_file <- file.path(ld_dir, paste0(gene, ".ld.meta.tsv"))
  vars_file <- file.path(ld_dir, paste0(gene, ".ld.vars"))
  bin_file <- file.path(ld_dir, paste0(gene, ".ld.bin"))
  if (!all(file.exists(c(meta_file, vars_file, bin_file)))) {
    stop(paste(gene, "LD files missing"))
  }
  meta <- read.delim(meta_file, stringsAsFactors=FALSE, check.names=FALSE)
  vars <- readLines(vars_file, warn=FALSE)
  vars <- vars[nzchar(vars)]
  if (nrow(meta) != length(vars) || any(meta$ref_id != vars)) {
    stop(paste(gene, "LD metadata/order mismatch"))
  }
  idx <- match(meta$snp, dat$snp)
  if (anyNA(idx)) stop(paste(gene, "LD SNP missing from summary input"))
  dat <- dat[idx, , drop=FALSE]
  n <- nrow(meta)
  con <- file(bin_file, "rb")
  on.exit(close(con), add=TRUE)
  vals <- readBin(con, what="numeric", n=n*n, size=4, endian="little")
  if (length(vals) != n*n) stop(paste(gene, "LD binary length mismatch"))
  LD <- matrix(vals, nrow=n, ncol=n, byrow=TRUE)
  # Matrix is symmetric, but explicitly symmetrize tiny binary/read artifacts.
  LD <- (LD + t(LD)) / 2
  signs <- as.numeric(meta$effect_vs_ld_major_sign)
  if (any(!signs %in% c(-1,1))) stop(paste(gene, "invalid LD orientation sign"))
  LD <- LD * tcrossprod(signs)
  diag(LD) <- 1
  dimnames(LD) <- list(dat$snp, dat$snp)
  if (any(!is.finite(LD))) stop(paste(gene, "non-finite LD values"))
  if (max(abs(LD - t(LD))) > 1e-5) stop(paste(gene, "LD is not symmetric"))
  list(dat=dat, LD=LD, meta=meta)
}

cs_count <- function(x) {
  if (is.null(x$sets) || is.null(x$sets$cs)) return(0L)
  length(x$sets$cs)
}

safe_susie <- function(d, suffix, maxit=1000L) {
  check_dataset(d, req="LD")
  tryCatch(
    {
      fit <- runsusie(
        d,
        suffix=suffix,
        maxit=maxit,
        repeat_until_convergence=FALSE,
        L=min(10L, max(1L, length(d$snp)-1L))
      )
      list(ok=TRUE, fit=fit, error="")
    },
    error=function(e) {
      list(ok=FALSE, fit=NULL, error=conditionMessage(e))
    }
  )
}

save_checkpoint <- function(all_rows, gene_rows, failures, output_dir) {
  if (length(all_rows)) {
    all_df <- do.call(rbind, all_rows)
    all_df <- all_df[order(all_df$gene_symbol, all_df$p12, -all_df$PP.H4), ]
    write.table(
      all_df,
      file=file.path(output_dir, "STAGE2C_SUSIE_ALL_PRIORS.tsv"),
      sep="\t", quote=FALSE, row.names=FALSE
    )
  }
  if (length(gene_rows)) {
    gene_df <- do.call(rbind, gene_rows)
    gene_df <- gene_df[order(-gene_df$max_PP.H4, na.last=TRUE), ]
    write.table(
      gene_df,
      file=file.path(output_dir, "STAGE2C_SUSIE_DEFAULT.tsv"),
      sep="\t", quote=FALSE, row.names=FALSE
    )
  }
  if (length(failures)) {
    fail_df <- do.call(rbind, failures)
    write.table(
      fail_df,
      file=file.path(output_dir, "STAGE2C_SUSIE_FAILURES.tsv"),
      sep="\t", quote=FALSE, row.names=FALSE
    )
  }
}

abf <- read.delim(abf_file, stringsAsFactors=FALSE, check.names=FALSE)
files <- sort(list.files(input_dir, pattern="\\.tsv\\.gz$", full.names=TRUE))
if (!length(files)) stop("no Stage 2C input files")

all_rows <- list()
gene_rows <- list()
failures <- list()

for (f in files) {
  gene <- sub("\\.tsv\\.gz$", "", basename(f))
  dat <- read.delim(gzfile(f), header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  required <- c(
    "snp","pos37","beta_pqtl","se_pqtl","maf_pqtl","n_pqtl",
    "beta_outcome_aligned","se_outcome","maf_outcome","n_outcome"
  )
  miss <- setdiff(required, names(dat))
  if (length(miss)) stop(paste(gene, "missing", paste(miss, collapse=",")))
  dat <- dat[
    is.finite(dat$beta_pqtl) & is.finite(dat$se_pqtl) & dat$se_pqtl > 0 &
    is.finite(dat$beta_outcome_aligned) & is.finite(dat$se_outcome) & dat$se_outcome > 0 &
    is.finite(dat$maf_pqtl) & dat$maf_pqtl > 0 & dat$maf_pqtl <= 0.5 &
    is.finite(dat$maf_outcome) & dat$maf_outcome > 0 & dat$maf_outcome <= 0.5,
    , drop=FALSE
  ]
  dat <- dat[!duplicated(dat$snp), , drop=FALSE]

  obj <- read_ld(gene, dat)
  dat <- obj$dat
  LD <- obj$LD
  n1 <- as.integer(round(median(dat$n_pqtl, na.rm=TRUE)))
  n2 <- as.integer(round(median(dat$n_outcome, na.rm=TRUE)))

  d1 <- list(
    beta=dat$beta_pqtl,
    varbeta=dat$se_pqtl^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    MAF=dat$maf_pqtl,
    N=n1,
    sdY=1,
    LD=LD,
    type="quant"
  )
  d2 <- list(
    beta=dat$beta_outcome_aligned,
    varbeta=dat$se_outcome^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    MAF=dat$maf_outcome,
    N=n2,
    LD=LD,
    type="quant"
  )

  cat("Stage2C SuSiE:", gene, "n=", nrow(dat), "\n")
  p1_rds <- file.path(output_dir, paste0(gene, "_pqtl_susie.rds"))
  p2_rds <- file.path(output_dir, paste0(gene, "_egfr_susie.rds"))

  if (file.exists(p1_rds)) {
    S1 <- readRDS(p1_rds)
    s1 <- list(ok=isTRUE(S1$converged), fit=S1,
               error=if (isTRUE(S1$converged)) "" else "cached pQTL fit not converged")
  } else {
    s1 <- safe_susie(d1, paste0(gene, "_pqtl"))
    if (s1$ok) saveRDS(s1$fit, file=p1_rds)
  }

  if (file.exists(p2_rds)) {
    S2 <- readRDS(p2_rds)
    s2 <- list(ok=isTRUE(S2$converged), fit=S2,
               error=if (isTRUE(S2$converged)) "" else "cached eGFR fit not converged")
  } else {
    s2 <- safe_susie(d2, paste0(gene, "_egfr"))
    if (s2$ok) saveRDS(s2$fit, file=p2_rds)
  }

  if (!s1$ok || !s2$ok) {
    failures[[length(failures)+1]] <- data.frame(
      gene_symbol=gene,
      nsnps=nrow(dat),
      pqtl_status=if (s1$ok) "converged" else "failed_or_nonconverged",
      egfr_status=if (s2$ok) "converged" else "failed_or_nonconverged",
      pqtl_message=s1$error,
      egfr_message=s2$error,
      stringsAsFactors=FALSE
    )
    a <- abf[abf$gene_symbol == gene, , drop=FALSE]
    abf_h4 <- if (nrow(a)) as.numeric(a$PP.H4[1]) else NA_real_
    abf_h3 <- if (nrow(a)) as.numeric(a$PP.H3[1]) else NA_real_
    gene_rows[[length(gene_rows)+1]] <- data.frame(
      gene_symbol=gene,
      nsnps=nrow(dat),
      pqtl_credible_sets=if (s1$ok) cs_count(s1$fit) else NA_integer_,
      egfr_credible_sets=if (s2$ok) cs_count(s2$fit) else NA_integer_,
      signal_pairs_default=0L,
      max_PP.H4=NA_real_,
      PP.H3_at_best_pair=NA_real_,
      best_pqtl_signal="",
      best_egfr_signal="",
      any_H4_ge_0.8=0L,
      any_H4_ge_0.9=0L,
      stage2b_ABF_PP.H3=abf_h3,
      stage2b_ABF_PP.H4=abf_h4,
      comparison_status="SuSiE_nonconverged_or_failed",
      N_pqtl=n1,
      N_egfr=n2,
      stringsAsFactors=FALSE
    )
    save_checkpoint(all_rows, gene_rows, failures, output_dir)
    cat("Stage2C SuSiE nonconvergence/failure:", gene,
        "pQTL=", s1$error, "eGFR=", s2$error, "\n")
    rm(LD, d1, d2)
    gc()
    next
  }

  S1 <- s1$fit
  S2 <- s2$fit
  default_res <- NULL
  for (p12 in c(1e-6, 1e-5, 1e-4)) {
    res <- coloc.susie(S1, S2, p1=1e-4, p2=1e-4, p12=p12)
    sm <- as.data.frame(res$summary, stringsAsFactors=FALSE)
    if (!nrow(sm)) {
      all_rows[[length(all_rows)+1]] <- data.frame(
        gene_symbol=gene,p1=1e-4,p2=1e-4,p12=p12,nsnps=nrow(dat),
        hit1="",hit2="",PP.H0=NA,PP.H1=NA,PP.H2=NA,PP.H3=NA,PP.H4=NA,
        idx1=NA,idx2=NA,stringsAsFactors=FALSE
      )
    } else {
      for (i in seq_len(nrow(sm))) {
        all_rows[[length(all_rows)+1]] <- data.frame(
          gene_symbol=gene,p1=1e-4,p2=1e-4,p12=p12,
          nsnps=as.integer(sm$nsnps[i]),
          hit1=as.character(sm$hit1[i]),hit2=as.character(sm$hit2[i]),
          PP.H0=as.numeric(sm$PP.H0.abf[i]),
          PP.H1=as.numeric(sm$PP.H1.abf[i]),
          PP.H2=as.numeric(sm$PP.H2.abf[i]),
          PP.H3=as.numeric(sm$PP.H3.abf[i]),
          PP.H4=as.numeric(sm$PP.H4.abf[i]),
          idx1=as.integer(sm$idx1[i]),idx2=as.integer(sm$idx2[i]),
          stringsAsFactors=FALSE
        )
      }
    }
    if (p12 == 1e-5) default_res <- res
  }

  default_rows <- Filter(function(x) x$gene_symbol[1] == gene && x$p12[1] == 1e-5, all_rows)
  dr <- if (length(default_rows)) do.call(rbind, default_rows) else data.frame()
  if (nrow(dr) && any(is.finite(dr$PP.H4))) {
    best <- dr[which.max(ifelse(is.finite(dr$PP.H4), dr$PP.H4, -Inf)), , drop=FALSE]
    max_h4 <- best$PP.H4
    best_h3 <- best$PP.H3
    best_hit1 <- best$hit1
    best_hit2 <- best$hit2
    pair_n <- nrow(dr)
  } else {
    max_h4 <- NA_real_; best_h3 <- NA_real_; best_hit1 <- ""; best_hit2 <- ""; pair_n <- 0L
  }
  a <- abf[abf$gene_symbol == gene, , drop=FALSE]
  abf_h4 <- if (nrow(a)) as.numeric(a$PP.H4[1]) else NA_real_
  abf_h3 <- if (nrow(a)) as.numeric(a$PP.H3[1]) else NA_real_

  status <- "unresolved"
  if (is.finite(max_h4) && max_h4 >= 0.8 && is.finite(abf_h4) && abf_h4 >= 0.8) status <- "ABF_and_SuSiE_shared_signal"
  if (is.finite(max_h4) && max_h4 >= 0.8 && is.finite(abf_h4) && abf_h4 < 0.2) status <- "SuSiE_multisignal_rescue"
  if (is.finite(max_h4) && max_h4 < 0.8 && is.finite(abf_h4) && abf_h4 >= 0.8) status <- "ABF_not_confirmed_by_SuSiE"

  gene_rows[[length(gene_rows)+1]] <- data.frame(
    gene_symbol=gene,
    nsnps=nrow(dat),
    pqtl_credible_sets=cs_count(S1),
    egfr_credible_sets=cs_count(S2),
    signal_pairs_default=pair_n,
    max_PP.H4=max_h4,
    PP.H3_at_best_pair=best_h3,
    best_pqtl_signal=best_hit1,
    best_egfr_signal=best_hit2,
    any_H4_ge_0.8=as.integer(is.finite(max_h4) && max_h4 >= 0.8),
    any_H4_ge_0.9=as.integer(is.finite(max_h4) && max_h4 >= 0.9),
    stage2b_ABF_PP.H3=abf_h3,
    stage2b_ABF_PP.H4=abf_h4,
    comparison_status=status,
    N_pqtl=n1,
    N_egfr=n2,
    stringsAsFactors=FALSE
  )
  save_checkpoint(all_rows, gene_rows, failures, output_dir)
  rm(LD, S1, S2, d1, d2)
  gc()
}

all_df <- if (length(all_rows)) do.call(rbind, all_rows) else data.frame()
gene_df <- if (length(gene_rows)) do.call(rbind, gene_rows) else data.frame()

if (nrow(all_df)) {
  all_df <- all_df[order(all_df$gene_symbol, all_df$p12, -all_df$PP.H4), ]
  write.table(
    all_df,
    file=file.path(output_dir, "STAGE2C_SUSIE_ALL_PRIORS.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
}
if (nrow(gene_df)) {
  gene_df <- gene_df[order(-gene_df$max_PP.H4, na.last=TRUE), ]
  write.table(
    gene_df,
    file=file.path(output_dir, "STAGE2C_SUSIE_DEFAULT.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
}
if (length(failures)) {
  fail_df <- do.call(rbind, failures)
  write.table(
    fail_df,
    file=file.path(output_dir, "STAGE2C_SUSIE_FAILURES.tsv"),
    sep="\t", quote=FALSE, row.names=FALSE
  )
  cat("CKD_STAGE2C_SUSIE_PARTIAL_PASS failures=", nrow(fail_df), "\n", sep="")
} else {
  cat("CKD_STAGE2C_SUSIE_PASS\n")
}
print(gene_df)
