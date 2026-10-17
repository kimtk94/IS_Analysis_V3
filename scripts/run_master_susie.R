#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 6 || length(args) > 7) {
  stop(paste(
    "usage: run_master_susie.R",
    "<input_dir> <ld_dir> <abf_default.tsv> <output_dir>",
    "<ancestry> <outcome_type=quant|cc|auto> [default_p12]"
  ))
}

input_dir <- args[[1]]
ld_dir <- args[[2]]
abf_file <- args[[3]]
output_dir <- args[[4]]
ancestry <- toupper(args[[5]])
outcome_type_arg <- tolower(args[[6]])
default_p12 <- if (length(args) >= 7) as.numeric(args[[7]]) else 1e-5

if (!ancestry %in% c("EUR","EAS","AFR","SAS","AMR")) stop("unsupported ancestry")
if (!outcome_type_arg %in% c("quant","cc","auto")) stop("outcome_type must be quant, cc, or auto")
if (!is.finite(default_p12) || default_p12 <= 0 || default_p12 >= 1) stop("invalid default_p12")

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages(library(coloc))
suppressPackageStartupMessages(library(susieR))

read_ld <- function(locus, dat) {
  meta_file <- file.path(ld_dir, paste0(locus, ".ld.meta.tsv"))
  vars_file <- file.path(ld_dir, paste0(locus, ".ld.vars"))
  bin_file <- file.path(ld_dir, paste0(locus, ".ld.bin"))
  if (!all(file.exists(c(meta_file, vars_file, bin_file)))) {
    stop(paste(locus, "LD files missing"))
  }

  meta <- read.delim(meta_file, stringsAsFactors=FALSE, check.names=FALSE)
  vars <- readLines(vars_file, warn=FALSE)
  vars <- vars[nzchar(vars)]

  required_meta <- c("ref_id","snp","effect_vs_ld_major_sign")
  miss_meta <- setdiff(required_meta, names(meta))
  if (length(miss_meta)) stop(paste(locus, "LD metadata missing", paste(miss_meta, collapse=",")))

  if ("ancestry" %in% names(meta)) {
    observed <- unique(toupper(as.character(meta$ancestry)))
    observed <- observed[nzchar(observed)]
    if (length(observed) != 1 || observed[[1]] != ancestry) {
      stop(paste(locus, "LD ancestry mismatch:", paste(observed, collapse=","), "!=", ancestry))
    }
  } else {
    warning(paste(locus, "LD metadata has no ancestry column; ancestry provenance cannot be verified"))
  }

  if (nrow(meta) != length(vars) || any(meta$ref_id != vars)) {
    stop(paste(locus, "LD metadata/order mismatch"))
  }

  idx <- match(meta$snp, dat$snp)
  if (anyNA(idx)) stop(paste(locus, "LD SNP missing from summary input"))
  dat <- dat[idx, , drop=FALSE]

  n <- nrow(meta)
  con <- file(bin_file, "rb")
  on.exit(close(con), add=TRUE)
  vals <- readBin(con, what="numeric", n=n*n, size=4, endian="little")
  if (length(vals) != n*n) stop(paste(locus, "LD binary length mismatch"))

  LD <- matrix(vals, nrow=n, ncol=n, byrow=TRUE)
  LD <- (LD + t(LD)) / 2
  signs <- as.numeric(meta$effect_vs_ld_major_sign)
  if (any(!signs %in% c(-1,1))) stop(paste(locus, "invalid LD orientation sign"))
  LD <- LD * tcrossprod(signs)
  diag(LD) <- 1
  dimnames(LD) <- list(dat$snp, dat$snp)

  if (any(!is.finite(LD))) stop(paste(locus, "non-finite LD values"))
  if (max(abs(LD-t(LD))) > 1e-5) stop(paste(locus, "LD is not symmetric"))

  list(dat=dat, LD=LD, meta=meta)
}

resolve_outcome_type <- function(dat) {
  if (outcome_type_arg != "auto") return(outcome_type_arg)
  if ("outcome_type" %in% names(dat)) {
    vals <- unique(tolower(as.character(dat$outcome_type)))
    vals <- vals[vals %in% c("quant","cc")]
    if (length(vals) == 1) return(vals[[1]])
  }
  "quant"
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
      list(ok=isTRUE(fit$converged), fit=fit,
           error=if (isTRUE(fit$converged)) "" else "SuSiE did not converge")
    },
    error=function(e) list(ok=FALSE, fit=NULL, error=conditionMessage(e))
  )
}

abf <- read.delim(abf_file, stringsAsFactors=FALSE, check.names=FALSE)
abf_key <- if ("locus" %in% names(abf)) "locus" else if ("gene_symbol" %in% names(abf)) "gene_symbol" else NA_character_
if (is.na(abf_key)) stop("ABF file requires locus or gene_symbol column")

files <- sort(list.files(input_dir, pattern="\\.tsv\\.gz$", full.names=TRUE))
if (!length(files)) stop("no SuSiE input files")

all_rows <- list()
locus_rows <- list()
failures <- list()

write_checkpoint <- function() {
  if (length(all_rows)) {
    x <- do.call(rbind, all_rows)
    x <- x[order(x$locus, x$p12, -x$PP.H4), ]
    write.table(x, file=file.path(output_dir,"MASTER_SUSIE_ALL_PRIORS.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)
  }
  if (length(locus_rows)) {
    x <- do.call(rbind, locus_rows)
    x <- x[order(-x$max_PP.H4, na.last=TRUE), ]
    write.table(x, file=file.path(output_dir,"MASTER_SUSIE_DEFAULT.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)
  }
  if (length(failures)) {
    x <- do.call(rbind, failures)
    write.table(x, file=file.path(output_dir,"MASTER_SUSIE_FAILURES.tsv"),
                sep="\t", quote=FALSE, row.names=FALSE)
  }
}

for (f in files) {
  locus <- sub("\\.tsv\\.gz$", "", basename(f))
  dat <- read.delim(gzfile(f), header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

  required <- c(
    "snp","pos37","beta_pqtl","se_pqtl","maf_pqtl","n_pqtl",
    "beta_outcome_aligned","se_outcome","maf_outcome","n_outcome"
  )
  miss <- setdiff(required, names(dat))
  if (length(miss)) stop(paste(locus, "missing", paste(miss, collapse=",")))

  dat <- dat[
    is.finite(dat$beta_pqtl) & is.finite(dat$se_pqtl) & dat$se_pqtl > 0 &
    is.finite(dat$beta_outcome_aligned) & is.finite(dat$se_outcome) & dat$se_outcome > 0 &
    is.finite(dat$maf_pqtl) & dat$maf_pqtl > 0 & dat$maf_pqtl <= 0.5 &
    is.finite(dat$maf_outcome) & dat$maf_outcome > 0 & dat$maf_outcome <= 0.5,
    , drop=FALSE
  ]
  dat <- dat[!duplicated(dat$snp), , drop=FALSE]
  if (nrow(dat) < 50) stop(paste(locus, "has too few shared SNPs:", nrow(dat)))

  obj <- read_ld(locus, dat)
  dat <- obj$dat
  LD <- obj$LD

  n1 <- as.integer(round(median(dat$n_pqtl, na.rm=TRUE)))
  n2 <- as.integer(round(median(dat$n_outcome, na.rm=TRUE)))
  outcome_type <- resolve_outcome_type(dat)

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
    type=outcome_type
  )

  if (outcome_type == "quant") {
    d2$sdY <- 1
  } else {
    if (!"case_fraction" %in% names(dat)) {
      stop(paste(locus, "case-control outcome requires case_fraction column"))
    }
    s <- median(dat$case_fraction, na.rm=TRUE)
    if (!is.finite(s) || s <= 0 || s >= 1) stop(paste(locus, "invalid case_fraction"))
    d2$s <- s
  }

  cat("MASTER SuSiE:", locus, "ancestry=", ancestry, "outcome_type=", outcome_type, "n=", nrow(dat), "\n")

  p1_rds <- file.path(output_dir, paste0(locus,"_pqtl_susie.rds"))
  p2_rds <- file.path(output_dir, paste0(locus,"_outcome_susie.rds"))

  if (file.exists(p1_rds)) {
    S1 <- readRDS(p1_rds)
    s1 <- list(ok=isTRUE(S1$converged), fit=S1,
               error=if (isTRUE(S1$converged)) "" else "cached pQTL fit not converged")
  } else {
    s1 <- safe_susie(d1, paste0(locus,"_pqtl"))
    if (!is.null(s1$fit)) saveRDS(s1$fit,p1_rds)
  }

  if (file.exists(p2_rds)) {
    S2 <- readRDS(p2_rds)
    s2 <- list(ok=isTRUE(S2$converged), fit=S2,
               error=if (isTRUE(S2$converged)) "" else "cached outcome fit not converged")
  } else {
    s2 <- safe_susie(d2, paste0(locus,"_outcome"))
    if (!is.null(s2$fit)) saveRDS(s2$fit,p2_rds)
  }

  a <- abf[abf[[abf_key]] == locus, , drop=FALSE]
  abf_h4 <- if (nrow(a) && "PP.H4" %in% names(a)) as.numeric(a$PP.H4[1]) else NA_real_
  abf_h3 <- if (nrow(a) && "PP.H3" %in% names(a)) as.numeric(a$PP.H3[1]) else NA_real_

  if (!s1$ok || !s2$ok) {
    failures[[length(failures)+1]] <- data.frame(
      locus=locus, ancestry=ancestry, outcome_type=outcome_type, nsnps=nrow(dat),
      pqtl_status=if (s1$ok) "converged" else "failed_or_nonconverged",
      outcome_status=if (s2$ok) "converged" else "failed_or_nonconverged",
      pqtl_message=s1$error, outcome_message=s2$error,
      stringsAsFactors=FALSE
    )
    locus_rows[[length(locus_rows)+1]] <- data.frame(
      locus=locus, ancestry=ancestry, outcome_type=outcome_type, nsnps=nrow(dat),
      pqtl_credible_sets=if (s1$ok) cs_count(s1$fit) else NA_integer_,
      outcome_credible_sets=if (s2$ok) cs_count(s2$fit) else NA_integer_,
      signal_pairs_default=0L, max_PP.H4=NA_real_, PP.H3_at_best_pair=NA_real_,
      best_pqtl_signal="", best_outcome_signal="",
      any_H4_ge_0.8=0L, any_H4_ge_0.9=0L,
      abf_PP.H3=abf_h3, abf_PP.H4=abf_h4,
      comparison_status="SuSiE_nonconverged_or_failed",
      N_pqtl=n1, N_outcome=n2,
      stringsAsFactors=FALSE
    )
    write_checkpoint()
    rm(LD,d1,d2); gc()
    next
  }

  S1 <- s1$fit
  S2 <- s2$fit
  priors <- sort(unique(c(1e-6, default_p12, 1e-4)))

  for (p12 in priors) {
    res <- coloc.susie(S1,S2,p1=1e-4,p2=1e-4,p12=p12)
    sm <- as.data.frame(res$summary, stringsAsFactors=FALSE)
    if (!nrow(sm)) {
      all_rows[[length(all_rows)+1]] <- data.frame(
        locus=locus, ancestry=ancestry, outcome_type=outcome_type,
        p1=1e-4,p2=1e-4,p12=p12,nsnps=nrow(dat),
        hit1="",hit2="",PP.H0=NA,PP.H1=NA,PP.H2=NA,PP.H3=NA,PP.H4=NA,
        idx1=NA,idx2=NA, stringsAsFactors=FALSE
      )
    } else {
      for (i in seq_len(nrow(sm))) {
        all_rows[[length(all_rows)+1]] <- data.frame(
          locus=locus, ancestry=ancestry, outcome_type=outcome_type,
          p1=1e-4,p2=1e-4,p12=p12,nsnps=as.integer(sm$nsnps[i]),
          hit1=as.character(sm$hit1[i]),hit2=as.character(sm$hit2[i]),
          PP.H0=as.numeric(sm$PP.H0.abf[i]),PP.H1=as.numeric(sm$PP.H1.abf[i]),
          PP.H2=as.numeric(sm$PP.H2.abf[i]),PP.H3=as.numeric(sm$PP.H3.abf[i]),
          PP.H4=as.numeric(sm$PP.H4.abf[i]),
          idx1=as.integer(sm$idx1[i]),idx2=as.integer(sm$idx2[i]),
          stringsAsFactors=FALSE
        )
      }
    }
  }

  default_rows <- Filter(function(x) {
    x$locus[1] == locus && abs(x$p12[1]-default_p12) < .Machine$double.eps^0.5
  }, all_rows)
  dr <- if (length(default_rows)) do.call(rbind,default_rows) else data.frame()

  if (nrow(dr) && any(is.finite(dr$PP.H4))) {
    best <- dr[which.max(ifelse(is.finite(dr$PP.H4),dr$PP.H4,-Inf)),,drop=FALSE]
    max_h4 <- best$PP.H4
    best_h3 <- best$PP.H3
    best_hit1 <- best$hit1
    best_hit2 <- best$hit2
    pair_n <- nrow(dr)
  } else {
    max_h4 <- NA_real_; best_h3 <- NA_real_
    best_hit1 <- ""; best_hit2 <- ""; pair_n <- 0L
  }

  status <- "unresolved"
  if (is.finite(max_h4) && max_h4 >= 0.8 && is.finite(abf_h4) && abf_h4 >= 0.8) {
    status <- "ABF_and_SuSiE_shared_signal"
  } else if (is.finite(max_h4) && max_h4 >= 0.8 && is.finite(abf_h4) && abf_h4 < 0.2) {
    status <- "SuSiE_multisignal_rescue"
  } else if (is.finite(max_h4) && max_h4 < 0.8 && is.finite(abf_h4) && abf_h4 >= 0.8) {
    status <- "ABF_not_confirmed_by_SuSiE"
  }

  locus_rows[[length(locus_rows)+1]] <- data.frame(
    locus=locus, ancestry=ancestry, outcome_type=outcome_type, nsnps=nrow(dat),
    pqtl_credible_sets=cs_count(S1), outcome_credible_sets=cs_count(S2),
    signal_pairs_default=pair_n, max_PP.H4=max_h4, PP.H3_at_best_pair=best_h3,
    best_pqtl_signal=best_hit1, best_outcome_signal=best_hit2,
    any_H4_ge_0.8=as.integer(is.finite(max_h4) && max_h4 >= 0.8),
    any_H4_ge_0.9=as.integer(is.finite(max_h4) && max_h4 >= 0.9),
    abf_PP.H3=abf_h3, abf_PP.H4=abf_h4,
    comparison_status=status, N_pqtl=n1, N_outcome=n2,
    stringsAsFactors=FALSE
  )

  write_checkpoint()
  rm(LD,S1,S2,d1,d2); gc()
}

write_checkpoint()
cat("MASTER_SUSIE_PASS ancestry=", ancestry, "\n", sep="")
