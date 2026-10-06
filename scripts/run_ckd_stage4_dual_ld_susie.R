#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop(
    "usage: run_ckd_stage4_dual_ld_susie.R <input_root> <ld_root> <abf_root> <output_dir>"
  )
}

input_root <- args[[1]]
ld_root <- args[[2]]
abf_root <- args[[3]]
output_dir <- args[[4]]

dir.create(
  output_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

suppressPackageStartupMessages(
  library(coloc)
)

suppressPackageStartupMessages(
  library(susieR)
)


read_ld <- function(gene, pop, dat) {

  meta_file <- file.path(
    ld_root,
    pop,
    paste0(gene, ".ld.meta.tsv")
  )

  bin_file <- file.path(
    ld_root,
    pop,
    paste0(gene, ".ld.bin")
  )

  meta <- read.delim(
    meta_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  n <- nrow(meta)

  con <- file(
    bin_file,
    "rb"
  )

  vals <- readBin(
    con,
    what="numeric",
    n=n*n,
    size=4,
    endian="little"
  )

  close(con)

  if (length(vals) != n*n) {
    stop(
      paste(
        gene,
        pop,
        "LD binary size mismatch"
      )
    )
  }

  LD <- matrix(
    vals,
    nrow=n,
    ncol=n,
    byrow=TRUE
  )

  LD <- (LD + t(LD)) / 2

  signs <- as.numeric(
    meta$effect_vs_ld_major_sign
  )

  if (any(!signs %in% c(-1, 1))) {
    stop(
      paste(
        gene,
        pop,
        "invalid LD orientation"
      )
    )
  }

  LD <- LD * tcrossprod(signs)

  diag(LD) <- 1

  idx <- match(
    dat$snp,
    meta$snp
  )

  if (anyNA(idx)) {
    stop(
      paste(
        gene,
        pop,
        "canonical SNP missing from LD"
      )
    )
  }

  LD <- LD[
    idx,
    idx,
    drop=FALSE
  ]

  meta <- meta[
    idx,
    ,
    drop=FALSE
  ]

  dimnames(LD) <- list(
    dat$snp,
    dat$snp
  )

  if (max(abs(LD - t(LD))) > 1e-5) {
    stop(
      paste(
        gene,
        pop,
        "LD not symmetric"
      )
    )
  }

  list(
    LD=LD,
    meta=meta
  )
}


cs_count <- function(x) {

  if (
    is.null(x$sets)
    ||
    is.null(x$sets$cs)
  ) {
    return(0L)
  }

  length(x$sets$cs)
}


safe_susie <- function(d, suffix) {

  tryCatch(
    {
      check_dataset(
        d,
        req="LD"
      )

      run_once <- function(maxit) {
        tryCatch(
          runsusie(
            d,
            suffix=suffix,
            maxit=maxit,
            repeat_until_convergence=FALSE,
            L=min(
              10L,
              max(
                1L,
                length(d$snp)-1L
              )
            )
          ),
          error=function(e) e
        )
      }

      fit <- run_once(1000L)

      if (
        inherits(fit, "error")
        && grepl(
          "did not converge",
          conditionMessage(fit),
          fixed=TRUE
        )
      ) {
        message(
          "SuSiE bounded rescue ",
          suffix,
          ": 1000 -> 2000"
        )
        fit <- run_once(2000L)
      }

      if (inherits(fit, "error")) {
        stop(
          paste0(
            "SuSiE bounded rescue failed: ",
            conditionMessage(fit)
          )
        )
      }

      list(
        ok=isTRUE(fit$converged),
        fit=fit,
        error=if (
          isTRUE(fit$converged)
        ) "" else "nonconverged"
      )
    },

    error=function(e) {
      list(
        ok=FALSE,
        fit=NULL,
        error=conditionMessage(e)
      )
    }
  )
}


read_abf <- function(phenotype) {

  path <- file.path(
    abf_root,
    phenotype,
    "coloc_results",
    "EAS_COLOC_DEFAULT.tsv"
  )

  read.delim(
    path,
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
}


all_rows <- list()
summary_rows <- list()
failures <- list()


checkpoint <- function() {

  if (length(all_rows)) {
    x <- do.call(
      rbind,
      all_rows
    )

    write.table(
      x,
      file=file.path(
        output_dir,
        "DUAL_LD_SUSIE_ALL_PRIORS.tsv"
      ),
      sep="\t",
      quote=FALSE,
      row.names=FALSE
    )
  }

  if (length(summary_rows)) {
    x <- do.call(
      rbind,
      summary_rows
    )

    write.table(
      x,
      file=file.path(
        output_dir,
        "DUAL_LD_SUSIE_DEFAULT.tsv"
      ),
      sep="\t",
      quote=FALSE,
      row.names=FALSE
    )
  }

  if (length(failures)) {
    x <- do.call(
      rbind,
      failures
    )

    write.table(
      x,
      file=file.path(
        output_dir,
        "DUAL_LD_SUSIE_FAILURES.tsv"
      ),
      sep="\t",
      quote=FALSE,
      row.names=FALSE
    )
  }
}


for (phenotype in c("eGFR", "BUN")) {

  input_dir <- file.path(
    input_root,
    phenotype
  )

  files <- sort(
    list.files(
      input_dir,
      pattern="\\.tsv\\.gz$",
      full.names=TRUE
    )
  )

  abf <- read_abf(
    phenotype
  )

  for (f in files) {

    gene <- sub(
      "\\.tsv\\.gz$",
      "",
      basename(f)
    )

    dat <- read.delim(
      gzfile(f),
      header=TRUE,
      sep="\t",
      stringsAsFactors=FALSE,
      check.names=FALSE
    )

    dat <- dat[
      is.finite(dat$beta_pqtl) &
      is.finite(dat$se_pqtl) &
      dat$se_pqtl > 0 &
      is.finite(dat$maf_pqtl) &
      dat$maf_pqtl > 0 &
      dat$maf_pqtl <= 0.5 &
      is.finite(dat$beta_outcome_aligned) &
      is.finite(dat$se_outcome) &
      dat$se_outcome > 0,
      ,
      drop=FALSE
    ]

    dat <- dat[
      !duplicated(dat$snp),
      ,
      drop=FALSE
    ]

    if (nrow(dat) < 100) {
      stop(
        paste(
          phenotype,
          gene,
          "too few SNPs"
        )
      )
    }

    EUR <- read_ld(
      gene,
      "EUR",
      dat
    )

    EAS <- read_ld(
      gene,
      "EAS",
      dat
    )

    eas_maf <- as.numeric(
      EAS$meta$reference_maf
    )

    keep <- which(
      is.finite(eas_maf) &
      eas_maf > 0 &
      eas_maf <= 0.5
    )

    dat <- dat[
      keep,
      ,
      drop=FALSE
    ]

    EUR$LD <- EUR$LD[
      keep,
      keep,
      drop=FALSE
    ]

    EAS$LD <- EAS$LD[
      keep,
      keep,
      drop=FALSE
    ]

    eas_maf <- eas_maf[
      keep
    ]

    n1 <- as.integer(
      round(
        median(
          dat$n_pqtl,
          na.rm=TRUE
        )
      )
    )

    n2 <- as.integer(
      round(
        median(
          dat$n_outcome,
          na.rm=TRUE
        )
      )
    )

    d1 <- list(
      beta=dat$beta_pqtl,
      varbeta=dat$se_pqtl^2,
      snp=as.character(dat$snp),
      position=as.integer(dat$pos37),
      MAF=dat$maf_pqtl,
      N=n1,
      sdY=1,
      LD=EUR$LD,
      type="quant"
    )

    d2 <- list(
      beta=dat$beta_outcome_aligned,
      varbeta=dat$se_outcome^2,
      snp=as.character(dat$snp),
      position=as.integer(dat$pos37),
      MAF=eas_maf,
      N=n2,
      sdY=1,
      LD=EAS$LD,
      type="quant"
    )

    cat(
      "\n============================================================\n"
    )

    cat(
      "Dual-LD SuSiE:",
      phenotype,
      gene,
      "n=",
      nrow(dat),
      "\n"
    )

    s1 <- safe_susie(
      d1,
      paste0(
        gene,
        "_pqtl_EUR"
      )
    )

    s2 <- safe_susie(
      d2,
      paste0(
        gene,
        "_",
        phenotype,
        "_EAS"
      )
    )

    if (!s1$ok || !s2$ok) {

      failures[[length(failures)+1]] <- data.frame(
        phenotype=phenotype,
        gene_symbol=gene,
        nsnps=nrow(dat),
        pqtl_status=if (s1$ok) "converged" else "failed",
        outcome_status=if (s2$ok) "converged" else "failed",
        pqtl_message=s1$error,
        outcome_message=s2$error,
        stringsAsFactors=FALSE
      )

      checkpoint()

      rm(EUR, EAS, d1, d2)
      gc()

      next
    }

    S1 <- s1$fit
    S2 <- s2$fit

    saveRDS(
      S1,
      file=file.path(
        output_dir,
        paste0(
          phenotype,
          "_",
          gene,
          "_pqtl_EUR.rds"
        )
      )
    )

    saveRDS(
      S2,
      file=file.path(
        output_dir,
        paste0(
          phenotype,
          "_",
          gene,
          "_outcome_EAS.rds"
        )
      )
    )

    default_rows <- list()

    for (p12 in c(1e-6, 1e-5, 1e-4)) {

      res <- tryCatch(
        coloc.susie(
          S1,
          S2,
          p1=1e-4,
          p2=1e-4,
          p12=p12
        ),
        error=function(e) NULL
      )

      if (
        is.null(res) ||
        is.null(res$summary)
      ) {
        next
      }

      sm <- as.data.frame(
        res$summary,
        stringsAsFactors=FALSE
      )

      if (!nrow(sm)) {
        next
      }

      for (i in seq_len(nrow(sm))) {

        row <- data.frame(
          phenotype=phenotype,
          gene_symbol=gene,
          p1=1e-4,
          p2=1e-4,
          p12=p12,
          nsnps=as.integer(sm$nsnps[i]),
          hit1=as.character(sm$hit1[i]),
          hit2=as.character(sm$hit2[i]),
          PP.H0=as.numeric(sm$PP.H0.abf[i]),
          PP.H1=as.numeric(sm$PP.H1.abf[i]),
          PP.H2=as.numeric(sm$PP.H2.abf[i]),
          PP.H3=as.numeric(sm$PP.H3.abf[i]),
          PP.H4=as.numeric(sm$PP.H4.abf[i]),
          stringsAsFactors=FALSE
        )

        all_rows[[length(all_rows)+1]] <- row

        if (p12 == 1e-5) {
          default_rows[[length(default_rows)+1]] <- row
        }
      }
    }

    if (length(default_rows)) {

      dr <- do.call(
        rbind,
        default_rows
      )

      best <- dr[
        which.max(
          ifelse(
            is.finite(dr$PP.H4),
            dr$PP.H4,
            -Inf
          )
        ),
        ,
        drop=FALSE
      ]

      max_h4 <- best$PP.H4
      best_h3 <- best$PP.H3
      hit1 <- best$hit1
      hit2 <- best$hit2
      pair_n <- nrow(dr)

    } else {

      max_h4 <- NA_real_
      best_h3 <- NA_real_
      hit1 <- ""
      hit2 <- ""
      pair_n <- 0L
    }

    a <- abf[
      abf$gene_symbol == gene,
      ,
      drop=FALSE
    ]

    abf_h3 <- if (nrow(a)) {
      as.numeric(a$PP.H3[1])
    } else {
      NA_real_
    }

    abf_h4 <- if (nrow(a)) {
      as.numeric(a$PP.H4[1])
    } else {
      NA_real_
    }

    status <- "unresolved"

    if (
      is.finite(max_h4) &&
      max_h4 >= 0.8
    ) {
      status <- "strong_shared_signal"

    } else if (
      is.finite(max_h4) &&
      max_h4 >= 0.5
    ) {
      status <- "suggestive_shared_signal"

    } else if (
      is.finite(best_h3) &&
      is.finite(max_h4) &&
      best_h3 > max_h4
    ) {
      status <- "distinct_signal_supported"
    }

    summary_rows[[length(summary_rows)+1]] <- data.frame(
      phenotype=phenotype,
      gene_symbol=gene,
      nsnps=nrow(dat),

      pqtl_credible_sets=
        cs_count(S1),

      outcome_credible_sets=
        cs_count(S2),

      signal_pairs_default=
        pair_n,

      max_PP.H4=
        max_h4,

      PP.H3_at_best_pair=
        best_h3,

      best_pqtl_signal=
        hit1,

      best_outcome_signal=
        hit2,

      any_H4_ge_0.5=
        as.integer(
          is.finite(max_h4) &&
          max_h4 >= 0.5
        ),

      any_H4_ge_0.8=
        as.integer(
          is.finite(max_h4) &&
          max_h4 >= 0.8
        ),

      ABF_PP.H3=
        abf_h3,

      ABF_PP.H4=
        abf_h4,

      comparison_status=
        status,

      N_pqtl=n1,
      N_outcome=n2,

      stringsAsFactors=FALSE
    )

    checkpoint()

    cat(
      "PASS",
      phenotype,
      gene,
      "pQTL_CS=",
      cs_count(S1),
      "outcome_CS=",
      cs_count(S2),
      "max_H4=",
      max_h4,
      "\n"
    )

    rm(
      EUR,
      EAS,
      S1,
      S2,
      d1,
      d2
    )

    gc()
  }
}


checkpoint()

if (length(failures)) {
  cat(
    "CKD_STAGE4_DUAL_LD_SUSIE_PARTIAL_PASS failures=",
    length(failures),
    "\n",
    sep=""
  )
} else {
  cat(
    "CKD_STAGE4_DUAL_LD_SUSIE_PASS\n"
  )
}

if (length(summary_rows)) {
  x <- do.call(
    rbind,
    summary_rows
  )

  print(
    x[
      order(
        x$phenotype,
        -x$max_PP.H4
      ),
    ]
  )
}
