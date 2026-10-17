#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 4) {
  stop("usage: run_master_coloc.R <input_dir> <output_dir> [outcome_type=quant|cc|auto] [default_p12]")
}
input_dir <- args[[1]]
output_dir <- args[[2]]
outcome_type_arg <- if (length(args) >= 3) args[[3]] else "auto"
default_p12 <- if (length(args) >= 4) as.numeric(args[[4]]) else 1e-5
if (!outcome_type_arg %in% c("quant","cc","auto")) stop("outcome_type must be quant, cc, or auto")
if (!is.finite(default_p12) || default_p12 <= 0 || default_p12 >= 1) stop("invalid default_p12")
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages(library(coloc))

files <- sort(list.files(input_dir, pattern="\\.tsv\\.gz$", full.names=TRUE))
if (!length(files)) stop("no coloc input files")

resolve_outcome_type <- function(dat) {
  if (outcome_type_arg != "auto") return(outcome_type_arg)
  if ("outcome_type" %in% names(dat)) {
    vals <- unique(tolower(as.character(dat$outcome_type)))
    vals <- vals[vals %in% c("quant","cc")]
    if (length(vals) == 1) return(vals[[1]])
  }
  "quant"
}

summary_rows <- list()
for (f in files) {
  locus <- sub("\\.tsv\\.gz$", "", basename(f))
  dat <- read.delim(gzfile(f), header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("snp","pos37","beta_pqtl","se_pqtl","maf_pqtl","n_pqtl",
                "beta_outcome_aligned","se_outcome","maf_outcome","n_outcome")
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
    type="quant"
  )

  d2 <- list(
    beta=dat$beta_outcome_aligned,
    varbeta=dat$se_outcome^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    MAF=dat$maf_outcome,
    N=n2,
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

  check_dataset(d1)
  check_dataset(d2)

  priors <- sort(unique(c(1e-6, default_p12, 1e-4)))
  for (p12 in priors) {
    res <- coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12)
    ss <- as.data.frame(as.list(res$summary), stringsAsFactors=FALSE)
    rr <- res$results
    top_i <- which.max(rr$SNP.PP.H4)
    row <- data.frame(
      locus=locus,
      outcome_type=outcome_type,
      p1=1e-4,p2=1e-4,p12=p12,
      nsnps=as.integer(ss$nsnps),
      PP.H0=as.numeric(ss$PP.H0.abf),
      PP.H1=as.numeric(ss$PP.H1.abf),
      PP.H2=as.numeric(ss$PP.H2.abf),
      PP.H3=as.numeric(ss$PP.H3.abf),
      PP.H4=as.numeric(ss$PP.H4.abf),
      H4_over_H3H4=as.numeric(ss$PP.H4.abf)/(as.numeric(ss$PP.H3.abf)+as.numeric(ss$PP.H4.abf)),
      top_snp=as.character(rr$snp[top_i]),
      top_snp_PP.H4=as.numeric(rr$SNP.PP.H4[top_i]),
      N_pqtl=n1,N_outcome=n2,
      stringsAsFactors=FALSE
    )
    summary_rows[[length(summary_rows)+1]] <- row
    if (isTRUE(all.equal(p12, default_p12))) {
      write.table(rr, file=file.path(output_dir,paste0(locus,"_coloc_snps.tsv")),
                  sep="\t",quote=FALSE,row.names=FALSE)
    }
  }
}

summary <- do.call(rbind,summary_rows)
summary <- summary[order(summary$p12, -summary$PP.H4),]
write.table(summary,file=file.path(output_dir,"MASTER_COLOC_ALL_PRIORS.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
default <- summary[abs(summary$p12-default_p12) < .Machine$double.eps^0.5,]
write.table(default,file=file.path(output_dir,"MASTER_COLOC_DEFAULT.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("MASTER_COLOC_PASS\n")
print(default[order(-default$PP.H4),])
