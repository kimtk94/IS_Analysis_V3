#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3) {
  stop("Usage: 52_stage3e2_coloc_abf.R <input_dir> <out.tsv> <R_lib>")
}

input_dir <- args[1]
outfile <- args[2]
rlib <- args[3]

.libPaths(c(rlib, .libPaths()))
if (!requireNamespace("coloc", quietly=TRUE)) {
  stop("R package 'coloc' missing. Run 50_stage3e0_coloc_preflight.sh with INSTALL_COLOC_PACKAGES=1")
}

files <- list.files(input_dir, pattern="\\.coloc\\.tsv\\.gz$", full.names=TRUE)

res <- list()
ix <- 0

for (f in files) {
  x <- read.delim(gzfile(f), stringsAsFactors=FALSE, check.names=FALSE)
  if (nrow(x) < 50) next

  gene <- unique(x$gene)
  trait <- unique(x$trait)
  otype <- unique(x$outcome_type)

  # T2D/case-control requires an externally verified case fraction s.
  # Never infer it here.
  if (length(otype) != 1 || otype == "cc") {
    ix <- ix + 1
    res[[ix]] <- data.frame(
      gene=gene[1], trait=trait[1], prior_p12=NA,
      nsnps=nrow(x), PP.H0=NA,PP.H1=NA,PP.H2=NA,PP.H3=NA,PP.H4=NA,
      status="HOLD_CASE_CONTROL_S_REQUIRED",
      stringsAsFactors=FALSE
    )
    next
  }

  ok <- is.finite(x$beta_exp) & is.finite(x$se_exp) &
        is.finite(x$beta_out) & is.finite(x$se_out) &
        is.finite(x$exp_eaf)

  x <- x[ok,]
  if (nrow(x) < 50) next

  # Use exposure EAF as ancestry-matched MAF proxy for both datasets when
  # outcome EAF is absent. When outcome EAF is present, use it.
  maf1 <- pmin(x$exp_eaf, 1-x$exp_eaf)
  out_eaf <- suppressWarnings(as.numeric(x$out_eaf))
  maf2 <- ifelse(is.finite(out_eaf), pmin(out_eaf,1-out_eaf), maf1)

  n1 <- suppressWarnings(as.numeric(x$n_exp))
  n2 <- suppressWarnings(as.numeric(x$n_out))

  N1 <- round(median(n1[is.finite(n1)], na.rm=TRUE))
  N2 <- round(median(n2[is.finite(n2)], na.rm=TRUE))

  if (!is.finite(N1) || !is.finite(N2)) {
    ix <- ix + 1
    res[[ix]] <- data.frame(
      gene=gene[1],trait=trait[1],prior_p12=NA,
      nsnps=nrow(x),PP.H0=NA,PP.H1=NA,PP.H2=NA,PP.H3=NA,PP.H4=NA,
      status="HOLD_N_REQUIRED",
      stringsAsFactors=FALSE
    )
    next
  }

  snp <- ifelse(x$rsid!="" & x$rsid!=".", x$rsid, paste0(x$chrom,":",x$pos))

  d1 <- list(
    beta=x$beta_exp,
    varbeta=x$se_exp^2,
    snp=snp,
    position=x$pos,
    type="quant",
    N=N1,
    MAF=maf1
  )

  d2 <- list(
    beta=x$beta_out,
    varbeta=x$se_out^2,
    snp=snp,
    position=x$pos,
    type="quant",
    N=N2,
    MAF=maf2
  )

  for (p12 in c(1e-6,1e-5,1e-4)) {
    fit <- coloc::coloc.abf(
      dataset1=d1,
      dataset2=d2,
      p1=1e-4,p2=1e-4,p12=p12
    )
    s <- fit$summary

    ix <- ix + 1
    res[[ix]] <- data.frame(
      gene=gene[1], trait=trait[1], prior_p12=p12,
      nsnps=unname(s["nsnps"]),
      PP.H0=unname(s["PP.H0.abf"]),
      PP.H1=unname(s["PP.H1.abf"]),
      PP.H2=unname(s["PP.H2.abf"]),
      PP.H3=unname(s["PP.H3.abf"]),
      PP.H4=unname(s["PP.H4.abf"]),
      status="PASS",
      stringsAsFactors=FALSE
    )
  }
}

out <- if (length(res)) do.call(rbind,res) else data.frame()
write.table(out,file=outfile,sep="\t",row.names=FALSE,quote=FALSE,na="NA")
cat("coloc rows =", nrow(out), "\n")
