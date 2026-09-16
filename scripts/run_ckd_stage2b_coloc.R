#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("usage: run_ckd_stage2b_coloc.R <input_dir> <output_dir>")
input_dir <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

suppressPackageStartupMessages(library(coloc))

files <- sort(list.files(input_dir, pattern="\\.tsv\\.gz$", full.names=TRUE))
if (!length(files)) stop("no coloc input files")

summary_rows <- list()
for (f in files) {
  gene <- sub("\\.tsv\\.gz$", "", basename(f))
  dat <- read.delim(gzfile(f), header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("snp","pos37","beta_pqtl","se_pqtl","maf_pqtl","n_pqtl",
                "beta_outcome_aligned","se_outcome","maf_outcome","n_outcome")
  miss <- setdiff(required, names(dat))
  if (length(miss)) stop(paste(gene, "missing", paste(miss, collapse=",")))
  dat <- dat[is.finite(dat$beta_pqtl) & is.finite(dat$se_pqtl) & dat$se_pqtl > 0 &
             is.finite(dat$beta_outcome_aligned) & is.finite(dat$se_outcome) & dat$se_outcome > 0 &
             is.finite(dat$maf_pqtl) & dat$maf_pqtl > 0 & dat$maf_pqtl <= 0.5 &
             is.finite(dat$maf_outcome) & dat$maf_outcome > 0 & dat$maf_outcome <= 0.5, ]
  dat <- dat[!duplicated(dat$snp), ]
  if (nrow(dat) < 50) stop(paste(gene, "has too few shared SNPs:", nrow(dat)))

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
    type="quant"
  )
  d2 <- list(
    beta=dat$beta_outcome_aligned,
    varbeta=dat$se_outcome^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    MAF=dat$maf_outcome,
    N=n2,
    type="quant"
  )
  check_dataset(d1)
  check_dataset(d2)

  for (p12 in c(1e-6,1e-5,1e-4)) {
    res <- coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=p12)
    ss <- as.data.frame(as.list(res$summary), stringsAsFactors=FALSE)
    rr <- res$results
    top_i <- which.max(rr$SNP.PP.H4)
    row <- data.frame(
      gene_symbol=gene,
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
      N_pqtl=n1,N_egfr=n2,
      stringsAsFactors=FALSE
    )
    summary_rows[[length(summary_rows)+1]] <- row
    if (p12 == 1e-5) {
      write.table(rr, file=file.path(output_dir,paste0(gene,"_coloc_snps.tsv")),
                  sep="\t",quote=FALSE,row.names=FALSE)
    }
  }
}
summary <- do.call(rbind,summary_rows)
summary <- summary[order(summary$p12, -summary$PP.H4),]
write.table(summary,file=file.path(output_dir,"STAGE2B_COLOC_SUMMARY.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
default <- summary[summary$p12==1e-5,]
write.table(default,file=file.path(output_dir,"STAGE2B_COLOC_DEFAULT.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("CKD_STAGE2B_COLOC_PASS\n")
print(default[order(-default$PP.H4),])
