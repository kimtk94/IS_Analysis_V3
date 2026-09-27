#!/usr/bin/env Rscript

args <- commandArgs(
  trailingOnly=TRUE
)

if (length(args) != 3) {
  stop(
    paste(
      "usage:",
      "run_ckd_stage4_eas_coloc.R",
      "<input_dir>",
      "<phenotype>",
      "<output_dir>"
    )
  )
}

input_dir <- args[[1]]
phenotype <- args[[2]]
output_dir <- args[[3]]

dir.create(
  output_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

suppressPackageStartupMessages(
  library(coloc)
)

files <- sort(
  list.files(
    input_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

if (!length(files)) {
  stop(
    paste(
      "no input files:",
      input_dir
    )
  )
}

all_rows <- list()

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

  required <- c(
    "snp",
    "pos37",
    "beta_pqtl",
    "se_pqtl",
    "n_pqtl",
    "beta_outcome_aligned",
    "se_outcome",
    "n_outcome"
  )

  miss <- setdiff(
    required,
    names(dat)
  )

  if (length(miss)) {
    stop(
      paste(
        gene,
        "missing:",
        paste(
          miss,
          collapse=","
        )
      )
    )
  }

  dat <- dat[
    is.finite(dat$beta_pqtl) &
    is.finite(dat$se_pqtl) &
    dat$se_pqtl > 0 &
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

  if (nrow(dat) < 50) {
    stop(
      paste(
        gene,
        "too few shared SNPs:",
        nrow(dat)
      )
    )
  }

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

  # UKB-PPP protein GWAS and Chen eGFR/BUN
  # are analysed on normalized quantitative scales.
  # Supplying sdY=1 means MAF is not required
  # for beta/varbeta coloc.abf.
  d1 <- list(
    beta=dat$beta_pqtl,
    varbeta=dat$se_pqtl^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    N=n1,
    sdY=1,
    type="quant"
  )

  d2 <- list(
    beta=dat$beta_outcome_aligned,
    varbeta=dat$se_outcome^2,
    snp=as.character(dat$snp),
    position=as.integer(dat$pos37),
    N=n2,
    sdY=1,
    type="quant"
  )

  check_dataset(d1)
  check_dataset(d2)

  cat(
    "EAS coloc:",
    phenotype,
    gene,
    "n=",
    nrow(dat),
    "\n"
  )

  for (
    p12 in c(
      1e-6,
      1e-5,
      1e-4
    )
  ) {

    res <- coloc.abf(
      dataset1=d1,
      dataset2=d2,
      p1=1e-4,
      p2=1e-4,
      p12=p12
    )

    ss <- as.data.frame(
      as.list(
        res$summary
      ),
      stringsAsFactors=FALSE
    )

    rr <- res$results

    top_i <- which.max(
      rr$SNP.PP.H4
    )

    h3 <- as.numeric(
      ss$PP.H3.abf
    )

    h4 <- as.numeric(
      ss$PP.H4.abf
    )

    all_rows[[
      length(all_rows)+1
    ]] <- data.frame(
      gene_symbol=gene,
      phenotype=phenotype,

      p1=1e-4,
      p2=1e-4,
      p12=p12,

      nsnps=
        as.integer(
          ss$nsnps
        ),

      PP.H0=
        as.numeric(
          ss$PP.H0.abf
        ),

      PP.H1=
        as.numeric(
          ss$PP.H1.abf
        ),

      PP.H2=
        as.numeric(
          ss$PP.H2.abf
        ),

      PP.H3=h3,
      PP.H4=h4,

      H4_over_H3H4=
        ifelse(
          (h3+h4) > 0,
          h4/(h3+h4),
          NA
        ),

      top_snp=
        as.character(
          rr$snp[top_i]
        ),

      top_snp_PP.H4=
        as.numeric(
          rr$SNP.PP.H4[
            top_i
          ]
        ),

      N_pqtl=n1,
      N_outcome=n2,

      stringsAsFactors=FALSE
    )

    if (p12 == 1e-5) {

      write.table(
        rr,
        file=file.path(
          output_dir,
          paste0(
            gene,
            "_coloc_snps.tsv"
          )
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
      )
    }
  }
}

out <- do.call(
  rbind,
  all_rows
)

out <- out[
  order(
    out$gene_symbol,
    out$p12,
    -out$PP.H4
  ),
]

write.table(
  out,
  file=file.path(
    output_dir,
    "EAS_COLOC_ALL_PRIORS.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

default <- out[
  out$p12 == 1e-5,
  ,
  drop=FALSE
]

default <- default[
  order(
    -default$PP.H4
  ),
]

write.table(
  default,
  file=file.path(
    output_dir,
    "EAS_COLOC_DEFAULT.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

cat(
  "CKD_EAS_COLOC_PASS\n"
)

print(default)
