#!/usr/bin/env Rscript

suppressWarnings({
  if (!requireNamespace("nlme", quietly=TRUE)) {
    stop("R package 'nlme' is required.")
  }
})

root <- Sys.getenv("IS_ANALYSIS_ROOT", unset="/srv/is-analysis")
input <- file.path(root, "results/metabolic_resilience/stage5_koges/STAGE5F_REPEATED_MBI_ANALYSIS_READY.tsv.gz")
outdir <- file.path(root, "results/metabolic_resilience/stage5_koges/models")
outfile <- file.path(outdir, "STAGE5F_REPEATED_MBI_LMM.tsv")
summaryfile <- file.path(outdir, "STAGE5F_REPEATED_MBI_LMM_SUMMARY.json")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(input)) {
  stop(
    paste0(
      "Missing analysis-ready input: ", input, "\n",
      "Required columns: participant_id, gene, protein_grs, time_years, mbi, ",
      "age, sex, and optional PC1..PCn."
    )
  )
}

d <- read.delim(gzfile(input), stringsAsFactors=FALSE, check.names=FALSE)

required <- c(
  "participant_id", "gene", "protein_grs", "time_years",
  "mbi", "age", "sex"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop(paste("Missing columns:", paste(missing, collapse=", ")))

pc_cols <- grep("^PC[0-9]+$", names(d), value=TRUE)
genes <- sort(unique(d$gene))
res <- list()

for (gene in genes) {
  x <- d[d$gene == gene, , drop=FALSE]
  vars <- c("participant_id","protein_grs","time_years","mbi","age","sex",pc_cols)
  x <- x[complete.cases(x[, vars, drop=FALSE]), , drop=FALSE]

  n_subjects <- length(unique(x$participant_id))
  if (n_subjects < 20 || nrow(x) < 40) {
    res[[length(res)+1]] <- data.frame(
      gene=gene, participants=n_subjects, observations=nrow(x),
      beta_grs_time=NA_real_, se=NA_real_, p=NA_real_,
      status="WARN_LOW_N", stringsAsFactors=FALSE
    )
    next
  }

  rhs <- c("time_years", "protein_grs", "protein_grs:time_years", "age", "sex", pc_cols)
  form <- as.formula(paste("mbi ~", paste(rhs, collapse=" + ")))

  fit <- try(
    nlme::lme(
      fixed=form,
      random=~ time_years | participant_id,
      data=x,
      method="REML",
      na.action=na.omit,
      control=nlme::lmeControl(returnObject=TRUE)
    ),
    silent=TRUE
  )

  if (inherits(fit, "try-error")) {
    res[[length(res)+1]] <- data.frame(
      gene=gene, participants=n_subjects, observations=nrow(x),
      beta_grs_time=NA_real_, se=NA_real_, p=NA_real_,
      status="FAIL_MODEL", stringsAsFactors=FALSE
    )
    next
  }

  tt <- summary(fit)$tTable
  term <- "time_years:protein_grs"
  if (!term %in% rownames(tt)) term <- "protein_grs:time_years"

  if (!term %in% rownames(tt)) {
    res[[length(res)+1]] <- data.frame(
      gene=gene, participants=n_subjects, observations=nrow(x),
      beta_grs_time=NA_real_, se=NA_real_, p=NA_real_,
      status="FAIL_INTERACTION_TERM", stringsAsFactors=FALSE
    )
    next
  }

  res[[length(res)+1]] <- data.frame(
    gene=gene,
    participants=n_subjects,
    observations=nrow(x),
    beta_grs_time=unname(tt[term,"Value"]),
    se=unname(tt[term,"Std.Error"]),
    p=unname(tt[term,"p-value"]),
    status="PASS",
    stringsAsFactors=FALSE
  )
}

out <- if (length(res)) do.call(rbind, res) else data.frame()
write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

json <- sprintf(
  '{"status":"%s","genes":%d,"pass":%d,"model":"repeated MBI linear mixed model","target_term":"protein_grs x time_years"}',
  ifelse(nrow(out) > 0, "PASS", "FAIL_ZERO"),
  nrow(out),
  ifelse(nrow(out) > 0, sum(out$status == "PASS"), 0)
)
writeLines(json, summaryfile)
cat(json, "\n")
