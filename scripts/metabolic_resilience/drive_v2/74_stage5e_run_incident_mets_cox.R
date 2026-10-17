#!/usr/bin/env Rscript

suppressWarnings({
  if (!requireNamespace("survival", quietly=TRUE)) {
    stop("R package 'survival' is required.")
  }
})

root <- Sys.getenv("IS_ANALYSIS_ROOT", unset="/srv/is-analysis")
input <- file.path(root, "results/metabolic_resilience/stage5_koges/STAGE5E_INCIDENT_METS_ANALYSIS_READY.tsv.gz")
outdir <- file.path(root, "results/metabolic_resilience/stage5_koges/models")
outfile <- file.path(outdir, "STAGE5E_INCIDENT_METS_COX.tsv")
summaryfile <- file.path(outdir, "STAGE5E_INCIDENT_METS_COX_SUMMARY.json")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(input)) {
  stop(
    paste0(
      "Missing analysis-ready input: ", input, "\n",
      "Required columns: participant_id, gene, protein_grs, baseline_healthy, ",
      "followup_years, incident_mets, age, sex, and optional PC1..PCn."
    )
  )
}

d <- read.delim(gzfile(input), stringsAsFactors=FALSE, check.names=FALSE)

required <- c(
  "participant_id", "gene", "protein_grs", "baseline_healthy",
  "followup_years", "incident_mets", "age", "sex"
)
missing <- setdiff(required, names(d))
if (length(missing)) {
  stop(paste("Missing columns:", paste(missing, collapse=", ")))
}

truth <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1","true","yes","y")
}

d <- d[truth(d$baseline_healthy), , drop=FALSE]
if (!nrow(d)) stop("No baseline-healthy participants in analysis-ready input.")

pc_cols <- grep("^PC[0-9]+$", names(d), value=TRUE)
genes <- sort(unique(d$gene))
res <- list()

for (gene in genes) {
  x <- d[d$gene == gene, , drop=FALSE]

  vars <- c(
    "protein_grs", "followup_years", "incident_mets",
    "age", "sex", pc_cols
  )
  x <- x[complete.cases(x[, vars, drop=FALSE]), , drop=FALSE]

  if (nrow(x) < 20 || length(unique(x$incident_mets)) < 2) {
    res[[length(res)+1]] <- data.frame(
      gene=gene,
      n=nrow(x),
      events=sum(x$incident_mets == 1, na.rm=TRUE),
      beta=NA_real_,
      se=NA_real_,
      HR=NA_real_,
      CI95_low=NA_real_,
      CI95_high=NA_real_,
      p=NA_real_,
      status="WARN_LOW_N_OR_NO_EVENT_VARIATION",
      stringsAsFactors=FALSE
    )
    next
  }

  rhs <- c("protein_grs", "age", "sex", pc_cols)
  form <- as.formula(
    paste0(
      "survival::Surv(followup_years, incident_mets) ~ ",
      paste(rhs, collapse=" + ")
    )
  )

  fit <- try(
    survival::coxph(form, data=x, ties="efron"),
    silent=TRUE
  )

  if (inherits(fit, "try-error")) {
    res[[length(res)+1]] <- data.frame(
      gene=gene,
      n=nrow(x),
      events=sum(x$incident_mets == 1, na.rm=TRUE),
      beta=NA_real_, se=NA_real_, HR=NA_real_,
      CI95_low=NA_real_, CI95_high=NA_real_, p=NA_real_,
      status="FAIL_MODEL",
      stringsAsFactors=FALSE
    )
    next
  }

  sm <- summary(fit)$coefficients
  if (!"protein_grs" %in% rownames(sm)) next

  beta <- unname(sm["protein_grs","coef"])
  se <- unname(sm["protein_grs","se(coef)"])
  p <- unname(sm["protein_grs","Pr(>|z|)"])

  res[[length(res)+1]] <- data.frame(
    gene=gene,
    n=nrow(x),
    events=sum(x$incident_mets == 1, na.rm=TRUE),
    beta=beta,
    se=se,
    HR=exp(beta),
    CI95_low=exp(beta - 1.96*se),
    CI95_high=exp(beta + 1.96*se),
    p=p,
    status="PASS",
    stringsAsFactors=FALSE
  )
}

out <- if (length(res)) do.call(rbind, res) else data.frame()
write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

json <- sprintf(
  '{"status":"%s","genes":%d,"pass":%d,"model":"Cox proportional hazards","primary_predictor":"protein_grs"}',
  ifelse(nrow(out) > 0, "PASS", "FAIL_ZERO"),
  nrow(out),
  ifelse(nrow(out) > 0, sum(out$status == "PASS"), 0)
)
writeLines(json, summaryfile)
cat(json, "\n")
