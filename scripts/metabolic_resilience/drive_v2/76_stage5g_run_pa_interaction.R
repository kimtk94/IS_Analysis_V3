#!/usr/bin/env Rscript

suppressWarnings({
  if (!requireNamespace("nlme", quietly=TRUE)) {
    stop("R package 'nlme' is required.")
  }
  if (!requireNamespace("jsonlite", quietly=TRUE)) {
    stop("R package 'jsonlite' is required.")
  }
})

root <- Sys.getenv("IS_ANALYSIS_ROOT", unset="/srv/is-analysis")
stage5 <- file.path(root, "results/metabolic_resilience/stage5_koges")
feasibility <- file.path(stage5, "STAGE5B_KOGES_FEASIBILITY_SUMMARY.json")
input <- file.path(stage5, "STAGE5G_PA_INTERACTION_ANALYSIS_READY.tsv.gz")
outdir <- file.path(stage5, "models")
outfile <- file.path(outdir, "STAGE5G_PA_INTERACTION.tsv")
summaryfile <- file.path(outdir, "STAGE5G_PA_INTERACTION_SUMMARY.json")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(feasibility)) {
  stop("Missing Stage5B feasibility summary.")
}

gate <- jsonlite::fromJSON(feasibility)
if (is.null(gate$physical_activity_usable_waves)) {
  stop("Stage5B feasibility summary lacks physical_activity_usable_waves.")
}
if (as.numeric(gate$physical_activity_usable_waves) < 3) {
  stop(
    "Repeated PA interaction is HOLD because physical activity is usable in <3 waves; use baseline PA stratification."
  )
}

if (!file.exists(input)) {
  stop(
    paste0(
      "Missing analysis-ready input: ", input, "\n",
      "Required columns: participant_id, gene, protein_grs, time_years, mbi, ",
      "physical_activity, age, sex, and optional PC1..PCn."
    )
  )
}

d <- read.delim(gzfile(input), stringsAsFactors=FALSE, check.names=FALSE)

required <- c(
  "participant_id","gene","protein_grs","time_years","mbi",
  "physical_activity","age","sex"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop(paste("Missing columns:", paste(missing, collapse=", ")))

pc_cols <- grep("^PC[0-9]+$", names(d), value=TRUE)
genes <- sort(unique(d$gene))
res <- list()

for (gene in genes) {
  x <- d[d$gene == gene, , drop=FALSE]
  vars <- c(
    "participant_id","protein_grs","time_years","mbi",
    "physical_activity","age","sex",pc_cols
  )
  x <- x[complete.cases(x[, vars, drop=FALSE]), , drop=FALSE]
  n_subjects <- length(unique(x$participant_id))

  if (n_subjects < 20 || nrow(x) < 40) {
    res[[length(res)+1]] <- data.frame(
      gene=gene, participants=n_subjects, observations=nrow(x),
      beta_grs_pa_time=NA_real_, se=NA_real_, p=NA_real_,
      status="WARN_LOW_N", stringsAsFactors=FALSE
    )
    next
  }

  rhs <- c(
    "time_years", "protein_grs", "physical_activity",
    "protein_grs:physical_activity",
    "protein_grs:time_years",
    "physical_activity:time_years",
    "protein_grs:physical_activity:time_years",
    "age", "sex", pc_cols
  )
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
      beta_grs_pa_time=NA_real_, se=NA_real_, p=NA_real_,
      status="FAIL_MODEL", stringsAsFactors=FALSE
    )
    next
  }

  tt <- summary(fit)$tTable
  candidates <- c(
    "time_years:protein_grs:physical_activity",
    "protein_grs:physical_activity:time_years",
    "protein_grs:time_years:physical_activity"
  )
  term <- candidates[candidates %in% rownames(tt)][1]

  if (is.na(term) || length(term) == 0) {
    res[[length(res)+1]] <- data.frame(
      gene=gene, participants=n_subjects, observations=nrow(x),
      beta_grs_pa_time=NA_real_, se=NA_real_, p=NA_real_,
      status="FAIL_THREE_WAY_TERM", stringsAsFactors=FALSE
    )
    next
  }

  res[[length(res)+1]] <- data.frame(
    gene=gene,
    participants=n_subjects,
    observations=nrow(x),
    beta_grs_pa_time=unname(tt[term,"Value"]),
    se=unname(tt[term,"Std.Error"]),
    p=unname(tt[term,"p-value"]),
    status="PASS",
    stringsAsFactors=FALSE
  )
}

out <- if (length(res)) do.call(rbind, res) else data.frame()
write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

payload <- list(
  status=ifelse(nrow(out) > 0, "PASS", "FAIL_ZERO"),
  genes=nrow(out),
  pass=ifelse(nrow(out) > 0, sum(out$status == "PASS"), 0),
  model="repeated MBI with ProteinGRS x physical activity x time",
  modifier="physical_activity"
)
writeLines(jsonlite::toJSON(payload, auto_unbox=TRUE, pretty=TRUE), summaryfile)
cat(jsonlite::toJSON(payload, auto_unbox=TRUE, pretty=TRUE), "\n")
