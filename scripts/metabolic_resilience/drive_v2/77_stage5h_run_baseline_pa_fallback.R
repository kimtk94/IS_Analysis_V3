#!/usr/bin/env Rscript

suppressWarnings({
  if (!requireNamespace("survival", quietly=TRUE)) {
    stop("R package 'survival' is required.")
  }
  if (!requireNamespace("jsonlite", quietly=TRUE)) {
    stop("R package 'jsonlite' is required.")
  }
})

root <- Sys.getenv("IS_ANALYSIS_ROOT", unset="/srv/is-analysis")
stage5 <- file.path(root, "results/metabolic_resilience/stage5_koges")
feasibility <- file.path(stage5, "STAGE5B_KOGES_FEASIBILITY_SUMMARY.json")
input <- file.path(stage5, "STAGE5H_BASELINE_PA_ANALYSIS_READY.tsv.gz")
outdir <- file.path(stage5, "models")
outfile <- file.path(outdir, "STAGE5H_BASELINE_PA_FALLBACK.tsv")
summaryfile <- file.path(outdir, "STAGE5H_BASELINE_PA_FALLBACK_SUMMARY.json")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

if (!file.exists(feasibility)) {
  stop("Missing Stage5B feasibility summary.")
}
gate <- jsonlite::fromJSON(feasibility)

if (is.null(gate$physical_activity_usable_waves)) {
  stop("Stage5B feasibility summary lacks physical_activity_usable_waves.")
}
if (as.numeric(gate$physical_activity_usable_waves) >= 3) {
  stop("Repeated PA is usable in >=3 waves; Stage5G is the prespecified primary PA interaction route.")
}

if (!file.exists(input)) {
  stop(
    paste0(
      "Missing analysis-ready fallback input: ", input, "\n",
      "Required columns: participant_id, gene, protein_grs, baseline_pa_group, ",
      "baseline_healthy, followup_years, incident_mets, age, sex, optional PC1..PCn."
    )
  )
}

d <- read.delim(gzfile(input), stringsAsFactors=FALSE, check.names=FALSE)
required <- c(
  "participant_id","gene","protein_grs","baseline_pa_group",
  "baseline_healthy","followup_years","incident_mets","age","sex"
)
missing <- setdiff(required, names(d))
if (length(missing)) stop(paste("Missing columns:", paste(missing, collapse=", ")))

truth <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1","true","yes","y")
}

d <- d[truth(d$baseline_healthy), , drop=FALSE]
d$baseline_pa_group <- as.factor(d$baseline_pa_group)

pc_cols <- grep("^PC[0-9]+$", names(d), value=TRUE)
genes <- sort(unique(d$gene))
res <- list()

for (gene in genes) {
  x <- d[d$gene == gene, , drop=FALSE]
  vars <- c(
    "protein_grs","baseline_pa_group","followup_years",
    "incident_mets","age","sex",pc_cols
  )
  x <- x[complete.cases(x[, vars, drop=FALSE]), , drop=FALSE]

  if (nrow(x) < 20 || length(unique(x$incident_mets)) < 2 || nlevels(droplevels(x$baseline_pa_group)) < 2) {
    res[[length(res)+1]] <- data.frame(
      gene=gene,n=nrow(x),term="",beta=NA_real_,se=NA_real_,HR=NA_real_,p=NA_real_,
      status="WARN_LOW_N_OR_VARIATION",stringsAsFactors=FALSE
    )
    next
  }

  rhs <- c("protein_grs * baseline_pa_group","age","sex",pc_cols)
  form <- as.formula(
    paste0(
      "survival::Surv(followup_years, incident_mets) ~ ",
      paste(rhs, collapse=" + ")
    )
  )

  fit <- try(survival::coxph(form, data=x, ties="efron"), silent=TRUE)
  if (inherits(fit, "try-error")) {
    res[[length(res)+1]] <- data.frame(
      gene=gene,n=nrow(x),term="",beta=NA_real_,se=NA_real_,HR=NA_real_,p=NA_real_,
      status="FAIL_MODEL",stringsAsFactors=FALSE
    )
    next
  }

  sm <- summary(fit)$coefficients
  terms <- rownames(sm)[grepl("protein_grs:baseline_pa_group|baseline_pa_group.*:protein_grs", rownames(sm))]

  if (!length(terms)) {
    res[[length(res)+1]] <- data.frame(
      gene=gene,n=nrow(x),term="",beta=NA_real_,se=NA_real_,HR=NA_real_,p=NA_real_,
      status="FAIL_INTERACTION_TERM",stringsAsFactors=FALSE
    )
    next
  }

  for (term in terms) {
    beta <- unname(sm[term,"coef"])
    se <- unname(sm[term,"se(coef)"])
    p <- unname(sm[term,"Pr(>|z|)"])
    res[[length(res)+1]] <- data.frame(
      gene=gene,n=nrow(x),term=term,beta=beta,se=se,HR=exp(beta),p=p,
      status="PASS",stringsAsFactors=FALSE
    )
  }
}

out <- if (length(res)) do.call(rbind,res) else data.frame()
write.table(out,outfile,sep="\t",quote=FALSE,row.names=FALSE,na="NA")

payload <- list(
  status=ifelse(nrow(out)>0,"PASS","FAIL_ZERO"),
  rows=nrow(out),
  pass=ifelse(nrow(out)>0,sum(out$status=="PASS"),0),
  route="baseline PA fallback",
  model="incident MetS Cox with ProteinGRS x baseline_pa_group"
)
writeLines(jsonlite::toJSON(payload,auto_unbox=TRUE,pretty=TRUE),summaryfile)
cat(jsonlite::toJSON(payload,auto_unbox=TRUE,pretty=TRUE),"\n")
