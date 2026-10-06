#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "/srv/is-analysis/IS_Analysis_V3"
outdir <- file.path(root, "data", "skin_beauty", "stage0_gwas")
auditdir <- file.path(root, "results", "skin_beauty", "stage0_audit")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(auditdir, recursive = TRUE, showWarnings = FALSE)

ids <- c(
  under_eye = "ebi-a-GCST90094903",
  crows_feet = "ebi-a-GCST90094904"
)

if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  install.packages("ieugwasr", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))
}

jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (!nzchar(jwt)) {
  cat("ERROR: OPENGWAS_JWT is not set.\n")
  cat("Create a token at https://api.opengwas.io/profile/ and export it before running.\n")
  quit(status = 2)
}

cat("===== OPENGWAS METADATA =====\n")
meta_list <- lapply(unname(ids), function(id) {
  x <- ieugwasr::gwasinfo(id, opengwas_jwt = jwt)
  print(x)
  x
})
meta <- do.call(rbind, meta_list)
write.table(meta, file.path(auditdir, "STAGE0_OPENGWAS_METADATA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n===== OPENGWAS FILE LINKS =====\n")
file_list <- lapply(unname(ids), function(id) {
  x <- ieugwasr::gwasinfo_files(id, opengwas_jwt = jwt)
  if (is.null(x) || NROW(x) == 0) {
    cat("NO FILE LINKS:", id, "\n")
    return(NULL)
  }
  x$study_id <- id
  print(x)
  x
})
files <- do.call(rbind, file_list)
if (!is.null(files) && NROW(files) > 0) {
  write.table(files, file.path(auditdir, "STAGE0_OPENGWAS_FILES.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Minimal data-path validation without burning allowance on millions of variants.
# Top hits establish that each outcome is queryable and expose expected alleles/effects.
cat("\n===== TOP HITS VALIDATION =====\n")
top_list <- lapply(unname(ids), function(id) {
  x <- ieugwasr::tophits(id = id, clump = 1, opengwas_jwt = jwt)
  if (is.null(x) || NROW(x) == 0) {
    cat("NO TOP HITS:", id, "\n")
    return(NULL)
  }
  x$study_id <- id
  print(utils::head(x, 20))
  x
})
tops <- do.call(rbind, top_list)
if (!is.null(tops) && NROW(tops) > 0) {
  write.table(tops, file.path(auditdir, "STAGE0_WRINKLE_TOPHITS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# Save exact download instructions rather than automatically downloading ~5M-variant
# files. Stage 1 MR only needs outcome associations for the exposure instrument SNPs;
# full files become useful for regional coloc / fine-mapping.
notes <- c(
  "Primary OpenGWAS IDs:",
  paste(names(ids), ids, sep = " = "),
  "",
  "Stage 1 strategy: query only cis-pQTL instrument SNPs against both outcomes using ieugwasr::associations(..., proxies=0).",
  "Stage 2 strategy: query genomic ranges around prioritized loci for coloc; use ancestry-matched EAS LD locally.",
  "gwasinfo_files() returns expiring signed URLs if full VCF/indices are available."
)
writeLines(notes, file.path(auditdir, "STAGE0_NEXT_STEPS.txt"))

cat("\nDONE\n")
cat("metadata =", file.path(auditdir, "STAGE0_OPENGWAS_METADATA.tsv"), "\n")
cat("top hits =", file.path(auditdir, "STAGE0_WRINKLE_TOPHITS.tsv"), "\n")
