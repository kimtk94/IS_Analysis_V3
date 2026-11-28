#!/usr/bin/env Rscript
# Evaluate the actual success-summary field without loading edgeR or count data.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 2L)
program <- parse(file=args[[1]], keep.source=FALSE)
is_call <- function(x, name) is.call(x) && identical(x[[1]], as.name(name))
fits <- Filter(function(x) is_call(x, "tryCatch"), as.list(program))
stopifnot(length(fits) == 1L, is_call(fits[[1]][[2]], "{"))
statements <- as.list(fits[[1]][[2]])[-1L]
success <- Filter(function(x) {
  is_call(x, "<-") && identical(x[[2]], as.name("summary")) &&
    is_call(x[[3]], "list") &&
    identical(as.list(x[[3]])$status, "GO_EDGER_MODEL_FIT_QC")
}, statements)
stopifnot(length(success) == 1L)
field <- as.list(success[[1]][[3]])$classification_score_threshold
stopifnot(!is.null(field))
cases <- read.delim(args[[2]], check.names=FALSE)
stopifnot(nrow(cases) == 2L)
for (i in seq_len(nrow(cases))) {
  actual <- eval(field, list(score_cutoff=cases$score_cutoff[[i]]))
  expected <- cases$expected_classification_score_threshold[[i]]
  if (!identical(actual, expected)) {
    stop("P5 success metadata does not preserve the selected score cutoff")
  }
}
cat("P5_SOURCE_PARSE_AND_CUTOFF_METADATA_PASS\n")
