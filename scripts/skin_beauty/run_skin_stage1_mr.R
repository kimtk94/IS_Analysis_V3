#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "/srv/is-analysis/IS_Analysis_V3"
exp_path <- if (length(args) >= 2) args[[2]] else Sys.getenv("SKIN_PQTL_EXPOSURE", unset = "")
outdir <- file.path(root, "results", "skin_beauty", "stage1_mr")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ids <- c(
  under_eye = "ebi-a-GCST90094903",
  crows_feet = "ebi-a-GCST90094904"
)

if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
if (!requireNamespace("ieugwasr", quietly = TRUE)) install.packages("ieugwasr", repos=c("https://mrcieu.r-universe.dev","https://cloud.r-project.org"))
if (!requireNamespace("TwoSampleMR", quietly = TRUE)) install.packages("TwoSampleMR", repos=c("https://mrcieu.r-universe.dev","https://cloud.r-project.org"))

jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (!nzchar(jwt)) {
  cat("ERROR: OPENGWAS_JWT is not set.\n")
  quit(status=2)
}

find_exposure <- function(root) {
  pats <- c("UKBPPP", "ST16", "cis", "pqtl")
  cmd <- sprintf("find %s -type f \\( -name '*.tsv' -o -name '*.tsv.gz' -o -name '*.txt' -o -name '*.txt.gz' \\) 2>/dev/null", shQuote(root))
  x <- tryCatch(system(cmd, intern=TRUE), error=function(e) character())
  if (!length(x)) return(NA_character_)
  score <- vapply(x, function(p) sum(vapply(pats, function(k) grepl(k, basename(p), ignore.case=TRUE), logical(1))), numeric(1))
  x <- x[order(score, decreasing=TRUE)]
  hit <- x[score[order(score, decreasing=TRUE)] >= 2]
  if (length(hit)) hit[[1]] else NA_character_
}

if (!nzchar(exp_path)) exp_path <- find_exposure(root)
if (is.na(exp_path) || !file.exists(exp_path)) {
  cat("ERROR: exposure cis-pQTL file not found.\n")
  cat("Set SKIN_PQTL_EXPOSURE=/path/to/cis_pqtl.tsv.gz or pass as arg2.\n")
  quit(status=3)
}

cat("EXPOSURE=", exp_path, "\n", sep="")
dt <- data.table::fread(exp_path, nThread=2, showProgress=FALSE)
orig <- names(dt)
norm <- tolower(gsub("[^a-z0-9]+", "_", orig))
lookup <- setNames(orig, norm)

pick <- function(cands, required=TRUE) {
  for (z in cands) if (z %in% names(lookup)) return(unname(lookup[[z]]))
  if (required) stop("Missing required column; candidates: ", paste(cands, collapse=","))
  NA_character_
}

csnp  <- pick(c("snp","rsid","rs_id","variant","variant_id"))
cbeta <- pick(c("beta","effect","b"))
cse   <- pick(c("se","standard_error","stderr"))
cea   <- pick(c("effect_allele","ea","a1","alt"))
coa   <- pick(c("other_allele","oa","a2","ref"))
ceaf  <- pick(c("effect_allele_frequency","eaf","a1freq","af","freq"), FALSE)
cp    <- pick(c("p","pval","p_value","pvalue"), FALSE)
cn    <- pick(c("n","samplesize","sample_size","n_total"), FALSE)
cprot <- pick(c("protein","protein_name","target","gene","gene_name","symbol","assay","aptamer"))

exp <- data.frame(
  SNP = as.character(dt[[csnp]]),
  beta.exposure = as.numeric(dt[[cbeta]]),
  se.exposure = as.numeric(dt[[cse]]),
  effect_allele.exposure = as.character(dt[[cea]]),
  other_allele.exposure = as.character(dt[[coa]]),
  eaf.exposure = if (!is.na(ceaf)) as.numeric(dt[[ceaf]]) else NA_real_,
  pval.exposure = if (!is.na(cp)) as.numeric(dt[[cp]]) else NA_real_,
  samplesize.exposure = if (!is.na(cn)) as.numeric(dt[[cn]]) else NA_real_,
  exposure = as.character(dt[[cprot]]),
  id.exposure = as.character(dt[[cprot]]),
  stringsAsFactors = FALSE
)
exp <- exp[!is.na(exp$SNP) & grepl("^rs[0-9]+$", exp$SNP) & is.finite(exp$beta.exposure) & is.finite(exp$se.exposure) & exp$se.exposure > 0, ]
exp$F_stat <- (exp$beta.exposure / exp$se.exposure)^2
exp <- exp[is.finite(exp$F_stat) & exp$F_stat > 10, ]
exp <- exp[!duplicated(paste(exp$exposure, exp$SNP)), ]

write.table(exp, file.path(outdir,"STAGE1_EXPOSURE_INSTRUMENTS.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("eligible instrument rows=", nrow(exp), " proteins=", length(unique(exp$exposure)), " SNPs=", length(unique(exp$SNP)), "\n", sep="")

snps <- unique(exp$SNP)
if (!length(snps)) quit(status=4)

query_one <- function(out_id, out_name) {
  cat("QUERY ", out_name, " ", out_id, " SNPs=", length(snps), "\n", sep="")
  # proxies=0 is deliberate: do not introduce EUR-proxy LD into an EAS outcome.
  chunks <- split(snps, ceiling(seq_along(snps) / 500))
  res <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    cat("  chunk ", i, "/", length(chunks), " n=", length(chunks[[i]]), "\n", sep="")
    x <- tryCatch(
      ieugwasr::associations(chunks[[i]], out_id, proxies=0, opengwas_jwt=jwt),
      error=function(e) {cat("  ERROR: ", conditionMessage(e), "\n", sep=""); NULL}
    )
    if (!is.null(x) && nrow(x)) res[[i]] <- x
  }
  y <- data.table::rbindlist(res, fill=TRUE)
  if (!nrow(y)) return(NULL)
  y$outcome_name <- out_name
  y$outcome_id <- out_id
  data.table::fwrite(y, file.path(outdir, paste0("STAGE1_OUTCOME_", toupper(out_name), ".tsv")), sep="\t")
  y
}

outs <- lapply(names(ids), function(nm) query_one(unname(ids[[nm]]), nm))
out <- data.table::rbindlist(outs, fill=TRUE)
if (!nrow(out)) {
  cat("ERROR: no outcome associations returned.\n")
  quit(status=5)
}

# Normalize current ieugwasr association output into TwoSampleMR outcome schema.
# OpenGWAS typically returns: rsid, beta, se, p, ea, nea, eaf, n, id, trait.
need <- function(df, cands) {
  z <- intersect(cands, names(df))
  if (!length(z)) return(rep(NA, nrow(df)))
  df[[z[[1]]]]
}
out2 <- data.frame(
  SNP = as.character(need(out, c("rsid","SNP"))),
  beta.outcome = as.numeric(need(out, c("beta","beta.outcome"))),
  se.outcome = as.numeric(need(out, c("se","se.outcome"))),
  effect_allele.outcome = as.character(need(out, c("ea","effect_allele","effect_allele.outcome"))),
  other_allele.outcome = as.character(need(out, c("nea","oa","other_allele","other_allele.outcome"))),
  eaf.outcome = as.numeric(need(out, c("eaf","eaf.outcome"))),
  pval.outcome = as.numeric(need(out, c("p","pval","pval.outcome"))),
  samplesize.outcome = as.numeric(need(out, c("n","samplesize","samplesize.outcome"))),
  outcome = as.character(out$outcome_name),
  id.outcome = as.character(out$outcome_id),
  stringsAsFactors = FALSE
)
out2 <- out2[!is.na(out2$SNP) & is.finite(out2$beta.outcome) & is.finite(out2$se.outcome), ]

harm <- TwoSampleMR::harmonise_data(exp, out2, action=2)
data.table::fwrite(harm, file.path(outdir,"STAGE1_HARMONISED.tsv"), sep="\t")

mr_input <- harm[harm$mr_keep %in% TRUE, ]
if (!nrow(mr_input)) {
  cat("ERROR: zero harmonised MR-eligible rows.\n")
  quit(status=6)
}

mrres <- TwoSampleMR::mr(mr_input)
data.table::fwrite(mrres, file.path(outdir,"STAGE1_MR_RESULTS.tsv"), sep="\t")

# Heterogeneity only where the selected method is estimable.
het <- tryCatch(TwoSampleMR::mr_heterogeneity(mr_input), error=function(e) NULL)
if (!is.null(het) && nrow(het)) data.table::fwrite(het, file.path(outdir,"STAGE1_HETEROGENEITY.tsv"), sep="\t")

pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(mr_input), error=function(e) NULL)
if (!is.null(pleio) && nrow(pleio)) data.table::fwrite(pleio, file.path(outdir,"STAGE1_EGGER_INTERCEPT.tsv"), sep="\t")

# Multiple-testing annotations across all reported MR rows.
mrres$fdr_bh <- p.adjust(mrres$pval, method="BH")
mrres$bonferroni <- p.adjust(mrres$pval, method="bonferroni")
data.table::fwrite(mrres, file.path(outdir,"STAGE1_MR_RESULTS_WITH_MULTIPLE_TESTING.tsv"), sep="\t")

cat("\n===== STAGE 1 COMPLETE =====\n")
cat("harmonised rows =", nrow(mr_input), "\n")
cat("MR rows =", nrow(mrres), "\n")
cat("FDR<0.05 rows =", sum(mrres$fdr_bh < 0.05, na.rm=TRUE), "\n")
cat("outdir =", outdir, "\n")
