#!/usr/bin/env Rscript
# Stream one delimited summary member, write a canonical table, then select cis instruments.

parse_args <- function(x) {
  out <- list()
  while (length(x) > 0) {
    if (!startsWith(x[[1]], "--") || length(x) < 2) stop("arguments must be --key value pairs")
    out[[substring(x[[1]], 3)]] <- x[[2]]
    x <- x[-c(1, 2)]
  }
  out
}
required <- c("archive", "gene", "ancestry", "coordinates", "standardized-dir",
              "instrument-dir", "legacy-dir", "limit-type", "limit")
a <- parse_args(commandArgs(trailingOnly=TRUE))
missing <- setdiff(required, names(a)); if (length(missing)) stop(paste("missing:", paste(missing, collapse=", ")))

dir.create(a[["standardized-dir"]], recursive=TRUE, showWarnings=FALSE)
dir.create(a[["instrument-dir"]], recursive=TRUE, showWarnings=FALSE)
dir.create(a[["legacy-dir"]], recursive=TRUE, showWarnings=FALSE)

# Archives may be a tar containing one summary member, or a plain fixture TSV.
archive <- a$archive
if (grepl("\\.tar(\\.gz|\\.tgz)?$", archive, ignore.case=TRUE)) {
  members <- utils::untar(archive, list=TRUE)
  members <- members[!grepl("/$", members)]
  if (!length(members)) stop("archive contains no regular member")
  command <- paste("tar -xOf", shQuote(archive), shQuote(members[[1]]))
  con <- pipe(command, open="r")
} else con <- file(archive, open="r")
on.exit(close(con), add=TRUE)

limit <- as.integer(a$limit)
if (a[["limit-type"]] == "lines") {
  lines <- readLines(con, n=limit, warn=FALSE)
} else {
  # readChar enforces an upper bound on the decompressed stream.
  text <- readChar(con, nchars=limit, useBytes=TRUE)
  lines <- strsplit(text, "\n", fixed=TRUE)[[1]]
  if (length(lines) && !endsWith(text, "\n")) lines <- lines[-length(lines)]
}
if (length(lines) < 2) stop("summary stream has no complete data rows")
input_text <- paste(lines, collapse="\n")
header <- lines[[1]]
if (grepl("\t", header, fixed=TRUE)) {
  tab <- read.delim(textConnection(input_text), check.names=FALSE,
                    stringsAsFactors=FALSE)
} else {
  tab <- read.table(textConnection(input_text), header=TRUE, sep="",
                    check.names=FALSE, stringsAsFactors=FALSE)
}

# UKB-PPP uses the BOLT-LMM summary-statistic convention: BETA is the effect
# per copy of ALLELE1, and A1FREQ is the frequency of that same allele.  Keep
# these source names explicit rather than relying on positional assumptions.
aliases <- list(chr=c("chr","chromosome","chrom"), pos=c("pos","position","bp","genpos"),
  effect_allele=c("effect_allele","ea","alt","allele1"),
  other_allele=c("other_allele","oa","ref","allele0"),
  beta=c("beta","effect"), se=c("se","stderr"), p_value=c("p_value","p","pval"),
  log10p=c("log10p"), eaf=c("eaf","effect_allele_frequency","a1freq"),
  rsid=c("rsid","snp","variant_id","id"))
lower <- tolower(names(tab))
has_column <- function(keys) any(lower %in% keys)
pick <- function(keys, required=TRUE) {
  hit <- which(lower %in% keys)
  if (!length(hit)) { if (required) stop(paste("missing column:", keys[[1]])); return(rep(NA, nrow(tab))) }
  tab[[hit[[1]]]]
}

# A direct p-value takes precedence.  With LOG10P, zero is valid (p=1), a
# finite positive value is transformed, and negative, infinite, or missing
# values become NA.  Very large finite values may underflow to zero, which is
# a valid probability and is retained.
direct_p <- has_column(aliases$p_value)
if (direct_p) {
  p_value <- suppressWarnings(as.numeric(pick(aliases$p_value)))
  log10p <- rep(NA_real_, nrow(tab))
} else {
  log10p <- suppressWarnings(as.numeric(pick(aliases$log10p)))
  p_value <- ifelse(is.finite(log10p) & log10p >= 0, 10^(-log10p), NA_real_)
}
pos_numeric <- suppressWarnings(as.numeric(pick(aliases$pos)))
eaf_present <- has_column(aliases$eaf)
canonical <- data.frame(gene=toupper(a$gene), ancestry=toupper(a$ancestry),
  chr=as.character(pick(aliases$chr)), pos=as.integer(pos_numeric), rsid=pick(aliases$rsid, FALSE),
  effect_allele=pick(aliases$effect_allele), other_allele=pick(aliases$other_allele),
  beta=suppressWarnings(as.numeric(pick(aliases$beta))),
  se=suppressWarnings(as.numeric(pick(aliases$se))), p_value=p_value,
  eaf=suppressWarnings(as.numeric(pick(aliases$eaf, FALSE))),
  stringsAsFactors=FALSE)
canonical$f_statistic <- (canonical$beta / canonical$se)^2

# Numeric QC is performed before instrument selection.  SE must be positive;
# positions must be positive whole numbers; frequencies and p-values are
# closed-interval probabilities.  Missing values are invalid when the source
# supplies that field (EAF remains optional for non-UKB inputs).
invalid_beta <- !is.finite(canonical$beta)
invalid_se <- !is.finite(canonical$se) | canonical$se <= 0
invalid_coordinate <- is.na(canonical$chr) | !nzchar(trimws(canonical$chr)) |
  !is.finite(pos_numeric) | pos_numeric <= 0 | pos_numeric != floor(pos_numeric)
invalid_eaf <- eaf_present & (!is.finite(canonical$eaf) | canonical$eaf < 0 | canonical$eaf > 1)
invalid_p <- !is.finite(canonical$p_value) | canonical$p_value < 0 | canonical$p_value > 1
invalid_f <- !is.finite(canonical$f_statistic) | canonical$f_statistic < 0
invalid_any <- invalid_beta | invalid_se | invalid_coordinate | invalid_eaf | invalid_p | invalid_f

batch <- "batch_001"
std_dir <- file.path(a[["standardized-dir"]], toupper(a$ancestry), batch)
dir.create(std_dir, recursive=TRUE, showWarnings=FALSE)
std_path <- file.path(std_dir, paste0(toupper(a$gene), ".tsv"))
qc <- data.frame(
  gene=toupper(a$gene), ancestry=toupper(a$ancestry), input_rows=nrow(canonical),
  invalid_beta=sum(invalid_beta), invalid_se=sum(invalid_se),
  invalid_coordinate=sum(invalid_coordinate), invalid_eaf=sum(invalid_eaf),
  invalid_p_value=sum(invalid_p), invalid_f_statistic=sum(invalid_f),
  invalid_any=sum(invalid_any), output_rows=sum(!invalid_any),
  p_value_source=ifelse(direct_p, "direct", "LOG10P"),
  log10p_zero=ifelse(direct_p, 0L, sum(is.finite(log10p) & log10p == 0)),
  log10p_negative=ifelse(direct_p, 0L, sum(is.finite(log10p) & log10p < 0)),
  log10p_nonfinite_or_missing=ifelse(direct_p, 0L, sum(!is.finite(log10p))),
  stringsAsFactors=FALSE)
canonical <- canonical[!invalid_any, , drop=FALSE]
write.table(canonical, std_path, sep="\t", quote=FALSE, row.names=FALSE)
write.table(qc, file.path(std_dir, paste0(toupper(a$gene), ".qc.tsv")),
            sep="\t", quote=FALSE, row.names=FALSE)

coord <- read.delim(a$coordinates, stringsAsFactors=FALSE)
coord <- coord[toupper(coord$gene) == toupper(a$gene), , drop=FALSE]
if (nrow(coord) != 1 || !all(c("chr","start","end","genome_build") %in% names(coord)))
  stop("coordinate table must contain one row with gene, chr, start, end, genome_build")
if (coord$genome_build[[1]] != "GRCh38") stop("only GRCh38 coordinates are accepted")
cis_window <- as.numeric(ifelse(is.null(a[["cis-window-bp"]]), 1e6, a[["cis-window-bp"]]))
p_threshold <- as.numeric(ifelse(is.null(a[["p-value-threshold"]]), 5e-8, a[["p-value-threshold"]]))
f_threshold <- as.numeric(ifelse(is.null(a[["f-statistic-threshold"]]), 10, a[["f-statistic-threshold"]]))
cis <- canonical$chr == coord$chr[[1]] & canonical$pos >= coord$start[[1]] - cis_window &
       canonical$pos <= coord$end[[1]] + cis_window
instruments <- canonical[cis & canonical$p_value <= p_threshold & canonical$f_statistic >= f_threshold, , drop=FALSE]
for (directory in c(a[["instrument-dir"]], file.path(a[["legacy-dir"]], toupper(a$ancestry)))) {
  dir.create(directory, recursive=TRUE, showWarnings=FALSE)
  write.table(instruments, file.path(directory, paste0(toupper(a$gene), ".tsv")),
              sep="\t", quote=FALSE, row.names=FALSE)
}
cat(sprintf("prepared %s %s: %d canonical, %d instruments; %s=%d\n",
            toupper(a$gene), toupper(a$ancestry), nrow(canonical), nrow(instruments),
            a[["limit-type"]], limit))
