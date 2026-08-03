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

aliases <- list(chr=c("chr","chromosome"), pos=c("pos","position","bp"),
  effect_allele=c("effect_allele","ea","alt"), other_allele=c("other_allele","oa","ref"),
  beta=c("beta","effect"), se=c("se","stderr"), p_value=c("p_value","p","pval"),
  eaf=c("eaf","effect_allele_frequency"), rsid=c("rsid","snp","variant_id"))
lower <- tolower(names(tab))
required_columns <- c("chr", "pos", "effect_allele", "other_allele", "beta", "se", "p_value")
missing_columns <- required_columns[!vapply(aliases[required_columns], function(keys) {
  any(lower %in% keys)
}, logical(1))]
if (ncol(tab) == 1 || length(missing_columns)) {
  details <- if (ncol(tab) == 1) "parsed only one column" else
    paste("missing required UKB-PPP columns:", paste(missing_columns, collapse=", "))
  stop(sprintf("failed to parse UKB-PPP summary '%s': %s; discovered header: %s",
               archive, details, paste(names(tab), collapse=" | ")))
}
pick <- function(keys, required=TRUE) {
  hit <- which(lower %in% keys)
  if (!length(hit)) { if (required) stop(paste("missing column:", keys[[1]])); return(rep(NA, nrow(tab))) }
  tab[[hit[[1]]]]
}
canonical <- data.frame(gene=toupper(a$gene), ancestry=toupper(a$ancestry),
  chr=pick(aliases$chr), pos=as.integer(pick(aliases$pos)), rsid=pick(aliases$rsid, FALSE),
  effect_allele=pick(aliases$effect_allele), other_allele=pick(aliases$other_allele),
  beta=as.numeric(pick(aliases$beta)), se=as.numeric(pick(aliases$se)),
  p_value=as.numeric(pick(aliases$p_value)), eaf=as.numeric(pick(aliases$eaf, FALSE)),
  stringsAsFactors=FALSE)
canonical$f_statistic <- (canonical$beta / canonical$se)^2
batch <- "batch_001"
std_dir <- file.path(a[["standardized-dir"]], toupper(a$ancestry), batch)
dir.create(std_dir, recursive=TRUE, showWarnings=FALSE)
std_path <- file.path(std_dir, paste0(toupper(a$gene), ".tsv"))
write.table(canonical, std_path, sep="\t", quote=FALSE, row.names=FALSE)

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
