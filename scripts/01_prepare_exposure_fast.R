#!/usr/bin/env Rscript
# Stream the gene's chromosome summary, normalize UKB-PPP fields, and select cis instruments.

parse_args <- function(x) {
  out <- list()
  while (length(x) > 0) {
    if (!startsWith(x[[1]], "--") || length(x) < 2) stop("arguments must be --key value pairs")
    out[[substring(x[[1]], 3)]] <- x[[2]]
    x <- x[-c(1, 2)]
  }
  out
}

normalize_chr <- function(x) {
  x <- toupper(sub("^chr", "", trimws(as.character(x)), ignore.case=TRUE))
  x[!is.na(x) & x == "23"] <- "X"
  x[!is.na(x) & x == "24"] <- "Y"
  x[!is.na(x) & x == "M"] <- "MT"
  x
}

path_component <- function(value, label) {
  if (length(value) != 1 || is.na(value) ||
      !grepl("^[A-Za-z0-9][A-Za-z0-9._-]*$", value)) {
    stop(paste(label, "must be one filesystem-safe identifier"))
  }
  value
}

read_coordinates <- function(path, gene) {
  coord <- read.delim(path, stringsAsFactors=FALSE, check.names=FALSE)
  required <- c("gene", "chr", "start", "end", "genome_build")
  if (!all(required %in% names(coord))) stop("coordinate table must contain gene, chr, start, end, genome_build")
  coord <- coord[!is.na(coord$gene) & toupper(coord$gene) == gene, , drop=FALSE]
  if (nrow(coord) != 1) stop("coordinate table must contain exactly one row for the requested gene")
  if (is.na(coord$genome_build[[1]]) || coord$genome_build[[1]] != "GRCh38") {
    stop("only GRCh38 coordinates are accepted")
  }
  coord$chr <- normalize_chr(coord$chr)
  coord$start <- suppressWarnings(as.numeric(coord$start))
  coord$end <- suppressWarnings(as.numeric(coord$end))
  if (!grepl("^([1-9][0-9]*|X|Y|MT)$", coord$chr[[1]]) ||
      !is.finite(coord$start[[1]]) || !is.finite(coord$end[[1]]) ||
      coord$start[[1]] < 1 || coord$end[[1]] < coord$start[[1]]) {
    stop("gene coordinates must have a valid chromosome and ordered positive positions")
  }
  coord
}

is_chromosome_member <- function(member) {
  grepl("(^|[_.-])chr([0-9]+|X|Y|MT?)([_.-]|$)", basename(member), ignore.case=TRUE)
}

select_archive_member <- function(archive, chromosome) {
  members <- utils::untar(archive, list=TRUE)
  members <- members[!grepl("/$", members)]
  if (!length(members)) stop("archive contains no summary member")
  # Do not accept chr1 when chr10 is requested, or silently take the first assay.
  basenames <- basename(members)
  pattern <- paste0("(^|[_.-])chr", chromosome, "([_.-]|$)")
  matched <- members[grepl(pattern, basenames, ignore.case=TRUE)]
  if (length(matched) == 1) return(matched[[1]])
  if (length(matched) > 1) stop(paste("ambiguous archive: multiple members for chromosome", chromosome))
  chromosome_named <- is_chromosome_member(members)
  if (any(chromosome_named) || length(members) != 1) {
    stop(paste("archive must contain exactly one summary member for chromosome", chromosome))
  }
  # A single generically named TSV is supported; its chromosome is checked below.
  members[[1]]
}

read_summary_lines <- function(archive, chromosome, limit_type, limit) {
  is_tar <- grepl("(\\.tar(\\.(gz|bz2|xz))?|\\.tgz|\\.tbz2|\\.txz)$",
                  archive, ignore.case=TRUE)
  member <- NULL
  if (is_tar) {
    member <- select_archive_member(archive, chromosome)
    command <- paste("tar -xOf", shQuote(archive), "--", shQuote(member))
    con <- pipe(command, open="rb")
    if (grepl("\\.gz$", member, ignore.case=TRUE)) con <- gzcon(con)
  } else {
    # gzfile also reads uncompressed text and detects compression from its header.
    con <- gzfile(archive, open="rb")
  }
  on.exit(close(con), add=TRUE)
  if (limit_type == "lines") {
    lines <- readLines(con, n=limit, warn=FALSE)
  } else {
    bytes <- readBin(con, what="raw", n=limit)
    newlines <- which(bytes == as.raw(10))
    if (!length(newlines)) stop("summary stream has no complete data rows")
    text <- rawToChar(bytes[seq_len(max(newlines))])
    lines <- strsplit(text, "\n", fixed=TRUE)[[1]]
  }
  if (length(lines) < 2) stop("summary stream has no complete data rows")
  list(lines=lines, member=member)
}

parse_summary <- function(lines) {
  if (grepl("\t", lines[[1]], fixed=TRUE)) {
    tab <- read.delim(text=paste(lines, collapse="\n"), check.names=FALSE,
                      stringsAsFactors=FALSE, colClasses="character",
                      quote="", comment.char="", row.names=NULL)
  } else {
    header <- strsplit(trimws(lines[[1]]), "[[:space:]]+")[[1]]
    values <- lapply(lines[-1], function(line) {
      fields <- strsplit(trimws(line), "[[:space:]]+")[[1]]
      if (length(fields) > length(header) && tolower(tail(header, 1)) == "extra") {
        fields <- c(head(fields, length(header) - 1),
                    paste(fields[length(header):length(fields)], collapse=" "))
      }
      if (length(fields) != length(header)) stop("summary row does not match its header")
      fields
    })
    tab <- as.data.frame(do.call(rbind, values), stringsAsFactors=FALSE)
    names(tab) <- header
  }
  if (anyDuplicated(tolower(names(tab)))) stop("summary contains duplicate column names")
  tab
}

canonical_summary <- function(tab, gene, ancestry) {
  aliases <- list(chr=c("chr", "chromosome", "chrom"),
    pos=c("pos", "position", "bp", "genpos"),
    effect_allele=c("effect_allele", "ea", "alt", "allele1"),
    other_allele=c("other_allele", "oa", "ref", "allele0"),
    beta=c("beta", "effect"), se=c("se", "stderr"),
    p_value=c("p_value", "p", "pval"),
    eaf=c("eaf", "effect_allele_frequency", "a1freq"),
    rsid=c("rsid", "snp", "variant_id", "id"))
  lower <- tolower(names(tab))
  required <- c("chr", "pos", "effect_allele", "other_allele", "beta", "se")
  missing <- required[!vapply(aliases[required], function(keys) any(lower %in% keys), logical(1))]
  if (!any(lower %in% c(aliases$p_value, "log10p"))) missing <- c(missing, "p_value or LOG10P")
  if (length(missing)) {
    stop(paste("missing required UKB-PPP columns:", paste(missing, collapse=", "),
               "; discovered header:", paste(names(tab), collapse=" | ")))
  }
  pick <- function(keys) {
    hit <- match(keys, lower, nomatch=0)
    hit <- hit[hit != 0]
    if (!length(hit)) return(rep(NA_character_, nrow(tab)))
    tab[[hit[[1]]]]
  }
  number <- function(x) suppressWarnings(as.numeric(x))
  direct_p <- any(lower %in% aliases$p_value)
  logp <- rep(NA_real_, nrow(tab))
  if (direct_p) {
    # A supplied direct probability takes precedence over LOG10P, including
    # when it is invalid. Do not silently repair it from a second field.
    p <- number(pick(aliases$p_value))
  } else {
    logp <- number(pick("log10p"))
    p <- rep(NA_real_, nrow(tab))
    valid_logp <- is.finite(logp) & logp >= 0
    # Large finite values can underflow to p=0. Inf is an invalid source value,
    # not evidence for an infinitely strong association.
    p[valid_logp] <- 10^(-logp[valid_logp])
  }
  canonical <- data.frame(gene=gene, ancestry=ancestry,
    chr=normalize_chr(pick(aliases$chr)), pos=number(pick(aliases$pos)),
    rsid=pick(aliases$rsid), effect_allele=toupper(pick(aliases$effect_allele)),
    other_allele=toupper(pick(aliases$other_allele)), beta=number(pick(aliases$beta)),
    se=number(pick(aliases$se)), p_value=p, eaf=number(pick(aliases$eaf)),
    stringsAsFactors=FALSE)
  # Retain the original row denominator for QC while selecting ADD downstream.
  additive <- rep(TRUE, nrow(tab))
  if ("test" %in% lower) {
    additive <- toupper(pick("test")) == "ADD"
    additive <- !is.na(additive) & additive
  }
  canonical$f_statistic <- (canonical$beta / canonical$se)^2
  list(data=canonical, additive=additive, log10p=logp,
       p_value_source=if (direct_p) "direct" else "LOG10P",
       eaf_present=any(lower %in% aliases$eaf))
}

main <- function() {
  a <- parse_args(commandArgs(trailingOnly=TRUE))
  required <- c("archive", "gene", "ancestry", "coordinates", "standardized-dir",
                "instrument-dir", "legacy-dir", "limit-type", "limit")
  missing <- setdiff(required, names(a))
  if (length(missing)) stop(paste("missing:", paste(missing, collapse=", ")))
  gene <- path_component(toupper(a$gene), "gene")
  ancestry <- path_component(toupper(a$ancestry), "ancestry")
  batch <- path_component(if (is.null(a[["batch-id"]])) "batch_001" else a[["batch-id"]], "batch-id")
  source_id <- if (is.null(a[["source-id"]])) NULL else path_component(a[["source-id"]], "source-id")
  limit <- suppressWarnings(as.integer(a$limit))
  if (!a[["limit-type"]] %in% c("lines", "bytes") || is.na(limit) || limit < 2) {
    stop("limit-type must be lines or bytes, with a positive limit including header and data")
  }
  coord <- read_coordinates(a$coordinates, gene)
  stream <- read_summary_lines(a$archive, coord$chr[[1]], a[["limit-type"]], limit)
  prepared <- canonical_summary(parse_summary(stream$lines), gene, ancestry)
  all_rows <- prepared$data
  on_chromosome <- !is.na(all_rows$chr) & all_rows$chr == coord$chr[[1]]
  selected <- prepared$additive & on_chromosome
  if (!any(selected)) stop("summary sample contains no additive variants on the gene chromosome")
  if (!is.null(stream$member) && is_chromosome_member(stream$member) &&
      any(prepared$additive & !on_chromosome)) {
    stop("chromosome-specific archive member contains a different chromosome")
  }
  canonical <- all_rows[selected, , drop=FALSE]

  threshold <- function(name, default) {
    value <- if (is.null(a[[name]])) default else suppressWarnings(as.numeric(a[[name]]))
    if (length(value) != 1 || !is.finite(value) || value < 0) stop(paste("invalid", name))
    value
  }
  cis_window <- threshold("cis-window-bp", 1e6)
  p_threshold <- threshold("p-value-threshold", 5e-8)
  f_threshold <- threshold("f-statistic-threshold", 10)
  if (p_threshold > 1) stop("p-value-threshold must be between zero and one")
  invalid <- list(
    invalid_beta=!is.finite(canonical$beta),
    invalid_se=!is.finite(canonical$se) | canonical$se <= 0,
    invalid_coordinate=!is.finite(canonical$pos) | canonical$pos < 1 |
      canonical$pos != floor(canonical$pos),
    invalid_eaf=prepared$eaf_present & (!is.finite(canonical$eaf) |
      canonical$eaf < 0 | canonical$eaf > 1),
    invalid_p_value=!is.finite(canonical$p_value) | canonical$p_value < 0 |
      canonical$p_value > 1,
    invalid_f_statistic=!is.finite(canonical$f_statistic) | canonical$f_statistic < 0,
    invalid_alleles=is.na(canonical$effect_allele) | is.na(canonical$other_allele) |
      canonical$effect_allele == canonical$other_allele
  )
  invalid_any <- Reduce(`|`, invalid)
  if (any(invalid_any)) warning(sprintf("excluding %d invalid summary rows", sum(invalid_any)))
  canonical <- canonical[!invalid_any, , drop=FALSE]
  cis <- canonical$pos >= coord$start[[1]] - cis_window & canonical$pos <= coord$end[[1]] + cis_window
  instruments <- canonical[cis & canonical$p_value <= p_threshold &
                           canonical$f_statistic >= f_threshold, , drop=FALSE]

  if (is.null(source_id)) {
    # Existing one-source CLI calls retain the canonical path; instruments gain
    # ancestry isolation. Batch-runner calls always provide a source identifier.
    relative <- file.path(ancestry, batch, paste0(gene, ".tsv"))
    legacy_relative <- if (is.null(a[["batch-id"]])) file.path(ancestry, paste0(gene, ".tsv")) else relative
  } else {
    relative <- file.path(ancestry, batch, gene, paste0(source_id, ".tsv"))
    legacy_relative <- relative
  }
  paths <- c(standardized=file.path(a[["standardized-dir"]], relative),
             instrument=file.path(a[["instrument-dir"]], relative),
             legacy=file.path(a[["legacy-dir"]], legacy_relative))
  qc_path <- sub("\\.tsv$", ".qc.tsv", paths[["standardized"]])
  logp <- prepared$log10p[selected]
  using_logp <- prepared$p_value_source == "LOG10P"
  qc <- data.frame(
    gene=gene, ancestry=ancestry, batch_id=batch,
    source_id=if (is.null(source_id)) "" else source_id,
    source_path=normalizePath(a$archive, mustWork=TRUE),
    archive_member=if (is.null(stream$member)) "" else stream$member,
    status=if (nrow(canonical)) "PREPARED" else "NO_VALID_ROWS",
    recorded_at_utc=format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz="UTC"),
    limit_type=a[["limit-type"]], limit=limit,
    input_rows=nrow(all_rows), non_additive_rows=sum(!prepared$additive),
    off_gene_chromosome_rows=sum(prepared$additive & !on_chromosome),
    evaluated_rows=sum(selected),
    # Invalid reason counts overlap and refer only to additive rows on the gene
    # chromosome; invalid_any is their union, not the sum of those counts.
    as.list(vapply(invalid, sum, integer(1))),
    invalid_any=sum(invalid_any), output_rows=nrow(canonical),
    instrument_rows=nrow(instruments),
    eaf_supplied=prepared$eaf_present, p_value_source=prepared$p_value_source,
    log10p_zero=if (using_logp) sum(is.finite(logp) & logp == 0) else 0L,
    log10p_negative=if (using_logp) sum(is.finite(logp) & logp < 0) else 0L,
    log10p_nonfinite_or_missing=if (using_logp) sum(!is.finite(logp)) else 0L,
    log10p_underflow=if (using_logp) sum(is.finite(logp) & logp >= 0 & 10^(-logp) == 0) else 0L,
    stringsAsFactors=FALSE
  )
  if (!nrow(canonical)) {
    # A failed rerun must not pair new failure QC with older successful TSVs.
    # Keep the last successful quartet intact and diagnose this attempt apart.
    failed_qc_path <- sub("\\.qc\\.tsv$", ".failed.qc.tsv", qc_path)
    dir.create(dirname(failed_qc_path), recursive=TRUE, showWarnings=FALSE)
    write.table(qc, failed_qc_path, sep="\t", quote=FALSE, row.names=FALSE)
    stop(paste("summary sample contains no valid additive variants on the gene chromosome; failed-attempt QC:",
               failed_qc_path))
  }
  dir.create(dirname(qc_path), recursive=TRUE, showWarnings=FALSE)
  write.table(qc, qc_path, sep="\t", quote=FALSE, row.names=FALSE)
  for (name in names(paths)) {
    dir.create(dirname(paths[[name]]), recursive=TRUE, showWarnings=FALSE)
    table <- if (name == "standardized") canonical else instruments
    write.table(table, paths[[name]], sep="\t", quote=FALSE, row.names=FALSE)
  }
  cat(sprintf("prepared %s %s %s: %d canonical, %d instruments; %s=%d\n",
              gene, ancestry, batch, nrow(canonical), nrow(instruments), a[["limit-type"]], limit))
  cat(paste(names(paths), paths, sep="="), sep="\n")
  cat(paste0("qc=", qc_path), sep="\n")
}

main()
