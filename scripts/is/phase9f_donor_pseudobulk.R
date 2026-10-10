# IS Phase9F-E R3_7 — donor-aware pseudobulk, descriptive reference evidence only.
# Expected call after loading Seurat RNA count layer. No package installation,
# no differential-expression P values and no cell-level pseudoreplication.
#
# Matrix columns MUST match annotation / patient / library_size exactly.
# Feature resolver is the existing R3_6 resolve_target_hits(layer_genes, symbol).

is_phase9f_donor_pseudobulk <- function(
  M, annotation, patient, library_size, genes, resolve_feature,
  min_cells=20L, min_paired_donors=4L
) {
  n <- ncol(M)
  if (length(annotation)!=n || length(patient)!=n || length(library_size)!=n) {
    stop("PSEUDOBULK_ALIGNMENT_FAIL: per-cell vectors do not align with matrix")
  }
  patient <- as.character(patient)
  annotation <- as.character(annotation)
  if (anyNA(patient) || any(!nzchar(patient)) || anyNA(annotation) ||
      any(!nzchar(annotation))) {
    stop("PSEUDOBULK_METADATA_FAIL: missing patient or author celltype labels")
  }
  if (any(!is.finite(library_size)) || any(library_size < 0)) {
    stop("PSEUDOBULK_LIBRARY_QC_FAIL: invalid library size")
  }
  if (length(unique(patient))<2L) stop("PSEUDOBULK_DONOR_QC_FAIL: <2 donors")
  min_cells <- as.integer(min_cells)
  min_paired_donors <- as.integer(min_paired_donors)
  if (min_cells < 1 || min_paired_donors < 2) stop("Invalid donor QC threshold")

  indices <- split(seq_len(n), interaction(patient,annotation,drop=TRUE,sep="|"))
  pairs <- unique(data.frame(
    patient=patient,celltype=annotation,stringsAsFactors=FALSE
  ))
  pairs <- pairs[order(pairs$patient,pairs$celltype),,drop=FALSE]
  # Indices from matched patient + author annotation avoid merging rare types.
  group_indices <- lapply(seq_len(nrow(pairs)),function(i) {
    which(patient==pairs$patient[i] & annotation==pairs$celltype[i])
  })
  group_n <- vapply(group_indices,length,integer(1))
  group_lib <- vapply(group_indices,function(ix) sum(library_size[ix]),numeric(1))
  base <- data.frame(
    donor=pairs$patient,celltype=pairs$celltype,
    n_cells=group_n,library_size=group_lib,stringsAsFactors=FALSE
  )
  base$coverage_status <- ifelse(base$n_cells>=min_cells & base$library_size>0,
                                 "ELIGIBLE","LOW_CELL_COUNT_OR_LIBRARY")
  rows <- vector("list",length(genes))
  all_feature_status <- vector("list",length(genes))
  for (j in seq_along(genes)) {
    g <- genes[j]
    lookup <- resolve_feature(rownames(M),g)
    ix <- lookup$hit
    absent <- length(ix)<1
    counts <- if (absent) NULL else as.numeric(Matrix::colSums(M[ix,,drop=FALSE]))
    r <- base
    r$gene <- g
    r$feature_status <- if (absent) "FEATURE_NOT_IN_MATRIX" else "FEATURE_PRESENT"
    r$target_sum_counts <- if (absent) NA_real_ else vapply(
      group_indices,function(q) sum(counts[q]),numeric(1))
    r$detected_cells <- if (absent) NA_integer_ else vapply(
      group_indices,function(q) as.integer(sum(counts[q]>0)),integer(1))
    r$detection_fraction <- if (absent) NA_real_ else r$detected_cells/r$n_cells
    r$pseudobulk_cpm <- if (absent) NA_real_ else
      ifelse(r$library_size>0,1e6*r$target_sum_counts/r$library_size,NA_real_)
    r$log2_cpm_plus1 <- log2(r$pseudobulk_cpm+1)
    r$coverage_status <- if (absent) "FEATURE_NOT_IN_MATRIX" else r$coverage_status
    rows[[j]] <- r[,c("gene","donor","celltype","feature_status",
                     "n_cells","library_size","target_sum_counts",
                     "detected_cells","detection_fraction","pseudobulk_cpm",
                     "log2_cpm_plus1","coverage_status")]
    all_feature_status[[j]] <- data.frame(
      gene=g,feature_status=if (absent) "FEATURE_NOT_IN_MATRIX" else "FEATURE_PRESENT",
      matched_features=if (absent) "" else paste(rownames(M)[ix],collapse=";"),
      stringsAsFactors=FALSE
    )
  }
  pb <- do.call(rbind,rows)
  rownames(pb) <- NULL
  fs <- do.call(rbind,all_feature_status)
  coverage <- aggregate(pb$coverage_status=="ELIGIBLE",
                        by=list(gene=pb$gene,celltype=pb$celltype),
                        FUN=sum)
  names(coverage)[3] <- "n_eligible_donors"
  all_n <- aggregate(pb$donor,by=list(gene=pb$gene,celltype=pb$celltype),
                     FUN=function(x) length(unique(x)))
  names(all_n)[3] <- "n_donors_with_cells"
  coverage <- merge(coverage,all_n,by=c("gene","celltype"),sort=TRUE)
  coverage$coverage_status <- ifelse(
    coverage$n_eligible_donors>=min_paired_donors,"DONOR_COVERAGE_PASS",
    ifelse(coverage$n_eligible_donors>=3,"LIMITED_DONOR_COVERAGE",
           "INSUFFICIENT_DONOR_COVERAGE")
  )
  absent <- coverage$gene %in% fs$gene[fs$feature_status=="FEATURE_NOT_IN_MATRIX"]
  coverage$coverage_status[absent] <- "FEATURE_NOT_IN_MATRIX"
  list(pseudobulk=pb,coverage=coverage,feature_status=fs,
       thresholds=list(min_cells=min_cells,min_paired_donors=min_paired_donors,
                       donors=sort(unique(patient)),n_cells=n))
}

is_phase9f_paired_comparisons <- function(result, comparisons,
                                          min_paired_donors=4L) {
  pb <- result$pseudobulk
  out <- vector("list",nrow(comparisons))
  for (i in seq_len(nrow(comparisons))) {
    spec <- comparisons[i,,drop=FALSE]
    a <- pb[pb$gene==spec$gene & pb$celltype==spec$celltype_a &
            pb$coverage_status=="ELIGIBLE",c("donor","log2_cpm_plus1")]
    b <- pb[pb$gene==spec$gene & pb$celltype==spec$celltype_b &
            pb$coverage_status=="ELIGIBLE",c("donor","log2_cpm_plus1")]
    joined <- merge(a,b,by="donor",suffixes=c("_a","_b"))
    eligible <- nrow(joined)
    delta <- if (eligible>0) joined$log2_cpm_plus1_a-joined$log2_cpm_plus1_b
    else numeric(0)
    out[[i]] <- data.frame(
      gene=spec$gene,celltype_a=spec$celltype_a,
      celltype_b=spec$celltype_b,
      n_paired_donors=eligible,
      paired_donors=if (eligible) paste(sort(joined$donor),collapse=";") else "",
      median_delta_log2_cpm_plus1=if (eligible) median(delta) else NA_real_,
      n_positive=if (eligible) sum(delta>0) else NA_integer_,
      n_negative=if (eligible) sum(delta<0) else NA_integer_,
      n_tied=if (eligible) sum(delta==0) else NA_integer_,
      positive_fraction=if (eligible) sum(delta>0)/eligible else NA_real_,
      qc_status=if (eligible>=min_paired_donors) "DESCRIPTIVE_PAIRED_QC_PASS"
                else "INSUFFICIENT_PAIRED_DONORS",
      inferential_p_value="NOT_ESTIMATED",
      stringsAsFactors=FALSE
    )
  }
  out <- do.call(rbind,out)
  rownames(out) <- NULL
  out
}

is_phase9f_predeclared_comparisons <- function() {
  data.frame(
    gene=c("SH3PXD2A","SH3PXD2A","SH3PXD2A","COL4A2",
           "COL4A2","COL4A2","COL4A1","COL4A1",
           "CALHM2","ALDH2"),
    celltype_a=c("Oligodendrocytes","Fibroblasts","Smooth muscle cells",
                 "Smooth muscle cells","Fibroblasts","Pericytes",
                 "Fibroblasts","Smooth muscle cells",
                 "Endothelial cells","Endothelial cells"),
    celltype_b=c("Microglia and Macrophages","Microglia and Macrophages",
                 "Microglia and Macrophages","Endothelial cells",
                 "Endothelial cells","Endothelial cells",
                 "Endothelial cells","Endothelial cells",
                 "Microglia and Macrophages","Microglia and Macrophages"),
    stringsAsFactors=FALSE
  )
}

is_phase9f_write_donor_outputs <- function(result,comparisons,out_dir) {
  if (!dir.exists(out_dir)) stop("Output directory missing")
  files <- list(
    HUMAN_DONOR_PSEUDOBULK=result$pseudobulk,
    HUMAN_DONOR_CELLTYPE_COVERAGE=result$coverage,
    HUMAN_DONOR_FEATURE_STATUS=result$feature_status,
    HUMAN_DONOR_PAIRED_COMPARISONS=comparisons
  )
  for (name in names(files)) {
    write.table(files[[name]],file.path(out_dir,paste0(name,".tsv")),
                sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  }
  pass <- sum(comparisons$qc_status=="DESCRIPTIVE_PAIRED_QC_PASS")
  report <- c(
    "# IS Phase9F-E human donor-aware pseudobulk (R3_7)",
    "",
    paste0("- Donors recorded: ",paste(result$thresholds$donors,collapse=";")),
    paste0("- RNA-assay cell count: ",result$thresholds$n_cells),
    paste0("- Min cells per donor/celltype: ",result$thresholds$min_cells),
    paste0("- Min paired donors for descriptive comparison: ",
           result$thresholds$min_paired_donors),
    paste0("- Predeclared comparisons with adequate coverage: ",
           pass, "/",nrow(comparisons)),
    "- Values are summed raw RNA counts / summed per-cell library counts (CPM).",
    "- Cell-level observations are not independent biological replicates.",
    "- No DE P values, cell-specific QTL, or genetic causal claims computed.",
    "- Normal adult temporal lobe reference is not an ischemic-stroke disease sample.",
    "- Microglia and Macrophages is the *original combined author class*.",
    "- FGF5 and C4orf22 absent features are NOT_ASSESSABLE; not expression zeros."
  )
  writeLines(report,file.path(out_dir,"HUMAN_DONOR_PSEUDOBULK_README.md"))
  invisible(files)
}
