#!/usr/bin/env Rscript
# Legacy coloc ABF posterior-prior reweighting ONLY, using coloc's own
# internal prior.adjust implementation (version-pinned 5.2.3).
# NOT a new SNP-level ABF analysis, SuSiE analysis, or independent QTL.
suppressPackageStartupMessages(library(coloc))

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) {
  stop("Usage: Rscript reweight_is_legacy_coloc_priors.R <ABF_MASTER.tsv> <NEW_EMPTY_OUTDIR>")
}
source <- normalizePath(args[[1]],mustWork=TRUE)
out <- args[[2]]
if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))>0) {
  stop("Output dir is not empty: refusing overwrite")
}
if (packageVersion("coloc") != "5.2.3") {
  stop("PINNED_VERSION_REQUIRED: coloc 5.2.3")
}

old_p1 <- 1e-4
old_p2 <- 1e-4
old_p12 <- 1e-5
prior_grid <- c(1e-6,3e-6,1e-5,3e-5,1e-4)
input <- read.delim(source,sep="\t",check.names=FALSE,
                    stringsAsFactors=FALSE,quote="",na.strings=c("","NA"))
required <- c("locus","dataset_key","gene_base","gene_symbol","status",
              "nsnps","qtl_n",paste0("PP.H",0:4))
if (!all(required %in% colnames(input))) {
  stop("Missing legacy master source columns")
}
if (nrow(input)!=646L || sum(input$status=="PASS")!=646L) {
  stop("Expected frozen original 646 PASS tests")
}
if (anyDuplicated(input[c("locus","dataset_key","gene_base")])) {
  stop("Duplicate gene/locus/dataset assay")
}
if (any(!is.finite(input$nsnps)) || any(input$nsnps<2L)) {
  stop("Invalid nsnps field")
}
posterior <- as.matrix(input[,paste0("PP.H",0:4)])
mode(posterior) <- "numeric"
if (any(!is.finite(posterior)) || any(posterior < 0) ||
    any(abs(rowSums(posterior)-1)>1e-7)) {
  stop("Invalid baseline H0-H4 posterior mass")
}
if (any(prior_grid < old_p1*old_p2) ||
    any(prior_grid > min(old_p1,old_p2))) {
  stop("Prior grid violates range assumed by coloc")
}

# Package function is internal and pinned by exact release.
# Retains model Bayes factors implicitly and rescales hypothesis priors.
prior_adjust <- get("prior.adjust",envir=asNamespace("coloc"))
results <- vector("list",nrow(input))
for (i in seq_len(nrow(input))) {
  row <- input[i,,drop=FALSE]
  summ <- c(as.numeric(row$nsnps),as.numeric(row[1,paste0("PP.H",0:4)]))
  names(summ) <- c("nsnps",paste0("PP.H",0:4,".abf"))
  res <- prior_adjust(summ,newp12=prior_grid,p1=old_p1,p2=old_p2,p12=old_p12)
  if (is.null(dim(res)) || nrow(res)!=length(prior_grid) ||
      any(!is.finite(res)) || any(res<0) ||
      any(abs(rowSums(res)-1)>1e-8)) {
    stop(paste("Invalid reweighted coloc hypothesis mass row",i))
  }
  b <- which(prior_grid==old_p12)
  if (length(b)!=1 || max(abs(res[b,]-unlist(row[1,paste0("PP.H",0:4)],use.names=FALSE)))>1e-9) {
    stop(paste("Baseline prior identity check failed row",i))
  }
  results[[i]] <- data.frame(
    source_row=i,
    locus=row$locus,gene_symbol=row$gene_symbol,
    gene_base=row$gene_base,dataset_key=row$dataset_key,
    nsnps=row$nsnps,qtl_n=row$qtl_n,
    old_p1=old_p1,old_p2=old_p2,old_p12=old_p12,
    conditional_p12=prior_grid,
    PP_H0=res[,1],PP_H1=res[,2],PP_H2=res[,3],
    PP_H3=res[,4],PP_H4=res[,5],
    H4_over_H3=ifelse(res[,4]>0,res[,5]/res[,4],NA_real_),
    H4_ge_0_8=res[,5]>=.8,
    H4_over_H3_ge_3=ifelse(res[,4]>0,res[,5]/res[,4]>=3,FALSE),
    model_note="SUMMARY_ONLY_CONDITIONAL_PRIOR_REWEIGHT",
    stringsAsFactors=FALSE
  )
}
all <- do.call(rbind,results)
if (nrow(all)!=646L*length(prior_grid)) stop("Output cardinality mismatch")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
write.table(all,file.path(out,"IS_LEGACY_646_ABF_P12_GRID.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE,na="NA")
focus <- c("FGF5","ALDH2","SH3PXD2A","COL4A2","CALHM2","COL4A1","NEURL1","C4orf22","INA")
focal <- all[all$gene_symbol %in% focus,,drop=FALSE]
if (nrow(focal)==0) stop("Focus genes not in master")
top <- list()
for (g in unique(focal$gene_symbol)) {
  baseline <- focal[focal$gene_symbol==g & focal$conditional_p12==old_p12,,drop=FALSE]
  pick <- baseline[which.max(baseline$PP_H4),,drop=FALSE]
  z <- focal[focal$gene_symbol==g & focal$dataset_key==pick$dataset_key &
             focal$locus==pick$locus,,drop=FALSE]
  top[[length(top)+1]] <- z
}
top <- do.call(rbind,top)
top <- top[order(match(top$gene_symbol,focus),top$conditional_p12),,drop=FALSE]
write.table(top,file.path(out,"IS_LEGACY_TOP_GENE_P12_GRID.tsv"),
            sep="\t",row.names=FALSE,quote=FALSE,na="NA")

# These summary-only reweighted diagnostics do NOT validate original
# source-prior provenance, balanced LD/QTL overlap or one-signal assumptions.
readme <- c(
  "# Historical 646 GTEx ABF tests — conditional prior sensitivity",
  "",
  "## Scope and claims",
  "- This is a **model-specific posterior-prior reweighting** using existing H0-H4 summaries, NOT an independently re-executed SNP-level coloc.",
  "- 646 original PASS gene-tissue rows on four legacy BBJ loci; expanded 80-window discovery remains UNTESTED by this input.",
  "- Assumed original priors: p1=1e-4, p2=1e-4, p12=1e-5. Original run metadata does NOT independently document these values; conditional analysis until verified.",
  "- P12 grid: 1e-6; 3e-6; 1e-5; 3e-5; 1e-4. The p1/p2 priors are held fixed.",
  "- Computed via coloc 5.2.3 internal prior.adjust and validated that baseline p12 reproduces original H0-H4 for every input row.",
  "- The best gene-tissue pair in the focus file was selected at the original presumed p12=1e-5; this is exploratory tissue selection across tests.",
  "- The 646 rows are NOT 646 independent associations and are NOT multiple-test adjusted.",
  "- Reweighted posterior robustness does NOT repair unknown ancestry mismatch, QTL alleles, missing SNPs, molecular LD, tissue selection, multi-signal structure, or independent replication.",
  "- SuSiE pair posterior is not to be reweighted as an ABF gene-tissue result.",
  "- Source model provenance and original variant-level inputs must be confirmed before manuscript-level inference."
)
writeLines(readme,file.path(out,"IS_LEGACY_ABF_P12_SENSITIVITY_README.md"))
cat("R3_8_PRIOR_REWEIGHT=PASS\n")
cat("SOURCE_ROWS=",nrow(input)," GRID_ROWS=",nrow(all),
    " GENE_ROWS=",nrow(top)," PACKAGE=",as.character(packageVersion("coloc")),"\n",sep="")
for (g in c("FGF5","SH3PXD2A","COL4A2","ALDH2","CALHM2")) {
  z <- top[top$gene_symbol==g,,drop=FALSE]
  if (nrow(z)) {
    cat("FOCAL",g,
        paste(paste0(format(z$conditional_p12,scientific=TRUE),":",sprintf("%.3f",z$PP_H4)),collapse=" "),"\n")
  }
}
