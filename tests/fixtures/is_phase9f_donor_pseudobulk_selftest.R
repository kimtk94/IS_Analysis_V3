#!/usr/bin/env Rscript
# Small deterministic, synthetic sanity test; never use raw GWAS or donor data.
args <- commandArgs(trailingOnly=TRUE)
source(args[[1]])
if (!requireNamespace("Matrix",quietly=TRUE)) {
  stop("Matrix required for this fixture")
}
donors <- c("D1","D2","D3","D4")
types <- c("Oligodendrocytes","Microglia and Macrophages","Endothelial cells")
meta <- do.call(rbind,lapply(donors,function(d) do.call(rbind,lapply(types,function(ct) {
  n <- if (d=="D4" && ct=="Endothelial cells") 5 else 25
  data.frame(donor=rep(d,n),celltype=rep(ct,n),stringsAsFactors=FALSE)
}))))
nn <- nrow(meta)
mat <- matrix(0,nrow=3,ncol=nn,dimnames=list(c("SH3PXD2A","COL4A2","HOUSE"),paste0("CELL",seq_len(nn))))
mat["HOUSE",] <- 10
mat["SH3PXD2A",] <- ifelse(meta$celltype=="Oligodendrocytes",5,1)
mat["COL4A2",] <- ifelse(meta$celltype=="Endothelial cells",5,1)
M <- Matrix::Matrix(mat,sparse=TRUE)
resolve <- function(feature_names,g) list(hit=which(feature_names==g),method="ROWNAME_SYMBOL")
res <- is_phase9f_donor_pseudobulk(M,meta$celltype,meta$donor,
                                  as.numeric(Matrix::colSums(M)),
                                  c("SH3PXD2A","COL4A2","FGF5"),
                                  resolve,min_cells=20,min_paired_donors=4)
stopifnot(length(unique(res$pseudobulk$donor))==4L)
stopifnot(nrow(res$pseudobulk)==36L)
stopifnot(res$thresholds$n_cells==nn)
absent <- res$pseudobulk[res$pseudobulk$gene=="FGF5",]
stopifnot(all(absent$feature_status=="FEATURE_NOT_IN_MATRIX"))
stopifnot(all(is.na(absent$pseudobulk_cpm)))
ok <- res$pseudobulk[res$pseudobulk$gene=="SH3PXD2A" &
                     res$pseudobulk$celltype=="Oligodendrocytes",]
stopifnot(nrow(ok)==4L)
stopifnot(all(ok$n_cells==25L))
stopifnot(all(ok$detected_cells==25L))
stopifnot(all(ok$target_sum_counts==125))
low <- res$pseudobulk[res$pseudobulk$gene=="COL4A2" &
                      res$pseudobulk$donor=="D4" &
                      res$pseudobulk$celltype=="Endothelial cells",]
stopifnot(nrow(low)==1L)
stopifnot(low$n_cells==5L)
stopifnot(low$coverage_status=="LOW_CELL_COUNT_OR_LIBRARY")
comparisons <- data.frame(
  gene=c("SH3PXD2A","COL4A2","FGF5"),
  celltype_a=c("Oligodendrocytes","Endothelial cells","Oligodendrocytes"),
  celltype_b=rep("Microglia and Macrophages",3),
  stringsAsFactors=FALSE)
res_comp <- is_phase9f_paired_comparisons(res,comparisons,min_paired_donors=4)
stopifnot(res_comp$n_paired_donors[1]==4L)
stopifnot(res_comp$qc_status[1]=="DESCRIPTIVE_PAIRED_QC_PASS")
stopifnot(res_comp$n_positive[1]==4L)
stopifnot(res_comp$n_paired_donors[2]==3L)
stopifnot(res_comp$qc_status[2]=="INSUFFICIENT_PAIRED_DONORS")
stopifnot(res_comp$n_paired_donors[3]==0L)
stopifnot(is.na(res_comp$median_delta_log2_cpm_plus1[3]))
stopifnot(nrow(is_phase9f_predeclared_comparisons())==10L)
out <- tempfile(pattern="donor_pb_fixture_")
dir.create(out)
is_phase9f_write_donor_outputs(res,res_comp,out)
stopifnot(file.exists(file.path(out,"HUMAN_DONOR_PSEUDOBULK.tsv")))
stopifnot(file.exists(file.path(out,"HUMAN_DONOR_PAIRED_COMPARISONS.tsv")))
cat("IS_PHASE9F_DONOR_PSEUDOBULK_SYNTHETIC_TEST_PASS n_cells=",nn,"\n",sep="")
