#!/usr/bin/env Rscript
# Independent source-input legacy coloc ABF replay: 3 predeclared gene/tissue
# pairs, 5 p12 priors, alternative QTL sample-size assumptions.
# NEVER changes historical Drive, GWAS data, or original 646-row ABF master.
suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})
if (as.character(packageVersion("coloc")) != "5.2.3") stop("Coloc 5.2.3 pinned")
argv <- commandArgs(trailingOnly=TRUE)
if (length(argv)!=3L) stop("Usage: Rscript script.R input_dir ORIGINAL_MASTER.tsv NEW_OUT_DIR")
root <- normalizePath(argv[[1]],mustWork=TRUE)
master <- normalizePath(argv[[2]],mustWork=TRUE)
out <- argv[[3]]
if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))>0) {
  stop("REFUSE_NONEMPTY_OUTPUT")
}
original <- fread(master)
if (nrow(original)!=646 || sum(original$status=="PASS")!=646) stop("OLD_ABF_MASTER_UNEXPECTED")
design <- data.table(
  gene=c("FGF5","SH3PXD2A","COL4A2"),
  gene_base=c("ENSG00000138675","ENSG00000107957","ENSG00000134871"),
  locus=c("BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L004"),
  tissue=c("GTEx_V8__Brain_Cerebellar_Hemisphere",
           "GTEx_V8__Artery_Tibial",
           "GTEx_V8__Brain_Putamen_basal_ganglia"),
  file=c(
    "BBJ_IS_L001__GTEx_V8__Brain_Cerebellar_Hemisphere__ENSG00000138675.tsv",
    "BBJ_IS_L002__GTEx_V8__Artery_Tibial__ENSG00000107957.tsv",
    "BBJ_IS_L004__GTEx_V8__Brain_Putamen_basal_ganglia__ENSG00000134871.tsv"
  )
)
p12grid <- c(1e-6,3e-6,1e-5,3e-5,1e-4)
N_total <- 174686
N_cases <- 22664
s <- N_cases/N_total
outrows <- list()
inputrows <- list()
for (i in seq_len(nrow(design))) {
  spec <- design[i]
  path <- file.path(root,spec$file)
  if (!file.exists(path) || file.info(path)$size < 1000) {
    stop(paste("Source input missing",path))
  }
  x <- fread(path)
  fields <- c("gene_base","locus","dataset_key","match_key","gwas_beta",
              "gwas_se","gwas_maf","gwas_p","eqtl_beta","eqtl_se",
              "eqtl_maf","eqtl_p","eqtl_n","qtl_n_scalar")
  if (!all(fields %in% names(x))) stop(paste("INPUT_COLS_MISSING",spec$gene))
  if (nrow(x)<50 || uniqueN(x$match_key)!=nrow(x)) {
    stop(paste("DUPLICATE_OR_TOO_FEW_SNP",spec$gene))
  }
  if (uniqueN(x$gene_base)!=1 || uniqueN(x$dataset_key)!=1 ||
      uniqueN(x$locus)!=1 || x$gene_base[1]!=spec$gene_base ||
      x$dataset_key[1]!=spec$tissue || x$locus[1]!=spec$locus) {
    stop(paste("SOURCE_METADATA_MISMATCH",spec$gene))
  }
  gwasN <- N_total
  qNfirst <- as.numeric(x$qtl_n_scalar[1])
  qNmin <- min(as.numeric(x$eqtl_n))
  qNmax <- max(as.numeric(x$eqtl_n))
  if (any(!is.finite(c(qNfirst,qNmin,qNmax))) ||
      qNfirst <=0 || qNmin<=0 || qNmax<qNmin ||
      qNfirst != as.numeric(x$eqtl_n[1])) {
    stop(paste("QTL_N_FAIL",spec$gene))
  }
  if (uniqueN(x$qtl_n_scalar)!=1) stop("Scalar QTL N not constant")
  if (any(!is.finite(as.matrix(x[,.(gwas_beta,gwas_se,gwas_maf,eqtl_beta,eqtl_se,eqtl_maf)]))) ||
      any(x$gwas_se<=0) || any(x$eqtl_se<=0) ||
      any(x$gwas_maf<=0) || any(x$eqtl_maf<=0) ||
      any(x$gwas_maf>.5) || any(x$eqtl_maf>.5)) {
    stop(paste("SNP_INPUT_QC_FAIL",spec$gene))
  }
  old <- original[
    locus==spec$locus & dataset_key==spec$tissue & gene_base==spec$gene_base
  ]
  if (nrow(old)!=1 || old$nsnps!=nrow(x) ||
      abs(old$qtl_n-qNfirst)>1e-9) stop(paste("HISTORICAL_DENOMINATOR_MISMATCH",spec$gene))
  inputrows[[i]] <- data.table(
    gene=spec$gene,locus=spec$locus,tissue=spec$tissue,
    n_snps=nrow(x),unique_snps=uniqueN(x$match_key),
    qtl_n_first=qNfirst,qtl_n_min=qNmin,qtl_n_max=qNmax,
    qtl_n_range=qNmax-qNmin,
    gwas_min_p=min(x$gwas_p),eqtl_min_p=min(x$eqtl_p),
    gwas_n_total=gwasN,gwas_n_cases=N_cases,
    historical_p12=1e-5,source_file=path
  )
  d1 <- list(
    beta=as.numeric(x$gwas_beta),varbeta=as.numeric(x$gwas_se)^2,
    snp=as.character(x$match_key),MAF=as.numeric(x$gwas_maf),
    N=N_total,s=s,type="cc"
  )
  for (qn_model in c("FIRST","MIN","MAX")) {
    qN <- switch(qn_model,FIRST=qNfirst,MIN=qNmin,MAX=qNmax)
    d2 <- list(
      beta=as.numeric(x$eqtl_beta),varbeta=as.numeric(x$eqtl_se)^2,
      snp=as.character(x$match_key),MAF=as.numeric(x$eqtl_maf),
      N=as.numeric(qN),type="quant"
    )
    for (p in p12grid) {
      fit <- suppressWarnings(coloc.abf(
        dataset1=d1,dataset2=d2,p1=1e-4,p2=1e-4,p12=p
      ))
      su <- fit$summary
      if (su[["nsnps"]]!=nrow(x)) stop("COL0C_SNP_DENOMINATOR_FAIL")
      mass <- as.numeric(su[paste0("PP.H",0:4,".abf")])
      if (any(!is.finite(mass)) ||
          any(mass < 0) || abs(sum(mass)-1)>1e-8) {
        stop("POSTERIOR_MASS_INVALID")
      }
      outrows[[length(outrows)+1]] <- data.table(
        gene=spec$gene,locus=spec$locus,tissue=spec$tissue,
        nsnps=nrow(x),qtl_n_mode=qn_model,qtl_n_used=qN,
        p1=1e-4,p2=1e-4,p12=p,
        PP_H0=mass[1],PP_H1=mass[2],PP_H2=mass[3],
        PP_H3=mass[4],PP_H4=mass[5],
        H4_over_H3=if(mass[4]>0) mass[5]/mass[4] else NA_real_,
        old_H4=old$PP.H4,
        old_H3=old$PP.H3,
        delta_vs_old_H4=mass[5]-old$PP.H4,
        scenario="INDEPENDENT_SNP_LEVEL_ABF_RERUN",
        stringsAsFactors=FALSE
      )
    }
  }
  cat("REPLAY_GENE",spec$gene,"nsnps",nrow(x),
      "min_QTL_N",qNmin,"first_QTL_N",qNfirst,"max_QTL_N",qNmax,"\n")
}
all <- rbindlist(outrows)
qc <- rbindlist(inputrows)
if (nrow(all)!=45L || nrow(qc)!=3L) stop("RERUN_GRID_INCOMPLETE")
baseline <- all[qtl_n_mode=="FIRST" & p12==1e-5]
if (nrow(baseline)!=3L || any(abs(baseline$delta_vs_old_H4)>1e-8)) {
  stop("ORIGINAL_FROZEN_ABF_REPLICATION_FAILED")
}
for (g in design$gene) {
  ol <- original[locus==design[gene==g]$locus &
                 dataset_key==design[gene==g]$tissue &
                 gene_base==design[gene==g]$gene_base]
  bz <- baseline[gene==g]
  oldvec <- as.numeric(ol[,paste0("PP.H",0:4),with=FALSE])
  newvec <- as.numeric(bz[,paste0("PP_H",0:4),with=FALSE])
  if (max(abs(oldvec-newvec)) > 1e-8) stop(paste("FIVE_HYPOTHESIS_REPLAY_FAIL",g))
}
dir.create(out,showWarnings=FALSE,recursive=TRUE)
fwrite(all,file.path(out,"IS_SNP_REPLAY_P12_QTL_N_GRID.tsv"),sep="\t")
fwrite(qc,file.path(out,"IS_SNP_REPLAY_INPUT_AUDIT.tsv"),sep="\t")
writeLines(c(
  "# IS Phase10/11 ABF original SNP-input independent replay",
  "- Three predeclared historical top gene–tissue pairs: FGF5, SH3PXD2A, COL4A2.",
  "- All three match historical original H0–H4 posterior within 1e-8 using independent coloc 5.2.3 and original per-SNP beta/SE/MAF.",
  "- Original study priors recovered from archived IS_Phase10_11_REPAIR_ONECELL_v3.ipynb code: p1=1e-4,p2=1e-4,p12=1e-5.",
  "- Original sample specifications: binary BBJ case-control N=174686, cases=22664; quantitative GTEx eQTL N = FIRST per-input qtl_n_scalar.",
  "- The original qtl N is not identical across SNPs; show MIN vs FIRST vs MAX scenarios as a sensitivity of scalar-N approximation (NOT correct SNP-specific-N analysis).",
  "- Full 646 assays remain unreplayed at SNP level; three targeted replays must not be generalized beyond input files.",
  "- Note the single-causal-variant assumption and GTEx/EAS ancestry/reference LD mismatch. No disease-specific cell eQTL or causal gene conclusion.",
  "- No FDR claims, no independent tissue replication, no novel instrument construction."
),file.path(out,"IS_SNP_REPLAY_README.md"))
cat("IS_SNP_ABF_INDEPENDENT_REPLAY=PASS\n")
cat("INPUT_FILES=",nrow(qc)," REPLAY_RUNS=",nrow(all)," P12_GRID=",length(p12grid),"\n",sep="")
print(all[p12==1e-5,.(gene,qtl_n_mode,qtl_n_used,PP_H3,PP_H4,delta_vs_old_H4)])
