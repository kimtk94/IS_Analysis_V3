#!/usr/bin/env Rscript
# Diagnostic: EAS GWAS vs GTEx QTL MAF agreement and coloc robustness.
# A minor-allele-frequency filter IS NOT harmonization, and ancestry MAF
# differences are not inherently errors. No SNP data is modified.
suppressPackageStartupMessages({library(data.table);library(coloc)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage Rscript script.R ORIGINAL_INPUT_DIR NEW_OUTPUT_DIR")
inp<-normalizePath(args[[1]],mustWork=TRUE)
out<-args[[2]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))>0)stop("REFUSE_NONEMPTY_OUTPUT")
config<-data.table(
 gene=c("FGF5","SH3PXD2A","COL4A2","CALHM2","NEURL1","C4orf22","INA","COL4A1"),
 file=c("BBJ_IS_L001__GTEx_V8__Brain_Cerebellar_Hemisphere__ENSG00000138675.tsv","BBJ_IS_L002__GTEx_V8__Artery_Tibial__ENSG00000107957.tsv","BBJ_IS_L004__GTEx_V8__Brain_Putamen_basal_ganglia__ENSG00000134871.tsv","BBJ_IS_L002__GTEx_V8__Brain_Cerebellar_Hemisphere__ENSG00000138172.tsv","BBJ_IS_L002__GTEx_V8__Brain_Cerebellar_Hemisphere__ENSG00000107954.tsv","BBJ_IS_L001__GTEx_V8__Brain_Cerebellar_Hemisphere__ENSG00000197826.tsv","BBJ_IS_L002__GTEx_V8__Brain_Anterior_cingulate_cortex_BA24__ENSG00000148798.tsv","BBJ_IS_L004__GTEx_V8__Brain_Amygdala__ENSG00000187498.tsv")
)
outrows<-list()
for(i in seq_len(nrow(config))){
 cfg<-config[i]
 d<-fread(file.path(inp,cfg$file))
 if(any(!is.finite(d$gwas_maf))||any(!is.finite(d$eqtl_maf)) ||
    any(d$gwas_maf<=0)||any(d$eqtl_maf<=0) ||
    any(d$gwas_maf>.5)||any(d$eqtl_maf>.5))stop("MAF_RANGE_FAIL")
 d[,maf_abs_delta:=abs(gwas_maf-eqtl_maf)]
 keyG<-d$match_key[which.min(d$gwas_p)]
 keyQ<-d$match_key[which.min(d$eqtl_p)]
 n0<-nrow(d)
 if(uniqueN(d$match_key)!=n0)stop("DUPLICATE_SHARED_SNP")
 for(thr in c(1,0.2,0.1)){
  z<-d[maf_abs_delta<=thr]
  n<-nrow(z)
  if(n<50)stop("INSUFFICIENT_FILTERED_SNPS")
  d1<-list(beta=as.numeric(z$gwas_beta),varbeta=as.numeric(z$gwas_se)^2,
      snp=as.character(z$match_key),MAF=as.numeric(z$gwas_maf),N=174686,
      s=22664/174686,type="cc")
  d2<-list(beta=as.numeric(z$eqtl_beta),varbeta=as.numeric(z$eqtl_se)^2,
      snp=as.character(z$match_key),MAF=as.numeric(z$eqtl_maf),
      N=as.numeric(z$qtl_n_scalar[1]),type="quant")
  su<-suppressWarnings(coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=1e-5)$summary)
  outrows[[length(outrows)+1]]<-data.table(
   gene=cfg$gene,maf_filter_max_diff=if(thr==1) NA_real_ else thr,
   scenario=if(thr==1) "UNFILTERED_HISTORICAL" else "MAF_DELTA_EXCLUSION_DIAGNOSTIC",
   n_original=n0,n_retained=n,n_excluded=n0-n,
   fraction_excluded=(n0-n)/n0,
   median_absolute_maf_delta=median(z$maf_abs_delta),
   gwas_lead_key=keyG,eqtl_lead_key=keyQ,
   lead_gwas_retained=keyG %in% z$match_key,
   lead_eqtl_retained=keyQ %in% z$match_key,
   min_gwas_p_retained=min(z$gwas_p),min_eqtl_p_retained=min(z$eqtl_p),
   PP_H3=as.numeric(su[["PP.H3.abf"]]),PP_H4=as.numeric(su[["PP.H4.abf"]]),
   H4_over_H3=as.numeric(su[["PP.H4.abf"]])/as.numeric(su[["PP.H3.abf"]])
  )
 }
}
res<-rbindlist(outrows)
if(nrow(res)!=24L)stop("INVALID_QC_GRID")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
fwrite(res,file.path(out,"IS_SNP_MAF_CROSS_ANCESTRY_DIAGNOSTIC.tsv"),sep="\t",na="NA")
writeLines(c(
 "# IS SNP MAF cross-population diagnostic (eight legacy gene–tissue input pairs)",
 "- Original SNPs are exact matched by GRCh38 variant coordinates and REF/ALT in the historical input.",
 "- GWAS sample ancestry (EAS/Japanese) differs from GTEx QTL cohort composition; large MAF differences are EXPECTED at some variants, not standalone allele harmonization failures.",
 "- This ad hoc post hoc MAF difference filtering is a stability SCREEN ONLY, not valid GWAS–QTL correction or a publication-ready coloc.",
 "- Allele matching, effect direction, genome build, sampled ancestry and regional coverage require independent source-level checks.",
 "- Baseline uses same p1=1e-4,p2=1e-4,p12=1e-5 and QTL N first scalar as recovered original notebook.",
 "- QTL and GWAS top p variants may differ due to LD and sample power; mismatch alone is not evidence against a shared causal allele.",
 "- Thresholds |MAF_GWAS - MAF_QTL| <=0.2 and <=0.1 are exploratory sensitivity choices; filtered posterior values cannot be used as causal evidence."
),file.path(out,"IS_SNP_MAF_DIAGNOSTIC_README.md"))
cat("IS_MAF_SENSITIVITY=PASS\n")
print(res[,.(gene,maf_filter_max_diff,n_retained,n_excluded,
 lead_gwas_retained,lead_eqtl_retained,PP_H3,PP_H4)])
