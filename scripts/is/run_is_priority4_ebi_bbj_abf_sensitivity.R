#!/usr/bin/env Rscript
# Separate, exploratory coloc.abf, EBI QTD000161 reprocessing × BBJ IS.
# GENOTYPE AN/2 does not independently attest the QTL expression-effective N.
# No SuSiE conditional inference or historical H4 revision.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2)stop("usage: Rscript runner.R INPUT.tsv OUTPUT_DIR")
inp <- args[[1]]; dest <- args[[2]]
if(!file.exists(inp))stop("missing input")
if(dir.exists(dest) && length(list.files(dest))>0)stop("refuse overwriting science output")
if(!requireNamespace("coloc",quietly=TRUE))stop("coloc R package unavailable")
set.seed(20261011)
d <- read.delim(inp,check.names=FALSE,stringsAsFactors=FALSE)
need <- c("gene_symbol","gene_id","locus","variant_GRCh37_REF_ALT","BBJ_beta_ALT",
 "BBJ_se","BBJ_p","BBJ_MAF","EBI_QTD_beta_ALT","EBI_QTD_se",
 "EBI_QTD_p","EBI_QTD_maf","EBI_QTD_AN")
if(!all(need %in% names(d)))stop("input source schema drift")
if(nrow(d)!=11050L || length(unique(d$gene_symbol))!=4L)stop("source cardinality drift")
for(v in c("BBJ_beta_ALT","BBJ_se","BBJ_p","BBJ_MAF","EBI_QTD_beta_ALT",
           "EBI_QTD_se","EBI_QTD_p","EBI_QTD_maf","EBI_QTD_AN")) {
 d[[v]] <- as.numeric(d[[v]])
 if(any(!is.finite(d[[v]])))stop(paste("nonfinite",v))
}
if(any(d$BBJ_se<=0|d$EBI_QTD_se<=0|d$BBJ_MAF<=0|d$BBJ_MAF>0.5|
       d$EBI_QTD_maf<=0|d$EBI_QTD_maf>0.5|
       d$BBJ_p<0|d$BBJ_p>1|d$EBI_QTD_p<0|d$EBI_QTD_p>1|
       d$EBI_QTD_AN<=0))stop("source invalid association, MAF or genotype AN")
if(anyDuplicated(paste(d$gene_symbol,d$variant_GRCh37_REF_ALT)))stop("duplicated gene-SNP")
if(any(!as.character(d$QTD_nominal_full_cis_tested) %in% c("True","TRUE")) ||
   any(!as.character(d$BBJ_GWAS_exact_FASTA_ref) %in% c("True","TRUE")))stop("upstream attestation flags absent")
if(any(d$harmonization!="ALT_AFTER_CHAIN_AND_FASTA_REF_EQUIVALENCE"))stop("not exact allele ALT")
genes <- c("FGF5","C4orf22","CALHM2","NEURL")
expected <- c(FGF5=2773,C4orf22=2772,CALHM2=2766,NEURL=2739)
ns <- c(125,150,175)
priors <- c(1e-6,3e-6,1e-5,3e-5,1e-4)
res <- list();qc <- list();rownum <- 0L
for(g in genes) {
 x <- d[d$gene_symbol==g,,drop=FALSE]
 if(nrow(x)!=expected[[g]])stop(paste("unexpected gene source intersection",g,nrow(x)))
 an <- x$EBI_QTD_AN/2
 if(any(an>175+1e-8) || any(an<50))stop("QTL genotype AN ceiling drift: cannot use prior sample proxy")
 qc[[g]] <- data.frame(gene_symbol=g,n_total=nrow(x),
   genomic_AN_over_2_min=min(an),genomic_AN_over_2_median=median(an),
   genomic_AN_over_2_max=max(an),
   BBJ_N_assumed=174686,BBJ_case_fraction=22664/174686,
   QTL_expression_effective_N_attested=FALSE,
   QTL_cohort_signed_LD_available=FALSE,
   QTL_study_original_coloc_version_proven=FALSE)
 for(n in ns) {
   d1 <- list(snp=x$variant_GRCh37_REF_ALT,
      beta=x$BBJ_beta_ALT,varbeta=x$BBJ_se^2,MAF=x$BBJ_MAF,
      N=174686,s=22664/174686,type="cc")
   d2 <- list(snp=x$variant_GRCh37_REF_ALT,
      beta=x$EBI_QTD_beta_ALT,varbeta=x$EBI_QTD_se^2,MAF=x$EBI_QTD_maf,
      N=n,type="quant")
   for(prior in priors) {
      fit <- tryCatch(
         suppressWarnings(suppressMessages(coloc::coloc.abf(
           dataset1=d1,dataset2=d2,p1=1e-4,p2=1e-4,p12=prior))),
         error=function(e)e)
      rownum <- rownum+1L
      if(inherits(fit,"error")) {
       res[[rownum]] <- data.frame(gene_symbol=g,gene_id=unique(x$gene_id),
          locus=unique(x$locus),snp_intersection=nrow(x),QTL_N_assumption=n,
          assumed_p1=1e-4,assumed_p2=1e-4,assumed_p12=prior,
          PP_H0=NA_real_,PP_H1=NA_real_,PP_H2=NA_real_,PP_H3=NA_real_,
          PP_H4=NA_real_,max_single_variant_SNP_PP_H4=NA_real_,
          status="FAILED_SENSITIVITY",diagnostic_message=substr(conditionMessage(fit),1,250),
          original_GTEx_v8_model_reproduced=FALSE,study_matched_QTL_LD_available=FALSE)
      } else {
       v<-as.numeric(fit$summary[c("PP.H0.abf","PP.H1.abf","PP.H2.abf","PP.H3.abf","PP.H4.abf")])
       if(any(!is.finite(v))||abs(sum(v)-1)>1e-6)stop("coloc summary posterior invalid")
       if(as.integer(fit$summary["nsnps"])!=nrow(x))stop("coloc output SNP set not same source")
       snppp <- as.numeric(fit$results$SNP.PP.H4)
       res[[rownum]] <- data.frame(gene_symbol=g,gene_id=unique(x$gene_id),
          locus=unique(x$locus),snp_intersection=nrow(x),QTL_N_assumption=n,
          assumed_p1=1e-4,assumed_p2=1e-4,assumed_p12=prior,
          PP_H0=v[[1]],PP_H1=v[[2]],PP_H2=v[[3]],PP_H3=v[[4]],
          PP_H4=v[[5]],max_single_variant_SNP_PP_H4=max(snppp,na.rm=TRUE),
          status="EXPLORATORY_ABF_SUCCESS",
          diagnostic_message="EBI reanalysis; inferred N scenario; single causal per trait assumption",
          original_GTEx_v8_model_reproduced=FALSE,study_matched_QTL_LD_available=FALSE)
      }
   }
 }
 cat("GENE_RUNS",g,"n",nrow(x),"completed",length(res),"\n",file=stderr())
}
if(!dir.exists(dest))dir.create(dest,recursive=TRUE)
r <- do.call(rbind,res)
if(nrow(r)!=60L)stop("Expected 4 genes x 3 N x 5 priors, never silently skip")
write.table(r,file=file.path(dest,"IS_PRIORITY4_EBI_BBJ_ABF_60_SCENARIO_SENSITIVITY.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE,na="")
write.table(do.call(rbind,qc),
 file=file.path(dest,"IS_PRIORITY4_EBI_BBJ_ABF_QTL_N_PROXY_QC.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE,na="")
if(any(r$status!="EXPLORATORY_ABF_SUCCESS"))stop("One or more sensitivity scenarios failed; inspect TSV")
cat("SENSITIVITY_PASS\n")
print(subset(r,QTL_N_assumption==175 & assumed_p12==1e-5,
 select=c(gene_symbol,snp_intersection,PP_H3,PP_H4)))
