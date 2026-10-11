#!/usr/bin/env Rscript
# Experimental coloc.abf - deliberately matched SNP sets in COL4A1/COL4A2,
# separated by QTL source/ AF difference/ p12 / inferred sdY-N perturbation.
# Association and posterior are NOT validation of shared causal variants.
suppressPackageStartupMessages(library(coloc))
args<-commandArgs(trailingOnly=TRUE)
value<-function(k) {
 i<-which(args==k)
 if(length(i)!=1L || i==length(args))stop(paste("Missing",k))
 args[[i+1L]]
}
inp<-value("--input")
source_qc<-value("--source-qc")
outdir<-value("--out-dir")
if(dir.exists(outdir) && length(list.files(outdir,all.files=TRUE,no..=TRUE))>0)
 stop("Refuse replacing science output")
x<-read.delim(inp,stringsAsFactors=FALSE,check.names=FALSE)
qc<-read.delim(source_qc,stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(qc)!=6L ||
   !setequal(unique(x$gene_symbol),c("COL4A1","COL4A2")) ||
   !setequal(unique(x$QTL_dataset),c("QTD000131","QTD000171","QTD000434")))
 stop("Missing exact paired cis dataset gene universe")
sample_n<-c(QTD000131=387,QTD000171=205,QTD000434=560)
prior_grid<-c(1e-6,3e-6,1e-5,3e-5,1e-4)
band_grid<-c("ALL_NONPAL_PAIRED","AF_DIFF_MAX_0P2_BOTH_GENES","AF_DIFF_MAX_0P1_BOTH_GENES")
Nscales<-c(0.8,1.0,1.2)
gwas_cases<-19032
gwas_ctrl<-237242
gwas_N<-256274
if(gwas_cases+gwas_ctrl!=gwas_N)stop("Original EAS AIS sample arithmetic")
result<-list()
diag<-list()
for(did in names(sample_n)){
 both<-x[x$QTL_dataset==did,,drop=FALSE]
 genes<-split(both,both$gene_symbol)
 if(length(genes)!=2L)stop("Two genes not present per tissue")
 ids_a<-genes$COL4A1$GRCh37_fasta_REF_ALT_variant
 ids_b<-genes$COL4A2$GRCh37_fasta_REF_ALT_variant
 if(anyDuplicated(ids_a) || anyDuplicated(ids_b) || !setequal(ids_a,ids_b))stop("Gene-paired SNP set changed")
 stopifnot(all(is.finite(as.numeric(both$QTL_n_variant_AN_div_2))))
 if(any(abs(both$QTL_n_variant_AN_div_2-sample_n[[did]])>22))
  stop("Per-SNP official QTL AN/2 far from study sample count")
 for(band in band_grid){
  ids<-ids_a
  if(band!="ALL_NONPAL_PAIRED"){
   maxdiff<-if(band=="AF_DIFF_MAX_0P1_BOTH_GENES")0.1 else 0.2
   sharedA<-genes$COL4A1$GRCh37_fasta_REF_ALT_variant[
    genes$COL4A1$QTL_EAS_AIS_AF_DELTA<=maxdiff]
   sharedB<-genes$COL4A2$GRCh37_fasta_REF_ALT_variant[
    genes$COL4A2$QTL_EAS_AIS_AF_DELTA<=maxdiff]
   ids<-intersect(sharedA,sharedB)
  }
  for(gene in c("COL4A1","COL4A2")){
   z<-genes[[gene]][match(ids,genes[[gene]]$GRCh37_fasta_REF_ALT_variant),,drop=FALSE]
   if(nrow(z)!=length(ids)||anyNA(z$QTL_beta_alt))stop("Missing gene-matched SNP")
   if(length(ids)<100L)stop("Insufficient matched SNPs for diagnostic")
   nums<-as.matrix(z[,c("QTL_beta_alt","QTL_se","QTL_p","QTL_ALT_AF",
     "EAS_AIS_beta_alt","EAS_AIS_se","EAS_AIS_p","EAS_AIS_ALT_AF",
     "BBJ_IS_p","QTL_n_variant_AN_div_2")])
   if(any(!is.finite(nums)) || any(z$QTL_se<=0) || any(z$EAS_AIS_se<=0) ||
      any(z$QTL_ALT_AF<=0|z$QTL_ALT_AF>=1) ||
      any(z$EAS_AIS_ALT_AF<=0|z$EAS_AIS_ALT_AF>=1))
    stop("Unusable beta SE p AF input for coloc")
   snps<-z$GRCh37_fasta_REF_ALT_variant
   # Use beta/varbeta estimates, N + genotype MAF to estimate sdY for quant.
   # Source AN/2 is independently checked against study record n=387/205/560.
   A<-list(snp=snps,type="quant",beta=z$QTL_beta_alt,
           varbeta=z$QTL_se^2,N=sample_n[[did]],
           MAF=pmin(z$QTL_ALT_AF,1-z$QTL_ALT_AF))
   B<-list(snp=snps,type="cc",beta=z$EAS_AIS_beta_alt,
           varbeta=z$EAS_AIS_se^2,N=gwas_N,
           s=gwas_cases/gwas_N,
           MAF=pmin(z$EAS_AIS_ALT_AF,1-z$EAS_AIS_ALT_AF))
   check_dataset(A);check_dataset(B)
   min_qtl<-min(z$QTL_p);min_gwas<-min(z$EAS_AIS_p)
   diag[[length(diag)+1]]<-data.frame(
    QTL_dataset=did,gene_symbol=gene,SNP_selection=band,nsnps=length(snps),
    original_cis_n_variants=qc$official_tested_GRCh38_variant_count[qc$QTL_dataset==did & qc$gene_symbol==gene],
    min_QTL_p=min_qtl,min_EAS_AIS_p=min_gwas,
    min_BBJ_IS_p=min(z$BBJ_IS_p),
    any_QTL_P_le_1e4=any(z$QTL_p<=1e-4),
    any_AIS_P_le_5e8=any(z$EAS_AIS_p<5e-8),
    sameSNP_QTL1e4_AIS5e8=sum(z$QTL_p<=1e-4 & z$EAS_AIS_p<5e-8),
    MAF_QTL_vs_EAS_AIS_abs_diff_median=median(z$QTL_EAS_AIS_AF_DELTA),
    stringsAsFactors=FALSE)
   for(factor in Nscales) {
    a<-A
    a$N<-round(sample_n[[did]]*factor)
    estimated_sdY<-suppressWarnings(tryCatch(
     coloc:::sdY.est(a$varbeta,a$MAF,a$N),
     error=function(e)NA_real_))
    if(!is.finite(estimated_sdY) || estimated_sdY<=0)
     stop("QTL sdY cannot be estimated from study N, MAF, varbeta")
    for(p12 in prior_grid){
     pp<-rep(NA_real_,5)
     names(pp)<-paste0("PP.H",0:4)
     success<-FALSE
     toplocus<-NA_character_
     topsnp<-NA_real_
     errored<-""
     abf<-tryCatch(suppressWarnings(coloc::coloc.abf(
       a,B,p1=1e-4,p2=1e-4,p12=p12)),
       error=function(e){errored<<-conditionMessage(e);NULL})
     if(!is.null(abf)) {
      for(j in 0:4){
       column<-paste0("PP.H",j,".abf")
       if(!(column%in%names(abf$summary)))
        stop("coloc package v5.2.3 posterior field changed")
       pp[[j+1L]]<-as.numeric(abf$summary[[column]])
      }
      if(any(!is.finite(pp)) || abs(sum(pp)-1)>1e-4)
       stop("Invalid coloc posterior mass")
      if(nrow(abf$results)!=length(snps))
       stop("coloc package merged SNP loss")
      topidx<-which.max(abf$results$SNP.PP.H4)
      toplocus<-as.character(abf$results$snp[[topidx]])
      topsnp<-as.numeric(abf$results$SNP.PP.H4[[topidx]])
      success<-TRUE
     }
     status<-if(!success)"ABF_FAILED_NO_BIOLOGICAL_STATUS" else
       if(!any(z$QTL_p<=1e-4))"DIAGNOSTIC_ONLY_QTL_GENE_SIGNAL_WEAK" else
       "EXPLORATORY_SINGLE_SIGNAL_PRIOR_AND_CROSS_ANCESTRY_SENSITIVE"
     result[[length(result)+1]]<-data.frame(
      QTL_dataset=did,gene_symbol=gene,SNP_selection=band,
      QTL_N_AN_over2_reference=sample_n[[did]],
      QTL_N_multiplier=factor,QTL_N_coloc_used=a$N,
      EAS_AIS_N_secondary_source=gwas_N,
      EAS_AIS_case_fraction_secondary=gwas_cases/gwas_N,
      nsnps=length(snps),QTL_min_p=min_qtl,
      EAS_AIS_min_p=min_gwas,p1=1e-4,p2=1e-4,p12=p12,
      estimated_QTL_sdY=estimated_sdY,
      PP_H0=pp[[1]],PP_H1=pp[[2]],PP_H2=pp[[3]],
      PP_H3=pp[[4]],PP_H4=pp[[5]],
      H4_over_H3=if(is.finite(pp[[4]]) && pp[[4]]>0)pp[[5]]/pp[[4]] else NA_real_,
      conditional_H4_given_H3_or_H4=if(all(is.finite(pp[4:5])) &&
        sum(pp[4:5])>0)pp[[5]]/sum(pp[4:5]) else NA_real_,
      top_H4_conditional_SNP=toplocus,
      top_H4_conditional_SNP_PP=topsnp,
      pipeline_gate=status,
      failure_reason=errored,
      co_localized_gene_claim="NOT_ESTABLISHED",
      EAS_BBJ_independent_replication_attested=FALSE,
      QTL_EAS_ancestry_matched_LD_attested=FALSE,
      stringsAsFactors=FALSE)
    }
   }
   cat("COL4A_PAIRED_COLOC_ABF_FINISHED",did,gene,band,
     "n",length(snps),"minQTL",format(min_qtl,digits=3),
     "minAIS",format(min_gwas,digits=3))
   cat("\n")
  }
 }
}
r<-do.call(rbind,result);dq<-do.call(rbind,diag)
if(nrow(r)!=270L||nrow(dq)!=18L)stop("Expected 270 posterior sensitivity rows / 18 input-QC records")
if(any(duplicated(dq[,c("QTL_dataset","gene_symbol","SNP_selection")])))
 stop("Duplicated paired source keys")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.table(r,file.path(outdir,"IS_COL4A1_COL4A2_THREE_QTL_P12_AF_N_ABF_SENSITIVITY.tsv"),
 row.names=FALSE,quote=FALSE,sep="\t")
write.table(dq,file.path(outdir,"IS_COL4A1_COL4A2_SNP_GENE_PAIRED_INPUT_CIS_QC.tsv"),
 row.names=FALSE,quote=FALSE,sep="\t")
cat("COL4A1_COL4A2_THREE_TISSUE_PAIRED_SINGLE_SIGNAL_ABF_COMPLETE",
 "rows",nrow(r),"PASS",sum(is.finite(r$PP_H4)),
 "new_verified_shared_causal_variant",0,"\n")
