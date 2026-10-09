#!/usr/bin/env Rscript
# GTEx (eQTL Catalogue QTD000131) nominal artery-aorta eQTL x EAS AIS coloc.abf
# DIAGNOSTIC ONLY: common-SNP set lacks full-locus GWS lead; 720kb truncation;
# assumes one causal signal per trait and uses study-wide total N.
suppressPackageStartupMessages(library(coloc))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/eqtl_catalogue_aorta_v1"
audit<-read.table(file.path(root,"G0022_EQTL_CATALOGUE_AORTA_INPUT_AUDIT.tsv"),
 sep="\t",header=TRUE,stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(audit)!=8L)stop("Expected exact eight prespecified genes")
src<-jsonlite::fromJSON(file.path(root,"G0022_EQTL_CATALOGUE_AORTA_INPUT_SUMMARY.json"))
if(src$dataset_id!="QTD000131" || src$study_n!=387L ||
   !isTRUE(src$GRCh37_reconverted_and_exact_allele_matched))
 stop("Invalid eQTL source metadata")
result<-list()
for(i in seq_len(nrow(audit))){
  gene<-audit$gene[[i]]
  name<-file.path(root,paste0("G0022_",gene,"_AIS_AORTA_ABF_INPUT.tsv"))
  if(!file.exists(name)){
    result[[length(result)+1L]]<-data.frame(
      gene=gene,nsnps=0L,min_gwas_p=NA_real_,min_qtl_p=NA_real_,
      PP.H0=NA_real_,PP.H1=NA_real_,PP.H2=NA_real_,
      PP.H3=NA_real_,PP.H4=NA_real_,primary_coloc_gate="NO_COMMON_VARIANTS",
      interpretation="NOT_TESTED_NOT_NEGATIVE",stringsAsFactors=FALSE)
    next
  }
  tab<-read.table(name,sep="\t",header=TRUE,stringsAsFactors=FALSE,
                  check.names=FALSE)
  if(length(unique(tab$variant_id))!=nrow(tab))stop("Duplicate SNP inputs")
  gbet<-as.numeric(tab$gwas_beta);gse<-as.numeric(tab$gwas_se)
  qbet<-as.numeric(tab$qtl_beta);qse<-as.numeric(tab$qtl_se)
  geaf<-as.numeric(tab$gwas_alt_eaf);qmaf<-as.numeric(tab$qtl_maf)
  np<-nrow(tab);gp<-min(as.numeric(tab$gwas_p));qp<-min(as.numeric(tab$qtl_p))
  numerically_valid<-np>=100&&all(is.finite(c(gbet,gse,qbet,qse,geaf,qmaf)))&&
    all(gse>0)&&all(qse>0)&&all(geaf>0&geaf<1)&&all(qmaf>0&qmaf<=.5)
  signal_gate<-numerically_valid && gp<=5e-8 && qp<=1e-4
  post<-setNames(rep(NA_real_,5),paste0("PP.H",0:4))
  if(numerically_valid){
    d_gwas<-list(snp=tab$variant_id,type="cc",beta=gbet,
       varbeta=gse^2,N=256274L,s=19032/256274,MAF=geaf)
    d_qtl<-list(snp=tab$variant_id,type="quant",beta=qbet,
       varbeta=qse^2,N=387L,MAF=qmaf)
    abf<-tryCatch(coloc.abf(d_gwas,d_qtl,p1=1e-4,p2=1e-4,p12=1e-5),
         error=function(e){message("coloc failure ",gene,": ",conditionMessage(e));NULL})
    if(!is.null(abf))for(j in 0:4){
      k<-paste0("PP.H",j,".abf")
      if(k%in%names(abf$summary))post[[j+1L]]<-as.numeric(abf$summary[[k]])
    }
  }
  status<-if(!numerically_valid)"FAILED_INPUT_QC" else
    if(!signal_gate)"DIAGNOSTIC_ONLY_MISSING_COMMON_GWAS_GWS_OR_QTL_SIGNAL" else
      "EXPLORATORY_ABF_ONLY_SINGLE_SIGNAL_REGION_TRUNCATED"
  result[[length(result)+1L]]<-data.frame(gene=gene,nsnps=np,
      min_gwas_p=gp,min_qtl_p=qp,
      PP.H0=post[[1]],PP.H1=post[[2]],PP.H2=post[[3]],
      PP.H3=post[[4]],PP.H4=post[[5]],
      primary_coloc_gate=status,
      interpretation="DO_NOT_CLAIM_SHARED_CAUSAL_VARIANT_WITHOUT_COMPLETE_SNP_UNIVERSE",
      stringsAsFactors=FALSE)
  cat("gene",gene,"n",np,"GWAS p",gp,"QTL p",qp,
      "PPH4",post[[5]],"GATE",status,"\n")
}
res<-do.call(rbind,result)
write.table(res,file.path(root,"G0022_EAS_AIS_GTEX_AORTA_COLOC_ABF_DIAGNOSTIC.tsv"),
 row.names=FALSE,quote=FALSE,sep="\t")
cat("TOTAL tested posterior distributions",
 sum(is.finite(res$PP.H4)),"PRIMARY VALIDATED COLOCALIZATION",0,"\n")
