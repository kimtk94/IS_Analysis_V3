#!/usr/bin/env Rscript
# Real QTD000131/136/141/171 GTEx gene-level nominal eQTL x EAS AIS.
# PP.H4 is DIAGNOSTIC ONLY until common GWS, complete cis universe, multi-signal LD.
suppressPackageStartupMessages(library(coloc))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>0)args[[1]] else
  "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/eqtl_catalogue_multitissue_v1"
studies <- c("QTD000131","QTD000136","QTD000141","QTD000171")
genes <- c("ACAD10","ALDH2","NAA25","HECTD4","BRAP","PTPN11","IFT81","ATP2A2")
out <- list()
for(study in studies){
  p <- file.path(root,study)
  if(!file.exists(file.path(p,"G0022_MULTITISSUE_QTL_INPUT_SUMMARY.json")))
    stop("Missing source-stamped original QTL extract: ",study)
  src <- jsonlite::fromJSON(file.path(p,"G0022_MULTITISSUE_QTL_INPUT_SUMMARY.json"))
  if(src$dataset_id!=study || src$status!="ALLELE_HARMONIZATION_COMPLETE_COLOC_UNVALIDATED")
    stop("Source manifest status mismatch")
  for(gene in genes){
    f <- file.path(p,paste0("G0022_",gene,"_AIS_",study,"_ABF_INPUT.tsv"))
    empty <- data.frame(
      dataset_id=study,tissue=src$tissue,gene=gene,matched_snps=0L,
      min_gwas_p=NA_real_,min_qtl_p=NA_real_,joint_GWS_present=FALSE,
      GWAS_CS_SNP_matches=0L,PP.H0=NA_real_,PP.H1=NA_real_,
      PP.H2=NA_real_,PP.H3=NA_real_,PP.H4=NA_real_,
      posterior_computed=FALSE,gate="NO_COMMON_ASSAYED_SNPS_NOT_NEGATIVE",
      conclusion="NOT_VALIDATED_COLOCALIZATION",stringsAsFactors=FALSE)
    if(!file.exists(f)){out[[length(out)+1L]]<-empty;next}
    x<-read.table(f,header=TRUE,sep="\t",stringsAsFactors=FALSE,
                   check.names=FALSE)
    if(anyDuplicated(x$variant_id_grch37))stop("Duplicated allele-oriented variant")
    nrowx<-nrow(x)
    gb<-as.numeric(x$gwas_beta_alt);gs<-as.numeric(x$gwas_se)
    qb<-as.numeric(x$qtl_beta_alt);qs<-as.numeric(x$qtl_se)
    gfreq<-as.numeric(x$gwas_alt_eaf)
    qmaf<-as.numeric(x$qtl_maf)
    gp<-min(as.numeric(x$gwas_p));qp<-min(as.numeric(x$qtl_p))
    cs<-sum(as.integer(x$gwas_95pct_cs_snp))
    gws<-gp<=5e-8
    base<-empty
    base$matched_snps<-nrowx
    base$min_gwas_p<-gp;base$min_qtl_p<-qp
    base$joint_GWS_present<-gws;base$GWAS_CS_SNP_matches<-cs
    basicok<-nrowx>=100L&&all(is.finite(c(gb,gs,qb,qs,gfreq,qmaf)))&&
      all(gs>0)&&all(qs>0)&&all(gfreq>0&gfreq<1)&&all(qmaf>0&qmaf<=.5)
    if(!basicok){base$gate<-"INPUT_QC_OR_MINIMUM_N_UNMET";out[[length(out)+1L]]<-base;next}
    # coloc.abf with beta/varbeta estimates and study-level sample sizes,
    # using allele freq -> MAF only for auxiliary prior sdY approximations.
    d1<-list(snp=x$variant_id_grch37,beta=gb,varbeta=gs^2,
             type="cc",N=256274L,s=19032/256274,MAF=pmin(gfreq,1-gfreq))
    d2<-list(snp=x$variant_id_grch37,beta=qb,varbeta=qs^2,
             type="quant",N=src$catalogue_sample_size,MAF=qmaf)
    withCallingHandlers({
      fit<-tryCatch(coloc::coloc.abf(d1,d2,p1=1e-4,p2=1e-4,p12=1e-5),
         error=function(e){message("FAIL ",study," ",gene,": ",conditionMessage(e));NULL})
      if(!is.null(fit)){
        for(i in 0:4){key<-paste0("PP.H",i,".abf")
          if(key%in%names(fit$summary))base[[paste0("PP.H",i)]]<-unname(fit$summary[[key]])}
        base$posterior_computed<-TRUE
      }
    },warning=function(w){invokeRestart("muffleWarning")})
    base$gate<-if(!base$posterior_computed)"COLOC_NUMERICAL_ERROR" else
      if(!gws)"BLOCKED_NO_SHARED_GENOMEWIDE_SIGNIFICANT_AIS_VARIANT" else
      if(cs==0L)"BLOCKED_NONE_OF_FOUR_AIS_CS_SNPS_ASSAYED" else
      if(qp>1e-4)"BLOCKED_WEAK_NOMINAL_QTL_IN_SHARED_SET" else
      "EXPLORATORY_ABF_REQUIRES_FULL_CIS_AND_MULTI_SIGNAL_VALIDATION"
    out[[length(out)+1L]]<-base
    cat(study,gene,"n",nrowx,"GWASp",gp,"QTLp",qp,
        "CS",cs,"diagnostic_pph4",base$PP.H4,"gate",base$gate,"\n")
    flush.console()
  }
}
r<-do.call(rbind,out)
write.table(r,file.path(root,"G0022_GTEX_MULTITISSUE_ABF_DIAGNOSTIC.tsv"),
 sep="\t",quote=FALSE,row.names=FALSE)
summary<-list(
  status="COMPUTATIONAL_QTL_COLOC_DIAGNOSTIC_ONLY",
  analyzed_studies=length(studies),gene_tissue_tests=nrow(r),
  computed_posteriors=sum(r$posterior_computed),
  datasets_with_any_AIS_GWS_in_matching_QTL_snps=
    length(unique(r$dataset_id[r$joint_GWS_present])),
  posterior_results_eligible_for_publication=0L,
  causal_genes_validated=0L,
  caveat="Common-SNP ascertainment, 720kb truncated QTL, ancestry mismatch, per-SNP N missing, GTEx QTL multi-signal uncertainty")
jsonlite::write_json(summary,file.path(root,"G0022_MULTITISSUE_ABF_SUMMARY.json"),
 pretty=TRUE,auto_unbox=TRUE)
print(summary)
