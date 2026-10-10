#!/usr/bin/env Rscript
# Diagnostic-only Figure: full original exposure SNPs and 1KG EAS504.
# Both corrected methods R_reliability_flag TRUE. Never publish causal claim.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/susie_rss_finite_ref_sandbox_v1"
dest<-if(length(args)>1)args[[2]] else
 "/srv/is-analysis/worktrees/is-broad-discovery-20261009/docs/figures/is"
x<-read.delim(file.path(root,"IS_ALCOHOL_FINITE_REF_LD_CORRECTION_MODEL_COMPARISON.tsv"),
 check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(x)!=4||!all(tolower(as.character(x$R_reliability_flag))=="true")||
   !all(tolower(as.character(x$converged))=="true") ||
   !identical(sort(unique(x$locus)),c("ADH1B","ALDH2")))
 stop("Missing sensitivity experiment QC or wrong model scope")
x$model_label<-ifelse(x$model=="FINITE_REF_504","Finite LD correction (B=504)",
                     "Finite LD + EB mismatch")
x$model_label<-factor(x$model_label,
 levels=c("Finite LD correction (B=504)","Finite LD + EB mismatch"))
x$locus<-factor(x$locus,levels=c("ADH1B","ALDH2"))
x$note<-paste0(x$exploratory_95pct_CS_n," candidate CS")
p<-ggplot(x,aes(x=model_label,y=exploratory_95pct_CS_n,fill=locus))+
 geom_col(width=.58,alpha=.88)+
 geom_text(aes(label=note),vjust=-.50,size=4)+
 facet_wrap(~locus,nrow=1)+
 scale_fill_manual(values=c("#2e7795","#c1884b"))+
 scale_y_continuous(limits=c(0,5),breaks=0:5)+
 labs(title="Corrected SuSiE-RSS converges but LD reliability still fails",
 subtitle="Original Japanese alcohol GWAS beta/SE; 1000G EAS 504 reference; isolated susieR v0.16.6",
      x=NULL,y="Exploratory 95% credible sets (NOT validated causal effects)",
      caption=paste("Both correction methods and BOTH regions have R_reliability_flag=TRUE.",
         "Median LD correction penalty: ADH1B ~1.16; ALDH2 ~66.4 (high).",
         "Baseline ALDH2 SuSiE was BLOCKED by GWAS-LD mismatch s=0.0908; these exploratory bars do not lift that gate.",
         "No independent AIS or alcohol-mediation causal inference was performed.",sep="\n"))+
 theme_minimal(base_size=11)+
 theme(legend.position="none",axis.text.x=element_text(angle=16,hjust=1),
       panel.grid.major.x=element_blank(),strip.text=element_text(face="bold",size=12),
       plot.title=element_text(face="bold",size=15),
       plot.caption=element_text(color="#52616f",hjust=0,size=8))
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
path<-file.path(dest,"IS_ALCOHOL_FINITE_LD_CORRECTION_RELIABILITY_FAIL.png")
ggsave(path,p,width=11.2,height=6.3,dpi=170,bg="white")
cat("SANDBOX_FIGURE_QC_PASS",path,file.info(path)$size,"bytes\n")
