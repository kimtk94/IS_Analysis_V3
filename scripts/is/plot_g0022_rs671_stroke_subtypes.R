#!/usr/bin/env Rscript
# Descriptive rs671-A association forest figure from six ancestry-source strata.
# Deliberately NOT a meta-analysis or statistical cross-subtype comparison.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
base<-if(length(args)>0)args[[1]] else
"/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1"
x<-read.delim(file.path(base,"G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv"),
  check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(x)!=6||anyDuplicated(x$study)||!all(x$ALT_beta<0))
 stop("Expected six verified ALT-A protective rs671 phenotype rows")
ordered<-c("Japanese_IS","EAS_AS","EAS_AIS","EAS_SVS","EAS_CES","EAS_LAS")
if(!setequal(x$trait,ordered))stop("Phenotype schema changed")
pretty<-c(Japanese_IS="BBJ Ischemic stroke",EAS_AS="GIGASTROKE Any stroke",
  EAS_AIS="GIGASTROKE Any ischemic",
  EAS_SVS="GIGASTROKE Small-vessel",
  EAS_CES="GIGASTROKE Cardioembolic",
  EAS_LAS="GIGASTROKE Large artery")
x$Label<-factor(unname(pretty[x$trait]),levels=rev(unname(pretty[ordered])))
x$label_or<-sprintf("%.3f (%.3f–%.3f)",x$ALT_OR,x$OR_95CI_lower,x$OR_95CI_upper)
x$label_p<-ifelse(x$p<.0001,format(x$p,scientific=TRUE,digits=2),
  sprintf("%.3f",x$p))
plot<-ggplot(x,aes(x=ALT_OR,y=Label))+
  geom_vline(xintercept=1,linetype="dashed",color="grey50",linewidth=.5)+
  geom_errorbar(aes(xmin=OR_95CI_lower,xmax=OR_95CI_upper),
     orientation="y",width=.15,color="#4d6476",linewidth=.65)+
  geom_point(size=2.75,color="#175a7e")+
  geom_text(aes(x=1.31,label=label_or),hjust=0,size=3.3,color="#263c49")+
  geom_text(aes(x=1.61,label=paste0("p=",label_p)),hjust=0,size=3.2,color="#263c49")+
  annotate("text",x=1.31,y=6.75,label="OR (95% CI)",
     fontface="bold",hjust=0,size=3.5)+
  annotate("text",x=1.61,y=6.75,label="Nominal p",
     fontface="bold",hjust=0,size=3.5)+
  scale_x_continuous(limits=c(.60,1.8),breaks=c(.6,.8,1.0,1.2))+
  labs(title="ALDH2 rs671-A and ischemic stroke phenotypes",
    subtitle="GWAS log-odds harmonized to rs671 ALT-A (GRCh37 chr12:112241766 G>A)",
    x="Odds ratio per A allele (association; not causal)",y=NULL,
    caption=paste("Descriptive only — BBJ/GIGASTROKE datasets may share cohorts;",
      "nested case definitions and overlapping controls prevent treating these rows as independent replication.",
      "Subtypes are not formally contrasted.",sep="\n"))+
  theme_minimal(base_size=11)+
  theme(axis.title.x=element_text(size=10),
     panel.grid.minor=element_blank(),
     panel.grid.major.y=element_blank(),
     plot.title=element_text(face="bold",size=15),
     plot.caption=element_text(hjust=0,size=8,color="#606d77"),
     plot.margin=margin(20,25,16,15),
     axis.text.y=element_text(color="#263c49",size=10))
out<-file.path(base,"G0022_RS671_STROKE_PHENOTYPE_FOREST.png")
ggsave(out,plot=plot,width=12.5,height=6.3,units="in",dpi=160,bg="white")
cat("FIGURE_OK",out,"bytes",file.info(out)$size,"\n")
