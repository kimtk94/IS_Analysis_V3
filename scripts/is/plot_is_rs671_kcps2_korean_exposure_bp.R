#!/usr/bin/env Rscript
# Fully-source-derived KCPS2 rs671 A alcohol + SBP/DBP association panel.
# All phenotypes are original inverse-normalized traits; neither mmHg nor
# drinking-gram reductions; do not interpret cross-trait relative bar length
# as a mediation or proportion explained.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>=1)args[[1]] else
"/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1"
alc<-read.delim(file.path(root,"G0022_RS671_KCPS2_KOREAN_ALCOHOL_AMOUNT.tsv"),check.names=FALSE)
sbp<-read.delim(file.path(root,"G0022_RS671_KCPS2_KOREAN_SBP_GWAS.tsv"),check.names=FALSE)
dbp<-read.delim(file.path(root,"G0022_RS671_KCPS2_KOREAN_DBP_GWAS.tsv"),check.names=FALSE)
if(nrow(alc)!=1||nrow(sbp)!=1||nrow(dbp)!=1)stop("Three source-validated rs671 rows required")
if(alc$effect_allele_harmonized!="A"||
   sbp$harmonized_effect_allele!="A"||
   dbp$harmonized_effect_allele!="A")stop("Only rs671-A harmonized effects permitted")
rows<-data.frame(
  label=c("Alcohol intake","Systolic blood pressure","Diastolic blood pressure"),
  beta=c(alc$KCPS2_alcohol_beta_A,sbp$KCPS2_beta_A,dbp$KCPS2_beta_A),
  se=c(alc$KCPS2_alcohol_se,sbp$KCPS2_se,dbp$KCPS2_se),
  n=c(alc$KCPS2_sample_size_at_variant,sbp$KCPS2_per_variant_N,dbp$KCPS2_per_variant_N)
)
if(any(!is.finite(rows$beta))||any(rows$se<=0)||any(rows$n<1000))stop("Nonfinite source effects")
rows$label<-factor(rows$label,levels=rev(rows$label))
rows$lower<-rows$beta-1.96*rows$se
rows$upper<-rows$beta+1.96*rows$se
rows$annotation<-sprintf("β = %.4f ± %.4f;  N = %s",rows$beta,rows$se,format(rows$n,big.mark=",",trim=TRUE,scientific=FALSE))
p<-ggplot(rows,aes(x=beta,y=label))+
 geom_vline(xintercept=0,linetype="dashed",color="grey55")+
 geom_errorbar(aes(xmin=lower,xmax=upper),width=.15,orientation="y",linewidth=.8,color="#30546e")+
 geom_point(size=3,color="#176087")+
 geom_text(aes(x=.025,label=annotation),hjust=0,size=3.4,color="#233646")+
 scale_x_continuous(breaks=seq(-.6,0,.2),limits=c(-.72,.60))+
 labs(title="Korean KCPS2: ALDH2 rs671-A associations",
      subtitle="Original verified GWAS source, 2025 KCPS2; rs671 GRCh37 chr12:112241766 G>A",
      x="Per-A-allele β (trait-specific inverse-rank normal scale)",y=NULL,
      caption=paste(
       "95% CI = beta ± 1.96 SE. Three separate quantitative phenotypes; scales",
       "are inverse-normalized within each trait. No causal mediation, drug effect",
       "or absolute mmHg/day drinking change is inferred.",sep="\n"))+
 theme_minimal(base_size=12)+
 theme(plot.title=element_text(face="bold",size=15),
       panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),
       plot.caption=element_text(hjust=0,color="#566773",size=9),
       axis.text.y=element_text(color="#233646",size=10),
       plot.margin=margin(20,25,20,18))
destination<-file.path(root,"G0022_RS671_KCPS2_KOREAN_ALCOHOL_BP_FOREST.png")
ggsave(destination,p,width=11,height=5,dpi=165,bg="white")
cat("FIGURE_KCPS2_KOREAN_RS671_PASS",destination,"bytes",file.info(destination)$size,"\n")
