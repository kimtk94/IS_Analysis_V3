#!/usr/bin/env Rscript
# Method sensitivity: GWAS/external 1KG EAS LD mismatch is a quality gate,
# exploratory SuSiE CS counts must never be equated to proven causal signals.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>0)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
out<-if(length(args)>1)args[[2]] else
 "/srv/is-analysis/worktrees/is-broad-discovery-20261009/docs/figures/is"
d<-read.delim(file.path(root,"ALCOHOL_ADH1B_ALDH2_GWAS_LD_MISMATCH_FILTER_SENSITIVITY.tsv"))
if(nrow(d)!=8)stop("Expected two regions x four SNP filter scenarios")
d$locus<-factor(d$locus,levels=c("ADH1B","ALDH2"))
levels<-c("SOURCE_QC_DEFAULT","REF_MAF_10_EAFDELTA_05",
 "REF_MAF_10_EAFDELTA_03","REF_MAF_05_EAFDELTA_03")
d$scenario<-factor(d$scenario,levels=levels,
 labels=c("Default MAF 5%, AF Δ10pp","MAF 10%, AF Δ5pp",
  "MAF 10%, AF Δ3pp","MAF 5%, AF Δ3pp"))
d$label<-sprintf("s=%.3f  (n=%d)",d$LD_mismatch_s,d$variant_n)
a<-ggplot(d,aes(x=LD_mismatch_s,y=scenario,color=locus))+
 geom_vline(xintercept=.05,linetype="dashed",color="#ca744e",linewidth=.7)+
 geom_point(size=3.2)+geom_text(aes(label=label),hjust=-.15,size=3.1,color="#384758")+
 facet_wrap(~locus,ncol=1,scales="free_y")+
 scale_x_continuous(limits=c(0,.154),breaks=c(0,.025,.05,.075,.10,.125,.15))+
 scale_color_manual(values=c("#286a8a","#b67b45"))+
 labs(title="Japanese alcohol GWAS vs external EAS LD: sensitivity",
 subtitle="SuSiE-RSS estimate_s_rss | 1000 Genomes EAS 504 participants",
 x="Estimated GWAS–LD mismatch s (lower is better)",y=NULL,
 caption=paste("Dashed s=.05 is an internal conservative research stop rule, NOT a universal methodological cutoff.",
 "Stricter SNP filters can lower s by deleting major parts of the locus; they cannot prove full-locus LD compatibility.",
 "No causal gene, stroke effect or mediation identified.",sep="\n"))+
 theme_minimal(base_size=11)+theme(legend.position="none",
 strip.text=element_text(face="bold"),
 panel.grid.major.y=element_blank(),plot.title=element_text(face="bold",size=14),
 plot.caption=element_text(hjust=0,size=8,color="#53606a"))
dir.create(out,showWarnings=FALSE,recursive=TRUE)
png1<-file.path(out,"IS_ALCOHOL_ALDH2_ADH1B_GWAS_LD_MISMATCH_SENSITIVITY.png")
ggsave(png1,a,width=12,height=7,units="in",dpi=160,bg="white")
m<-read.delim(file.path(root,"ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.tsv"))
if(nrow(m)!=6)stop("Expected six ADH1B SuSiE sensitivity models")
m$filter<-factor(m$scenario,levels=c("BASELINE_MAF05_DELTA10",
 "RESTRICTED_MAF10_DELTA03"),labels=c("637 SNPs, default QC",
 "261 SNPs, stringent AF matching"))
b<-ggplot(m,aes(x=factor(L_maxcausal_assumed),y=exploratory_95pct_credible_sets,fill=filter))+
 geom_col(position=position_dodge(.7),width=.6)+
 geom_text(aes(label=exploratory_95pct_credible_sets),
  position=position_dodge(.7),vjust=-.35,size=4)+
 scale_fill_manual(values=c("#286a8a","#b67b45"))+
 scale_y_continuous(limits=c(0,7),breaks=0:7)+
 labs(title="ADH1B SuSiE credible sets are model-dependent",
 subtitle="All six exploratory models converged; external 1KG EAS n=504",
 x="SuSiE maximum effects parameter L",
 y="Exploratory 95% credible sets",fill=NULL,
 caption=paste("Counts change with SNP filtering and L; top-PIP variant also changes.",
 "These are NOT independent causal variants and do not establish ADH1B-to-stroke mediation.",sep="\n"))+
 theme_minimal(base_size=11)+theme(legend.position="top",
 plot.title=element_text(face="bold",size=14),
 plot.caption=element_text(hjust=0,size=8,color="#53606a"))
png2<-file.path(out,"IS_ALCOHOL_ADH1B_SUSIE_MODEL_STABILITY.png")
ggsave(png2,b,width=10,height=5.8,units="in",dpi=160,bg="white")
cat("ALCOHOL_CONDITIONAL_FIGURES_PASS",png1,round(file.info(png1)$size/1024),"KiB;",
    png2,round(file.info(png2)$size/1024),"KiB\n")
