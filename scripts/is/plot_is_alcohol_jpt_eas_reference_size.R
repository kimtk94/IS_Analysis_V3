#!/usr/bin/env Rscript
# Original EAS504 versus Japanese JPT104 versus 12 random EAS104 panels.
# Displays reference sample-size / LD mismatch s as QC only.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
root<-if(length(args)>=1)args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1"
out<-if(length(args)>=2)args[[2]] else
 "/srv/is-analysis/worktrees/is-broad-discovery-20261009/docs/figures/is"
x<-read.delim(file.path(root,"ALCOHOL_EAS104_RANDOM_DOWNSAMPLING_LD_MISMATCH.tsv"))
y<-read.delim(file.path(root,"ALCOHOL_EAS504_JPT104_FIXED_IDENTICAL_SNP_LD_MISMATCH.tsv"))
if(nrow(x)!=24||nrow(y)!=4)stop("Original source diagnostic row count mismatch")
if(nrow(y)!=4||any(y$source_snp_count!=y$used_snp_count))
 stop("Exact same-SNP source LD comparison missing")
category<-c("EAS504"="All EAS (n=504)",
 "JPT104"="Japanese JPT (n=104)","RANDOM_EAS104"="Random EAS (n=104)")
full<-rbind(
 data.frame(locus=x$locus,reference="RANDOM_EAS104",s=x$diagnostic_s),
 data.frame(locus=y$locus,reference=y$reference,s=y$LD_mismatch_s))
full$reference<-factor(category[full$reference],
 levels=c("All EAS (n=504)","Japanese JPT (n=104)","Random EAS (n=104)"))
if(any(!is.finite(full$s)))stop("Missing LD mismatch values")
med<-aggregate(s~locus+reference, data=full[full$reference=="Random EAS (n=104)",],FUN=median)
p<-ggplot(full,aes(x=reference,y=s,fill=reference))+
 geom_jitter(data=full[full$reference=="Random EAS (n=104)",],
    width=.10,height=0,size=2.3,shape=21,alpha=.85)+
 geom_point(data=full[full$reference!="Random EAS (n=104)",],
    size=3.6,shape=23,color="#344657")+
 geom_point(data=med,aes(x=reference,y=s),inherit.aes=FALSE,
    size=7,shape=95,color="#3b4e62")+
 facet_wrap(~locus,ncol=2)+
 scale_fill_manual(values=c("#6fa1c2","#deb080","#86a99a"))+
 scale_y_continuous(limits=c(0,.33),breaks=seq(0,.3,.05))+
 labs(title="LD reference size matters for Japanese alcohol GWAS consistency",
 subtitle="Identical source β/SE and SNP sets; EAS n=504, JPT n=104, 12 fixed-seed random EAS draws n=104",
 x=NULL,y="GWAS–LD mismatch s (smaller is better)",
 caption=paste(
  "Points: source EAS504 and JPT104 estimates; 12 random EAS104 repeats (black bar: median).",
  "Repeat panels overlap; this is an exploratory reference-size sensitivity experiment, not a causal test.",
  "Reference sample size/ancestry and study GWAS LD mismatch remain unresolved, especially at ALDH2.",
  sep="\n"))+
 theme_minimal(base_size=11)+
 theme(legend.position="none",panel.grid.major.x=element_blank(),
 axis.text.x=element_text(angle=18,hjust=1,size=9),
 strip.text=element_text(face="bold",size=12),
 plot.title=element_text(face="bold",size=14),
 plot.caption=element_text(hjust=0,size=8,color="#657080"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
png<-file.path(out,"IS_ALCOHOL_JPT_EAS_LD_REFERENCE_SIZE_DIAGNOSTIC.png")
ggsave(png,p,width=12.2,height=6.2,dpi=170,bg="white")
cat("LD_REF_SIZE_FIGURE_PASS",png,file.info(png)$size,"bytes\n")
