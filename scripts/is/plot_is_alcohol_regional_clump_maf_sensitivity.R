#!/usr/bin/env Rscript
# Regional Koyanagi alcohol GWAS clump counts. Not independent signals.
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1) args[[1]] else
 "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_gws_regional_clumping"
out <- if(length(args)>=2) args[[2]] else
 "/srv/is-analysis/worktrees/is-broad-discovery-20261009/docs/figures/is"
src <- file.path(root,"IS_ALCOHOL_7REGION_CLUMP_SENSITIVITY.tsv")
x <- read.delim(src,check.names=FALSE)
if(nrow(x)!=21)stop("Expected exactly 7 regions x 3 MAF thresholds")
x$MAF <- factor(ifelse(x$MAF_cutoff==0,"No MAF filter",
                     ifelse(x$MAF_cutoff==.01,"MAF >= 1%","MAF >= 5%")),
 levels=c("No MAF filter","MAF >= 1%","MAF >= 5%"))
x$locus <- factor(x$locus,levels=rev(c("GCKR","KLB","ADH1B",
  "ALDH1B1","ALDH1A1","chr12_106Mb","ALDH2")))
x$clump_index_count <- as.integer(x$clump_index_count)
p<-ggplot(x,aes(x=clump_index_count,y=locus,fill=MAF))+
 geom_col(position=position_dodge(width=.72),width=.65)+
 geom_text(aes(label=clump_index_count),
           position=position_dodge(width=.72),hjust=-.2,size=3)+
 scale_fill_manual(values=c("#586f8d","#2e8796","#d39b56"))+
 scale_x_continuous(limits=c(0,67),breaks=c(0,10,20,30,40,50,60))+
 labs(title="EAS reference clumping is highly MAF-sensitive at ALDH2",
 subtitle="Japanese Koyanagi 2024 alcohol GWAS  |  1000 Genomes EAS n=504",
 x="PLINK2 clump count (regional positional clusters; not causal signals)",
 y=NULL,fill=NULL,
 caption=paste("7 predefined +/- 500 kb regions; clump p1<5e-8, p2<1e-4, r2=0.1, radius 1 Mb.",
 "Numerical p=0 was floored at 1e-300 (ties); some index SNPs unstable.",
 "Do not interpret 57 or 20 ALDH2 clumps as independent biological causal variants.",
 sep="\n"))+
 theme_minimal(base_size=11)+
 theme(legend.position="top",panel.grid.major.y=element_blank(),
 plot.title=element_text(face="bold",size=14),
 plot.caption=element_text(color="#5b6470",size=8,hjust=0),
 plot.margin=margin(15,20,15,10))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
png<-file.path(out,"IS_ALCOHOL_7REGION_EAS_CLUMP_MAF_SENSITIVITY.png")
ggsave(png,p,width=12.2,height=6,units="in",dpi=170,bg="white")
cat("REGIONAL_CLUMP_FIGURE_PASS",png,file.info(png)$size,"bytes\n")
