#!/usr/bin/env Rscript
# Real 1000 Genomes Phase3 GRCh37 n504 EAS and n104 JPT ALT-dosage r^2
# Only within-chromosome sentinel pairs displayed (the other 18 include
# different chromosomes). This is an LD sensitivity graphic, NOT MR validity.
suppressPackageStartupMessages(library(ggplot2))
args<-commandArgs(trailingOnly=TRUE)
source<-if(length(args)>=1)args[[1]] else
  "/srv/is-analysis/data/is/ld_reference/nonaldh2_alcohol_eas_20261010"
dest<-if(length(args)>=2)args[[2]] else
  "/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1"
data<-read.delim(file.path(source,"IS_RS671_7SNP_EAS504_JPT104_LD_SENSITIVITY.tsv"),
  check.names=FALSE)
if(nrow(data)!=21||sum(data$same_chr==1)!=3)
  stop("Expected 21 complete LD pairs (3 same chromosome)")
x<-data[data$same_chr==1,]
x$pair<-c("chr4 KLB – ADH1B",
          "chr9 ALDH1B1 – ALDH1A1",
          "chr12 intergenic – ALDH2 rs671")
x$pair<-factor(x$pair,levels=rev(x$pair))
e<-data.frame(pair=x$pair,Population="East Asian panel (504)",
 r2=x$EAS_504_ALT_dosage_r2)
j<-data.frame(pair=x$pair,Population="Japanese JPT (104)",
 r2=x$JPT_104_ALT_dosage_r2)
plotdata<-rbind(e,j)
if(any(!is.finite(plotdata$r2)))stop("Unmeasured Japanese LD pair")
plotdata$Population<-factor(plotdata$Population,
  levels=c("East Asian panel (504)","Japanese JPT (104)"))
plotdata$annotation<-sprintf("%.4f",plotdata$r2)
pos<-position_dodge(width=.55)
p<-ggplot(plotdata,aes(x=r2,y=pair,color=Population))+
 geom_segment(aes(x=0,xend=r2,y=pair,yend=pair),position=pos,
              linewidth=.5,alpha=.45)+
 geom_point(size=3,position=pos)+
 geom_text(aes(label=annotation),position=pos,hjust=-.2,size=3.3,
           color="#243746")+
 scale_x_continuous(limits=c(0,.061),breaks=c(0,.01,.02,.03,.04,.05),
  labels=scales::label_number(accuracy=.01))+
 scale_color_manual(values=c("#23688f","#c07839"))+
 labs(title="Alcohol-associated sentinels: observed 1000G LD",
 subtitle="Exact GRCh37 allele-aligned EAS n=504 and JPT n=104 dosages",
 x="Observed genotype dosage r² between positional sentinels",y=NULL,color=NULL,
 caption=paste("Three same-chromosome comparisons shown; all 21 pairs checked.",
  "Reference n=104 JPT gives noisy LD for rare alleles (minor allele count as low as 5).",
  "Low pairwise LD does not demonstrate no horizontal pleiotropy or valid MR instruments.",
  sep="\n"))+
 theme_minimal(base_size=11)+
 theme(legend.position="top",panel.grid.major.y=element_blank(),
 plot.title=element_text(face="bold",size=15),
 plot.caption=element_text(hjust=0,color="#52606d",size=8.5),
 axis.text.y=element_text(size=10,color="#243746"),
 plot.margin=margin(15,30,15,12))
dir.create(dest,showWarnings=FALSE,recursive=TRUE)
output<-file.path(dest,"IS_ALCOHOL_7SNP_EAS_JPT_REAL_LD_COMPARISON.png")
ggsave(output,p,width=11.8,height=5.2,units="in",dpi=170,bg="white")
cat("LD_EAS_JPT_FIGURE_PASS",output,"bytes",file.info(output)$size,"\n")
