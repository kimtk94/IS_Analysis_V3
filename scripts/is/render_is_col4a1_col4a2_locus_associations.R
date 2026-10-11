#!/usr/bin/env Rscript
# GRCh37 locus chart: tag-SNP associations, no LD or shared causal variant.
args<-commandArgs(trailingOnly=TRUE)
v<-function(a){i<-which(args==a);if(length(i)!=1||i==length(args))stop(paste("Missing",a));args[[i+1L]]}
f<-v("--input");out<-v("--out")
if(file.exists(out))stop("Refuse plot overwrite")
d<-read.delim(f,stringsAsFactors=FALSE,check.names=FALSE)
d<-d[d$QTL_dataset=="QTD000434",]
if(nrow(d)!=7892 || !setequal(unique(d$gene_symbol),c("COL4A1","COL4A2")))
 stop("ROSMAP COL4A1/COL4A2 matched 3,946 SNVs each expected")
tags<-c("13:110954237:A:G","13:110954510:T:C","13:111040798:A:C")
if(!all(tags%in%d$GRCh37_fasta_REF_ALT_variant))
 stop("Both local signal SNPs and established literature SNP must be present")
gwas<-d[d$gene_symbol=="COL4A2",]
positions<-as.integer(sub("^13:([0-9]+):.*$","\\1",gwas$GRCh37_fasta_REF_ALT_variant))
lo<-min(positions);hi<-max(positions)
lim<-c(lo,hi)
svg(out,width=11.7,height=8.2,pointsize=10.5)
par(mfrow=c(2,1),mar=c(3.5,5.2,3,1.1),oma=c(1,0,3.1,0),family="sans")
for(gene in c("COL4A1","COL4A2")){
 z<-d[d$gene_symbol==gene,,drop=FALSE]
 x<-as.integer(sub("^13:([0-9]+):.*$","\\1",z$GRCh37_fasta_REF_ALT_variant))/1e6
 qp<-pmax(as.numeric(z$QTL_p),1e-250)
 gp<-pmax(as.numeric(z$EAS_AIS_p),1e-250)
 ylim<-c(0,max(-log10(qp),-log10(gp),9)*1.10)
 plot(x,-log10(gp),pch=16,cex=.52,col=adjustcolor("#9BA9B5",alpha.f=.42),
  xlim=lim/1e6,ylim=ylim,
  xlab="GRCh37 chr13 position (Mb)",
  ylab=expression(-log[10](p)),main=paste(gene," — ROSMAP DLPFC eQTL vs EAS AIS"),
  cex.main=1.05)
 points(x,-log10(qp),pch=16,cex=.56,col=adjustcolor("#1E7F7C",alpha.f=.43))
 abline(h=-log10(5e-8),col="#8D98A3",lty=3)
 abline(h=-log10(1e-4),col="#1E7F7C",lty=3)
 for(tag in tags){
  q<-z[z$GRCh37_fasta_REF_ALT_variant==tag,]
  px<-as.numeric(sub("^13:([0-9]+):.*$","\\1",tag))/1e6
  clr<-if(tag=="13:111040798:A:C")"#8C4A77" else "#C45A45"
  abline(v=px,col=adjustcolor(clr,alpha.f=.55),lty=2)
  points(px,-log10(as.numeric(q$QTL_p)),pch=21,cex=1.6,bg=clr,col="white",lwd=.8)
  points(px,-log10(as.numeric(q$EAS_AIS_p)),pch=24,cex=1.4,bg=clr,col="white",lwd=.8)
 }
 legend("topright",bty="n",cex=.79,
  legend=c("GIGASTROKE EAS AIS","ROSMAP DLPFC cis-eQTL",
           "Our two same-SNP tag loci","Published rs9515201 (lacunar IS / SVD)"),
  col=c("#8C99A8","#1E7F7C","#C45A45","#8C4A77"),pch=c(16,16,21,21))
}
mtext("COL4A1 / COL4A2 paired full-cis genomic overlap (no LD-based conditioning)",side=3,outer=TRUE,line=1.2,cex=1.2,font=2)
mtext("rs9515201 is a published locus in primarily European lacunar-IS / ICH data, not established replication in East Asian AIS.",
 side=1,outer=TRUE,line=-.8,cex=.76,col="#555555")
dev.off()
cat("COL4A1_COL4A2_ROSMAP_CIS_VS_EAS_GWAS_LOCUS_FIGURE_PASS",out,"\n")
