#!/usr/bin/env Rscript
# Frozen historical coloc.abf: baseline vs high-prior posterior.
# Visualization is sensitivity only, NOT gene causality evidence.
getopt<-function(k) {
 a<-commandArgs(trailingOnly=TRUE);i<-which(a==k)
 if(length(i)!=1||i==length(a))stop(paste("Missing",k))
 a[[i+1]]
}
inp<-getopt("--assays");out<-getopt("--out-dir")
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))>0)stop("No overwrite")
d<-read.delim(inp,check.names=FALSE)
if(nrow(d)!=646||sum(d$original_H4_gt_H3)!=40||
   sum(d$high_H4_gt_H3)!=565||sum(d$high_H4_ge_0p8)!=12)
 stop("Expected historical 646 input and posterior counts")
dir.create(out,showWarnings=FALSE,recursive=TRUE)
cols<-c("#31739A","#6964A0","#178E75","#C37C38")
loci<-c("BBJ_IS_L001","BBJ_IS_L002","BBJ_IS_L003","BBJ_IS_L004")
cl<-adjustcolor(cols[match(d$locus,loci)],alpha.f=.28)
svg(file.path(out,"IS_646_ABF_H3_H4_BASELINE_VS_HIGHPRIOR.svg"),
    width=12.8,height=7.1,pointsize=11)
par(mfrow=c(1,2),mar=c(4.8,4.9,4.1,1),oma=c(2,0,1.6,0),family="sans")
for(which in c("original","high")){
 if(which=="original"){
  x<-d$original_p12_1e5_H3;y<-d$original_p12_1e5_H4
  tt<-"Original prior: p12=1e-5";count<-40;hi<-0
 }else{
  x<-d$high_p12_1e4_H3;y<-d$high_p12_1e4_H4
  tt<-"Sensitivity only: p12=1e-4";count<-565;hi<-12
 }
 plot(x,y,xlim=c(0,1),ylim=c(0,1),asp=1,pch=16,cex=.7,col=cl,
      xlab="PP.H3: distinct causal variants",
      ylab="PP.H4: shared causal variant",main=tt)
 abline(0,1,col="#686868",lty=2)
 abline(h=.8,col="#A94B42",lty=3)
 if(which=="high"){
  k<-which(d$high_H4_ge_0p8)
  points(x[k],y[k],pch=21,cex=1.12,bg="#FAEBC6",col="#A14F3C",lwd=1.2)
 }
 mtext(sprintf("H4 > H3: %d/646 | H4 >= 0.8: %d/646",count,hi),
       side=3,line=.4,cex=.85,font=2)
}
mtext("646 original gene-tissue assays: changing p12 does not add replication or confirm a causal gene",
      side=1,outer=TRUE,line=.7,cex=.82,col="#555555")
mtext("BBJ IS + GTEx V8 eQTL: sensitivity of coloc.abf posterior hypotheses",
      side=3,outer=TRUE,line=.6,font=2,cex=.95)
dev.off()
cat("IS_646_H3_H4_FIGURE_COMPLETE\n")
