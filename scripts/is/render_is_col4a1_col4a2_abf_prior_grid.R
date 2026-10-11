#!/usr/bin/env Rscript
# Paired COL4A1-v-COL4A2 single-causal model p12 posterior heatmap;
# no shared causal inference or independent GWAS cohort replication implied.
args<-commandArgs(trailingOnly=TRUE)
opt<-function(k){ i<-which(args==k);if(length(i)!=1||i==length(args))stop(paste("Missing",k));args[[i+1]] }
x<-read.delim(opt("--input"),check.names=FALSE,stringsAsFactors=FALSE)
out<-opt("--out")
if(file.exists(out))stop("Refuse figure overwrite")
if(nrow(x)!=270L || sum(is.finite(x$PP_H4))!=270L)stop("Incomplete ABF result grid")
p12s<-c(1e-6,3e-6,1e-5,3e-5,1e-4)
ids<-rbind(
 c("QTD000131","COL4A1"),c("QTD000131","COL4A2"),
 c("QTD000171","COL4A1"),c("QTD000171","COL4A2"),
 c("QTD000434","COL4A1"),c("QTD000434","COL4A2"))
labels<-c("GTEx Aorta   COL4A1","GTEx Aorta   COL4A2",
          "GTEx Cortex  COL4A1","GTEx Cortex  COL4A2",
          "ROSMAP DLPFC COL4A1","ROSMAP DLPFC COL4A2")
m<-matrix(NA_real_,nrow=6,ncol=length(p12s))
mH3<-m
for(i in seq_len(6)) for(j in seq_along(p12s)){
 z<-x[x$QTL_dataset==ids[i,1] & x$gene_symbol==ids[i,2] &
   x$SNP_selection=="ALL_NONPAL_PAIRED" &
   abs(x$QTL_N_multiplier-1)<1e-9 &
   abs(x$p12-p12s[j])<1e-15,,drop=FALSE]
 if(nrow(z)!=1)stop("Missing exact primary posterior record")
 m[i,j]<-z$PP_H4
 mH3[i,j]<-z$PP_H3
}
if(anyNA(m)||anyNA(mH3)||any(m<0|m>1))stop("Invalid posterior map")
svg(out,width=12.4,height=7.1,pointsize=11)
par(mar=c(6,3,4.3,1),family="sans")
plot.new()
plot.window(xlim=c(-.7,5.45),ylim=c(-.7,6.45),xaxs="i",yaxs="i")
pal<-colorRampPalette(c("#F0F2F4","#D3E5E8","#82B9BC","#217E86","#1B4B5A"))(101)
for(i in seq_len(6)){
 y<-6-i
 for(j in seq_along(p12s)){
  val<-m[i,j]
  rect(j-.5,y-.43,j+.45,y+.41,col=pal[[min(101,max(1,floor(val*100)+1))]],
       border="white",lwd=3)
  fg<-if(val>.65)"white" else "#253D45"
  text(j-.025,y+.115,sprintf("H4 %.3f",val),col=fg,cex=.87,font=2)
  text(j-.025,y-.18,sprintf("H3 %.3f",mH3[i,j]),col=fg,cex=.77)
 }
}
axis(1,at=1:5,labels=c("1e-6","3e-6","1e-5","3e-5","1e-4"),tick=FALSE,line=.45,cex.axis=.95)
axis(2,at=5:0,labels=labels,las=1,tick=FALSE,line=-.25,cex.axis=.91)
box(bty="n")
title("COL4A1 vs COL4A2: shared-causal H4 and distinct-signal H3",
 cex.main=1.25,line=2)
mtext("coloc.abf single-signal model · 3 QTL cohorts × 2 genes · same common SNPs",side=3,line=.35,cex=.87,col="#586168")
mtext("Prior p12 (p1=p2=1e-4); all nonpalindromic SNPs, nominal QTL source N, EAS AIS sample case fraction external.",
 side=1,line=3,cex=.86)
mtext("Exploratory / ancestry-unmatched eQTL and AIS; no verified multi-signal conditioning or causal transcript assignment.",
 side=1,line=4.35,cex=.80,col="#666666")
dev.off()
cat("COL4A1_COL4A2_SINGLE_SIGNAL_H3_H4_PRIOR_SENSITIVITY_FIGURE_PASS",out,"\n")
