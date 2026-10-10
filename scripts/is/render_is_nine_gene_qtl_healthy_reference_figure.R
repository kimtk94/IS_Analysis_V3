#!/usr/bin/env Rscript
# Descriptive figure: nine preselected historical ischemic stroke candidates.
# Same archived best-baseline gene/tissue assay at three p12 values;
# independent healthy adult brain reference localization is side metadata.
# The two evidence domains must NEVER be combined into a causal score.
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(k) {
 i<-which(args==k)
 if(length(i)!=1 || i>=length(args)) stop(paste("Missing",k))
 args[[i+1]]
}
p <- getarg("--nine-gene-table")
outdir <- getarg("--out-dir")
if(!file.exists(p)) stop("Nine-gene matrix missing")
if(dir.exists(outdir) && length(list.files(outdir,all.files=TRUE,no..=TRUE))>0)
 stop("Do not overwrite previous figures")
d<-read.delim(p,check.names=FALSE,stringsAsFactors=FALSE)
expected<-c("FGF5","ALDH2","SH3PXD2A","COL4A2","COL4A1",
            "CALHM2","NEURL1","C4orf22","INA")
if(nrow(d)!=9 || !setequal(d$historic_symbol,expected))
 stop("Nine-gene scientific comparisons changed")
if(any(!is.finite(d$best_H4_original_p12_1e5)) ||
   any(!is.finite(d$same_assay_H4_strict_p12_1e6)) ||
   any(!is.finite(d$same_assay_H4_high_p12_1e4)))
 stop("Prior sensitivity data incomplete")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d<-d[order(d$best_H4_original_p12_1e5,decreasing=TRUE),]
svg(file.path(outdir,"IS_NINE_GENE_PRIOR_H4_HEALTHY_CONTROL_REFERENCE.svg"),
    width=13,height=6.9,pointsize=11)
par(mar=c(5.5,8,4.7,2.3),family="sans")
n<-nrow(d)
top<-n:1
plot(NA,xlim=c(0,1.72),ylim=c(.4,n+.55),axes=FALSE,
     xlab="",ylab="",main="")
axis(1,at=seq(0,1,.2),labels=seq(0,1,.2))
axis(2,at=top,labels=d$historic_symbol,las=1,tick=FALSE)
abline(v=seq(0,1,.2),col="#e8e8e8")
abline(v=.8,col="#909090",lty=3,lwd=1)
segments(d$same_assay_H4_strict_p12_1e6,top,
         d$same_assay_H4_high_p12_1e4,top,
         col="#a7a7a7",lwd=2)
points(d$same_assay_H4_strict_p12_1e6,top,pch=16,cex=.9,col="#617A99")
points(d$best_H4_original_p12_1e5,top,pch=18,cex=1.2,col="#176957")
points(d$same_assay_H4_high_p12_1e4,top,pch=17,cex=1,col="#BB7836")
mtext("PP.H4 at the same preselected gene-tissue assay across p12 scenarios",
      side=1,line=3.1,cex=.89)
title(main="IS candidate interpretation: coloc prior sensitivity and healthy-brain context",
      cex.main=1.04,line=2.3)
legend("topleft",
       legend=c("p12 = 1e-6 (strict)","p12 = 1e-5 (archived baseline)","p12 = 1e-4 (high)"),
       pch=c(16,18,17),col=c("#617A99","#176957","#BB7836"),
       ncol=3,cex=.8,bty="n",inset=c(0,-.015),xpd=TRUE)
text(1.08,n+.46,"Healthy adult brain reference: top annotated cell type",
     adj=0,font=2,cex=.79)
for (i in seq_len(n)) {
 label<-d$healthy_reference_top_celltype[[i]]
 if(label=="" || label=="NOT_ASSESSABLE") label="Not assessable (feature absent)"
 text(1.08,top[[i]],label,adj=0,cex=.73,col=if(label=="Not assessable (feature absent)")"#777777" else "#202020")
}
mtext("Six healthy donors; tissue localization is descriptive, not ischemic-stroke DE or gene causality.",
      side=1,line=4.7,cex=.76,col="#555555")
dev.off()
cat("IS_NINE_GENE_FIGURE_COMPLETE",
    file.path(outdir,"IS_NINE_GENE_PRIOR_H4_HEALTHY_CONTROL_REFERENCE.svg"),"\n")
