#!/usr/bin/env Rscript
# Genome-wide support and untested QTL footprint in 30 PROVISIONAL EAS regions.
# This plot is a genomic interval triage, not causal gene prioritization.
arg<-commandArgs(trailingOnly=TRUE)
getarg<-function(k) {
 i<-which(arg==k);if(length(i)!=1||i==length(arg))stop(paste("Missing",k))
 arg[[i+1]]
}
src<-getarg("--queue");out<-getarg("--out")
if(file.exists(out))stop("Refuse figure overwrite")
d<-read.delim(src,check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(d)!=30||sum(d$n_EAS_positional_genes_without_archived_coloc)!=826||
   sum(d$n_archived_QTL_tested_gene_ids)!=43)stop("30 region original QTL gate changed")
tier<-substring(d$primary_GWAS_EAS_BBJ_AIS_tier,1,2)
if(!identical(as.integer(table(factor(tier,levels=c("T1","T2","T3","T4")))),
              c(4L,1L,2L,23L)))stop("Original GWAS support tier changed")
cols<-c(T1="#206D72",T2="#B77837",T3="#4F81A5",T4="#ACB6BE")
y<-rev(seq_len(nrow(d)))
lab<-sub("^IS_XDATA_","",d$region_id)
total<-d$n_positional_gene_ids
tested<-d$n_archived_QTL_tested_gene_ids
missing<-d$n_EAS_positional_genes_without_archived_coloc
svg(out,width=13,height=12,pointsize=10)
layout(matrix(c(1,2),1,2),widths=c(1.18,1.4))
par(mar=c(5.1,5.3,4.3,1.0),family="sans")
plot(NA,xlim=c(0,max(total)*1.16),ylim=c(0,31),axes=FALSE,xlab="",ylab="")
abline(v=c(0,25,50,75,100),col="#ECECEC")
rect(0,y-.31,missing,y+.31,col=cols[tier],border=NA)
rect(missing,y-.31,total,y+.31,col="#DDDDDD",border=NA)
axis(1,at=c(0,25,50,75,100))
axis(2,at=y,labels=lab,las=1,tick=FALSE,cex.axis=.8)
mtext("Number of positional genes",side=1,line=3.3)
title("QTL coverage gaps",line=1.5,cex.main=1.15)
legend("bottomright",legend=c("Untested positional genes","Archived QTL-tested genes"),
       fill=c(cols["T1"],"#DDDDDD"),bty="n",cex=.75)
par(mar=c(5.1,1.1,4.3,2.2))
to_log<-function(x) pmin(-log10(pmax(x,1e-300)),20)
maxp<-20
plot(NA,xlim=c(0,maxp),ylim=c(0,31),axes=FALSE,xlab="",ylab="")
abline(v=c(0,5,10,15,20),col="#ECECEC")
abline(v=-log10(5e-8),lty=2,lwd=1.3,col="#707070")
points(to_log(d$BBJ_500kb_min_p),y+.14,pch=16,cex=.88,col="#2C6C99")
points(to_log(d$GWS_EAS_AIS_500kb_min_p),y-.14,pch=17,cex=.85,col="#BC7036")
axis(1,at=c(0,5,10,15,20))
mtext(expression(-log[10](p)),side=1,line=3.3)
title("GWAS signal in each ±500kb context",line=1.5,cex.main=1.07)
legend("topright",legend=c("BBJ Japanese IS","GIGASTROKE EAS AIS","P = 5e-8 threshold"),
       col=c("#2C6C99","#BC7036","#707070"),pch=c(16,17,NA),
       lty=c(NA,NA,2),bty="n",cex=.78)
mtext("T1: AIS + BBJ associated | T2: CES-specific | T3: AS-only | T4: exploratory",
      side=3,outer=TRUE,line=.4,cex=.85)
dev.off()
cat("IS_EAS30_REGION_COVERAGE_FIGURE_COMPLETE ",out,"\n")
