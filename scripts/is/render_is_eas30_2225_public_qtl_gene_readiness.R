#!/usr/bin/env Rscript
# Positional gene molecular-QTL readiness, not an inferred causal-gene heatmap.
arg<-commandArgs(trailingOnly=TRUE)
getarg<-function(k) {
 i<-which(arg==k);if(length(i)!=1||i>=length(arg))stop(paste("Missing",k))
 arg[[i+1]]
}
inp<-getarg("--gene-readiness");queue<-getarg("--queue");out<-getarg("--out")
if(file.exists(out))stop("No plot overwrite")
g<-read.delim(inp,check.names=FALSE,stringsAsFactors=FALSE)
q<-read.delim(queue,check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(g)!=2225||nrow(q)!=30)stop("Original EAS positional universe incomplete")
eas<-g[nzchar(g$region_id),]
if(nrow(eas)!=869)stop("Expected 869 East Asian positional genes")
old<-eas$archived_646_coloc_assays>0
new<-eas$n_new_QTL_tissues_with_any_nominal_rows>0 & !old
not_new<-!old & !new
if(sum(old)!=43||sum(new)!=124||sum(not_new)!=702)stop("Input source state drift")
regions<-q$region_id
n<-length(regions);values<-matrix(0,nrow=3,ncol=n)
rownames(values)<-c("Archived 646-assay ABF tested",
                    "New raw nominal QTL, no coloc",
                    "No archived ABF or pilot source data")
for(i in seq_along(regions)){
 x<-eas[eas$region_id==regions[i],]
 values[,i]<-c(sum(x$archived_646_coloc_assays>0),
               sum(x$n_new_QTL_tissues_with_any_nominal_rows>0 & x$archived_646_coloc_assays==0),
               sum(x$n_new_QTL_tissues_with_any_nominal_rows==0 & x$archived_646_coloc_assays==0))
 if(sum(values[,i])!=q$n_positional_gene_ids[i])stop("Per-region source states do not partition genes")
}
svg(out,width=11.6,height=12.3,pointsize=11)
par(mar=c(4.8,7.2,5.3,1.9),family="sans")
bp<-barplot(values,horiz=TRUE,beside=FALSE,las=1,
  names.arg=sub("^IS_XDATA_","",regions),cex.names=.77,
  col=c("#36759A","#32927B","#D5DADF"),border=NA,
  xlim=c(0,max(colSums(values))*1.18),
  xlab="Distinct positional genes in each provisional EAS GWAS interval",
  main="IS QTL evidence-readiness: 43 historic, 124 newly sourced, 702 unassessed")
axis(1,at=c(0,20,40,60,80,100))
for(i in seq_len(n)){
 text(sum(values[,i])+2,bp[i],labels=q$n_positional_gene_ids[i],adj=0,cex=.76,col="#444444")
}
legend("topright",inset=c(.01,.01),bty="n",cex=.83,
       fill=c("#36759A","#32927B","#D5DADF"),
       legend=c("Archived BBJ/GTEx ABF-tested (43)",
                "Public nominal QTL rows; GWAS matching pending (124)",
                "Not examined in these QTL sources (702)"))
mtext("30 positional groups are NOT independent causal loci; missing QTL rows are NOT negative results.",
      side=1,line=3.6,cex=.8,col="#555555")
dev.off()
cat("IS_EAS30_2225_PUBLIC_QTL_READINESS_FIGURE_PASS ",out,"\n")
