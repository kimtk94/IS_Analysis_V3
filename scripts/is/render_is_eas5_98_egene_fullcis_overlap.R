#!/usr/bin/env Rscript
# Public eGene source SNP overlap before vs after full cis QTL extraction.
# Technical coverage only; not coloc and not independent BBJ/EAS replication.
args<-commandArgs(trailingOnly=TRUE)
opt<-function(k){
  i<-which(args==k)
  if(length(i)!=1||i==length(args))stop(paste("argument required",k))
  args[[i+1]]
}
f<-opt("--input");out<-opt("--out")
if(file.exists(out))stop("Refuse replacing existing figure")
x<-read.delim(f,stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(x)!=98||length(unique(x$gene_id))!=66)stop("Expected 98 source BH eGene assays / 66 genes")
regions<-c("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007",
           "IS_XDATA_G0025","IS_XDATA_G0018")
old<-tapply(as.numeric(x$OLD_2MB_BBJ_AIS_exact_nonpal_n),x$region_id,sum)
new<-tapply(as.numeric(x$NEW_fullcis_BOTH_nonpal_SNP_n),x$region_id,sum)
values<-rbind(
 "Original narrow QTL region"=vapply(regions,function(s)if(is.na(old[s]))0 else old[s],numeric(1)),
 "Rebuilt full cis QTL tested region"=vapply(regions,function(s)if(is.na(new[s]))0 else new[s],numeric(1)))
if(any(values[2,]<values[1,]))stop("Expanded verified subset must preserve original match counts")
if(sum(values)<=0)stop("No data")
svg(out,width=11.75,height=7.5,pointsize=11)
par(mar=c(5.6,5.5,4.4,1.5),family="sans")
cols<-c("#859AA8","#1D7E86")
mat<-barplot(values,beside=TRUE,col=cols,border=NA,
 names.arg=sub("IS_XDATA_","",regions),cex.names=.95,
 ylab="Gene × QTL dataset × SNP exact matches (both GWAS)",
 xlab="Provisional EAS ischemic-stroke region",ylim=c(0,max(values)*1.18),
 main="Exact GRCh37 BBJ + EAS AIS SNP overlap in public FDR-significant eGenes")
abline(h=pretty(c(0,max(values))),lty=3,col="#E2E7E9")
barplot(values,beside=TRUE,col=cols,border=NA,add=TRUE,axes=FALSE,ann=FALSE)
labels<-ifelse(values>1000,format(round(values/1000,1),trim=TRUE,nsmall=1),values)
labels<-ifelse(values>1000,paste0(labels,"k"),as.character(labels))
text(mat,values+max(values)*.018,labels=labels,cex=.82)
legend("topleft",fill=cols,legend=rownames(values),bty="n",cex=.84)
mtext("Data-unit duplicates across genes/datasets possible; no molecular QTL donor LD or coloc posteriors assessed.",
 side=1,line=4.2,cex=.78,col="#586069")
dev.off()
cat("IS_EAS5_FULLCIS_BBJ_GIGA_QTL_SNP_OVERLAP_COMPARISON_FIGURE_PASS",out,"\n")
