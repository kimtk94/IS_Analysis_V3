#!/usr/bin/env Rscript
# Three distinct input/QC states: historical nominal interval, full cis with
# historically misnamed GIGA variant_id interpreted strictly as REF:ALT,
# and full cis with GIGASTROKE OA:EA independently FASTA-anchored.
# These are coverage counts, not independent SNPs or colocalization results.
args<-commandArgs(trailingOnly=TRUE)
opt<-function(k){
 i<-which(args==k)
 if(length(i)!=1||i==length(args))stop(paste("Required",k))
 args[[i+1]]
}
a<-read.delim(opt("--strict"),stringsAsFactors=FALSE,check.names=FALSE)
b<-read.delim(opt("--anchored"),stringsAsFactors=FALSE,check.names=FALSE)
out<-opt("--out")
if(file.exists(out))stop("Refuse replacing an existing figure")
keys<-function(z)paste(z$region_id,z$QTL_dataset,z$gene_id,sep="|")
if(nrow(a)!=98||nrow(b)!=98||length(unique(a$gene_id))!=66||
   length(unique(b$gene_id))!=66||!setequal(keys(a),keys(b)))
 stop("Assays, gene set, or input universe drift")
a<-a[match(keys(b),keys(a)),]
if(any(a$OLD_2MB_BBJ_AIS_exact_nonpal_n!=b$OLD_2MB_BBJ_AIS_exact_nonpal_n)||
   any(as.numeric(b$NEW_fullcis_BOTH_nonpal_SNP_n)<as.numeric(a$NEW_fullcis_BOTH_nonpal_SNP_n)))
 stop("FASTA source allele repair must not discard ref-verified matches")
regions<-c("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0025","IS_XDATA_G0018")
byreg<-function(column,table)
 vapply(regions,function(r)sum(as.numeric(table[table$region_id==r,column])),numeric(1))
mat<-rbind(
 "Original 2 Mb public nominal QTL slice"=byreg("OLD_2MB_BBJ_AIS_exact_nonpal_n",a),
 "Full cis, strict historical variant ID"=byreg("NEW_fullcis_BOTH_nonpal_SNP_n",a),
 "Full cis, GRCh37 FASTA-anchored GWAS alleles"=byreg("NEW_fullcis_BOTH_nonpal_SNP_n",b))
if(sum(mat[1,])!=86965||sum(mat[2,])!=145267||sum(mat[3,])!=255753)
 stop("Scientific SNP overlap coverage invariant drift")
svg(out,width=12.3,height=7.4,pointsize=10.5)
par(mar=c(5.9,6.2,5.15,1.4),family="sans")
color<-c("#9EAFBA","#E0B16D","#237E85")
at<-barplot(mat,beside=TRUE,col=color,border=NA,ylim=c(0,max(mat)*1.16),
 axes=FALSE,names.arg=sub("IS_XDATA_","",regions),cex.names=.98)
axis(2,las=1,at=pretty(c(0,max(mat))))
abline(h=pretty(c(0,max(mat))),lty=3,col="#E7EAED")
barplot(mat,beside=TRUE,col=color,border=NA,axes=FALSE,add=TRUE)
for(i in seq_len(nrow(mat))) for(j in seq_len(ncol(mat))){
 val<-mat[i,j]
 if(val>0) text(at[i,j],val+max(mat)*.019,
   labels=if(val>1000) paste0(format(round(val/1000,1),trim=TRUE,nsmall=1),"k") else as.character(val),
   cex=.82,font=2)
}
title("FASTA-anchored GWAS allele recoding restores public cis-eQTL overlap",
      cex.main=1.14,line=1.65)
legend("topleft",fill=color,legend=rownames(mat),bty="n",cex=.83)
mtext("Gene × QTL dataset × non-palindromic SNP records shared with BBJ IS and EAS AIS GWAS",
      side=2,line=4.55,cex=.93)
mtext("Source: 66 eGenes / 98 study-QTL assays; GIGASTROKE raw MD5 + EAS build + FASTA checked.",
      side=3,line=.25,cex=.88,col="#555555")
mtext("No new colocalization posteriors; overlapping stroke cohorts and molecular-QTL LD remain unresolved.",
      side=1,line=4.6,cex=.81,col="#555555")
dev.off()
cat("IS_EAS5_HISTORICAL_OA_EA_TO_FASTA_REFALT_COVERAGE_FIGURE_PASS",out,"\n")
