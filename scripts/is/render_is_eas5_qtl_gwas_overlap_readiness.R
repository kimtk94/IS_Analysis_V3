#!/usr/bin/env Rscript
# Static scientific QC: public nominal QTL source vs exact GRCh37
# non-palindromic variant overlap in both BBJ IS and EAS AIS.
# No colocalization or causal gene labels in this figure.
cmd<-commandArgs(trailingOnly=TRUE)
arg<-function(key) {
 ix<-which(cmd==key);if(length(ix)!=1||ix==length(cmd))stop(paste("Missing",key))
 cmd[[ix+1]]
}
p1<-arg("--source-gene-readiness")
p2<-arg("--allele-gene-readiness")
out<-arg("--out")
if(file.exists(out))stop("Refuse overwrite of scientific figure")
all<-read.delim(p1,check.names=FALSE,stringsAsFactors=FALSE)
al<-read.delim(p2,check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(all)!=2225)stop("All original 2225 gene IDs required")
regions<-c("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007",
           "IS_XDATA_G0025","IS_XDATA_G0018")
orig<-c(71,14,10,92,53)
sourced<-c(46,9,4,61,29)
M<-matrix(0,nrow=4,ncol=5,dimnames=list(
 c("Both GWAS >=100 nonpal SNPs","Both GWAS 10-99 nonpal SNPs",
 "Both GWAS 1-9 nonpal SNPs","No both-GWAS nonpal SNPs"),regions))
for(i in seq_along(regions)) {
 rid<-regions[i]
 group<-all[all$region_id==rid,]
 sourced_genes<-group$gene_id_stable[group$n_new_QTL_tissues_with_any_nominal_rows>0]
 if(nrow(group)!=orig[i]||length(sourced_genes)!=sourced[i])
  stop("Original 149 nominal-source gene denominator changed")
 aln<-al[al$region_id==rid,]
 mx<-setNames(rep(0,length(sourced_genes)),sourced_genes)
 for(j in seq_len(nrow(aln))) {
  gene<-aln$gene_id[j]
  if(!gene %in% names(mx))stop("QTL/GWAS matched gene missing in source-union")
  mx[gene]<-max(mx[gene],aln$exact_NONPALINDROMIC_SNP_in_BOTH_BBJ_AIS[j])
 }
 M[,i]<-c(sum(mx>=100),sum(mx>=10&mx<100),sum(mx>0&mx<10),sum(mx==0))
 if(sum(M[,i])!=sourced[i])stop("Locus denominator mismatch")
}
cols<-c("#256B8A","#509A8C","#D4A354","#CDD5D9")
svg(out,width=11.4,height=7.2,pointsize=11)
par(mar=c(5,6.1,5,2),family="sans")
bp<-barplot(M,col=cols,beside=FALSE,border=NA,
 names.arg=sub("IS_XDATA_","",regions),
 ylab="Distinct source-tested positional genes",
 xlab="Provisional GRCh37 IS cross-study genomic group",
 main="Public nominal QTL gene coverage vs allele-exact GWAS variant overlap",
 ylim=c(0,1.14*max(colSums(M))),cex.names=.9)
text(bp,colSums(M)+2,labels=sourced,cex=.95,font=2,col="#333333")
legend("topright",legend=rownames(M),fill=cols,bty="n",cex=.86)
mtext("Source: 3 public nominal QTL datasets, 5 regions | only original QTL gene-SNP tests",
 side=3,line=.35,cex=.8,col="#555555")
mtext("Variant overlap != colocalization; molecular-QTL LD, case/control N and independence still under review.",
 side=1,line=3.8,cex=.82,col="#555555")
dev.off()
cat("IS_5REGION_GWAS_QTL_SOURCE_EXACT_OVERLAP_FIGURE_PASS",out,"\n")
