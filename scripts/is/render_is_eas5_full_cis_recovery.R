#!/usr/bin/env Rscript
# Variant-count completeness vs public EBI permutation gene-specific cis test universes.
# "Complete" = equal counts PLUS separately verified permutation-lead identity.
# Not a coloc/fine-mapping result or a genetic causal assignment.
cmd<-commandArgs(trailingOnly=TRUE)
getopt<-function(k) {
 i<-which(cmd==k);if(length(i)!=1||i==length(cmd))stop(paste("Missing",k))
 cmd[[i+1]]
}
inp<-getopt("--merged");out<-getopt("--out")
if(file.exists(out))stop("Refuse overwriting figure")
d<-read.delim(inp,check.names=FALSE,stringsAsFactors=FALSE)
if(nrow(d)!=720)stop("Expected full five-region 720 gene-dataset matrix")
regions<-c("IS_XDATA_G0015","IS_XDATA_G0023","IS_XDATA_G0007",
           "IS_XDATA_G0025","IS_XDATA_G0018")
state<-matrix(0,nrow=3,ncol=10)
rownames(state)<-c("Source QTL cis variant counts recovered",
                   "Regional slice truncated",
                   "No EBI permutation gene test")
labels<-character(10)
for(i in seq_along(regions)){
 group<-d[d$region_id==regions[i],]
 for(j in seq_len(2)) {
  k<-2*(i-1)+j
  complete<-if(j==1)
   sum(group$original_cis_completeness=="FULL_CIS_EQUIVALENCE_VS_EBI_NVARIANTS_COUNT_ONLY")
  else sum(group$extended_cis_completeness=="MATCHES_PERMUTED_FULL_CIS_VARIANT_COUNT")
  missing<-if(j==1)
   sum(group$original_cis_completeness=="NO_EBI_PERMUTED_GENE_RECORD_CANNOT_ASSESS")
  else sum(group$extended_cis_completeness=="NO_PERMUTED_GENE_TEST_ON_OFFICIAL_SOURCE")
  incomplete<-nrow(group)-complete-missing
  state[,k]<-c(complete,incomplete,missing)
  labels[k]<-if(j==1)"Initial" else "Extended"
 }
}
if(sum(state[1,seq(1,10,by=2)])!=10||
   sum(state[1,seq(2,10,by=2)])!=393||
   sum(state[3,seq(1,10,by=2)])!=327||
   sum(state[3,seq(2,10,by=2)])!=327)
 stop("Case-study cis universe state altered")
svg(out,width=12.7,height=7.45,pointsize=11)
par(mar=c(6.3,4.7,4.7,1.2),family="sans")
cols<-c("#287E72","#DCA251","#D3D9E0")
at<-barplot(state,beside=FALSE,space=c(.12,.8),ylim=c(0,320),
 col=cols,border=NA,axes=FALSE,names.arg=labels,cex.names=.9,
 ylab="Gene × QTL dataset comparisons")
axis(2,las=1,at=seq(0,300,50))
abline(h=seq(50,300,50),col="#E8EAEC",lty=3)
for(i in seq_along(at)) {
 q=sum(state[,i])
 text(at[i],q+10,labels=as.character(q),font=2,cex=.85)
 if(state[1,i]>0)text(at[i],state[1,i]/2,labels=state[1,i],
   cex=.84,col="white",font=2)
}
for(i in seq_along(regions)){
 x<-mean(at[c(2*i-1,2*i)])
 mtext(sub("IS_XDATA_","",regions[i]),side=1,at=x,line=3.8,font=2,cex=.9)
}
title("Full cis-tested QTL source recovery in five provisional IS intervals",
      cex.main=1.18,line=1.5)
legend("topright",bty="n",cex=.88,fill=cols,legend=rownames(state))
mtext("Official eQTL Catalogue n_variants; matched counts + independent permutation lead-ID check",
      side=3,line=.2,cex=.85,col="#4D5559")
mtext("Not a shared causal-variant test; original 646 coloc results unchanged.",
      side=1,line=5.5,cex=.82,col="#5C6065")
dev.off()
cat("IS_EAS5_ALL_CIS_COUNT_COMPLETENESS_FIGURE_PASS ",out,"\n")
