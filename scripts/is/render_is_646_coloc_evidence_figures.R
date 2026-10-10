#!/usr/bin/env Rscript
# Reproducible publication-style descriptive figures, using only base R.
# Inputs are the full 646 direct ABF prior grid, eight fixed-priority assays
# and the fully retained 2,225-gene matrix. Zero causal inference is claimed.
args <- commandArgs(trailingOnly=TRUE)
getarg <- function(k) {
 i <- which(args==k)
 if(length(i)!=1 || i>=length(args))stop(paste("Missing",k))
 args[[i+1]]
}
grid <- read.delim(getarg("--prior-grid"),check.names=FALSE,stringsAsFactors=FALSE)
priority <- read.delim(getarg("--priority8"),check.names=FALSE,stringsAsFactors=FALSE)
genes <- read.delim(getarg("--genes"),check.names=FALSE,stringsAsFactors=FALSE)
outdir <- getarg("--out-dir")
expected <- c(1e-6,3e-6,1e-5,3e-5,1e-4)
if(nrow(grid)!=3230 || nrow(priority)!=8 || nrow(genes)!=2225)
 stop("Input evidence universe changed")
if(!setequal(unique(grid$p12),expected))
 stop("Five predeclared priors were altered")
if(dir.exists(outdir) && length(list.files(outdir,all.files=TRUE,no..=TRUE))>0)
 stop("Do not overwrite existing figures")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

ordered <- grid[order(grid$p12,grid$locus,grid$dataset_key,grid$gene_base),]
s <- lapply(expected,function(x) {
 d <- ordered[abs(ordered$p12-x)<1e-13,,drop=FALSE]
 if(nrow(d)!=646)stop("Nonunique/missing assays in prior grid")
 c(H4_GT_H3=sum(d$PP_H4>d$PP_H3),H4_GE_0P8=sum(d$PP_H4>=0.8),
   H4_GE_0P5=sum(d$PP_H4>=0.5))
})
summary <- do.call(rbind,s)
labels <- c("1e-6","3e-6","1e-5","3e-5","1e-4")
draw_prior_counts <- function(path) {
 svg(path,width=8.8,height=5.4,pointsize=11)
 on.exit(dev.off())
 par(mar=c(4.6,5.4,3.3,1),family="sans")
 barplot(t(summary[,c("H4_GT_H3","H4_GE_0P8")]),beside=TRUE,
         col=c("#5975A4","#C08047"),
         names.arg=labels,ylab="Gene-tissue coloc assays / 646",
         xlab="Prior p12 (p1 = p2 = 1e-4)",
         ylim=c(0,620),main="Ischemic stroke x GTEx: prior sensitivity")
 legend("topleft",legend=c("PP.H4 > PP.H3","PP.H4 >= 0.8"),
        fill=c("#5975A4","#C08047"),bty="n",cex=.84)
 mtext("Direct 3,230 SNP-level ABF models | exploratory, not causality",
       side=1,line=3.4,cex=.77,col="#555555")
}
draw_prior_counts(file.path(outdir,"IS_646_DIRECT_P12_THRESHOLD_COUNTS.svg"))

key <- function(d)paste(d$locus,d$dataset_key,d$gene_base,sep="::")
keys <- key(grid)
selected <- vector("list",8)
for(i in seq_len(8)){
 p <- priority[i,]
 m <- grid$locus==p$locus & grid$dataset_key==p$dataset_key
 # Known priority8 source row provides symbol but not gene ID. Choose by
 # matching baseline posterior PP.H4 to source, which fails closed on ties.
 m <- m & abs(grid$p12-1e-5)<1e-13 &
   abs(grid$PP_H4-as.numeric(p$pp_h4_baseline))<1e-8
 candidates <- unique(key(grid[m,]))
 if(length(candidates)!=1)stop(paste("Priority source not uniquely aligned:",p$gene))
 all <- grid[keys==candidates[[1]],]
 all <- all[order(all$p12),]
 if(nrow(all)!=5 || !setequal(all$p12,expected))stop("Incomplete priority assay curve")
 selected[[i]] <- list(gene=p$gene,source=candidates[[1]],h4=all$PP_H4)
}
draw_priority8 <- function(path){
 svg(path,width=10.3,height=6.2,pointsize=11)
 on.exit(dev.off())
 par(mar=c(4.7,5.0,3.4,2.4),family="sans")
 x <- log10(expected)
 cols <- grDevices::hcl.colors(8,"Dark 3")
 plot(NA,xlim=range(x),ylim=c(0,1),
      xaxt="n",xlab="Predeclared coloc prior p12",
      ylab="Posterior PP.H4 (shared-signal hypothesis)",
      main="Eight preselected gene-tissue assays: direct SNP ABF")
 axis(1,at=x,labels=labels)
 abline(h=0.8,lty=3,col="#707070")
 abline(v=log10(1e-5),lty=2,col="#777777")
 for(i in seq_len(8)) {
  lines(x,selected[[i]]$h4,col=cols[i],lwd=2,lty=1+(i>4))
  points(x,selected[[i]]$h4,col=cols[i],pch=15+i%%5,cex=.75)
 }
 legend("topleft",legend=vapply(selected,function(z)z$gene,character(1)),
   col=cols,lty=c(rep(1,4),rep(2,4)),lwd=2,
   pch=15+(1:8)%%5,ncol=2,cex=.82,bty="n")
 mtext("H4 = 0.8 reference line. Strong p12 dependence; matched QTL LD not verified.",
       side=1,line=3.4,cex=.77,col="#555555")
}
draw_priority8(file.path(outdir,"IS_PRIORITY8_DIRECT_P12_H4_SENSITIVITY.svg"))

draw_universe <- function(path) {
 assessed <- sum(as.integer(genes$new_DIRECT_prior_grid_assays)>0)
 if(assessed!=43)stop("Expected 43 historically assessed positional gene IDs")
 svg(path,width=8.4,height=4.6,pointsize=11)
 on.exit(dev.off())
 par(mar=c(4.0,7.5,3.3,1.2),family="sans")
 amounts <- c(assessed,2225-assessed)
 barplot(amounts,horiz=TRUE,col=c("#517D9C","#B6BEC6"),
         names.arg=c("GTEx ABF: 646 assays","No QTL assay in this branch"),
         las=1,xlim=c(0,2400),xlab="Positional gene IDs",
         main="Broad ischemic stroke gene universe retained (n = 2,225)")
 text(amounts+27,2:1,labels=amounts,pos=4,xpd=TRUE)
 mtext("Missing in this coloc branch != no biological association",
       side=1,line=2.8,cex=.79,col="#555555")
}
draw_universe(file.path(outdir,"IS_2225_GENE_QTL_COVERAGE.svg"))

cat("IS_FIGURES_COMPLETE",
    paste(basename(list.files(outdir,full.names=FALSE)),collapse=","),
    "\n")
