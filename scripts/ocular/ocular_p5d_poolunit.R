#!/usr/bin/env Rscript
# P5d deliberately conservative experiment-unit sensitivity:
# 9 pools * 2 conditions = 18 distinct 10x libraries, paired by pool.
# Not independent drug-vs-batch validation; same experiment/time.
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript ocular_p5d_poolunit.R <runner-private-dir> <primary-gene-level.tsv.gz> <public-output-dir>")
private<-args[1];prior_file<-args[2];out<-args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
result<-list(status="HOLD",study="GSE313736",
  caution="Pool-level libraries are independent units, but each treatment arm has a distinct 10x library; treatment and library prep remain linked.")
persist<-function()write_json(result,file.path(out,"p5d_poolunit_model_summary.json"),auto_unbox=TRUE,pretty=TRUE,digits=8,na="null")
tryCatch({
  m<-read.delim(file.path(private,"samples.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  tab<-read.delim(gzfile(file.path(private,"paired_raw_counts.tsv.gz")),check.names=FALSE,stringsAsFactors=FALSE)
  prior<-read.delim(gzfile(prior_file),check.names=FALSE,stringsAsFactors=FALSE)
  if(!identical(names(tab)[-1],m$sample_id))stop("COUNT_METADATA_ORDER_INCONSISTENT")
  if(nrow(m)!=114L||nrow(tab)!=36601L)stop("SOURCE_DIMENSIONS_CHANGED")
  mat<-as.matrix(tab[,-1,drop=FALSE]);storage.mode(mat)<-"integer"
  rownames(mat)<-as.character(tab$gene_id)
  if(anyNA(mat)||any(mat<0L))stop("INVALID_UMI")
  m$pool<-factor(m$pool,levels=LETTERS[1:9])
  m$treatment<-factor(m$treatment,levels=c("Vehicle","DEX"))
  if(!all(table(m$pool,m$treatment)[,1]==table(m$pool,m$treatment)[,2]))
    stop("INCONSISTENT_DONOR_COUNTS_ACROSS_CONDITIONS")
  pool_arms<-expand.grid(pool=LETTERS[1:9],treatment=c("Vehicle","DEX"),
                         stringsAsFactors=FALSE)
  pool_counts<-sapply(seq_len(nrow(pool_arms)),function(i){
    pool<-pool_arms$pool[i];arm<-pool_arms$treatment[i]
    ix<-which(as.character(m$pool)==pool&as.character(m$treatment)==arm)
    if(length(ix)<1L)stop("EMPTY_POOL_ARM")
    rowSums(mat[,ix,drop=FALSE])
  })
  stopifnot(ncol(pool_counts)==18L)
  pool_arms$pool<-factor(pool_arms$pool,levels=LETTERS[1:9])
  pool_arms$treatment<-factor(pool_arms$treatment,levels=c("Vehicle","DEX"))
  dmat<-model.matrix(~pool+treatment,data=pool_arms)
  if(qr(dmat)$rank!=ncol(dmat)||ncol(dmat)!=10L)stop("POOL_DESIGN_BAD_RANK")
  d<-DGEList(counts=pool_counts)
  keep<-filterByExpr(d,group=pool_arms$treatment,min.count=10,min.total.count=15)
  d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
  d<-estimateDisp(d,dmat,robust=TRUE)
  fit<-glmQLFit(d,dmat,robust=TRUE)
  test<-glmQLFTest(fit,coef="treatmentDEX")
  tt<-topTags(test,n=Inf,sort.by="none")$table
  tt$gene_id<-rownames(tt)
  prior_sub<-prior[!duplicated(prior$gene_id),c("gene_id","logFC","PValue","FDR","gene_symbol")]
  cidx<-match(prior_sub$gene_id,tt$gene_id)
  common<-which(!is.na(cidx))
  if(length(common)<10000L)stop("TOO_FEW_COMMON_TESTED_GENES")
  fold_cor<-cor(prior_sub$logFC[common],tt$logFC[cidx[common]],method="pearson")
  spearman<-cor(prior_sub$logFC[common],tt$logFC[cidx[common]],method="spearman")
  both_sig<-sum(prior_sub$FDR[common]<.05&tt$FDR[cidx[common]]<.05)
  base_sig<-sum(prior_sub$FDR[common]<.05)
  # Gene-level results summarize *pooled library* model only, not donor profiles.
  names(tt)[names(tt)=="logFC"]<-"pool_logFC"
  names(tt)[names(tt)=="FDR"]<-"pool_FDR"
  names(tt)[names(tt)=="PValue"]<-"pool_PValue"
  filtered<-tt[,c("gene_id","pool_logFC","logCPM","F","pool_PValue","pool_FDR")]
  filtered$gene_symbol<-prior_sub$gene_symbol[match(filtered$gene_id,prior_sub$gene_id)]
  filtered<-filtered[order(filtered$pool_PValue),]
  out_f<-gzfile(file.path(out,"p5d_poolunit_all_gene_results.tsv.gz"),"wt")
  write.table(filtered,out_f,sep="\t",quote=FALSE,row.names=FALSE)
  close(out_f)
  primary_top<-head(prior_sub[order(prior_sub$PValue),],100)
  top_ix<-match(primary_top$gene_id,filtered$gene_id)
  valid<-which(!is.na(top_ix))
  directional<-sum(sign(primary_top$logFC[valid])==
     sign(filtered$pool_logFC[top_ix[valid]]))
  # Separate 9 library-pool directions for top 100 using library-level TMM logCPM
  pool_cpm<-cpm(d,log=TRUE,prior.count=2)
  detail<-list()
  for(g in primary_top$gene_id){
    if(!g %in% rownames(pool_cpm))next
    diffs<-numeric()
    for(pool in LETTERS[1:9]){
      no<-which(pool_arms$pool==pool & pool_arms$treatment=="Vehicle")
      de<-which(pool_arms$pool==pool & pool_arms$treatment=="DEX")
      diffs<-c(diffs,as.numeric(pool_cpm[g,de]-pool_cpm[g,no]))
    }
    orig<-primary_top$logFC[match(g,primary_top$gene_id)]
    detail[[g]]<-data.frame(gene_id=g,
         pools_agree=sum(sign(diffs)==sign(orig)),
         pooled_mean_delta=mean(diffs),
         pool_delta_sd=sd(diffs))
  }
  tab100<-do.call(rbind,detail)
  write.table(tab100,file.path(out,"p5d_top100_pool_consistency.tsv"),
    sep="\t",row.names=FALSE,quote=FALSE)
  # MDS comparison at *18 pool-level count libraries*.
  md<-plotMDS(d,plot=FALSE)
  png(file.path(out,"p5d_18_pool_library_MDS.png"),width=1000,height=800,res=140)
  plot(md$x,md$y,pch=16,cex=1.2,
    col=ifelse(pool_arms$treatment=="DEX","#b53d36","#236ca2"),
    xlab="MDS dimension 1",ylab="MDS dimension 2",
    main="9 paired pools (18 pooled 10x libraries)")
  for(pool in LETTERS[1:9]){
    a<-which(pool_arms$pool==pool)
    segments(md$x[a[1]],md$y[a[1]],md$x[a[2]],md$y[a[2]],
             col="grey60",lty=2,lwd=1)
  }
  legend("topright",c("Vehicle","DEX"),pch=16,
    col=c("#236ca2","#b53d36"),bty="n")
  dev.off()
  result<-list(status="PASS_9_POOL_PAIRED_EDGE_R_QC",
      primary_biological_donors=57L,
      pooling_unit="9 separate source donor pools x 2 treatments",
      pooled_libraries=18L,design_rank=qr(dmat)$rank,
      residual_df=18L-qr(dmat)$rank,
      raw_genes=nrow(mat),genes_tested=nrow(tt),
      fdr_lt_0_05=sum(filtered$pool_FDR<.05),
      fdr_and_abs_logFC_ge_0_25=sum(filtered$pool_FDR<.05&abs(filtered$pool_logFC)>=.25),
      pearson_logFC_vs_donor_model=unname(fold_cor),
      spearman_logFC_vs_donor_model=unname(spearman),
      common_genes=length(common),
      primary_FDR_genes_common=base_sig,
      FDR_genes_shared=both_sig,
      primary_FDR_retain_frac=both_sig/base_sig,
      primary_top100_direction_agree=directional,
      primary_top100_considered=length(valid),
      top100_all9_pool_signs=sum(tab100$pools_agree==9),
      per_library_median_umi=median(colSums(pool_counts)),
      common_dispersion=d$common.dispersion,
      batch_confounded=TRUE,
      clinical_delta_IOP_unobserved=TRUE,
      notes=c("Pooled technical library model is a sensitivity analysis, not a causal identification of DEX vs library effect",
          "9 pools are not 9 independent laboratories. Cell-level pseudoreplication is avoided.",
          "All 18 libraries are derived from same protocol; reference experiment remains external and unmapped."))
  persist()
  cat("P5D_POOLUNIT_PASS tested=",nrow(tt),
      " FDR=",sum(filtered$pool_FDR<.05),
      " pearson=",round(fold_cor,5),
      " overlap=",both_sig,"/",base_sig,
      " top100_direction=",directional,"/",length(valid),"\n")
},error=function(e){
  result$status<-"HOLD_POOLUNIT_ERROR"
  result$error<-substr(conditionMessage(e),1,240)
  persist()
  stop(e)
})
