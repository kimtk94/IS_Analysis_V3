#!/usr/bin/env Rscript
# P5o: targeted 96-gene paired-donor / pool-level descriptive sensitivity.
# NOT 9 independent edgeR QL fits; source is confidential ephemeral donor counts.
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript ocular_p5o_targeted_paired_pool.R <private-dir> <p5l-96.tsv> <aggregate-dir>")
priv<-args[1];conflicts<-args[2];out<-args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
summary_path<-file.path(out,"P5o_targeted_pool_sensitivity_summary.json")
result<-list(status="HOLD")
persist<-function()write_json(result,summary_path,pretty=TRUE,auto_unbox=TRUE,digits=8,na="null")
tryCatch({
  meta<-read.delim(file.path(priv,"samples.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  raw<-read.delim(gzfile(file.path(priv,"paired_raw_counts.tsv.gz")),
      check.names=FALSE,stringsAsFactors=FALSE)
  dis<-read.delim(conflicts,check.names=FALSE,stringsAsFactors=FALSE)
  need<-c("gene_id","gene_symbol","primary_log2fc","external_2026_TM_log2fc",
      "poolunit_log2fc","poolunit_FDR","triangulated_conflict")
  if(!all(need %in% names(dis)) || nrow(dis)!=96L || anyDuplicated(dis$gene_id))
      stop("SOURCE_96_GENE_IDS_CHANGED")
  if(!identical(names(raw)[-1],meta$sample_id))
      stop("PRIVATE_SAMPLE_ORDER_MISMATCH")
  mat<-as.matrix(raw[,-1,drop=FALSE]);storage.mode(mat)<-"integer"
  rownames(mat)<-as.character(raw$gene_id)
  if(nrow(mat)!=36601L || ncol(mat)!=114L || anyNA(mat) || any(mat<0L))
      stop("SOURCE_INTEGER_COUNT_MATRIX_CHANGED")
  if(sum(mat)!=1927677307)stop("RAW_UMI_COUNT_CONSERVATION_CHANGED")
  meta$pair_id<-factor(meta$pair_id)
  meta$treatment<-factor(meta$treatment,levels=c("Vehicle","DEX"))
  meta$pool<-factor(meta$pool,levels=LETTERS[1:9])
  if(!all(table(meta$pair_id,meta$treatment)==1L))
      stop("PAIRED_DONOR_BALANCE_INVALID")
  pairnames<-levels(meta$pair_id)
  aidx<-match(paste0(pairnames,"_Vehicle"),meta$sample_id)
  bidx<-match(paste0(pairnames,"_DEX"),meta$sample_id)
  if(anyNA(aidx)||anyNA(bidx))stop("PAIRED_SAMPLE_NAMES_INVALID")
  pair_pool<-as.character(meta$pool[aidx])
  if(any(pair_pool!=as.character(meta$pool[bidx])))
      stop("POOL_ID_NOT_PAIRED")
  if(!identical(sort(unique(pair_pool)),LETTERS[1:9]))
      stop("EXPECTED_9_DONOR_POOLS")
  d<-DGEList(counts=mat)
  keep<-filterByExpr(d,group=meta$treatment,min.count=10,min.total.count=15)
  d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
  genes<-as.character(dis$gene_id)
  if(!all(genes %in% rownames(d$counts)))
      stop("P5L_CONFLICT_GENE_DROPPED_IN_TMM_CPM")
  selected<-cpm(d,log=TRUE,prior.count=2)[genes,,drop=FALSE]
  donor_deltas<-selected[,bidx,drop=FALSE]-selected[,aidx,drop=FALSE]
  if(ncol(donor_deltas)!=57L || nrow(donor_deltas)!=96L)
      stop("PAIR_DELTA_MATRIX_DIMENSION_ERROR")
  pool_delta<-sapply(LETTERS[1:9],function(p)rowMeans(
      donor_deltas[,which(pair_pool==p),drop=FALSE]))
  stopifnot(ncol(pool_delta)==9L,nrow(pool_delta)==96L)
  # Resample independent *pool* means, not donor counts (fixed RNG seed).
  set.seed(314159)
  b_ind<-matrix(sample.int(9,9L*10000L,replace=TRUE),nrow=9L)
  rows<-list()
  strong<-c("ERFE","CHAC1","PARM1")
  for(i in seq_len(nrow(dis))){
    v<-as.numeric(donor_deltas[i,])
    pools<-as.numeric(pool_delta[i,])
    original_sign<-sign(dis$primary_log2fc[i])
    loo<-vapply(seq_len(9L),function(k)mean(pools[-k]),numeric(1))
    # Equal weight to each source pool to avoid oversized pool dominance.
    pool_boot<-matrix(pools[b_ind],nrow=9L)
    boot_means<-colMeans(pool_boot)
    ci<-quantile(boot_means,c(.025,.975),names=FALSE,type=7)
    same_loo<-sum(sign(loo)==original_sign)
    rows[[i]]<-data.frame(
      gene_id=genes[i],gene_symbol=as.character(dis$gene_symbol[i]),
      primary_log2fc=dis$primary_log2fc[i],
      external_100nM_48h_log2fc=dis$external_2026_TM_log2fc[i],
      poolunit_edger_log2fc=dis$poolunit_log2fc[i],
      poolunit_edger_FDR=dis$poolunit_FDR[i],
      paired_57_mean_logCPM_delta=mean(v),
      paired_57_median_logCPM_delta=median(v),
      paired_57_sign_same=sum(sign(v)==original_sign),
      paired_57_sign_opposite=sum(sign(v)==-original_sign),
      independent_pool_mean_logCPM_delta=mean(pools),
      pool_sign_same=sum(sign(pools)==original_sign),
      leave_one_pool_out_mean_min=min(loo),
      leave_one_pool_out_mean_max=max(loo),
      leave_one_pool_out_mean_sign_same=same_loo,
      pool_bootstrap_95_lo=ci[1],
      pool_bootstrap_95_hi=ci[2],
      pool_bootstrap_95_excludes_zero=ci[1]*ci[2]>0,
      test_class="DESCRIPTIVE_TMM_LOGCPM_NOT_EDGE_R_QL")
  }
  stats<-do.call(rbind,rows)
  stats<-stats[order(-abs(stats$primary_log2fc)),]
  write.table(stats,file.path(out,"P5o_96_targeted_paired_pool_sensitivity.tsv"),
     row.names=FALSE,quote=FALSE,sep="\t")
  select<-stats[stats$gene_symbol %in% strong,,drop=FALSE]
  if(nrow(select)!=3L)stop("MISSING_TIER_A_GENES")
  # Scientific diagnostic, at pool-aggregate level; no donor-specific values shown.
  png(file.path(out,"P5o_three_targets_pool_mean_directions.png"),
      width=1100,height=760,res=145)
  plot(c(1,9),range(c(-2,2,pool_delta[dis$gene_symbol %in% strong,])),
       type="n",xlab="Original 10x donor pool (A–I)",
       ylab="Paired mean log2 CPM delta (DEX - Vehicle)",
       xaxt="n",main="Tier-A targets: 9 source-pool signed effects")
  axis(1,at=1:9,labels=LETTERS[1:9])
  abline(h=0,col="grey60",lty=2)
  cols<-c(ERFE="#7b4ab3",CHAC1="#c76c22",PARM1="#207e82")
  for(g in strong){
    ix<-which(as.character(dis$gene_symbol)==g)
    lines(1:9,pool_delta[ix,],col=cols[[g]],lwd=2)
    points(1:9,pool_delta[ix,],col=cols[[g]],pch=16,cex=1.1)
  }
  legend("topleft",strong,col=unname(cols[strong]),lty=1,pch=16,bty="n")
  dev.off()
  result<-list(
      status="PASS_TARGETED_96_PAIRED_POOL_DESCRIPTIVE",
      data="GSE313736 18 original 10x matrices",
      source_discordant_run=41158660082,
      pairs=57L,pools=9L,genes=96L,
      method="TMM-normalized log2 CPM donor-paired differences, equal-weight pool mean, leave-one-pool-out pool-average and pool cluster bootstrap",
      timepoint="nominal 500nM DEX 48h",
      not_method="No robust QL fitting or FDR recalculation at the 9 leave-pool-out subsets",
      targets=lapply(seq_len(nrow(select)),function(i)as.list(select[i,])),
      all_9_pool_sign_same=sum(stats$pool_sign_same==9),
      all_9_leave_pool_out_sign_same=sum(stats$leave_one_pool_out_mean_sign_same==9),
      bootstrap_95_excludes_zero_count=sum(stats$pool_bootstrap_95_excludes_zero),
      caveats=c(
       "Bootstrap intervals condition on nine non-independent source-pool technical-library units; they are descriptive",
       "P5b formal nine edgeR leave-pool-out fits did not include these 96 genes",
       "TMM normalization and post-treatment high-confidence BeamA cell selection can still bias measured contrast",
       "A stable signed pool effect cannot separate DEX pharmacology from correlated treatment-library preparation",
       "External official 100nM bulk TM can only supply published gene-level evidence, not matched 500nM dose causality",
       "Clinical steroid-induced IOP effects unobserved"))
  persist()
  cat("P5O_PASS 96 genes; full_9_pool_signs=",result$all_9_pool_sign_same,
      " full_9_LOO_signs=",result$all_9_leave_pool_out_sign_same,
      " CI_excludes0=",result$bootstrap_95_excludes_zero_count,"\n")
},error=function(e){
  result$status<-"HOLD_TARGETED_PAIRED_POOL_ERROR"
  result$error<-substr(conditionMessage(e),1,350)
  persist()
  stop(e)
})
