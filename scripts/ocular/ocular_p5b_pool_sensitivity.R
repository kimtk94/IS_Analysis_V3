#!/usr/bin/env Rscript
# P5b: formal leave-one-donor-pool-out edgeR fits and treatment-vs-pool diagnostics.
# All donor-level sample and raw counts remain in transient CI runner storage.
suppressPackageStartupMessages({library(edgeR);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Usage: Rscript ocular_p5b_pool_sensitivity.R <private_dir> <primary_all_gene_results.tsv.gz> <output_dir>")
input_dir<-args[1];baseline_path<-args[2];out_dir<-args[3]
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
result<-list(dataset="GSE313736",status="HOLD",method="leave_one_pool_out_edgeR_QL",
  scope="in_vitro_BeamA_57_paired_donors",
  caveat="Pool exclusion cannot remove library-vs-treatment confounding; all arms share one study protocol")
save_result<-function()write_json(result,file.path(out_dir,"p5b_pool_sensitivity.json"),
   pretty=TRUE,auto_unbox=TRUE,digits=8,na="null")
tryCatch({
  m<-read.delim(file.path(input_dir,"samples.tsv"),stringsAsFactors=FALSE,check.names=FALSE)
  tab<-read.delim(gzfile(file.path(input_dir,"paired_raw_counts.tsv.gz")),
                  stringsAsFactors=FALSE,check.names=FALSE)
  baseline<-read.delim(gzfile(baseline_path),stringsAsFactors=FALSE,check.names=FALSE)
  if(!all(c("gene_id","logFC","FDR") %in% names(baseline)))stop("BASELINE_SCHEMA_INVALID")
  if(!identical(names(tab)[-1],m$sample_id))stop("PRIVATE_COUNTS_SAMPLE_ORDER_MISMATCH")
  ids<-as.character(tab$gene_id)
  if(anyDuplicated(ids))stop("DUPLICATED_GENE_IDS")
  counts<-as.matrix(tab[,-1,drop=FALSE]);storage.mode(counts)<-"integer"
  rownames(counts)<-ids
  if(ncol(counts)!=114L || nrow(counts)!=36601L || anyNA(counts))stop("BAD_SOURCE_DIMENSIONS")
  m$pair_id<-factor(m$pair_id)
  m$treatment<-factor(m$treatment,levels=c("Vehicle","DEX"))
  m$pool<-factor(m$pool,levels=LETTERS[1:9])
  if(!all(table(m$pair_id,m$treatment)==1L))stop("PAIRING_IMBALANCE")
  if(any(table(m$pool,m$treatment)[,1]!=table(m$pool,m$treatment)[,2]))stop("POOL_TREATMENT_UNBALANCED")

  # Technical diagnostics, not proof of a causally isolated drug effect.
  d0<-DGEList(counts=counts)
  keep0<-filterByExpr(d0,group=m$treatment,min.count=10,min.total.count=15)
  d0<-calcNormFactors(d0[keep0,,keep.lib.sizes=FALSE],method="TMM")
  # Anonymous MDS; all points visible and paired by line segments, no donor labels.
  md<-plotMDS(d0,plot=FALSE)
  if(length(md$x)!=nrow(m) || length(md$y)!=nrow(m))stop("MDS_WRONG_SAMPLE_COUNT")
  png(file.path(out_dir,"p5b_MDS_annotated.png"),width=1050,height=820,res=140)
  plot(md$x,md$y,type="n",xlab="MDS dimension 1",ylab="MDS dimension 2",
       main="BeamA donor-paired pseudobulk: treatment and pool diagnostics")
  for(p in levels(m$pair_id)){
    a<-which(m$pair_id==p)
    if(length(a)==2L)segments(md$x[a[1]],md$y[a[1]],md$x[a[2]],md$y[a[2]],
                            col=adjustcolor("grey70",alpha.f=.65),lwd=.8)
  }
  points(md$x,md$y,pch=16,cex=.8,
         col=ifelse(m$treatment=="DEX","#b63d37","#236ca2"))
  legend("topleft",c("Vehicle","DEX"),pch=16,col=c("#236ca2","#b63d37"),bty="n")
  dev.off()
  mds_table<-data.frame(pool=m$pool,treatment=m$treatment,x=md$x,y=md$y)
  # Only summary R2 metrics are emitted; no coordinates / per-donor trajectories.
  r2<-function(z,term)summary(lm(z~term))$r.squared
  result$mds<-list(n_samples=nrow(m),n_paired_donors=nlevels(m$pair_id),
    first_axis_treatment_r2=unname(r2(md$x,m$treatment)),
    first_axis_pool_r2=unname(r2(md$x,m$pool)),
    second_axis_treatment_r2=unname(r2(md$y,m$treatment)),
    second_axis_pool_r2=unname(r2(md$y,m$pool)))
  result$library_depth<-list(min=min(d0$samples$lib.size),
    max=max(d0$samples$lib.size),
    norm_factor_min=min(d0$samples$norm.factors),
    norm_factor_max=max(d0$samples$norm.factors))
  baseline<-baseline[!duplicated(baseline$gene_id),]
  baseline<-baseline[!is.na(baseline$logFC)&!is.na(baseline$FDR),]
  top100<-head(baseline[order(baseline$PValue),],100)
  fixed<-unique(c(top100$gene_id,baseline$gene_id[
    baseline$gene_symbol %in% c("FKBP5","FAM107A","MYOC","ANGPTL7","TSC22D3","UBE2E3")]))
  rows<-list();bygene<-list()
  for(pool in LETTERS[1:9]){
    ix<-which(m$pool!=pool)
    sub<-droplevels(m[ix,,drop=FALSE])
    sub$pair_id<-droplevels(sub$pair_id)
    sub$treatment<-factor(sub$treatment,levels=c("Vehicle","DEX"))
    design<-model.matrix(~pair_id+treatment,data=sub)
    if(qr(design)$rank!=ncol(design))stop(paste("LOO_RANK_DEFICIENT",pool))
    d<-DGEList(counts=counts[,ix,drop=FALSE])
    keep<-filterByExpr(d,group=sub$treatment,min.count=10,min.total.count=15)
    d<-calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
    if(nrow(d)<5000L)stop(paste("INSUFFICIENT_GENES",pool))
    # Sensitivity: estimate dispersion with robust=FALSE for bounded compute.
    # The primary 57-donor model used robust=TRUE. Effects and FDR cannot be
    # interpreted as identical estimators; biological sign stability is primary.
    d<-estimateDisp(d,design,robust=FALSE)
    fit<-glmQLFit(d,design,robust=FALSE)
    qlf<-glmQLFTest(fit,coef="treatmentDEX")
    tt<-topTags(qlf,n=Inf,sort.by="none")$table
    tt$gene_id<-rownames(tt)
    if(!all(c("logFC","FDR","gene_id")%in%names(tt)))stop("LOO_FIELDS_MISSING")
    j<-match(baseline$gene_id,tt$gene_id)
    both<-which(!is.na(j))
    xs<-baseline$logFC[both];ys<-tt$logFC[j[both]]
    r<-cor(xs,ys,method="pearson")
    retained<-sum(baseline$FDR[both]<.05 & tt$FDR[j[both]]<.05)
    baseline_sig<-sum(baseline$FDR[both]<.05)
    topindex<-match(top100$gene_id,tt$gene_id)
    eligible<-which(!is.na(topindex))
    agreement<-sum(sign(top100$logFC[eligible])==sign(tt$logFC[topindex[eligible]]))
    rows[[pool]]<-data.frame(excluded_pool=pool,paired_donors=nlevels(sub$pair_id),
       tested_genes=nrow(tt),fdr_lt_0_05=sum(tt$FDR<.05),
       baseline_genes_in_common=length(both),
       baseline_fdr_lt_0_05_in_common=baseline_sig,
       baseline_fdr_retained=retained,
       baseline_fdr_retention_rate=if(baseline_sig)retained/baseline_sig else NA_real_,
       pearson_all_common_logFC=r,
       top100_present=length(eligible),top100_sign_agree=agreement,
       top100_direction_concordance=agreement/length(eligible))
    ids_fixed<-intersect(fixed,tt$gene_id)
    if(length(ids_fixed)>0){
      q<-tt[match(ids_fixed,tt$gene_id),c("gene_id","logFC","FDR"),drop=FALSE]
      names(q)[2:3]<-c("sensitivity_logFC","sensitivity_FDR")
      q$excluded_pool<-pool
      bygene[[pool]]<-q
    }
    cat("P5B_LOO_POOL",pool,"N=",nlevels(sub$pair_id),
        "DE=",sum(tt$FDR<.05),"r=",round(r,6),
        "top100_direction=",agreement,"/",length(eligible),"\n",flush=TRUE)
  }
  stats<-do.call(rbind,rows)
  if(nrow(stats)!=9L)stop("EXPECTED_NINE_POOL_REFITS")
  # Only gene-level aggregate statistics; no donor IDs, pool-specific samples.
  write.table(stats,file.path(out_dir,"p5b_leave_one_pool_out.tsv"),
      row.names=FALSE,quote=FALSE,sep="\t")
  fixed_stats<-do.call(rbind,bygene)
  write.table(fixed_stats,file.path(out_dir,"p5b_top_gene_leave_pool.tsv"),
      row.names=FALSE,quote=FALSE,sep="\t")
  png(file.path(out_dir,"p5b_pool_robustness.png"),width=1150,height=850,res=150)
  plot(seq_len(nrow(stats)),stats$pearson_all_common_logFC,ylim=c(0,1.02),
       pch=19,xaxt="n",xlab="Omitted donor pool",ylab="log2FC Pearson vs primary",
       main="Leave-one-pool-out concordance: 9 secondary edgeR fits")
  axis(1,seq_len(nrow(stats)),stats$excluded_pool)
  abline(h=c(.95,.98),lty=2,col=c("grey65","grey50"))
  dev.off()
  result$pool_refits<-9L
  result$secondary_fit_method<-"edgeR_NB_QL_robust_FALSE_for_sensitivity_ONLY"
  result$primary_method<-"edgeR_NB_QL_robust_TRUE"
  result$primary_significant_fdr_0_05<-sum(baseline$FDR<.05)
  result$loo_logFC_pearson_min<-min(stats$pearson_all_common_logFC)
  result$loo_logFC_pearson_median<-median(stats$pearson_all_common_logFC)
  result$loo_top100_direction_agreement_min<-min(stats$top100_direction_concordance)
  result$loo_baseline_fdr_retention_min<-min(stats$baseline_fdr_retention_rate)
  result$loo_fdr_genes_min<-min(stats$fdr_lt_0_05)
  result$loo_fdr_genes_max<-max(stats$fdr_lt_0_05)
  result$per_pool<-lapply(seq_len(nrow(stats)),function(i)as.list(stats[i,]))
  result$batch_confounding_unresolved<-TRUE
  result$clinical_steroid_delta_IOP_unobserved<-TRUE
  result$external_independent_replication_unverified<-TRUE
  result$status<-"PASS_NINE_LPO_EDGE_R_SENSITIVITY_NUMERIC"
  save_result()
  cat("P5B_SUMMARY",result$status,
      "r_min=",round(result$loo_logFC_pearson_min,5),
      "top100_min=",round(result$loo_top100_direction_agreement_min,3),
      "fdr_retention_min=",round(result$loo_baseline_fdr_retention_min,3),"\n")
},error=function(e){
  result$status<-"HOLD_P5B_ERROR"
  result$error<-substr(conditionMessage(e),1,280)
  save_result()
  stop(e)
})
