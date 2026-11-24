#!/usr/bin/env Rscript
# Gate B: same donor-blocked edgeR QL model refit at BeamA >=0.50, raw integer pseudobulk counts.
# The paired count matrix and pseudonymous sample IDs remain runner-local and private.
suppressPackageStartupMessages({
  library(edgeR)
  library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Rscript ocular_gate_b_050_edger.R <runner_input_dir> <summary_output_dir>")
input_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
summary_path <- file.path(out_dir,"p5_edger_summary.json")
summary <- list(status="HOLD",dataset="GSE313736",
  phenotype="in_vitro_trabecular_meshwork_dexamethasone_not_clinical_IOP",
  primary_cell_state="BeamA_score_ge_0.50",
  method="edgeR_paired_NB_QL_BH_FDR",
  caveat="DEX and vehicle preparations sequenced in separate 10x libraries")
persist <- function() write_json(summary,summary_path,pretty=TRUE,auto_unbox=TRUE,digits=8,na="null")
tryCatch({
  count_path <- file.path(input_dir,"paired_raw_counts.tsv.gz")
  sample_path <- file.path(input_dir,"samples.tsv")
  feature_path <- file.path(input_dir,"gene_features.tsv")
  if (!all(file.exists(c(count_path,sample_path,feature_path)))) stop("RUNNER_LOCAL_FILES_MISSING")
  meta <- read.delim(sample_path,check.names=FALSE,stringsAsFactors=FALSE)
  counts_df <- read.delim(gzfile(count_path),check.names=FALSE,stringsAsFactors=FALSE)
  feats <- read.delim(feature_path,check.names=FALSE,stringsAsFactors=FALSE)
  stopifnot(!anyDuplicated(counts_df$gene_id),!anyDuplicated(meta$sample_id))
  if(!identical(names(counts_df)[-1],meta$sample_id)) stop("SAMPLE_ORDER_MISMATCH")
  if(!identical(as.character(counts_df$gene_id),as.character(feats$gene_id))) stop("GENE_ORDER_MISMATCH")
  counts <- as.matrix(counts_df[,-1,drop=FALSE])
  storage.mode(counts) <- "integer"
  rownames(counts) <- as.character(counts_df$gene_id)
  if(anyNA(counts) || any(counts<0L)) stop("NEGATIVE_OR_MISSING_COUNTS")
  if(ncol(counts)!=114L || nrow(counts)!=36601L)stop("GATE_B_SCORE050_57_PAIRS_36601_GENES_REQUIRED")
  meta$pair_id <- factor(meta$pair_id)
  meta$treatment <- factor(meta$treatment,levels=c("Vehicle","DEX"))
  meta$pool <- factor(meta$pool)
  tab <- table(meta$pair_id,meta$treatment)
  if(nrow(tab)!=57L || !all(tab==1L))stop("GATE_B_050_UNPAIRED_OR_DUPLICATED_DONOR")
  if(any(tapply(as.character(meta$pool),meta$pair_id,function(x)length(unique(x)))!=1L))
    stop("POOL_CROSSED_BETWEEN_PAIRED_ARMS")
  design <- model.matrix(~pair_id+treatment,data=meta)
  if(qr(design)$rank!=ncol(design) || !"treatmentDEX" %in% colnames(design))
    stop("DESIGN_RANK_DEFICIENT")
  d <- DGEList(counts=counts,genes=feats)
  # For paired designs filter on the treatment groups, NOT the high-rank donor block.
  # edgeR User Guide section 2.7 explicitly recommends group=Treatment here.
  keep <- filterByExpr(d,group=meta$treatment,min.count=10,min.total.count=15)
  if(sum(keep)<5000L)stop("UNEXPECTED_LOW_EXPRESSED_GENE_NUMBER")
  d <- d[keep,,keep.lib.sizes=FALSE]
  d <- calcNormFactors(d,method="TMM")
  # Same 57 donor eligibility verified via metadata-only Gate B; expression effects
  # remain observational DEX-associated and confounded with condition-specific libraries.
  # paired donor is the biological replicate; each pair has Vehicle and DEX
  d <- estimateDisp(d,design,robust=TRUE)
  fit <- glmQLFit(d,design,robust=TRUE)
  test <- glmQLFTest(fit,coef="treatmentDEX")
  tab <- topTags(test,n=Inf,sort.by="none")$table
  if (!all(c("logFC","logCPM","PValue","FDR")%in%names(tab)))stop("EXPECTED_EDGER_FIELDS_MISSING")
  if(anyNA(tab$FDR) || any(tab$FDR<0 | tab$FDR>1))stop("BAD_FDR")
  tab$gene_id <- rownames(tab)
  if (!"gene_symbol"%in%names(tab)) {
    mapper <- setNames(as.character(feats$gene_symbol),as.character(feats$gene_id))
    tab$gene_symbol <- mapper[tab$gene_id]
  }
  tab <- tab[,c("gene_id","gene_symbol","logFC","logCPM","F","PValue","FDR")]
  tab <- tab[order(tab$PValue),]
  q <- tab$FDR < 0.05
  priority <- q & abs(tab$logFC)>=0.25
  top <- head(tab,100L)
  write.table(top,file.path(out_dir,"p5_top100_genes.tsv"),
              sep="\t",row.names=FALSE,quote=FALSE)
  all_out <- gzfile(file.path(out_dir,"p5_all_gene_results.tsv.gz"),open="wt")
  write.table(tab,all_out,sep="\t",row.names=FALSE,quote=FALSE)
  close(all_out)
  # Pool-specific descriptive concordance, NOT leave-one-pool-out QL re-fitting.
  cpm_mat <- cpm(d,log=TRUE,prior.count=2)
  pool_sign <- data.frame(gene_id=top$gene_id,gene_symbol=top$gene_symbol,
       global_logFC=top$logFC,pools_agree=integer(nrow(top)),
       pools_tested=integer(nrow(top)))
  for(i in seq_len(nrow(top))) {
    id <- top$gene_id[i]
    if(!id %in% rownames(cpm_mat))next
    diffs <- numeric()
    for(p in levels(meta$pool)){
      ix <- which(meta$pool==p)
      sam <- meta[ix,]
      vals <- cpm_mat[id,ix]
      by_pair <- split(seq_along(ix),as.character(sam$pair_id))
      within <- vapply(by_pair,function(inds){
        a <- inds[sam$treatment[inds]=="Vehicle"]
        b <- inds[sam$treatment[inds]=="DEX"]
        if(length(a)!=1L||length(b)!=1L) return(NA_real_)
        as.numeric(vals[b]-vals[a])
      },numeric(1))
      if(any(is.finite(within)))diffs <- c(diffs,mean(within,na.rm=TRUE))
    }
    pool_sign$pools_tested[i] <- length(diffs)
    pool_sign$pools_agree[i] <- sum(sign(diffs)==sign(top$logFC[i]))
  }
  write.table(pool_sign,file.path(out_dir,"p5_top100_pool_direction.tsv"),
              sep="\t",row.names=FALSE,quote=FALSE)
  # Diagnostic figures: do not label donors. These are not clinical biomarkers.
  png(file.path(out_dir,"p5_MDS_aggregate.png"),width=1100,height=850,res=145)
  # Blank labels suppress plotting in some edgeR/limma versions, so draw
  # anonymous points explicitly from MDS coordinates (never donor labels).
  mds <- plotMDS(d,plot=FALSE)
  plot(mds$x,mds$y,pch=16,cex=.7,
       col=ifelse(meta$treatment=="DEX","firebrick","steelblue"),
       xlab="MDS dimension 1",ylab="MDS dimension 2",
       main="BeamA donor-paired 10x pseudobulk: anonymous samples")
  legend("topright",legend=c("Vehicle","DEX"),
         col=c("steelblue","firebrick"),pch=16,bty="n")
  dev.off()
  png(file.path(out_dir,"p5_MA_aggregate.png"),width=1100,height=850,res=145)
  plot(tab$logCPM,tab$logFC,pch=16,cex=.28,
       col=ifelse(priority,"firebrick","grey70"),xlab="Average logCPM",
       ylab="DEX vs Vehicle log2FC",main="Paired donor BeamA edgeR QL")
  abline(h=c(-.25,.25),col="grey40",lty=2)
  dev.off()
  controls <- c("FKBP5","FAM107A","TSC22D3","MYOC","ANGPTL7","ZBTB16","UBE2E3")
  result_controls <- tab[tab$gene_symbol %in% controls,
                 c("gene_symbol","logFC","PValue","FDR")]
  summary <- list(status="GO_EDGER_MODEL_FIT_QC_SCORE050",dataset="GSE313736",
    model="~ pair_id + treatment; treatmentDEX versus Vehicle",
    cell_state="celltypist_chen2025_fine/BeamA",classification_score_threshold=.50,
    min_cells_per_arm=30L,biological_pairs=57L,samples=ncol(counts),
    total_features=nrow(counts),tested_genes=nrow(d),excluded_genes=sum(!keep),
    model_rank=qr(design)$rank,design_columns=ncol(design),
    residual_df=ncol(counts)-qr(design)$rank,
    edgeR_version=as.character(packageVersion("edgeR")),
    R_version=R.version.string,
    fdr_lt_0_05=sum(q),fdr_lt_0_05_and_abs_logFC_ge_0_25=sum(priority),
    dex_higher_fdr=sum(q&tab$logFC>0),dex_lower_fdr=sum(q&tab$logFC<0),
    tmm_factor_min=min(d$samples$norm.factors),
    tmm_factor_max=max(d$samples$norm.factors),
    library_size_min=min(d$samples$lib.size),
    library_size_max=max(d$samples$lib.size),
    common_dispersion=d$common.dispersion,
    positive_control_genes=lapply(seq_len(nrow(result_controls)),function(i)as.list(result_controls[i,])),
    top_10_genes=lapply(seq_len(min(10L,nrow(top))),function(i)as.list(top[i,c("gene_id","gene_symbol","logFC","FDR")])),
    warnings=c("FDR significance is not clinical risk prediction",
      "Gene expression differences do not establish mechanism or independent replication",
      "Treatment is confounded with separate 10x libraries; 9 donor-pool consistency is descriptive",
      "No individual donor-level count data or metadata uploaded"))
  persist()
  cat("P5_EDGE_R_PASS donors=57 tested=",nrow(d)," FDR=",sum(q),
      " prioritized=",sum(priority)," edgeR=",as.character(packageVersion("edgeR")),"\n")
},error=function(e){
  summary$status <- "HOLD_EDGER_MODEL"
  summary$error <- paste(class(e)[1],substr(conditionMessage(e),1,350))
  persist()
  stop(e)
})
