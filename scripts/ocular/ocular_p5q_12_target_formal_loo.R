#!/usr/bin/env Rscript
# P5q: 12 predefined direction-discordant targets, formal 9x LOO edgeR QL.
# Two analyses: paired donor (secondary) and pooled 10x library (conservative).
# 57 donor pairs are not 57 independent treatment libraries. Public outputs
# are only gene-level inference, no donor metadata or per-donor counts.
suppressPackageStartupMessages({library(edgeR); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)
  stop("Usage: Rscript ocular_p5q_12_target_formal_loo.R <runner-private-dir> <P5l_96.tsv> <public-dir>")
priv <- args[1]; conflicts_path <- args[2]; out <- args[3]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
result <- list(status="HOLD",study="GSE313736",
  design="9x paired donor LOO plus 9x paired pooled-library LOO",
  caveat="Within-study sensitivity only: DEX vs preparation-specific 10x libraries remains confounded")
persist <- function() write_json(result,file.path(out,"P5q_12_target_formal_loo_summary.json"),
   auto_unbox=TRUE,pretty=TRUE,digits=8,na="null")
strong <- c("ERFE","CHAC1","PARM1")
adequate <- c("ADAMTS15","EGR3","LPAR1","ABHD3","DDIT4L","TFPI","STK17B","RAB32","SDC1")
targets <- c(strong,adequate)
tryCatch({
  m <- read.delim(file.path(priv,"samples.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
  raw <- read.delim(gzfile(file.path(priv,"paired_raw_counts.tsv.gz")),
     check.names=FALSE,stringsAsFactors=FALSE)
  dis <- read.delim(conflicts_path,check.names=FALSE,stringsAsFactors=FALSE)
  needed <- c("gene_id","gene_symbol","primary_log2fc","external_2026_TM_log2fc",
              "poolunit_log2fc","poolunit_FDR","triangulated_conflict")
  if(!all(needed %in% names(dis)) || nrow(dis)!=96L || anyDuplicated(dis$gene_id))
     stop("P5L_96_SCHEMA_OR_UNIVERSE_DRIFT")
  target_rows <- dis[dis$gene_symbol %in% targets,,drop=FALSE]
  if(nrow(target_rows)!=12L || anyDuplicated(target_rows$gene_symbol) ||
      !setequal(target_rows$gene_symbol,targets))
     stop("P5Q_EXACT_12_TARGET_SYMBOLS_REQUIRED")
  if(any(!as.character(target_rows$triangulated_conflict) %in% c("True","TRUE","1")))
     stop("TARGETS_NOT_ALL_POOLUNIT_SUPPORTED")
  if(any(!is.finite(target_rows$primary_log2fc)) ||
     any(target_rows$primary_log2fc*target_rows$external_2026_TM_log2fc>=0))
     stop("TARGET_ORIGINAL_OPPOSITE_SIGN_INVALID")
  if(nrow(m)!=114L || nrow(raw)!=36601L ||
     !identical(names(raw)[-1],m$sample_id) || anyDuplicated(raw$gene_id))
     stop("PRIVATE_SOURCE_SCHEMA_CHANGED")
  counts <- as.matrix(raw[,-1,drop=FALSE]); storage.mode(counts)<-"integer"
  rownames(counts) <- as.character(raw$gene_id)
  if(anyNA(counts) || any(counts<0L) || sum(counts)!=1927677307)
     stop("SOURCE_UMI_CONSERVATION_FAILED")
  if(!all(target_rows$gene_id %in% rownames(counts)))
     stop("TARGET_GENE_MISSING_FROM_RAW")
  rm(raw); gc(verbose=FALSE)
  m$pair_id <- factor(m$pair_id)
  m$pool <- factor(m$pool,levels=LETTERS[1:9])
  m$treatment <- factor(m$treatment,levels=c("Vehicle","DEX"))
  if(nlevels(m$pair_id)!=57L || !all(table(m$pair_id,m$treatment)==1L) ||
     any(table(m$pool,m$treatment)[,1]!=table(m$pool,m$treatment)[,2]))
     stop("PAIR_OR_POOL_IMBALANCE")
  if(any(tapply(as.character(m$pool),m$pair_id,
                function(x) length(unique(x)))!=1L))
     stop("PAIRED_ARMS_CROSS_POOL")
  arms <- expand.grid(pool=LETTERS[1:9],treatment=c("Vehicle","DEX"),
                      stringsAsFactors=FALSE)
  pooled <- sapply(seq_len(nrow(arms)),function(i){
    ix <- which(as.character(m$pool)==arms$pool[i] &
                as.character(m$treatment)==arms$treatment[i])
    if(length(ix)==0L)stop("EMPTY_POOL_TREATMENT_ARM")
    rowSums(counts[,ix,drop=FALSE])
  })
  rownames(pooled) <- rownames(counts)
  arms$pool <- factor(arms$pool,levels=LETTERS[1:9])
  arms$treatment <- factor(arms$treatment,levels=c("Vehicle","DEX"))
  if(sum(pooled)!=sum(counts))stop("POOL_AGGREGATION_NOT_CONSERVED")
  fit_subset <- function(mat,meta,design_formula,kind,omitted,robust_fit){
    meta <- droplevels(meta)
    meta$treatment <- factor(meta$treatment,levels=c("Vehicle","DEX"))
    design <- model.matrix(design_formula,data=meta)
    if(qr(design)$rank!=ncol(design) || !"treatmentDEX" %in% colnames(design))
       stop(paste("INVALID_LOO_DESIGN",kind,omitted))
    d <- DGEList(counts=mat)
    keep <- filterByExpr(d,group=meta$treatment,min.count=10,min.total.count=15)
    d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
    if(nrow(d)<5000L)stop(paste("INSUFFICIENT_GENES",kind,omitted))
    d <- estimateDisp(d,design,robust=robust_fit)
    fit <- glmQLFit(d,design,robust=robust_fit)
    qlf <- glmQLFTest(fit,coef="treatmentDEX")
    all <- topTags(qlf,n=Inf,sort.by="none")$table
    all$gene_id <- rownames(all)
    j <- match(target_rows$gene_id,all$gene_id)
    # Explicitly report filterByExpr dropouts: never assume 9/9 coverage.
    outrow <- data.frame(
      gene_id=as.character(target_rows$gene_id),
      gene_symbol=as.character(target_rows$gene_symbol),
      tier=ifelse(target_rows$gene_symbol %in% strong,"LOW_EXPRESSION_LARGE_EFFECT","ADEQUATE_EXPRESSION"),
      model=kind,omitted_pool=omitted,
      tested=!is.na(j),
      log2fc=as.numeric(all$logFC[j]),
      PValue=as.numeric(all$PValue[j]),
      FDR=as.numeric(all$FDR[j]),
      primary_log2fc=as.numeric(target_rows$primary_log2fc),
      primary_sign_agree=ifelse(is.na(j),FALSE,
         sign(all$logFC[j])==sign(target_rows$primary_log2fc)),
      tested_genes=nrow(all),model_residual_df=nrow(meta)-ncol(design),
      robust_ql=robust_fit,stringsAsFactors=FALSE)
    outrow
  }
  rows <- list()
  for(p in LETTERS[1:9]){
    dix <- as.character(m$pool)!=p
    pix <- as.character(arms$pool)!=p
    rows[[paste0("donor_",p)]] <- fit_subset(
       counts[,dix,drop=FALSE],m[dix,,drop=FALSE],
       ~pair_id+treatment,"paired_donor",p,FALSE)
    rows[[paste0("library_",p)]] <- fit_subset(
       pooled[,pix,drop=FALSE],arms[pix,,drop=FALSE],
       ~pool+treatment,"paired_pool_library",p,TRUE)
    cat("P5Q_OMITTED",p,"donor_n=",sum(dix),
        "pool_libraries=",sum(pix),"\n",flush=TRUE)
  }
  detail <- do.call(rbind,rows); rownames(detail)<-NULL
  if(nrow(detail)!=216L || nrow(unique(detail[c("gene_id","model","omitted_pool")]))!=216L)
     stop("EXPECTED_12_TARGETS_X_2_MODELS_X_9_POOLS")
  write.table(detail,file.path(out,"P5q_12_target_18_loo_refits.tsv"),
     sep="\t",row.names=FALSE,quote=FALSE,na="NA")
  summary_rows <- list()
  for(k in seq_len(nrow(target_rows))){
    for(kind in c("paired_donor","paired_pool_library")){
      sub <- detail[detail$gene_id==target_rows$gene_id[k] & detail$model==kind,,drop=FALSE]
      valid <- sub$tested & is.finite(sub$FDR) & is.finite(sub$log2fc)
      summary_rows[[length(summary_rows)+1L]] <- data.frame(
        gene_id=as.character(target_rows$gene_id[k]),
        gene_symbol=as.character(target_rows$gene_symbol[k]),
        tier=ifelse(target_rows$gene_symbol[k] %in% strong,
                    "LOW_EXPRESSION_LARGE_EFFECT","ADEQUATE_EXPRESSION"),
        model=kind,n_refits=9L,n_tested=sum(valid),
        same_sign=sum(sub$primary_sign_agree & valid),
        fdr_below_0_05=sum(valid & sub$FDR<.05),
        significant_and_same_sign=sum(valid & sub$FDR<.05 & sub$primary_sign_agree),
        max_fdr=if(any(valid)) max(sub$FDR[valid]) else NA_real_,
        min_abs_log2fc=if(any(valid)) min(abs(sub$log2fc[valid])) else NA_real_,
        all_nine_sign_stable=(sum(valid)==9L && all(sub$primary_sign_agree[valid])),
        all_nine_fdr_and_sign=(sum(valid)==9L &&
              all(sub$primary_sign_agree[valid] & sub$FDR[valid]<.05)),
        stringsAsFactors=FALSE)
    }
  }
  stats <- do.call(rbind,summary_rows)
  write.table(stats,file.path(out,"P5q_12_target_formal_loo_by_model.tsv"),
      sep="\t",row.names=FALSE,quote=FALSE,na="NA")
  result <- list(status=if(all(detail$tested))
      "PASS_18_FORMAL_EDGE_R_QL_REFITS" else "PASS_18_REFITS_WITH_TARGET_FILTER_DROPOUT",
    source="GSE313736 18 original 10x matrices; 57 paired BeamA donors",
    source_crosswalk="P5l 96 published direction-discordant gene-level results",
    targets=12L,low_expression_targets=strong,adequate_expression_targets=adequate,
    refits=18L,omitted_pools=LETTERS[1:9],
    donor_refits=9L,conservative_pooled_library_refits=9L,
    donor_model="~pair_id+treatment (robust=FALSE; sensitivity estimator)",
    pooled_library_model="~pool+treatment (robust=TRUE; eight pooled-library pairs)",
    refits_all_targets_tested=all(detail$tested),
    primary_fdr_basis="Per-refit BH across all retained genes, not 12-only BH",
    donor_model_full_9_sig_sign=sum(stats$model=="paired_donor" & stats$all_nine_fdr_and_sign),
    pooled_model_full_9_sig_sign=sum(stats$model=="paired_pool_library" & stats$all_nine_fdr_and_sign),
    requires_independent_dose_randomized_biological_replication=TRUE,
    no_mechanistic_or_clinical_claim=TRUE,
    caveats=c("Both models are within the same original source and cannot distinguish drug from condition-specific library preparation",
       "FDR stability is post-hoc sensitivity, not an independent validation or replication p-value",
       "Cells were assigned to high-confidence BeamA state after exposure; cell selection may bias contrast",
       "A pooled library contains cells from multiple donors and is not a separately randomized biological culture",
       "External 100nM TM 48h is not matched to our 500nM BeamA; dose-reversal causality unavailable",
       "No patient-level steroid delta-IOP or adverse-response endpoint"))
  persist()
  cat("P5Q_COMPLETED",result$status,"pool_all9_sign_fdr=",
       result$pooled_model_full_9_sig_sign,"/12\n",flush=TRUE)
},error=function(e){
  result$status <- "HOLD_P5Q_FORMAL_REFIT_ERROR"
  result$error <- substr(conditionMessage(e),1,320)
  persist()
  stop(e)
})
