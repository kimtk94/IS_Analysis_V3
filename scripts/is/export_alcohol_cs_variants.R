#!/usr/bin/env Rscript
# Exploratory variant-level SuSiE credible-set export; no canonical promotion.
suppressPackageStartupMessages(library(jsonlite))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop("Usage: Rscript export_alcohol_cs_variants.R <sandbox> <original_alcohol_root>")
sandbox <- normalizePath(args[[1]],mustWork=TRUE)
source_root <- normalizePath(args[[2]],mustWork=TRUE)
output <- list()
for(locus in c("ADH1B","ALDH2")){
  variants <- read.delim(file.path(source_root,locus,"variants.tsv"),check.names=FALSE)
  if(!"ID" %in% names(variants)) stop(paste("Missing ID column",locus))
  locus_rows <- list()
  for(mode in c("FINITE_REF_504","FINITE_REF_504_EB_MISMATCH")){
    fit <- readRDS(file.path(sandbox,"models",locus,paste0(mode,".rds")))
    if(length(fit$pip)!=nrow(variants) || ncol(fit$alpha)!=nrow(variants))
      stop(paste("Variant order/dimension mismatch",locus,mode))
    sets <- fit$sets$cs
    if(is.null(sets)) sets <- list()
    for(csname in names(sets)){
      ix <- as.integer(sets[[csname]])
      if(anyNA(ix)||any(ix<1|ix>nrow(variants))) stop("Invalid credible set indices")
      component <- suppressWarnings(as.integer(sub("^L","",csname)))
      if(is.na(component)||component<1||component>nrow(fit$alpha))
        stop(paste("Unmapped SER component",csname))
      for(i in ix){
        locus_rows[[length(locus_rows)+1L]] <- data.frame(
          locus=locus,model=mode,cs=csname,variant_index=i,
          variant_id=as.character(variants$ID[i]),
          pip=as.numeric(fit$pip[i]),
          ser_alpha=as.numeric(fit$alpha[component,i]),
          stringsAsFactors=FALSE)
      }
    }
  }
  if(length(locus_rows)) output[[locus]] <- do.call(rbind,locus_rows)
}
if(!length(output)) stop("No credible-set variants found")
tab <- do.call(rbind,output)
tab <- tab[order(tab$locus,tab$model,tab$cs,-tab$pip),]
out <- file.path(sandbox,"ALCOHOL_CS_VARIANTS_EXPLORATORY.tsv")
write.table(tab,out,sep="\t",quote=FALSE,row.names=FALSE)
write_json(list(status="EXPLORATORY_NOT_CAUSAL",rows=nrow(tab),
  loci=c("ADH1B","ALDH2"),source="original_alcohol_variants.tsv",
  ld_reliability_not_cleared=TRUE),
  file.path(sandbox,"ALCOHOL_CS_VARIANTS_EXPLORATORY_MANIFEST.json"),
  pretty=TRUE,auto_unbox=TRUE)
cat("CS_VARIANT_EXPORT_OK",out,"rows",nrow(tab),"\n")
