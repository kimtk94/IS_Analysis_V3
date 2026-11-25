# Fixed acute exercise signature. Read only pinned published aggregate objects.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- args[1]; out <- args[2]
dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_object <- function(file, name) {
  e <- new.env(parent=baseenv())
  got <- load(file.path(root,file), envir=e)
  stopifnot(identical(got,name))
  x <- e[[name]]
  stopifnot(inherits(x,"data.frame"))
  class(x) <- "data.frame"
  x
}
da <- read_object("MUSCLE_TRNSCRPT_DA.rda","MUSCLE_TRNSCRPT_DA")
required <- c("tissue","assay","full_model","contrast","contrast_short","contrast_type","contrast_category","feature_id","logFC","p_value","adj_p_value")
stopifnot(all(required %in% names(da)))
print(dim(da)); print(names(da)); print(table(da$contrast_short))
short <- "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)"
i <- which(da$tissue=="muscle" & da$assay=="transcript-rna-seq" &
           da$contrast_type=="exercise_with_controls" & da$contrast_category=="RE-CON" &
           da$contrast_short==short)
x <- da[i, required, drop=FALSE]
stopifnot(nrow(x)>1000, length(unique(x$full_model))==1L,
          length(unique(x$contrast))==1L,!anyDuplicated(x$feature_id))
print(unique(x$full_model));print(unique(x$contrast))
x$feature_id <- as.character(x$feature_id)
stopifnot(nrow(x)==16509)
canonical <- grepl("^ENSG[0-9]+([.][0-9]+)?$",x$feature_id)
stopifnot(sum(canonical)==16493,
          all(grepl("^ENSG[0-9]+[.][0-9]+_PAR_Y$",x$feature_id[!canonical])))
source_rows <- nrow(x)
source_sig <- x$adj_p_value<0.05 & x$logFC!=0
source_up <- sum(source_sig & x$logFC>0); source_down <- sum(source_sig & x$logFC<0)
cat("FEATURE_ID_GATE",source_rows,sum(canonical),sum(!canonical),"\n")
x <- x[canonical,,drop=FALSE]
rm(da); invisible(gc())
stopifnot(all(is.finite(x$logFC)),all(is.finite(x$p_value)),all(is.finite(x$adj_p_value)),
          all(x$p_value>=0 & x$p_value<=1),all(x$adj_p_value>=0 & x$adj_p_value<=1))
cat("Loading publisher feature map using native R\n");flush.console()
m <- read_object("HUMAN_FEATURE_TO_GENE.rda","HUMAN_FEATURE_TO_GENE")
stopifnot(nrow(m)==1920618L)
cat("MAPPING_SOURCE",nrow(m),paste(names(m),collapse="|"),"\n");flush.console()
needed <- c("assay","feature_id","entrez_gene","ensembl_gene","gene_symbol")
stopifnot(all(needed %in% names(m)))
rna <- which(m$assay=="transcript-rna-seq")
rna <- rna[as.character(m$feature_id[rna]) %in% x$feature_id]
m <- m[rna,needed[-1],drop=FALSE]
cat("EXACT_SOURCE_MAP_ROWS",nrow(m),"\n");flush.console()
m[] <- lapply(m,as.character)
m <- unique(m)
m$ensembl_gene <- sub("[.][0-9]+$","",m$ensembl_gene)
m <- m[!is.na(m$entrez_gene) & grepl("^[0-9]+$",m$entrez_gene) &
       !is.na(m$ensembl_gene) & grepl("^ENSG[0-9]+$",m$ensembl_gene) &
       m$ensembl_gene==sub("[.][0-9]+$","",m$feature_id),,drop=FALSE]
# Reject ambiguous feature-to-gene identities and duplicate stable gene/Entrez keys.
m <- m[!duplicated(m$feature_id) & !duplicated(m$feature_id,fromLast=TRUE),,drop=FALSE]
m <- m[!duplicated(m$entrez_gene) & !duplicated(m$entrez_gene,fromLast=TRUE) &
       !duplicated(m$ensembl_gene) & !duplicated(m$ensembl_gene,fromLast=TRUE),,drop=FALSE]
j <- match(x$feature_id,m$feature_id)
x$entrez <- m$entrez_gene[j];x$ensembl <- m$ensembl_gene[j];x$symbol <- m$gene_symbol[j]
sig <- x$adj_p_value<0.05 & x$logFC!=0
cat("SOURCE_COUNTS",nrow(x),sum(sig & x$logFC>0),sum(sig & x$logFC<0),
    sum(!is.na(x$entrez)),sum(sig & x$logFC>0 & !is.na(x$entrez)),
    sum(sig & x$logFC<0 & !is.na(x$entrez)),"\n")
y <- x[sig & !is.na(x$entrez),,drop=FALSE]
select <- function(up) {
  a <- y[if(up) y$logFC>0 else y$logFC<0,,drop=FALSE]
  a <- a[order(a$p_value,-abs(a$logFC),a$ensembl),,drop=FALSE]
  stopifnot(nrow(a)>=20)
  head(a,100)
}
panel <- rbind(select(TRUE),select(FALSE))
panel$sign <- ifelse(panel$logFC>0,1,-1)
stopifnot(!anyDuplicated(panel$entrez),!anyDuplicated(panel$ensembl))
write.table(panel[,c("entrez","ensembl","symbol","sign","logFC","p_value","adj_p_value")],
            file.path(out,"panel.tsv"), sep="\t",quote=FALSE,row.names=FALSE,na="")
cat("PANEL_COUNTS",sum(panel$sign==1),sum(panel$sign==-1),"\n")
writeLines(c(paste("source_rows",source_rows,sep="\t"),
             paste("canonical_rows",nrow(x),sep="\t"),
             paste("source_fdr_up",source_up,sep="\t"),
             paste("source_fdr_down",source_down,sep="\t"),
             paste("eligible_fdr_up",sum(sig & x$logFC>0 & !is.na(x$entrez)),sep="\t"),
             paste("eligible_fdr_down",sum(sig & x$logFC<0 & !is.na(x$entrez)),sep="\t"),
             paste("mapped_1to1",sum(!is.na(x$entrez)),sep="\t"),
             paste("panel_up",sum(panel$sign==1),sep="\t"),
             paste("panel_down",sum(panel$sign==-1),sep="\t"),
             paste("contrast_short",short,sep="\t"),
             paste("model",as.character(unique(x$full_model)),sep="\t"),
             paste("contrast",as.character(unique(x$contrast)),sep="\t")),file.path(out,"source.tsv"))
