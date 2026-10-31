#!/usr/bin/env Rscript
# Public, source-pinned Ng/Suter 2022 HCC RDS quantitative layer audit.
# Only aggregate dimensions, column-label overlap counts and numeric coverage.
# Never emit sample IDs, patient identifiers, feature labels, clinical values,
# processed matrices, or raw RDS files in the public GitHub Actions artifacts.

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1) stop("usage: Rscript audit_ng_public_multiomics_rds.R OUTPUT_JSON")
out <- args[[1]]
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
ref <- "438f52271d0a86a6701baa7a96e80c516ba575c2"
host <- paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
files <- c("HCCfull.rds","HCCfullz.rds","HCCdata.rds","info.rds")
source_bytes_max <- 16000000
source_files <- list()
objects <- list()
redact <- function(x) {
  x <- as.character(x)
  ifelse(is.na(x),NA_character_,x)
}
sha256_file <- function(path) {
  x <- suppressWarnings(system2("sha256sum",path,stdout=TRUE,stderr=FALSE))
  if (length(x)!=1 || !grepl("^[0-9a-f]{64} ",x)) return(NA_character_)
  sub(" .*","",x)
}
audit_layer <- function(x,sample_axis="columns") {
  if (!(is.matrix(x)||is.data.frame(x))) return(list(type="NOT_A_MATRIX",numeric_matrix=FALSE))
  nc<-ncol(x); nr<-nrow(x)
  sample_count<-if(sample_axis=="rows")nr else nc
  feature_count<-if(sample_axis=="rows")nc else nr
  sample_ids<-if(sample_axis=="rows")rownames(x) else colnames(x)
  feature_ids<-if(sample_axis=="rows")colnames(x) else rownames(x)
  sample_present<-!is.null(sample_ids)&&length(sample_ids)==sample_count&&all(nzchar(sample_ids))
  feature_present<-!is.null(feature_ids)&&length(feature_ids)==feature_count&&all(nzchar(feature_ids))
  mat<-as.matrix(x)
  numeric_type<-is.numeric(mat)||is.integer(mat)
  values<-if(numeric_type)as.numeric(mat) else numeric(0)
  finite_values<-values[is.finite(values)]
  value_quantiles<-if(length(finite_values))as.list(setNames(
    round(as.numeric(quantile(finite_values,c(.05,.25,.5,.75,.95),names=FALSE)),4),
    c("p05","p25","p50","p75","p95"))) else list()
  feature_patterns<-if(feature_present)c(
    ensembl_gene=sum(grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$",feature_ids)),
    uniprot_swissprot=sum(grepl("^[A-NR-Z][0-9][A-Z0-9]{3}[0-9]$",feature_ids)),
    nonempty=sum(nzchar(feature_ids))
  ) else c(ensembl_gene=0,uniprot_swissprot=0,nonempty=0)
  list(
    type=if(is.matrix(x))"matrix" else "data.frame",
    orientation=if(sample_axis=="rows")"samples_by_features" else "features_by_samples",
    native_rows=as.integer(nr), native_cols=as.integer(nc),
    feature_count=as.integer(feature_count),sample_count=as.integer(sample_count),
    numeric_matrix=as.logical(numeric_type),
    finite_numeric_count=if(numeric_type)as.integer(sum(is.finite(values))) else NA_integer_,
    missing_or_nonfinite_count=if(numeric_type)as.integer(sum(!is.finite(values))) else NA_integer_,
    negative_numeric_count=if(numeric_type)as.integer(sum(finite_values<0)) else NA_integer_,
    zero_numeric_count=if(numeric_type)as.integer(sum(finite_values==0)) else NA_integer_,
    value_quantiles_global=value_quantiles,
    feature_label_style_counts=as.list(feature_patterns),
    feature_labels_nonempty=as.logical(feature_present),
    sample_labels_nonempty=as.logical(sample_present),
    unique_sample_labels=if(sample_present)as.integer(length(unique(sample_ids))) else NA_integer_,
    duplicate_sample_label_count=if(sample_present)as.integer(sample_count-length(unique(sample_ids))) else NA_integer_,
    raw_labels_or_values_exported=FALSE
  )
}
audit_object <- function(obj, name) {
  if(name=="info.rds") {
    return(list(status="CLINICAL_OBJECT_STRUCTURE_ONLY",container_class=class(obj)[1],
      rows=if(!is.null(dim(obj)))as.integer(nrow(obj)) else NA_integer_,
      columns=if(!is.null(dim(obj)))as.integer(ncol(obj)) else NA_integer_,
      event_coding_verified=FALSE,raw_case_or_clinical_values_exported=FALSE))
  }
  wanted<-c("M","CN","T","P","PP")
  if(!is.list(obj)||!all(wanted %in% names(obj))) {
    return(list(status="UNRECOGNIZED_MULTIMODAL_CONTAINER",multimodal_names_verified=FALSE))
  }
  # Suter clustering HCCdata.rds is samples x reduced features (50 x [24,292,188,116,158]).
  # HCCfull/HCCfullz are features x samples. Never confuse feature columns with cases.
  sample_axis<-if(name=="HCCdata.rds")"rows" else "columns"
  layers<-setNames(lapply(wanted,function(k)audit_layer(obj[[k]],sample_axis)),wanted)
  getsamples<-function(k) {
    x<-obj[[k]]
    if(!(is.matrix(x)||is.data.frame(x)))return(NULL)
    if(sample_axis=="rows")rownames(x) else colnames(x)
  }
  tnames<-getsamples("T");pnames<-getsamples("P")
  identical_tp<-!is.null(tnames)&&!is.null(pnames)&&identical(tnames,pnames)
  intersect_tp<-if(!is.null(tnames)&&!is.null(pnames))length(intersect(tnames,pnames)) else NA_integer_
  same_all<-all(vapply(wanted,function(k) {
    v<-getsamples(k)
    !is.null(v)&&!is.null(tnames)&&identical(v,tnames)
  },logical(1)))
  tf<-if(sample_axis=="rows")colnames(obj[["T"]]) else rownames(obj[["T"]])
  pf<-if(sample_axis=="rows")colnames(obj[["P"]]) else rownames(obj[["P"]])
  shared_feat<-if(!is.null(tf)&&!is.null(pf))length(intersect(tf,pf)) else NA_integer_
  list(
    status="PUBLIC_MULTIMODAL_RDS_SHAPES_AUDITED",
    sample_axis=sample_axis,
    layers=layers,
    transcript_protein_identical_sample_labels=identical_tp,
    transcript_protein_shared_sample_labels=as.integer(intersect_tp),
    all_modality_exact_sample_order=same_all,
    transcript_protein_shared_feature_label_count=as.integer(shared_feat),
    biological_patient_sample_provenance_verified=FALSE,
    quantitative_assay_normalization_approved=FALSE,
    same_biopsy_rna_protein_scientific_join_verified=FALSE
  )
}

for (fname in files) {
  tmp <- tempfile(fileext=".rds")
  item <- list(source_public_git_commit=ref,filename=fname,
               status="UNVERIFIED",original_rds_uploaded=FALSE)
  tryCatch({
    suppressWarnings(download.file(paste0(host,fname),tmp,
                                   quiet=TRUE,mode="wb",method="libcurl"))
    bytes <- file.info(tmp)$size
    if(is.na(bytes)||bytes<100||bytes>source_bytes_max)
      stop("unsafe RDS bytes")
    item$bytes <- as.integer(bytes)
    item$sha256 <- sha256_file(tmp)
    x <- readRDS(tmp)
    # Optional *ephemeral only* sample keys for within-run exact original
    # supplementary clinical matching. The file is never uploaded or logged.
    private_keys <- Sys.getenv("NG_PRIVATE_KEYS_OUT", unset="")
    if(fname=="HCCfull.rds" && nzchar(private_keys)) {
      ids <- if(is.list(x) && !is.null(x[["T"]])) colnames(x[["T"]]) else NULL
      if(!is.null(ids) && length(ids)==50 && all(nzchar(ids))) {
        writeLines(ids,private_keys,useBytes=TRUE)
        Sys.chmod(private_keys,mode="0600")
      }
    }
    objects[[fname]] <- audit_object(x,fname)
    item$status <- "PUBLIC_PINNED_RDS_READ_IN_MEMORY"
    rm(x)
  },error=function(e){item$status <- "SOURCE_NOT_VALIDATED";item$error_class <- "DOWNLOAD_OR_RDS_PARSE_FAILURE"})
  if(file.exists(tmp)) unlink(tmp)
  source_files[[fname]] <- item
}
# Independent clinical metadata in public repository is a structural reference
# only: never inspect patient/outcome values or claim event interpretation.
required <- c("HCCfull.rds","HCCfullz.rds","HCCdata.rds")
passed <- all(vapply(required,function(f) source_files[[f]]$status=="PUBLIC_PINNED_RDS_READ_IN_MEMORY",logical(1)))
report <- list(
  schema_version="hcc-ng-public-rds-audit-v1",
  source_repository="https://github.com/cbg-ethz/HCC",
  source_commit=ref,
  source_publication_doi="10.1371/journal.pcbi.1009767",
  original_cohort_doi="10.1038/s41467-022-29960-8",
  files=source_files,
  objects=objects,
  required_three_multimodal_objects_read=passed,
  restricted_ega_bypass_attempted=FALSE,
  exact_patients_from_matched_clinical_codebook_verified=FALSE,
  true_patient_by_gene_and_protein_quantitative_join_scientifically_validated=FALSE,
  V02="BLOCKED_SOURCE_PROVENANCE_AND_SAMPLE_CODEBOOK_PENDING",
  raw_identifiers_or_numeric_matrices_in_artifact=FALSE
)
# Only base R: no JSON serialization dependencies needed.
json_escape <- function(x) {
  gsub("\n","\\n",gsub('"','\\\\"',gsub("\\\\","\\\\\\\\",x,fixed=TRUE),fixed=TRUE),fixed=TRUE)
}
encode <- function(x) {
 if(is.null(x)) return("null")
 if(is.list(x)) {
   nm<-names(x)
   if(is.null(nm))return(paste0("[",paste(vapply(x,encode,character(1)),collapse=","),"]"))
   return(paste0("{",paste(vapply(seq_along(x),function(i)
      paste0('"',json_escape(nm[[i]]),'":',encode(x[[i]])),character(1)),collapse=","),"}"))
 }
 if(is.logical(x))return(if(length(x)==1){if(is.na(x))"null" else if(x)"true" else "false"}else paste0("[",paste(ifelse(is.na(x),"null",ifelse(x,"true","false")),collapse=","),"]"))
 if(is.numeric(x))return(if(length(x)==1){if(is.na(x)||!is.finite(x))"null" else as.character(x)}else paste0("[",paste(ifelse(is.na(x)|!is.finite(x),"null",as.character(x)),collapse=","),"]"))
 if(is.character(x))return(if(length(x)==1){if(is.na(x))"null" else paste0('"',json_escape(x),'"')}else paste0("[",paste(vapply(x,encode,character(1)),collapse=","),"]"))
 "null"
}
writeLines(encode(report),out,useBytes=TRUE)
for(fname in required){
  cat(fname,source_files[[fname]]$status,"\n")
}
if(!passed)quit(status=1)
