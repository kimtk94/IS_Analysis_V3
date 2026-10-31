#!/usr/bin/env Rscript
# Stage22: audit original full-resolution RNA/proteome identifier namespaces.
# Pinned public RDS files. NO patient IDs, feature IDs, matrix values or
# restricted EGA records in logs/artifact. Do not claim biological effects.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop("usage Rscript audit_ng_full_feature_namespace.R OUTPUT")
out<-args[[1]]
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
host<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
fetch<-function(name,max_bytes) {
  file<-tempfile(fileext=".rds")
  on.exit(if(file.exists(file))unlink(file),add=TRUE)
  suppressWarnings(download.file(paste0(host,name),file,mode="wb",quiet=TRUE,method="libcurl"))
  size<-file.info(file)$size
  if(is.na(size)||size<100||size>max_bytes)stop("RDS size invalid")
  sha<-suppressWarnings(system2("sha256sum",file,stdout=TRUE,stderr=FALSE))
  if(length(sha)!=1||!grepl("^[a-f0-9]{64} ",sha))stop("SHA256 missing")
  list(object=readRDS(file),sha=sub(" .*","",sha))
}
obj<-fetch("HCCfull.rds",16000000)
meta<-fetch("HCCinfo.rds",1000000)
info<-meta$object
full<-obj$object
if(!is.list(full)||!all(c("T","P","PP")%in%names(full)))stop("missing T/P/PP")
tids<-rownames(full[["T"]]);pids<-rownames(full[["P"]])
if(length(tids)!=20356||length(pids)!=5585)stop("unexpected full matrix")
if(ncol(full[["T"]])!=50||ncol(full[["P"]])!=50||
   !identical(colnames(full[["T"]]),colnames(full[["P"]])))stop("source sample axis mismatch")
normal<-function(z)sub("\\.[0-9]+$","",as.character(z))
ensg<-grepl("^ENSG[0-9]{11}(\\.[0-9]+)?$",tids)
# UniProt canonical 6- or 10-character accessions and known common forms.
uniprot6<-"^[OPQ][0-9][A-Z0-9]{3}[0-9](-[0-9]+)?$"
uniprot10<-"^[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){2}(-[0-9]+)?$"
# Vectorized UniProt-like label screening (no identifier export).
uniprot<-grepl(uniprot6,pids)|grepl(uniprot10,pids)
uniprot_pipe<-grepl("^(sp|tr)\\|[A-Z0-9]+\\|",pids)
counts<-list(
  transcript_features=length(tids),proteome_features=length(pids),
  transcript_Ensembl_gene_like=sum(ensg),proteome_uniprot_accession_like=sum(uniprot),
  proteome_pipe_accession_like=sum(uniprot_pipe),
  transcript_unique_row_labels=length(unique(tids)),
  protein_unique_row_labels=length(unique(pids)),
  identical_raw_feature_name_overlap=length(intersect(tids,pids)),
  transcript_Ensembl_version_suffix_n=sum(grepl("^ENSG[0-9]{11}\\.[0-9]+$",tids))
)
sourceP<-if(is.list(info)&&is.list(info$allnames))info$allnames$P else NULL
sourceT<-if(is.list(info)&&is.list(info$allnames))info$allnames$T else NULL
crosswalk_counts<-list(status="SOURCE_SELECTED_SUBSET_ONLY",
                       full_ensembl_to_uniprot_join_verified=FALSE)
if(!is.null(sourceP)&&!is.null(sourceT)) {
  sp<-as.data.frame(sourceP,stringsAsFactors=FALSE)
  st<-as.data.frame(sourceT,stringsAsFactors=FALSE)
  if(all(c("name","gene")%in%colnames(sp))&&all(c("name","gene")%in%colnames(st))) {
    crosswalk_counts<-list(
      status="SOURCE_REDUCED_LOOKUP_SCREENED",
      selected_T_rows=nrow(st),selected_P_rows=nrow(sp),
      selected_P_names_identical_to_full_protein_label=sum(as.character(sp$name)%in%pids),
      selected_T_names_identical_to_full_RNA_label=sum(as.character(st$name)%in%tids),
      full_ensembl_to_uniprot_join_verified=FALSE)
  }
}
# fixed keys only, numeric aggregate only; avoid exposing any feature identifiers.
fmt<-paste0('{',
 '"schema":"hcc-ng-full-namespace-v1",',
 '"source_commit":"%s","full_sha256":"%s","info_sha256":"%s",',
 '"n_full_RNA":%d,"n_full_protein":%d,',
 '"RNA_ENSG_like_n":%d,"protein_UniProt_like_n":%d,"protein_piped_UniProt_like_n":%d,',
 '"RNA_unique_feature_labels_n":%d,"protein_unique_feature_labels_n":%d,',
 '"raw_feature_names_direct_overlap_n":%d,"RNA_ENSG_version_suffix_n":%d,',
 '"source_selected_RNA_lookup_n":%d,"source_selected_protein_lookup_n":%d,',
 '"selected_RNA_lookup_matches_full_n":%d,"selected_protein_lookup_matches_full_n":%d,',
 '"source_selected_mapping_status":"%s",',
 '"Ng2022_paper_comparable_genes_benchmark":5490,',
 '"Ng2022_paper_mRNA_protein_correlation_benchmark":5481,',
 '"full_ensembl_uniprot_mapping_verified":false,',
 '"same_biopsy_rna_protein_abundance_effect_tested":false,',
 '"feature_or_patient_ids_exported":false}')
get<-function(key)if(is.null(crosswalk_counts[[key]]))0L else as.integer(crosswalk_counts[[key]])
line<-sprintf(fmt,ref,obj$sha,meta$sha,
  counts$transcript_features,counts$proteome_features,
  counts$transcript_Ensembl_gene_like,counts$proteome_uniprot_accession_like,
  counts$proteome_pipe_accession_like,counts$transcript_unique_row_labels,
  counts$protein_unique_row_labels,counts$identical_raw_feature_name_overlap,
  counts$transcript_Ensembl_version_suffix_n,
  get("selected_T_rows"),get("selected_P_rows"),
  get("selected_T_names_identical_to_full_RNA_label"),
  get("selected_P_names_identical_to_full_protein_label"),
  crosswalk_counts$status)
writeLines(line,out,useBytes=TRUE)
cat(line,"\n")
if(counts$transcript_Ensembl_gene_like==0)quit(status=1)
