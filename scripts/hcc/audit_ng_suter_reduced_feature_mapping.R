#!/usr/bin/env Rscript
# Pinned-source public Suter HCCinfo gene mappings for reduced cross-omics.
# All identifiers/individual values stay in-memory; output counts only.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop("usage: Rscript audit_ng_suter_reduced_feature_mapping.R OUT_JSON")
out<-args[[1]];dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
src<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/")
fetch<-function(fname,max_bytes) {
 tmp<-tempfile(fileext=".rds")
 on.exit(if(file.exists(tmp))unlink(tmp),add=TRUE)
 suppressWarnings(download.file(paste0(src,fname),tmp,quiet=TRUE,mode="wb",method="libcurl"))
 bytes<-file.info(tmp)$size
 if(is.na(bytes)||bytes<100||bytes>max_bytes)stop("invalid immutable public RDS byte count")
 sha<-suppressWarnings(system2("sha256sum",tmp,stdout=TRUE,stderr=FALSE))
 if(length(sha)!=1||!grepl("^[a-f0-9]{64} ",sha))stop("no source digest")
 list(object=readRDS(tmp),sha256=sub(" .*","",sha),bytes=as.integer(bytes))
}
inp<-fetch("HCCdata.rds",500000)
inf<-fetch("HCCinfo.rds",500000)
data<-inp$object; info<-inf$object
if(!is.list(data)||!all(c("T","P","PP")%in%names(data))||!is.list(info)||!is.list(info$allnames))
 stop("required public mapping structures absent")
if(!all(c("T","P")%in%names(info$allnames)))stop("T/P gene annotation tables missing")
id_map<-function(omic) {
 x<-as.matrix(data[[omic]])
 colids<-colnames(x)
 tbl<-as.data.frame(info$allnames[[omic]],stringsAsFactors=FALSE)
 if(!all(c("name","gene")%in%colnames(tbl)))stop("mapping source lacks name/gene")
 ids<-as.character(tbl[["name"]])
 if(anyDuplicated(ids))stop("ambiguous source feature key")
 index<-match(colids,ids)
 genes<-as.character(tbl[["gene"]][index])
 genes<-toupper(trimws(genes))
 genes[is.na(genes)|genes%in%c("","NA","N/A","NULL")]<-NA_character_
 list(feature_total=ncol(x),sample_n=nrow(x),n_feature_key_matched=sum(!is.na(index)),
      genes=genes,source_feature_key_duplicates=anyDuplicated(ids)>0)
}
a<-id_map("T");b<-id_map("P")
if(a$sample_n!=50||b$sample_n!=50||a$feature_total!=188||b$feature_total!=116)
 stop("unexpected source reduced geometry")
genesT<-a$genes[!is.na(a$genes)];genesP<-b$genes[!is.na(b$genes)]
shared<-intersect(unique(genesT),unique(genesP))
cntT<-table(genesT);cntP<-table(genesP)
unique_pairs<-sum(cntT[shared]==1&cntP[shared]==1)
status<-if(length(shared)>0)"SOURCE_REDUCED_GENE_CROSSWALK_PRESENT" else "SOURCE_REDUCED_GENE_CROSSWALK_ABSENT"
# JSON is composed from fixed, well-controlled keys and numeric/boolean values;
# raw accession IDs or gene names cannot reach this report or job stdout.
fmt<-paste0('{"status":"%s","source_commit":"%s",',
    '"HCCdata_sha256":"%s","HCCinfo_sha256":"%s",',
    '"source_mapping_object_class":"%s","samples_per_modality":50,',
    '"transcript_feature_total":188,"protein_feature_total":116,',
    '"T_source_feature_keys_matched_n":%d,"P_source_feature_keys_matched_n":%d,',
    '"T_unique_gene_symbols_mapped_n":%d,"P_unique_gene_symbols_mapped_n":%d,',
    '"T_P_shared_gene_symbols_n":%d,"T_P_unambiguous_one_to_one_genes_n":%d,',
    '"full_feature_Ensembl_protein_mapping_verified":false,',
    '"source_gene_mapping_independently_curated":false,',
    '"clinical_effects_verified":false,',
    '"patient_ids_feature_names_expression_values_exported":false}')
klass<-if(inherits(info,"bnInfo"))"bnInfo" else "bnInfo_like_allnames"
# Original RDS may retain a plain list instead of the bnInfo class attribute;
# structure and immutable source SHA256, not class spelling, anchor provenance.
line<-sprintf(fmt,status,ref,inp$sha256,inf$sha256,klass,
    a$n_feature_key_matched,b$n_feature_key_matched,
    length(unique(genesT)),length(unique(genesP)),length(shared),unique_pairs)
writeLines(line,out,useBytes=TRUE)
cat(line,"\n")
if(length(shared)==0)quit(status=1)
