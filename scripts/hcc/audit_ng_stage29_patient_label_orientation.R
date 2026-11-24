#!/usr/bin/env Rscript
# Stage29: using SHA-pinned official Ng tumor-biopsy to unique-patient metadata
# index in private runner memory, evaluate the full↔fullz inferred mapping's
# patient-group consistency. Cannot prove which RDS is canonically labelled.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("usage: Rscript audit_ng_stage29_patient_label_orientation.R PRIVATE_50_PATIENT_INDEX OUTPUT_JSON")
patients<-suppressWarnings(as.integer(readLines(args[[1]],warn=FALSE)))
if(length(patients)!=50L||anyNA(patients)||
  !identical(sort(unique(patients)),seq_len(48L))||
  sum(tabulate(patients,nbins=48)-1L)!=2L)stop("official-source patient cardinality failed")
out<-args[[2]]
base<-"https://raw.githubusercontent.com/cbg-ethz/HCC/438f52271d0a86a6701baa7a96e80c516ba575c2/HCCinputs/"
sha<-c(HCCfull.rds="461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a",
 HCCfullz.rds="9c1e4e05472a0aa477b81aa7b4c4e1b27809df9a2b71c1f57e247bb031263007")
get<-function(name){
 f<-tempfile(fileext=".rds");on.exit(unlink(f),add=TRUE)
 suppressWarnings(download.file(paste0(base,name),f,mode="wb",method="libcurl",quiet=TRUE))
 h<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
 if(length(h)!=1L||sub(" .*","",h)!=sha[[name]])stop("source drift")
 readRDS(f)
}
full<-get("HCCfull.rds");z<-get("HCCfullz.rds")
mods<-c("T","P","PP")
original_names<-colnames(full$T)
if(length(original_names)!=50||anyDuplicated(original_names)||
   !all(vapply(c(full[c("T","P","PP")],z[c("T","P","PP")]),
      function(x)length(colnames(x))==50L &&
        !anyDuplicated(colnames(x))&&setequal(original_names,colnames(x)),
      logical(1))))stop("source RDS biopsy label SET differs from official 50")
full_order_identical<-identical(original_names,colnames(full$P))&&
 identical(original_names,colnames(full$PP))
fullz_order_identical<-identical(colnames(z$T),colnames(z$P))&&
 identical(colnames(z$T),colnames(z$PP))
cross_source_T_order_identical<-identical(original_names,colnames(z$T))
order_map<-function(k){
 a<-as.matrix(full[[k]])
 b<-as.matrix(z[[k]])
 rows<-intersect(rownames(a),rownames(b));cols<-original_names
 if(length(rows)<1000||!all(cols%in%colnames(a))||!all(cols%in%colnames(b)))
    stop("source axes mismatch")
 a<-a[rows,cols,drop=FALSE];b<-b[rows,cols,drop=FALSE]
 if(any(!is.finite(a))||any(!is.finite(b)))stop("invalid source values")
 ix<-order(rownames(a))
 a<-a[ix,,drop=FALSE];b<-b[ix,,drop=FALSE]
 train<-seq(1,nrow(a),by=2L)
 a<-a[train,,drop=FALSE];b<-b[train,,drop=FALSE]
 sigmaa<-apply(a,1,sd);sigmab<-apply(b,1,sd)
 ok<-is.finite(sigmaa)&is.finite(sigmab)&sigmaa>1e-12&sigmab>1e-12
 if(sum(ok)<1000)stop("too few informative rows")
 a<-a[ok,,drop=FALSE];b<-b[ok,,drop=FALSE]
 a<-sweep(sweep(a,1,rowMeans(a),"-"),1,sigmaa[ok],"/")
 b<-sweep(sweep(b,1,rowMeans(b),"-"),1,sigmab[ok],"/")
 max.col(crossprod(a,b)/nrow(a),ties.method="first")
}
p<-lapply(mods,order_map);names(p)<-mods
if(!all(vapply(p,function(i)identical(sort(i),seq_len(50)),logical(1))))
 stop("continuously normalized data not a bijective permutation")
if(!identical(p$T,p$P)||!identical(p$T,p$PP))
 stop("independent omic sample permutations disagree")
# The patient-index vector from official Ng supplement is ordered by the
# source HCCfull column order. Use it ONLY as labels for hypothesized
# full->fullz correspondence, not a claim that fullz has correct metadata.
same_label<-sum(p$T==seq_len(50))
same_patient<-sum(patients[p$T]==patients)
same_patient_swaps<-sum(patients[p$T]==patients & p$T!=seq_len(50))
different_patient<-sum(patients[p$T]!=patients)
if(same_patient+different_patient!=50)stop("sample count not partitioned")
payload<-paste0('{"schema":"hcc-stage29-patient-group-label-orientation-v1",',
 '"source_original_supplement_patient_groups_verified":true,',
 '"official_source_biopsies_n":50,"official_unique_patients_n":48,',
 '"source_full_sha256":"',sha[["HCCfull.rds"]],'",',
 '"source_fullz_sha256":"',sha[["HCCfullz.rds"]],'",',
 '"original_label_sets_identical_across_all_six_T_P_PP_layers":true,',
 '"original_full_T_P_PP_same_order":',tolower(as.character(full_order_identical)),',',
 '"original_fullz_T_P_PP_same_order":',tolower(as.character(fullz_order_identical)),',',
 '"full_T_vs_fullz_T_sample_name_order_equal":',tolower(as.character(cross_source_T_order_identical)),',',
 '"crossmodal_inferred_permutation_identical":true,',
 '"same_biopsy_label_position_best_match_n":',same_label,',',
 '"same_official_patient_group_under_candidate_crossfile_mapping_n":',same_patient,',',
 '"different_official_patient_group_under_candidate_crossfile_mapping_n":',different_patient,',',
 '"within_same_official_patient_group_but_different_biopsy_position_n":',same_patient_swaps,',',
 '"authoritative_orientation_of_source_biopsy_labels_verified":false,',
 '"original_preprocessing_creator_reproduced":false,',
 '"source_sample_labels_or_patient_IDs_exported":false,',
 '"biological_effects_inferred":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(payload,out,useBytes=TRUE)
cat(payload,"\n")
