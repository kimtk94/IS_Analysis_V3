#!/usr/bin/env Rscript
# Stage24: descriptive RNA-protein cross-patient Spearman on 1:1 mapped genes.
# PUBLIC Ng/Suter HCCfull processed arrays. Not Ng original raw-count inference.
# NO feature names, patient IDs, expression values or per-gene statistics
# may appear in stdout or output artifact. The original 50 biopsies represent
# only 48 biological patients. Multiple biopsies are explicitly sensitivity QC.

args<-commandArgs(trailingOnly=TRUE)
if(!(length(args)%in%c(3L,4L)))stop("usage: Rscript audit_ng_stage24_patient_spearman.R PRIVATE_PAIRS PRIVATE_PATIENT_INDEX OUTPUT_JSON [STAGE25_INFLUENCE_JSON]")
pairpath<-args[[1]];patientpath<-args[[2]];out<-args[[3]]
if(!file.exists(pairpath)||!file.exists(patientpath))stop("missing ephemeral source mapping")
pairs<-read.delim(pairpath,sep="\t",header=FALSE,quote="",stringsAsFactors=FALSE,
                  col.names=c("gene","protein"))
if(nrow(pairs)<5000 || nrow(pairs)>5500 || ncol(pairs)!=2L ||
   anyDuplicated(pairs$gene)||anyDuplicated(pairs$protein)||
   any(!grepl("^ENSG[0-9]{11}$",pairs$gene)))
  stop("not a strict feature-unique one-to-one pair list")
patients<-suppressWarnings(as.integer(readLines(patientpath,warn=FALSE)))
if(length(patients)!=50L||anyNA(patients)||any(patients<1L)||
   length(unique(patients))!=48L ||
   !identical(sort(unique(patients)),seq_len(48L)))
  stop("invalid source-verified patient assignment")
counts<-tabulate(patients,nbins=48L)
if(sum(counts-1L)!=2L)stop("biopsy-patient cardinality contradiction")

ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
expected<-"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
url<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/HCCfull.rds")
f<-tempfile(fileext=".rds")

suppressWarnings(download.file(url,f,mode="wb",quiet=TRUE,method="libcurl"))
if(is.na(file.info(f)$size)||file.info(f)$size<100L||file.info(f)$size>16000000L)
  stop("unexpected source RDS size")
sha<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
if(length(sha)!=1L||!identical(sub(" .*","",sha),expected))stop("source checksum drift")
full<-readRDS(f)
unlink(f)
if(!is.list(full)||!all(c("T","P")%in%names(full)))stop("missing RDS layer")
T<-full$T;P<-full$P
# The source uses numeric data.frames for some layers; preserve exact axes
# while converting to matrices. Never require is.matrix() as an input type.
if(!(is.matrix(T)||is.data.frame(T))||
   !(is.matrix(P)||is.data.frame(P))||
   nrow(T)!=20356L||nrow(P)!=5585L||
   ncol(T)!=50L||ncol(P)!=50L||!identical(colnames(T),colnames(P)))
  stop("source quantitative layer geometry or sample order changed")
T<-as.matrix(T);P<-as.matrix(P)
if(!is.numeric(T)||!is.numeric(P))stop("source data not numeric")
if(!all(pairs$gene%in%rownames(T))||!all(pairs$protein%in%rownames(P)))
  stop("mapped ID not present in source array")
rna<-T[pairs$gene,,drop=FALSE]
prot<-P[pairs$protein,,drop=FALSE]
if(!is.numeric(rna)||!is.numeric(prot)||any(!is.finite(rna))||any(!is.finite(prot)))
  stop("unexpected nonfinite numeric value in processed data")
# Mean is a *sensitivity assumption* for the two extra biopsy measurements,
# not an original quantitation method. Include 46 singleton-only control.
patient_average<-function(m){
  do.call(cbind,lapply(seq_len(48L),function(j)rowMeans(m[,patients==j,drop=FALSE])))
}
rna48<-patient_average(rna);prot48<-patient_average(prot)
singleton<-which(counts==1L)
if(length(singleton)<45L||length(singleton)>47L)stop("unexpected duplicate distribution")
rna_single<-rna48[,singleton,drop=FALSE]
prot_single<-prot48[,singleton,drop=FALSE]
row_spearman<-function(a,b){
  vapply(seq_len(nrow(a)),function(i){
    x<-a[i,];y<-b[i,]
    if(length(unique(x))<2L||length(unique(y))<2L)return(NA_real_)
    suppressWarnings(cor(x,y,method="spearman"))
  },numeric(1))
}
rho48<-row_spearman(rna48,prot48)
rho_single<-row_spearman(rna_single,prot_single)
rho50<-row_spearman(rna,prot)
q<-function(v,p)if(length(v)==0L)NA_real_ else unname(quantile(v,p,na.rm=TRUE))
summarize<-function(v){
  z<-v[is.finite(v)]
  list(finite_feature_n=length(z),undefined_constant_feature_n=length(v)-length(z),
       median=q(z,.5),q10=q(z,.1),q90=q(z,.9),
       fraction_positive=if(length(z))mean(z>0) else NA_real_,
       fraction_abs_ge_0_5=if(length(z))mean(abs(z)>=.5) else NA_real_)
}
s48<-summarize(rho48);s46<-summarize(rho_single);s50<-summarize(rho50)
num<-function(x)if(length(x)!=1L||is.na(x)||!is.finite(x))"null" else sprintf("%.6f",x)
obj<-function(z)paste0("{",
  '"finite_feature_n":',z$finite_feature_n,",",
  '"undefined_constant_feature_n":',z$undefined_constant_feature_n,",",
  '"median":',num(z$median),",",
  '"q10":',num(z$q10),",",
  '"q90":',num(z$q90),",",
  '"fraction_positive":',num(z$fraction_positive),",",
  '"fraction_abs_ge_0_5":',num(z$fraction_abs_ge_0_5),"}")
line<-paste0("{",
  '"schema":"hcc-ng-stage24-descriptive-spearman-v1",',
  '"source_public_git_commit":"',ref,'",',
  '"source_RDS_sha256":"',expected,'",',
  '"source_type":"Suter_preprocessed_normalized_imputed_batch_adjusted_not_raw_counts",',
  '"mapping_scope":"source_RNA_x_current_UniProt_reviewed_strict_one_to_one",',
  '"crosswalk_protein_gene_pair_n":',nrow(pairs),",",
  '"source_biopsy_n":50,',
  '"independent_patient_n":48,',
  '"patient_repeated_biopsy_groups_n":',sum(counts>1L),",",
  '"singleton_only_patient_n":',length(singleton),",",
  '"patient_mean_48":',obj(s48),",",
  '"singleton_only":',obj(s46),",",
  '"biopsy_50_nonindependent_sensitivity":',obj(s50),",",
  '"median_rho_patient_mean_minus_singleton_only":',num(s48$median-s46$median),",",
  '"full_gene_universe_original_Ng_reproduced":false,',
  '"original_Ng_5481_gene_correlation_reproduced":false,',
  '"multiple_testing_correction_or_gene_level_pvalues_calculated":false,',
  '"direct_etiology_or_recurrence_tested":false,',
  '"scientific_significance_or_causality_proven":false,',
  '"feature_patient_identifiers_or_expression_exported":false',
  "}")
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(line,out,useBytes=TRUE)
cat(line,"\n")
if(s48$finite_feature_n<5000L||s46$finite_feature_n<5000L)
  stop("too few informative mapped pairs for descriptive QC")

# Stage25: optional leave-one-independent-patient-out and cluster/bootstrap
# uncertainty. This is technical stability evidence, not a clinical CI or
# claim about individual genes. All resampling is by *patient*, never biopsy.
if(length(args)==4L){
  influence_out<-args[[4]]
  baseline<-s48$median
  loo<-vapply(seq_len(48L),function(j){
    keep<-seq_len(48L)!=j
    v<-row_spearman(rna48[,keep,drop=FALSE],prot48[,keep,drop=FALSE])
    if(sum(is.finite(v))<5000L)stop("too few defined LOO rho")
    median(v,na.rm=TRUE)
  },numeric(1))
  if(anyNA(loo)||any(!is.finite(loo)))stop("nonfinite leave-one-out rho")
  set.seed(20261031L)
  B<-80L
  boots<-vapply(seq_len(B),function(b){
    indexes<-sample.int(48L,48L,replace=TRUE)
    v<-row_spearman(rna48[,indexes,drop=FALSE],prot48[,indexes,drop=FALSE])
    if(sum(is.finite(v))<5000L)stop("too few defined bootstrap rho")
    median(v,na.rm=TRUE)
  },numeric(1))
  if(anyNA(boots)||any(!is.finite(boots)))stop("nonfinite bootstrap median")
  qlo<-unname(quantile(boots,.025));qhi<-unname(quantile(boots,.975))
  sorted<-sort(loo)
  ldelta<-max(abs(loo-baseline))
  # Flag is descriptive risk threshold; never infer gene-level confidence
  # from a stable cohort-wide median.
  influence_line<-paste0("{",
    '"schema":"hcc-ng-stage25-patient-influence-v1",',
    '"source_RDS_sha256":"',expected,'",',
    '"mapping_scope":"source_RNA_x_current_UniProt_reviewed_strict_one_to_one",',
    '"strict_gene_protein_pair_n":',nrow(pairs),",",
    '"patient_n":48,',
    '"leave_one_patient_out_iterations":48,',
    '"primary_median_spearman_rho":',num(baseline),",",
    '"leave_one_out_median_min":',num(min(loo)),",",
    '"leave_one_out_median_max":',num(max(loo)),",",
    '"leave_one_out_median_q05":',num(unname(quantile(loo,.05))),",",
    '"leave_one_out_median_q95":',num(unname(quantile(loo,.95))),",",
    '"leave_one_out_max_absolute_shift":',num(ldelta),",",
    '"leave_one_out_max_shift_gt_0_05_review_flag":',tolower(as.character(ldelta>0.05)),",",
    '"patient_bootstrap_resamples":',B,",",
    '"patient_bootstrap_seed":20261031,',
    '"patient_bootstrap_median_rho_q025":',num(qlo),",",
    '"patient_bootstrap_median_rho_q975":',num(qhi),",",
    '"bootstrap_interval_interpretation":"empirical_patient_resampling_descriptive_only",',
    '"technical_preprocessing_robustness_proven":false,',
    '"original_ng_correlation_reproduced":false,',
    '"gene_level_FDR_or_per_gene_CI_calculated":false,',
    '"clinical_or_causal_validity_proven":false,',
    '"patient_gene_rho_or_expression_values_exported":false',
    "}")
  dir.create(dirname(influence_out),recursive=TRUE,showWarnings=FALSE)
  writeLines(influence_line,influence_out,useBytes=TRUE)
  cat(influence_line,"\n")
}
