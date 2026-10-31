#!/usr/bin/env Rscript
# Stage29: original Ng Figure2b per-gene rhos versus Suter HCCfull on
# strict shared 1:1 Ensembl/UniProt RNA-Protein pairs, 48 patient vs 50 biopsy.
# No patient identifiers, feature-level scores or private files exported/logged.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=4L)stop("usage Rscript audit_ng_stage29_fig2b_gene_effects.R ORIGINAL_RHO_TSV STRICT_PAIRS_TSV PATIENT_INDEX OUTPUT_JSON")
pub_path<-args[[1]];pair_path<-args[[2]];index_path<-args[[3]];out<-args[[4]]
original<-read.delim(pub_path,header=FALSE,sep="\t",col.names=c("gene","published_rho"),stringsAsFactors=FALSE)
pairs<-read.delim(pair_path,header=FALSE,sep="\t",col.names=c("gene","protein"),stringsAsFactors=FALSE)
pi<-suppressWarnings(as.integer(readLines(index_path,warn=FALSE)))
if(nrow(original)!=5481L||anyDuplicated(original$gene)||
   any(!grepl("^ENSG[0-9]{11}$",original$gene))||
   anyNA(original$published_rho)||any(!is.finite(original$published_rho))||
   any(abs(original$published_rho)>1))stop("incorrect published Fig2b original rhos")
if(nrow(pairs)!=5343L||anyDuplicated(pairs$gene)||anyDuplicated(pairs$protein)||
   any(!grepl("^ENSG[0-9]{11}$",pairs$gene)))stop("non-strict 5343 source map")
if(length(pi)!=50L||anyNA(pi)||length(unique(pi))!=48L||
   !identical(sort(unique(pi)),seq_len(48L))||
   sum(tabulate(pi,nbins=48L)-1L)!=2L)stop("wrong patient group original source")
shared<-intersect(pairs$gene,original$gene)
if(length(shared)!=5253L)stop("original-current exact ENSEMBL intersection drift")
pairs<-pairs[match(shared,pairs$gene),,drop=FALSE]
pub_rho<-original$published_rho[match(shared,original$gene)]
ref<-"438f52271d0a86a6701baa7a96e80c516ba575c2"
expected<-"461bd5aad6024a3dfef2abffe9d830ed2222c69e0b159bcbfdd35fe779c3869a"
source<-paste0("https://raw.githubusercontent.com/cbg-ethz/HCC/",ref,"/HCCinputs/HCCfull.rds")
f<-tempfile(fileext=".rds")
suppressWarnings(download.file(source,f,mode="wb",method="libcurl",quiet=TRUE))
if(is.na(file.info(f)$size)||file.info(f)$size<100L||file.info(f)$size>16000000L)stop("source RDS bytes changed")
digest<-suppressWarnings(system2("sha256sum",f,stdout=TRUE,stderr=FALSE))
if(length(digest)!=1L||sub(" .*","",digest)!=expected)stop("source RDS SHA mismatch")
obj<-readRDS(f);unlink(f)
if(!is.list(obj)||!all(c("T","P")%in%names(obj)))stop("source T/P omics missing")
a<-as.matrix(obj$T);b<-as.matrix(obj$P)
if(ncol(a)!=50L||ncol(b)!=50L||
   anyDuplicated(colnames(a))||anyDuplicated(colnames(b))||
   !setequal(colnames(a),colnames(b))||
   !all(pairs$gene%in%rownames(a))||!all(pairs$protein%in%rownames(b)))
  stop("original source sample SET or RNA/protein feature alignment mismatch")
# Source assay column order may differ: use verified exact label SET and
# explicitly align P to T, never silently equate column positions.
ordered_biopsies<-colnames(a)
if(is.null(ordered_biopsies)||length(ordered_biopsies)!=50L)
  stop("source biopsies missing")
rna<-a[pairs$gene,ordered_biopsies,drop=FALSE]
prot<-b[pairs$protein,ordered_biopsies,drop=FALSE]
if(any(!is.finite(rna))||any(!is.finite(prot)))stop("nonfinite source quantity")
rna48<-vapply(seq_len(48L),function(k)rowMeans(rna[,pi==k,drop=FALSE]),numeric(nrow(rna)))
prot48<-vapply(seq_len(48L),function(k)rowMeans(prot[,pi==k,drop=FALSE]),numeric(nrow(prot)))
if(!identical(dim(rna48),c(length(shared),48L)))stop("patient average shape changed")
counts<-tabulate(pi,nbins=48L)
singleton<-which(counts==1L)
calc<-function(a,b){
  vapply(seq_len(nrow(a)),function(i){
    if(sd(a[i,])<1e-12||sd(b[i,])<1e-12)return(NA_real_)
    suppressWarnings(cor(a[i,],b[i,],method="spearman"))
  },numeric(1))
}
rho50<-calc(rna,prot);rho48<-calc(rna48,prot48)
rho46<-calc(rna48[,singleton,drop=FALSE],prot48[,singleton,drop=FALSE])
qnum<-function(x)if(length(x)==1L&&is.finite(x))sprintf("%.9f",x)else "null"
metrics<-function(v,baseline){
 idx<-which(is.finite(v)&is.finite(baseline))
 if(length(idx)<5000L)stop("too few matched source gene-wise effect estimates")
 x<-v[idx];y<-baseline[idx]
 agree<-mean(sign(x)==sign(y))
 list(n=length(idx),median=median(x),published=median(y),
    pearson=suppressWarnings(cor(x,y)),
    spearman=suppressWarnings(cor(x,y,method="spearman")),
    sign=agree,
    absdelta=median(abs(x-y)),
    delta=median(x-y),
    q10=unname(quantile(x-y,.1)),
    q90=unname(quantile(x-y,.9)),
    fracDeltaAbove02=mean(abs(x-y)>0.2),
    fracSignOpp=mean(sign(x)!=sign(y)))
}
m50<-metrics(rho50,pub_rho)
m48<-metrics(rho48,pub_rho)
m46<-metrics(rho46,pub_rho)
item<-function(k,x)paste0('{"cohort":"',k,'",',
 '"matched_effect_gene_n":',x$n,',',
 '"source_median_rho":',qnum(x$median),',',
 '"published_same_gene_median_rho":',qnum(x$published),',',
 '"between_studies_gene_rho_pearson":',qnum(x$pearson),',',
 '"between_studies_gene_rho_spearman":',qnum(x$spearman),',',
 '"gene_effect_sign_concordance_fraction":',qnum(x$sign),',',
 '"median_absolute_gene_rho_difference":',qnum(x$absdelta),',',
 '"median_signed_gene_rho_difference":',qnum(x$delta),',',
 '"gene_rho_difference_q10":',qnum(x$q10),',',
 '"gene_rho_difference_q90":',qnum(x$q90),',',
 '"fraction_abs_gene_rho_difference_gt_0_2":',qnum(x$fracDeltaAbove02),',',
 '"gene_effect_sign_discordance_fraction":',qnum(x$fracSignOpp),'}')
msg<-paste0('{"schema":"hcc-stage29-original-fig2b-per-gene-effect-benchmark-v1",',
 '"original_publication_source":"Ng_2022_Figure2b_real_Source_Data",',
 '"official_original_published_gene_n":5481,',
 '"current_source_strict_RNA_protein_pair_n":5343,',
 '"matching_original_published_ENSG_n":5253,',
 '"source_git_commit":"',ref,'",',
 '"source_full_RDS_sha256":"',expected,'",',
 '"source_biopsy_n":50,',
 '"independent_patient_n":48,',
 '"singleton_patient_n":',length(singleton),',',
 '"cohort_metrics":[',paste(c(item("naive_50_biopsies",m50),
                         item("patient_mean_48",m48),
                         item("singleton_46",m46)),collapse=","),'],',
 '"original_2022_quantitation_recreated":false,',
 '"original_2022_2026_gene_annotation_identical":false,',
 '"source_original_patient_tumor_normal_and_batch_controls_verified":false,',
 '"independent_cohort_replication_completed":false,',
 '"per_gene_rho_and_patient_identifiers_exported":false}')
dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
writeLines(msg,out,useBytes=TRUE)
cat(msg,"\n")
