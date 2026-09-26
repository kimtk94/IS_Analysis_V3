#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) stop('Usage: 74_stage5e_koges_longitudinal_models.R <analysis_table.tsv> <outdir>')
infile <- args[1]; outdir <- args[2]; dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- read.delim(infile,stringsAsFactors=FALSE,check.names=FALSE)
req <- c('participant_id','time_years','MBI','incident_mets','event_time','age','sex')
miss <- setdiff(req,names(d)); if(length(miss)) stop(paste('Missing:',paste(miss,collapse=', ')))
grs_cols <- grep('^GRS_',names(d),value=TRUE)
if(!length(grs_cols)) stop('No GRS_* columns found')
if(!requireNamespace('survival',quietly=TRUE)) stop('survival package required')
if(!requireNamespace('nlme',quietly=TRUE)) stop('nlme package required')
res <- list(); ix <- 0
base <- d[!duplicated(d$participant_id),]
for(g in grs_cols){
  covs <- c('age','sex')
  pc <- grep('^PC[0-9]+$',names(base),value=TRUE)
  covs <- c(covs,pc)
  form <- as.formula(paste0('survival::Surv(event_time, incident_mets) ~ ',g,' + ',paste(covs,collapse=' + ')))
  fit <- survival::coxph(form,data=base,x=TRUE)
  co <- summary(fit)$coefficients
  if(g %in% rownames(co)){
    ix<-ix+1; res[[ix]]<-data.frame(model='COX_INCIDENT_METS',grs=g,beta=co[g,'coef'],se=co[g,'se(coef)'],p=co[g,'Pr(>|z|)'],stringsAsFactors=FALSE)
  }
  dd <- d[is.finite(d$MBI) & is.finite(d$time_years),]
  if(nrow(dd)>0){
    f2 <- as.formula(paste0('MBI ~ time_years + ',g,' + ',g,':time_years'))
    lmm <- nlme::lme(fixed=f2,random=~time_years|participant_id,data=dd,na.action=na.omit,control=nlme::lmeControl(opt='optim'))
    tt <- summary(lmm)$tTable
    term <- paste0(g,':time_years'); if(!(term %in% rownames(tt))) term <- paste0('time_years:',g)
    if(term %in% rownames(tt)){
      ix<-ix+1; res[[ix]]<-data.frame(model='LMM_MBI_SLOPE',grs=g,beta=tt[term,'Value'],se=tt[term,'Std.Error'],p=tt[term,'p-value'],stringsAsFactors=FALSE)
    }
  }
}
out <- if(length(res)) do.call(rbind,res) else data.frame()
write.table(out,file=file.path(outdir,'STAGE5E_KOGES_LONGITUDINAL_MODELS.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
cat('rows=',nrow(out),'\n')
