#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop('Usage: 77_stage5h_model_diagnostics.R <analysis_table.tsv> <outdir>')
infile<-args[1]; outdir<-args[2]; dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d<-read.delim(infile,stringsAsFactors=FALSE,check.names=FALSE)
if(!requireNamespace('survival',quietly=TRUE)) stop('survival package required')
grs<-grep('^GRS_',names(d),value=TRUE); if(!length(grs)) stop('No GRS_* columns')
base<-d[!duplicated(d$participant_id),]
res<-list(); ix<-0
for(g in grs){
  covs<-c('age','sex',grep('^PC[0-9]+$',names(base),value=TRUE))
  need<-c('event_time','incident_mets',g,covs); if(any(!need %in% names(base))) next
  fit<-survival::coxph(as.formula(paste0('survival::Surv(event_time,incident_mets) ~ ',g,' + ',paste(covs,collapse=' + '))),data=base,x=TRUE)
  z<-survival::cox.zph(fit)
  tab<-as.data.frame(z$table)
  tab$term<-rownames(tab)
  for(i in seq_len(nrow(tab))){ix<-ix+1;res[[ix]]<-data.frame(grs=g,term=tab$term[i],rho=tab$rho[i],chisq=tab$chisq[i],p=tab$p[i],stringsAsFactors=FALSE)}
}
out<-if(length(res))do.call(rbind,res)else data.frame()
write.table(out,file=file.path(outdir,'STAGE5H_COX_PH_DIAGNOSTICS.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
cat('diagnostic rows=',nrow(out),'\n')
