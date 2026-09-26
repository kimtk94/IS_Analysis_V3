#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<2) stop('Usage: 75_stage5f_koges_pa_interaction.R <analysis_table.tsv> <outdir>')
infile<-args[1]; outdir<-args[2]; dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d<-read.delim(infile,stringsAsFactors=FALSE,check.names=FALSE)
req<-c('participant_id','time_years','MBI','physical_activity')
miss<-setdiff(req,names(d)); if(length(miss)) stop(paste('Missing:',paste(miss,collapse=', ')))
if(!requireNamespace('nlme',quietly=TRUE)) stop('nlme package required')
grs_cols<-grep('^GRS_',names(d),value=TRUE); if(!length(grs_cols)) stop('No GRS_* columns found')
res<-list(); ix<-0
for(g in grs_cols){
  dd<-d[is.finite(d$MBI)&is.finite(d$time_years)&is.finite(d$physical_activity),]
  if(nrow(dd)==0) next
  form<-as.formula(paste0('MBI ~ time_years + ',g,' + physical_activity + ',g,':physical_activity + ',g,':physical_activity:time_years'))
  fit<-nlme::lme(fixed=form,random=~time_years|participant_id,data=dd,na.action=na.omit,control=nlme::lmeControl(opt='optim'))
  tt<-summary(fit)$tTable
  terms<-rownames(tt)[grepl(g,rownames(tt),fixed=TRUE)&grepl('physical_activity',rownames(tt),fixed=TRUE)]
  for(term in terms){ix<-ix+1;res[[ix]]<-data.frame(grs=g,term=term,beta=tt[term,'Value'],se=tt[term,'Std.Error'],p=tt[term,'p-value'],stringsAsFactors=FALSE)}
}
out<-if(length(res)) do.call(rbind,res) else data.frame()
write.table(out,file=file.path(outdir,'STAGE5F_KOGES_PA_INTERACTION.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
cat('rows=',nrow(out),'\n')
