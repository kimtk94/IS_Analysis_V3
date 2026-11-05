args <- commandArgs(trailingOnly=TRUE)
mode <- args[[1]]
need <- function(packages) {
  for (p in packages) if (!requireNamespace(p,quietly=TRUE)) stop(paste('Install dependency before run:',p))
}
write_tsv <- function(x,p) write.table(x,p,sep='\t',row.names=FALSE,quote=FALSE,na='NA')
if (mode == 'mr') {
  need('TwoSampleMR')
  d <- read.delim(args[[2]],check.names=FALSE)
  instruments <- read.delim(args[[3]],check.names=FALSE)
  d <- d[d$key %in% instruments$key,,drop=FALSE]
  k <- nrow(d); if (!k) stop('No harmonized instruments')
  pars <- TwoSampleMR::default_parameters()
  rows <- list()
  for (method in c('mr_weighted_median','mr_egger_regression')) {
    if (k < 3) {
      rows[[method]] <- data.frame(method=method,nsnp=k,beta=NA,se=NA,pval=NA,status='NOT_RUN_INSUFFICIENT_INSTRUMENTS')
    } else {
      set.seed(20261003)
      res <- get(method,asNamespace('TwoSampleMR'))(d$beta_x,d$beta_y,d$se_x,d$se_y,pars)
      if (!is.finite(res$b) || !is.finite(res$se)) stop(paste('Invalid MR fit',method))
      rows[[method]] <- data.frame(method=method,nsnp=k,beta=res$b,se=res$se,pval=res$pval,status='SUCCESS',
         intercept=if(method=='mr_egger_regression') res$b_i else NA,
         intercept_pval=if(method=='mr_egger_regression') res$pval_i else NA)
    }
  }
  # Align schemas for insufficient-method rows.
  rows <- lapply(rows,function(x) {if(!'intercept'%in%names(x)){x$intercept<-NA;x$intercept_pval<-NA};x})
  write_tsv(do.call(rbind,rows),args[[4]])
} else if (mode %in% c('coloc','susie')) {
  need(c('jsonlite','coloc',if(mode=='susie') 'susieR' else character()))
  d <- read.delim(args[[2]],check.names=FALSE)
  meta <- jsonlite::fromJSON(args[[3]])
  if(nrow(d)<meta$min_regional_snps) stop('Insufficient regional SNPs')
  dataset <- function(suffix,spec) {
    n <- unique(d[[paste0('n_',suffix)]])
    if(length(n)!=1 || n<3) stop('SuSiE/coloc requires explicit constant per-region sample size')
    ans <- list(beta=d[[paste0('beta_',suffix)]],varbeta=d[[paste0('se_',suffix)]]^2,
                snp=d$key,position=d$pos_x,N=n,type=spec$type,MAF=pmin(d[[paste0('eaf_',suffix)]],1-d[[paste0('eaf_',suffix)]]))
    if(spec$type=='quant') {
      if(identical(spec$sdY_estimation,'coloc')) {
        # Explicit legacy-compatible estimate from varbeta, MAF and N.
        # Never enable this by default or replace a known phenotype SD.
      } else {
        if(is.null(spec$sdY) || !is.finite(spec$sdY) || spec$sdY<=0) stop('Explicit phenotype sdY required')
        ans$sdY <- spec$sdY
      }
    } else if(spec$type=='cc') {
      if(is.null(spec$case_fraction) || spec$case_fraction<=0 || spec$case_fraction>=1) stop('Case fraction required')
      ans$s <- spec$case_fraction
    } else stop('Unsupported trait type')
    coloc::check_dataset(ans); ans
  }
  x <- dataset('x',meta$exposure); y <- dataset('y',meta$outcome)
  if(mode=='coloc') {
    result <- lapply(meta$p12_grid,function(prior){
      fit <- coloc::coloc.abf(x,y,p1=meta$p1,p2=meta$p2,p12=prior)
      data.frame(p12=prior,as.list(fit$summary),check.names=FALSE)
    })
    write_tsv(do.call(rbind,result),args[[4]])
  } else {
    rx <- as.matrix(read.delim(args[[4]],row.names=1,check.names=FALSE))
    ry <- as.matrix(read.delim(args[[5]],row.names=1,check.names=FALSE))
    if(!identical(rownames(rx),d$key)||!identical(rownames(ry),d$key)) stop('LD order mismatch')
    x$LD <- rx; y$LD <- ry
    bounded_runsusie <- function(dat,label) {
      base_maxit <- as.integer(meta$maxit)
      if(!is.finite(base_maxit) || base_maxit < 1L) stop('Invalid SuSiE maxit')
      rescue_maxit <- if(!is.null(meta$maxit_rescue)) {
        as.integer(meta$maxit_rescue)
      } else {
        min(2000L,max(1500L,base_maxit*2L))
      }
      run_once <- function(iter) {
        tryCatch(
          coloc::runsusie(dat,maxit=iter,repeat_until_convergence=FALSE),
          error=function(e)e
        )
      }
      fit <- run_once(base_maxit)
      if(!inherits(fit,'error')) return(fit)
      msg <- conditionMessage(fit)
      if(!grepl('did not converge',msg,fixed=TRUE) || rescue_maxit <= base_maxit) {
        stop(paste('SuSiE',label,'failed:',msg))
      }
      message('SuSiE bounded rescue ',label,': ',base_maxit,' -> ',rescue_maxit)
      fit <- run_once(rescue_maxit)
      if(inherits(fit,'error')) {
        stop(paste('SuSiE',label,'bounded rescue failed at',rescue_maxit,'iterations:',conditionMessage(fit)))
      }
      fit
    }
    sx <- bounded_runsusie(x,'exposure')
    sy <- bounded_runsusie(y,'outcome')
    if(!isTRUE(sx$converged)||!isTRUE(sy$converged)) stop('SuSiE not converged; unresolved, never shared')
    if(!length(sx$sets$cs)||!length(sy$sets$cs)) {
      write_tsv(data.frame(status='UNRESOLVED_NO_CREDIBLE_SET',PP.H4.abf=NA),args[[6]])
    } else {
      fit <- coloc::coloc.susie(sx,sy,p1=meta$p1,p2=meta$p2,p12=meta$p12)
      if(!nrow(fit$summary)) stop('No SuSiE comparisons')
      result<-fit$summary;result$status<-'SUCCESS'
      result$exposure_credible_sets<-length(sx$sets$cs);result$outcome_credible_sets<-length(sy$sets$cs)
      write_tsv(result,args[[6]])
    }
  }
} else if (mode == 'cohort') {
  need(c('lme4','survival'))
  d<-read.delim(args[[2]])
  required<-c('id','time','egfr','score','age','sex')
  if(!all(required%in%names(d))) stop('Cohort input schema missing')
  if(anyDuplicated(d[c('id','time')])) stop('Repeated id/time')
  if(any(!complete.cases(d[required]))) stop('Missing cohort covariates: define policy upstream')
  covars<-strsplit(args[[4]],',',fixed=TRUE)[[1]];covars<-covars[nzchar(covars)]
  if(!all(covars%in%names(d)) || !length(covars)) stop('Explicit ancestry PC covariates required')
  if(any(!complete.cases(d[covars]))) stop('Missing PC covariates')
  counts<-table(d$id);d<-d[d$id%in%names(counts[counts>=3]),]
  if(length(unique(d$id))<30) stop('Too few longitudinal participants')
  d$sex<-factor(d$sex)
  formula<-as.formula(paste('egfr ~ time*score + age + sex +',paste(covars,collapse=' + '),'+ (1 + time | id)'))
  fit<-lme4::lmer(formula,data=d,REML=FALSE)
  if(lme4::isSingular(fit) || length(fit@optinfo$conv$lme4$messages)) stop('Mixed model singular/nonconverged')
  a<-coef(summary(fit));out<-data.frame(term=rownames(a),beta=a[,1],se=a[,2],model='egfr_mixed',pval=NA_real_)
  # No unqualified lmer t-to-p conversion; confidence intervals reported.
  out$ci_low<-out$beta-1.96*out$se;out$ci_high<-out$beta+1.96*out$se
  write_tsv(out,args[[3]])
} else if (mode == 'incident') {
  need('survival')
  d<-read.delim(args[[2]])
  pcs<-strsplit(args[[4]],',',fixed=TRUE)[[1]]
  required<-c('id','followup','event','baseline_ckd','score','age','sex',pcs)
  if(!all(required%in%names(d)) || any(!complete.cases(d[required])) || anyDuplicated(d$id)) stop('Incident cohort input invalid')
  if(!all(d$baseline_ckd%in%c(0,1)) || !all(d$event%in%c(0,1)) || any(d$followup<=0)) stop('Invalid survival endpoints')
  d<-d[d$baseline_ckd==0,,drop=FALSE];d$sex<-factor(d$sex)
  if(sum(d$event)<20) stop('Too few incident events')
  fit<-survival::coxph(as.formula(paste('survival::Surv(followup,event) ~ score + age + sex +',paste(pcs,collapse=' + '))),data=d,x=TRUE)
  a<-summary(fit)$coefficients
  if(any(!is.finite(a))) stop('Cox coefficients not finite')
  ph<-survival::cox.zph(fit)$table
  write_tsv(data.frame(term=rownames(a),beta=a[,1],se=a[,3],pval=a[,5],hr=a[,2],ph_pval=ph[match(rownames(a),rownames(ph)),'p'],global_ph_pval=ph['GLOBAL','p']),args[[3]])
} else stop('Unknown mode')
