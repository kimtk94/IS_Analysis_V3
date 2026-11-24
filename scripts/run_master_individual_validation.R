#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

usage <- paste(
  "usage: run_master_individual_validation.R",
  "--subject-file <tsv> --predictor <col> --output <tsv>",
  "[--baseline-endpoint <col>] [--slope-endpoint <col>]",
  "[--incident-endpoint <col>] [--covariates a,b,c]",
  "[--time-to-event <col> --event <col>]",
  "[--longitudinal-file <tsv> --id-col <col> --time-col <col> --egfr-col <col>]"
)

get_arg <- function(flag, default=NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("missing value for", flag))
  args[[i+1]]
}

subject_file <- get_arg("--subject-file")
predictor <- get_arg("--predictor")
output_file <- get_arg("--output")
if (is.null(subject_file) || is.null(predictor) || is.null(output_file)) stop(usage)

baseline_endpoint <- get_arg("--baseline-endpoint")
slope_endpoint <- get_arg("--slope-endpoint")
incident_endpoint <- get_arg("--incident-endpoint")
covariate_arg <- get_arg("--covariates", "")
time_to_event <- get_arg("--time-to-event")
event_col <- get_arg("--event")
longitudinal_file <- get_arg("--longitudinal-file")
id_col <- get_arg("--id-col", "participant_id")
time_col <- get_arg("--time-col", "time_years")
egfr_col <- get_arg("--egfr-col", "egfr")

covariates <- trimws(strsplit(covariate_arg, ",", fixed=TRUE)[[1]])
covariates <- covariates[nzchar(covariates)]

dat <- read.delim(subject_file, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
if (!predictor %in% names(dat)) stop(paste("predictor column missing:", predictor))

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

fit_lm <- function(endpoint, model_name) {
  if (is.null(endpoint) || !endpoint %in% names(dat)) {
    return(data.frame(
      model=model_name, endpoint=ifelse(is.null(endpoint),"",endpoint),
      method="linear_regression", status="endpoint_missing",
      n=NA, events=NA, beta=NA, se=NA, statistic=NA, p=NA,
      stringsAsFactors=FALSE
    ))
  }
  vars <- unique(c(endpoint,predictor,covariates))
  miss <- setdiff(vars,names(dat))
  if (length(miss)) stop(paste(model_name,"missing columns:",paste(miss,collapse=",")))
  d <- dat[,vars,drop=FALSE]
  for (v in vars) d[[v]] <- safe_numeric(d[[v]])
  d <- d[complete.cases(d),,drop=FALSE]
  if (nrow(d) < max(20, length(vars)+5)) {
    return(data.frame(
      model=model_name, endpoint=endpoint, method="linear_regression",
      status="insufficient_n", n=nrow(d), events=NA,
      beta=NA,se=NA,statistic=NA,p=NA, stringsAsFactors=FALSE
    ))
  }
  rhs <- paste(c(predictor,covariates),collapse=" + ")
  f <- as.formula(paste(endpoint,"~",rhs))
  fit <- lm(f,data=d)
  co <- summary(fit)$coefficients
  if (!predictor %in% rownames(co)) stop(paste("predictor coefficient unavailable in",model_name))
  x <- co[predictor,]
  data.frame(
    model=model_name,endpoint=endpoint,method="linear_regression",status="ok",
    n=nrow(d),events=NA,beta=unname(x[1]),se=unname(x[2]),
    statistic=unname(x[3]),p=unname(x[4]),stringsAsFactors=FALSE
  )
}

fit_logistic <- function(endpoint, model_name) {
  if (is.null(endpoint) || !endpoint %in% names(dat)) {
    return(data.frame(
      model=model_name, endpoint=ifelse(is.null(endpoint),"",endpoint),
      method="logistic_regression", status="endpoint_missing",
      n=NA,events=NA,beta=NA,se=NA,statistic=NA,p=NA,
      stringsAsFactors=FALSE
    ))
  }
  vars <- unique(c(endpoint,predictor,covariates))
  miss <- setdiff(vars,names(dat))
  if (length(miss)) stop(paste(model_name,"missing columns:",paste(miss,collapse=",")))
  d <- dat[,vars,drop=FALSE]
  for (v in vars) d[[v]] <- safe_numeric(d[[v]])
  d <- d[complete.cases(d),,drop=FALSE]
  events <- sum(d[[endpoint]]==1,na.rm=TRUE)
  if (nrow(d) < max(30,length(vars)+10) || events < 5 || events >= nrow(d)) {
    return(data.frame(
      model=model_name,endpoint=endpoint,method="logistic_regression",
      status="insufficient_events_or_n",n=nrow(d),events=events,
      beta=NA,se=NA,statistic=NA,p=NA,stringsAsFactors=FALSE
    ))
  }
  rhs <- paste(c(predictor,covariates),collapse=" + ")
  f <- as.formula(paste(endpoint,"~",rhs))
  fit <- glm(f,data=d,family=binomial())
  co <- summary(fit)$coefficients
  if (!predictor %in% rownames(co)) stop(paste("predictor coefficient unavailable in",model_name))
  x <- co[predictor,]
  data.frame(
    model=model_name,endpoint=endpoint,method="logistic_regression",status="ok",
    n=nrow(d),events=events,beta=unname(x[1]),se=unname(x[2]),
    statistic=unname(x[3]),p=unname(x[4]),stringsAsFactors=FALSE
  )
}

fit_cox <- function() {
  if (is.null(time_to_event) || is.null(event_col)) return(NULL)
  if (!requireNamespace("survival",quietly=TRUE)) {
    return(data.frame(
      model="incident_ckd_cox",endpoint=event_col,method="cox_ph",
      status="survival_package_missing",n=NA,events=NA,
      beta=NA,se=NA,statistic=NA,p=NA,stringsAsFactors=FALSE
    ))
  }
  vars <- unique(c(time_to_event,event_col,predictor,covariates))
  miss <- setdiff(vars,names(dat))
  if (length(miss)) stop(paste("cox missing columns:",paste(miss,collapse=",")))
  d <- dat[,vars,drop=FALSE]
  for (v in vars) d[[v]] <- safe_numeric(d[[v]])
  d <- d[complete.cases(d),,drop=FALSE]
  events <- sum(d[[event_col]]==1,na.rm=TRUE)
  if (nrow(d)<30 || events<5) {
    return(data.frame(
      model="incident_ckd_cox",endpoint=event_col,method="cox_ph",
      status="insufficient_events_or_n",n=nrow(d),events=events,
      beta=NA,se=NA,statistic=NA,p=NA,stringsAsFactors=FALSE
    ))
  }
  rhs <- paste(c(predictor,covariates),collapse=" + ")
  f <- as.formula(paste0("survival::Surv(",time_to_event,",",event_col,") ~ ",rhs))
  fit <- survival::coxph(f,data=d)
  co <- summary(fit)$coefficients
  if (!predictor %in% rownames(co)) stop("predictor coefficient unavailable in cox")
  x <- co[predictor,]
  data.frame(
    model="incident_ckd_cox",endpoint=event_col,method="cox_ph",status="ok",
    n=nrow(d),events=events,beta=unname(x["coef"]),se=unname(x["se(coef)"]),
    statistic=unname(x["z"]),p=unname(x["Pr(>|z|)"]),stringsAsFactors=FALSE
  )
}

fit_lme <- function() {
  if (is.null(longitudinal_file)) return(NULL)
  if (!requireNamespace("lme4",quietly=TRUE)) {
    return(data.frame(
      model="longitudinal_egfr_lme",endpoint=egfr_col,method="linear_mixed_model",
      status="lme4_missing",n=NA,events=NA,beta=NA,se=NA,statistic=NA,p=NA,
      stringsAsFactors=FALSE
    ))
  }
  ld <- read.delim(longitudinal_file,header=TRUE,sep="\t",stringsAsFactors=FALSE,check.names=FALSE)
  vars <- unique(c(id_col,time_col,egfr_col,predictor,covariates))
  miss <- setdiff(vars,names(ld))
  if (length(miss)) stop(paste("LME missing columns:",paste(miss,collapse=",")))
  d <- ld[,vars,drop=FALSE]
  for (v in setdiff(vars,id_col)) d[[v]] <- safe_numeric(d[[v]])
  d <- d[complete.cases(d),,drop=FALSE]
  if (nrow(d)<50 || length(unique(d[[id_col]]))<20) {
    return(data.frame(
      model="longitudinal_egfr_lme",endpoint=egfr_col,method="linear_mixed_model",
      status="insufficient_n",n=nrow(d),events=NA,
      beta=NA,se=NA,statistic=NA,p=NA,stringsAsFactors=FALSE
    ))
  }
  interaction_term <- paste0(time_col,":",predictor)
  rhs <- paste(c(time_col,predictor,paste0(time_col,"*",predictor),covariates,
                 paste0("(1|",id_col,")")),collapse=" + ")
  f <- as.formula(paste(egfr_col,"~",rhs))
  fit <- lme4::lmer(f,data=d,REML=FALSE)
  co <- summary(fit)$coefficients
  rn <- rownames(co)
  hit <- rn[rn %in% c(interaction_term,paste0(predictor,":",time_col))]
  if (!length(hit)) stop("time-by-predictor interaction unavailable in LME")
  x <- co[hit[[1]],]
  z <- unname(x["Estimate"]/x["Std. Error"])
  data.frame(
    model="longitudinal_egfr_lme",endpoint=egfr_col,method="linear_mixed_model",
    status="ok",n=nrow(d),events=NA,beta=unname(x["Estimate"]),
    se=unname(x["Std. Error"]),statistic=z,p=2*pnorm(abs(z),lower.tail=FALSE),
    stringsAsFactors=FALSE
  )
}

rows <- list(
  fit_lm(baseline_endpoint,"baseline_egfr"),
  fit_lm(slope_endpoint,"egfr_slope"),
  fit_logistic(incident_endpoint,"incident_ckd_logistic")
)
cx <- fit_cox(); if (!is.null(cx)) rows[[length(rows)+1]] <- cx
lmix <- fit_lme(); if (!is.null(lmix)) rows[[length(rows)+1]] <- lmix

res <- do.call(rbind,rows)
res$predictor <- predictor
res$effect_scale <- ifelse(res$method %in% c("logistic_regression","cox_ph"),"log_odds_or_log_hazard","linear_beta")
res$effect_ratio <- ifelse(res$method %in% c("logistic_regression","cox_ph") & is.finite(res$beta), exp(res$beta), NA_real_)
res$ci95_low <- ifelse(is.finite(res$beta) & is.finite(res$se), res$beta-1.96*res$se, NA_real_)
res$ci95_high <- ifelse(is.finite(res$beta) & is.finite(res$se), res$beta+1.96*res$se, NA_real_)
res$ratio_ci95_low <- ifelse(res$method %in% c("logistic_regression","cox_ph") & is.finite(res$ci95_low), exp(res$ci95_low), NA_real_)
res$ratio_ci95_high <- ifelse(res$method %in% c("logistic_regression","cox_ph") & is.finite(res$ci95_high), exp(res$ci95_high), NA_real_)

dir.create(dirname(output_file),recursive=TRUE,showWarnings=FALSE)
write.table(res,file=output_file,sep="\t",quote=FALSE,row.names=FALSE,na="")

cat("MASTER_INDIVIDUAL_VALIDATION_PASS\n")
print(res)
