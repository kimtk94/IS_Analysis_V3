args <- commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
    stop(
        "Usage: Rscript run_ckd_koges_final_models.R ",
        "<phenotype.tsv> <genotype.tsv> <output_dir>"
    )
}

phenotype_file <- args[1]
genotype_file <- args[2]
output_dir <- args[3]

dir.create(
    output_dir,
    recursive=TRUE,
    showWarnings=FALSE
)

suppressPackageStartupMessages({
    library(lme4)
    library(survival)
})

pheno <- read.delim(
    phenotype_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

geno <- read.delim(
    genotype_file,
    stringsAsFactors=FALSE,
    check.names=FALSE
)

required_pheno <- c(
    "participant_id",
    "time_years",
    "egfr",
    "sex",
    "baseline_ckd"
)

required_geno <- c(
    "participant_id",
    "gene_symbol",
    "dosage"
)

missing_pheno <- setdiff(
    required_pheno,
    names(pheno)
)

missing_geno <- setdiff(
    required_geno,
    names(geno)
)

if (length(missing_pheno) > 0) {
    stop(
        "Missing phenotype columns: ",
        paste(
            missing_pheno,
            collapse=", "
        )
    )
}

if (length(missing_geno) > 0) {
    stop(
        "Missing genotype columns: ",
        paste(
            missing_geno,
            collapse=", "
        )
    )
}


# ---------------------------------------------------------
# baseline age
# ---------------------------------------------------------

if ("baseline_age" %in% names(pheno)) {

    pheno$baseline_age_final <-
        as.numeric(pheno$baseline_age)

} else if ("age" %in% names(pheno)) {

    pheno$age <- as.numeric(pheno$age)
    pheno$time_years <- as.numeric(
        pheno$time_years
    )

    baseline_age_map <- aggregate(
        age ~ participant_id,
        data=pheno[
            order(
                pheno$participant_id,
                pheno$time_years
            ),
        ],
        FUN=function(x) x[1]
    )

    names(baseline_age_map)[2] <-
        "baseline_age_final"

    pheno <- merge(
        pheno,
        baseline_age_map,
        by="participant_id",
        all.x=TRUE,
        sort=FALSE
    )

} else {

    stop(
        "Neither baseline_age nor age available"
    )
}


pheno$time_years <- as.numeric(
    pheno$time_years
)

pheno$egfr <- as.numeric(
    pheno$egfr
)

geno$dosage <- as.numeric(
    geno$dosage
)


# ---------------------------------------------------------
# covariates
# ---------------------------------------------------------

pc_cols <- intersect(
    paste0("PC", 1:10),
    names(pheno)
)

optional_covars <- intersect(
    c(
        "site",
        "batch"
    ),
    names(pheno)
)

fixed_covars <- c(
    "baseline_age_final",
    "sex",
    pc_cols,
    optional_covars
)


# ---------------------------------------------------------
# helper
# ---------------------------------------------------------

extract_coef <- function(
    fit,
    possible_names
) {

    cf <- summary(fit)$coefficients

    rn <- rownames(cf)

    hit <- intersect(
        possible_names,
        rn
    )

    if (length(hit) == 0) {
        return(NULL)
    }

    x <- cf[hit[1], ]

    list(
        term=hit[1],
        beta=unname(x[1]),
        se=unname(x[2]),
        stat=unname(x[3]),
        p=ifelse(
            ncol(cf) >= 4,
            unname(x[4]),
            2 * pnorm(
                abs(unname(x[3])),
                lower.tail=FALSE
            )
        )
    )
}


fit_lmer_safe <- function(
    dat,
    fixed_formula
) {

    slope_formula <- as.formula(
        paste0(
            fixed_formula,
            " + (1 + time_years | participant_id)"
        )
    )

    intercept_formula <- as.formula(
        paste0(
            fixed_formula,
            " + (1 | participant_id)"
        )
    )

    warnings <- character()

    fit1 <- withCallingHandlers(

        tryCatch(
            lmer(
                slope_formula,
                data=dat,
                REML=FALSE,
                control=lmerControl(
                    optimizer="bobyqa",
                    optCtrl=list(
                        maxfun=200000
                    )
                )
            ),
            error=function(e) e
        ),

        warning=function(w) {
            warnings <<- c(
                warnings,
                conditionMessage(w)
            )
            invokeRestart("muffleWarning")
        }
    )

    fallback_reason <- ""

    if (inherits(fit1, "error")) {

        fallback_reason <- paste0(
            "random_slope_error:",
            conditionMessage(fit1)
        )

    } else {

        singular <- isSingular(
            fit1,
            tol=1e-4
        )

        conv_messages <- fit1@optinfo$conv$lme4$messages

        has_conv <- (
            length(conv_messages) > 0
            || any(
                grepl(
                    "failed to converge|max\\|grad\\|",
                    warnings,
                    ignore.case=TRUE
                )
            )
        )

        if (!singular && !has_conv) {

            return(list(
                fit=fit1,
                structure="random_intercept_slope",
                fallback=FALSE,
                reason="",
                warnings=paste(
                    unique(warnings),
                    collapse=" | "
                )
            ))
        }

        fallback_reason <- paste(
            c(
                if (singular)
                    "singular_fit",
                if (has_conv)
                    "convergence_warning"
            ),
            collapse=";"
        )
    }


    warnings2 <- character()

    fit2 <- withCallingHandlers(

        tryCatch(
            lmer(
                intercept_formula,
                data=dat,
                REML=FALSE,
                control=lmerControl(
                    optimizer="bobyqa",
                    optCtrl=list(
                        maxfun=200000
                    )
                )
            ),
            error=function(e) e
        ),

        warning=function(w) {
            warnings2 <<- c(
                warnings2,
                conditionMessage(w)
            )
            invokeRestart("muffleWarning")
        }
    )


    if (inherits(fit2, "error")) {

        return(list(
            fit=NULL,
            structure="FAILED",
            fallback=TRUE,
            reason=paste0(
                fallback_reason,
                ";random_intercept_error:",
                conditionMessage(fit2)
            ),
            warnings=paste(
                unique(
                    c(
                        warnings,
                        warnings2
                    )
                ),
                collapse=" | "
            )
        ))
    }


    list(
        fit=fit2,
        structure="random_intercept",
        fallback=TRUE,
        reason=fallback_reason,
        warnings=paste(
            unique(
                c(
                    warnings,
                    warnings2
                )
            ),
            collapse=" | "
        )
    )
}


genes <- sort(
    unique(
        geno$gene_symbol
    )
)

longitudinal_results <- list()
baseline_results <- list()
incident_results <- list()
cox_results <- list()


for (gene in genes) {

    message(
        "===== ",
        gene,
        " ====="
    )

    g <- geno[
        geno$gene_symbol == gene,
        c(
            "participant_id",
            "dosage"
        )
    ]

    d <- merge(
        pheno,
        g,
        by="participant_id",
        all=FALSE
    )


    # =====================================================
    # PRIMARY LONGITUDINAL
    # all longitudinally eligible participants
    # =====================================================

    d_long <- d[
        !is.na(d$egfr)
        & !is.na(d$time_years)
        & !is.na(d$dosage),
    ]

    fixed_terms <- c(
        "time_years",
        "dosage",
        "time_years:dosage",
        fixed_covars
    )

    fixed_formula <- paste(
        "egfr ~",
        paste(
            fixed_terms,
            collapse=" + "
        )
    )

    lf <- fit_lmer_safe(
        d_long,
        fixed_formula
    )


    if (is.null(lf$fit)) {

        longitudinal_results[[gene]] <-
            data.frame(
                gene_symbol=gene,
                status="FAIL",
                n_participants=
                    length(
                        unique(
                            d_long$participant_id
                        )
                    ),
                n_observations=
                    nrow(d_long),
                beta_genotype_time=NA,
                se_genotype_time=NA,
                z_genotype_time=NA,
                p_genotype_time=NA,
                random_effect_structure=
                    lf$structure,
                fallback_used=
                    lf$fallback,
                diagnostic=
                    lf$reason,
                warning=
                    lf$warnings,
                stringsAsFactors=FALSE
            )

    } else {

        cc <- extract_coef(
            lf$fit,
            c(
                "time_years:dosage",
                "dosage:time_years"
            )
        )

        longitudinal_results[[gene]] <-
            data.frame(
                gene_symbol=gene,
                status=
                    ifelse(
                        is.null(cc),
                        "FAIL",
                        "PASS"
                    ),
                n_participants=
                    length(
                        unique(
                            d_long$participant_id
                        )
                    ),
                n_observations=
                    nrow(d_long),
                beta_genotype_time=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$beta
                    ),
                se_genotype_time=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$se
                    ),
                z_genotype_time=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$stat
                    ),
                p_genotype_time=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$p
                    ),
                random_effect_structure=
                    lf$structure,
                fallback_used=
                    lf$fallback,
                diagnostic=
                    lf$reason,
                warning=
                    lf$warnings,
                stringsAsFactors=FALSE
            )
    }


    # =====================================================
    # BASELINE eGFR SENSITIVITY
    # =====================================================

    db <- d[
        order(
            d$participant_id,
            d$time_years
        ),
    ]

    db <- db[
        !duplicated(
            db$participant_id
        ),
    ]

    base_terms <- c(
        "dosage",
        fixed_covars
    )

    base_formula <- as.formula(
        paste(
            "egfr ~",
            paste(
                base_terms,
                collapse=" + "
            )
        )
    )

    bf <- tryCatch(
        lm(
            base_formula,
            data=db
        ),
        error=function(e) e
    )


    if (inherits(bf, "error")) {

        baseline_results[[gene]] <-
            data.frame(
                gene_symbol=gene,
                status="FAIL",
                n=nrow(db),
                beta=NA,
                se=NA,
                p=NA,
                message=
                    conditionMessage(bf),
                stringsAsFactors=FALSE
            )

    } else {

        cc <- extract_coef(
            bf,
            "dosage"
        )

        baseline_results[[gene]] <-
            data.frame(
                gene_symbol=gene,
                status=
                    ifelse(
                        is.null(cc),
                        "FAIL",
                        "PASS"
                    ),
                n=nrow(db),
                beta=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$beta
                    ),
                se=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$se
                    ),
                p=
                    ifelse(
                        is.null(cc),
                        NA,
                        cc$p
                    ),
                message="",
                stringsAsFactors=FALSE
            )
    }


    # =====================================================
    # INCIDENT CKD
    # =====================================================

    if (
        "incident_ckd" %in% names(d)
        && "interval_index" %in% names(d)
    ) {

        di <- d[
            d$baseline_ckd == 0
            & !is.na(d$incident_ckd)
            & !is.na(d$dosage),
        ]

        incident_terms <- c(
            "dosage",
            fixed_covars,
            "factor(interval_index)"
        )

        incident_formula <- as.formula(
            paste(
                "incident_ckd ~",
                paste(
                    incident_terms,
                    collapse=" + "
                )
            )
        )

        gf <- tryCatch(
            glm(
                incident_formula,
                data=di,
                family=binomial(
                    link="cloglog"
                )
            ),
            error=function(e) e
        )


        if (inherits(gf, "error")) {

            incident_results[[gene]] <-
                data.frame(
                    gene_symbol=gene,
                    status="FAIL",
                    n_participants=
                        length(
                            unique(
                                di$participant_id
                            )
                        ),
                    n_intervals=nrow(di),
                    events=sum(
                        di$incident_ckd == 1,
                        na.rm=TRUE
                    ),
                    beta=NA,
                    se=NA,
                    hazard_ratio=NA,
                    p=NA,
                    message=
                        conditionMessage(gf),
                    stringsAsFactors=FALSE
                )

        } else {

            cc <- extract_coef(
                gf,
                "dosage"
            )

            incident_results[[gene]] <-
                data.frame(
                    gene_symbol=gene,
                    status=
                        ifelse(
                            is.null(cc),
                            "FAIL",
                            "PASS"
                        ),
                    n_participants=
                        length(
                            unique(
                                di$participant_id
                            )
                        ),
                    n_intervals=nrow(di),
                    events=sum(
                        di$incident_ckd == 1,
                        na.rm=TRUE
                    ),
                    beta=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$beta
                        ),
                    se=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$se
                        ),
                    hazard_ratio=
                        ifelse(
                            is.null(cc),
                            NA,
                            exp(cc$beta)
                        ),
                    p=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$p
                        ),
                    message="",
                    stringsAsFactors=FALSE
                )
        }
    }


    # =====================================================
    # COX SENSITIVITY
    # only when participant-level event time is available
    # =====================================================

    if (
        "event_time_years" %in% names(d)
        && "incident_event" %in% names(d)
    ) {

        dc <- d[
            order(
                d$participant_id,
                d$time_years
            ),
        ]

        dc <- dc[
            !duplicated(
                dc$participant_id
            ),
        ]

        dc <- dc[
            dc$baseline_ckd == 0
            & !is.na(dc$event_time_years)
            & !is.na(dc$incident_event),
        ]

        cox_terms <- c(
            "dosage",
            fixed_covars
        )

        cox_formula <- as.formula(
            paste(
                "Surv(event_time_years, incident_event) ~",
                paste(
                    cox_terms,
                    collapse=" + "
                )
            )
        )

        cf <- tryCatch(
            coxph(
                cox_formula,
                data=dc
            ),
            error=function(e) e
        )


        if (!inherits(cf, "error")) {

            cc <- extract_coef(
                cf,
                "dosage"
            )

            cox_results[[gene]] <-
                data.frame(
                    gene_symbol=gene,
                    status=
                        ifelse(
                            is.null(cc),
                            "FAIL",
                            "PASS"
                        ),
                    n=nrow(dc),
                    events=sum(
                        dc$incident_event == 1,
                        na.rm=TRUE
                    ),
                    beta=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$beta
                        ),
                    se=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$se
                        ),
                    hazard_ratio=
                        ifelse(
                            is.null(cc),
                            NA,
                            exp(cc$beta)
                        ),
                    p=
                        ifelse(
                            is.null(cc),
                            NA,
                            cc$p
                        ),
                    stringsAsFactors=FALSE
                )
        }
    }
}


# ---------------------------------------------------------
# bind + multiplicity
# ---------------------------------------------------------

longitudinal_results <- do.call(
    rbind,
    longitudinal_results
)

baseline_results <- do.call(
    rbind,
    baseline_results
)


longitudinal_results$p_bonferroni <-
    p.adjust(
        longitudinal_results$p_genotype_time,
        method="bonferroni"
    )

longitudinal_results$p_fdr <-
    p.adjust(
        longitudinal_results$p_genotype_time,
        method="BH"
    )


baseline_results$p_bonferroni <-
    p.adjust(
        baseline_results$p,
        method="bonferroni"
    )

baseline_results$p_fdr <-
    p.adjust(
        baseline_results$p,
        method="BH"
    )


write.table(
    longitudinal_results,
    file=file.path(
        output_dir,
        "KOGES_LONGITUDINAL_EGFR_RESULTS.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

write.table(
    baseline_results,
    file=file.path(
        output_dir,
        "KOGES_BASELINE_EGFR_RESULTS.tsv"
    ),
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)


if (length(incident_results) > 0) {

    incident_results <- do.call(
        rbind,
        incident_results
    )

    incident_results$p_bonferroni <-
        p.adjust(
            incident_results$p,
            method="bonferroni"
        )

    incident_results$p_fdr <-
        p.adjust(
            incident_results$p,
            method="BH"
        )

    write.table(
        incident_results,
        file=file.path(
            output_dir,
            "KOGES_INCIDENT_CKD_RESULTS.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
}


if (length(cox_results) > 0) {

    cox_results <- do.call(
        rbind,
        cox_results
    )

    cox_results$p_bonferroni <-
        p.adjust(
            cox_results$p,
            method="bonferroni"
        )

    cox_results$p_fdr <-
        p.adjust(
            cox_results$p,
            method="BH"
        )

    write.table(
        cox_results,
        file=file.path(
            output_dir,
            "KOGES_COX_SENSITIVITY_RESULTS.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )
}


message(
    "KOGES FINAL MODEL ENGINE COMPLETE"
)
