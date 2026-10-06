#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop(
        paste0(
            "Usage: Rscript run_ckd_koges_longitudinal_models.R ",
            "<phenotype.tsv> <dosage.tsv> <output_dir>"
        )
    )
}

phenotype_file <- args[[1]]
dosage_file <- args[[2]]
output_dir <- args[[3]]

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


required_packages <- c(
    "lme4",
    "survival"
)

missing_packages <- required_packages[
    !vapply(
        required_packages,
        requireNamespace,
        logical(1),
        quietly = TRUE
    )
]

if (length(missing_packages) > 0) {
    stop(
        paste(
            "Missing R packages:",
            paste(missing_packages, collapse = ", ")
        )
    )
}


# ------------------------------------------------------------
# READ
# ------------------------------------------------------------

pheno <- read.delim(
    phenotype_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)

geno <- read.delim(
    dosage_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_pheno <- c(
    "participant_id",
    "visit_id",
    "time_years",
    "age",
    "sex",
    "egfr",
    "baseline_egfr",
    "baseline_ckd"
)

required_geno <- c(
    "participant_id",
    "gene_symbol",
    "target_rsid",
    "analysis_rsid",
    "effect_allele",
    "other_allele",
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
        paste(
            "Missing phenotype columns:",
            paste(missing_pheno, collapse = ", ")
        )
    )
}

if (length(missing_geno) > 0) {
    stop(
        paste(
            "Missing genotype columns:",
            paste(missing_geno, collapse = ", ")
        )
    )
}


# ------------------------------------------------------------
# BASIC QC
# ------------------------------------------------------------

pheno$participant_id <- as.character(
    pheno$participant_id
)

geno$participant_id <- as.character(
    geno$participant_id
)

genes <- sort(
    unique(geno$gene_symbol)
)

message(
    "Phenotype rows: ",
    nrow(pheno)
)

message(
    "Participants: ",
    length(unique(pheno$participant_id))
)

message(
    "Genes: ",
    paste(genes, collapse = ", ")
)


# ------------------------------------------------------------
# COVARIATES
# ------------------------------------------------------------

pc_candidates <- paste0(
    "PC",
    1:10
)

pcs <- intersect(
    pc_candidates,
    names(pheno)
)

optional_covariates <- c(
    pcs,
    intersect("site", names(pheno)),
    intersect("batch", names(pheno))
)


# ------------------------------------------------------------
# LONGITUDINAL eGFR
#
# Primary coefficient:
# dosage:time_years
#
# Negative beta:
# effect allele associated with faster decline.
# ------------------------------------------------------------

long_results <- list()


for (gene in genes) {

    g <- geno[
        geno$gene_symbol == gene,
        ,
        drop = FALSE
    ]

    # One dosage per participant/gene expected.
    if (anyDuplicated(g$participant_id)) {
        stop(
            paste(
                "Duplicate participant dosage:",
                gene
            )
        )
    }

    d <- merge(
        pheno,
        g,
        by = "participant_id",
        all = FALSE
    )

    d <- d[
        d$baseline_ckd == 0 &
        !is.na(d$egfr) &
        !is.na(d$time_years) &
        !is.na(d$dosage),
        ,
        drop = FALSE
    ]

    # Primary adjustment:
    # age + sex + PCs + optional site/batch.
    fixed_covars <- c(
        "age",
        "sex",
        optional_covariates
    )

    fixed_covars <- unique(
        fixed_covars[
            fixed_covars %in% names(d)
        ]
    )

    rhs <- paste(
        c(
            "time_years",
            "dosage",
            "time_years:dosage",
            fixed_covars
        ),
        collapse = " + "
    )

    formula_text <- paste0(
        "egfr ~ ",
        rhs,
        " + (1 + time_years | participant_id)"
    )

    fit <- tryCatch(
        lme4::lmer(
            as.formula(formula_text),
            data = d,
            REML = FALSE
        ),
        error = function(e) e
    )

    random_effect_structure <- "random_intercept_slope"

    # Sparse or weakly varying real datasets can produce a singular
    # random-slope fit. Fall back to random intercept rather than
    # treating this automatically as biological model failure.
    if (
        !inherits(fit, "error") &&
        lme4::isSingular(fit, tol = 1e-4)
    ) {

        formula_text_fallback <- paste0(
            "egfr ~ ",
            rhs,
            " + (1 | participant_id)"
        )

        fit2 <- tryCatch(
            lme4::lmer(
                as.formula(formula_text_fallback),
                data = d,
                REML = FALSE
            ),
            error = function(e) e
        )

        if (!inherits(fit2, "error")) {
            fit <- fit2
            formula_text <- formula_text_fallback
            random_effect_structure <- "random_intercept_fallback"
        }
    }

    if (inherits(fit, "error")) {

        long_results[[gene]] <- data.frame(
            gene_symbol = gene,
            status = "MODEL_FAILED",
            n_participants =
                length(unique(d$participant_id)),
            n_observations = nrow(d),
            beta_genotype_time = NA,
            se_genotype_time = NA,
            z_genotype_time = NA,
            p_genotype_time = NA,
            model = formula_text,
            message = conditionMessage(fit),
            stringsAsFactors = FALSE
        )

        next
    }

    cf <- summary(fit)$coefficients

    term_candidates <- c(
        "time_years:dosage",
        "dosage:time_years"
    )

    term <- term_candidates[
        term_candidates %in% rownames(cf)
    ]

    if (length(term) != 1) {

        long_results[[gene]] <- data.frame(
            gene_symbol = gene,
            status = "INTERACTION_TERM_MISSING",
            n_participants =
                length(unique(d$participant_id)),
            n_observations = nrow(d),
            beta_genotype_time = NA,
            se_genotype_time = NA,
            z_genotype_time = NA,
            p_genotype_time = NA,
            model = formula_text,
            message = "",
            stringsAsFactors = FALSE
        )

        next
    }

    term <- term[[1]]

    beta <- cf[term, "Estimate"]
    se <- cf[term, "Std. Error"]
    z <- beta / se
    p <- 2 * pnorm(
        abs(z),
        lower.tail = FALSE
    )

    long_results[[gene]] <- data.frame(
        gene_symbol = gene,
        status = "PASS",
        n_participants =
            length(unique(d$participant_id)),
        n_observations = nrow(d),
        beta_genotype_time = beta,
        se_genotype_time = se,
        z_genotype_time = z,
        p_genotype_time = p,
        model = formula_text,
        random_effect_structure = random_effect_structure,
        message = "",
        stringsAsFactors = FALSE
    )
}


long_results <- do.call(
    rbind,
    long_results
)

write.table(
    long_results,
    file.path(
        output_dir,
        "KOGES_LONGITUDINAL_EGFR_RESULTS.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ------------------------------------------------------------
# INCIDENT CKD
#
# Requires preconstructed event rows:
#
# incident_ckd
# interval_index
#
# Each row represents an at-risk interval.
#
# cloglog link approximates discrete-time hazard.
# ------------------------------------------------------------

incident_results <- list()


if (
    "incident_ckd" %in% names(pheno) &&
    "interval_index" %in% names(pheno)
) {

    for (gene in genes) {

        g <- geno[
            geno$gene_symbol == gene,
            ,
            drop = FALSE
        ]

        d <- merge(
            pheno,
            g,
            by = "participant_id",
            all = FALSE
        )

        d <- d[
            d$baseline_ckd == 0 &
            !is.na(d$incident_ckd) &
            !is.na(d$interval_index) &
            !is.na(d$dosage),
            ,
            drop = FALSE
        ]

        covars <- c(
            "dosage",
            "age",
            "sex",
            optional_covariates,
            "factor(interval_index)"
        )

        covars <- unique(covars)

        formula_text <- paste(
            "incident_ckd ~",
            paste(
                covars,
                collapse = " + "
            )
        )

        fit <- tryCatch(
            glm(
                as.formula(formula_text),
                family = binomial(
                    link = "cloglog"
                ),
                data = d
            ),
            error = function(e) e
        )

        if (inherits(fit, "error")) {

            incident_results[[gene]] <- data.frame(
                gene_symbol = gene,
                status = "MODEL_FAILED",
                n_participants =
                    length(unique(d$participant_id)),
                n_intervals = nrow(d),
                events =
                    sum(d$incident_ckd == 1),
                beta = NA,
                se = NA,
                hazard_ratio = NA,
                p = NA,
                model = formula_text,
                message = conditionMessage(fit),
                stringsAsFactors = FALSE
            )

            next
        }

        cf <- summary(fit)$coefficients

        if (!"dosage" %in% rownames(cf)) {
            next
        }

        beta <- cf[
            "dosage",
            "Estimate"
        ]

        se <- cf[
            "dosage",
            "Std. Error"
        ]

        p <- cf[
            "dosage",
            "Pr(>|z|)"
        ]

        incident_results[[gene]] <- data.frame(
            gene_symbol = gene,
            status = "PASS",
            n_participants =
                length(unique(d$participant_id)),
            n_intervals = nrow(d),
            events =
                sum(d$incident_ckd == 1),
            beta = beta,
            se = se,
            hazard_ratio = exp(beta),
            p = p,
            model = formula_text,
            message = "",
            stringsAsFactors = FALSE
        )
    }
}


if (length(incident_results) > 0) {

    incident_results <- do.call(
        rbind,
        incident_results
    )

    write.table(
        incident_results,
        file.path(
            output_dir,
            "KOGES_INCIDENT_CKD_RESULTS.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


message("KOGES MODEL ENGINE COMPLETE")
