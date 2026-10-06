#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase8b/$RUN_ID"
STATUS="$LOGROOT/STATUS.tsv"

OUTROOT="$ROOT/results/is/stage4_cross_eas/phase8_qc"

mkdir -p "$LOGROOT" "$OUTROOT"

cd "$REPO" || exit 1

printf "step\tstatus\tmessage\tlog\n" > "$STATUS"

record() {
    printf "%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" >> "$STATUS"
}

run_step() {

    STEP="$1"
    shift

    LOG="$LOGROOT/${STEP}.log"

    echo
    echo "===================================================="
    echo "STEP=$STEP"
    echo "START=$(date)"
    echo "LOG=$LOG"
    echo "===================================================="

    (
        "$@"
    ) > >(tee "$LOG") 2>&1

    RC=$?

    if [ "$RC" -eq 0 ]; then

        echo "[PASS] $STEP"

        record \
            "$STEP" \
            "PASS" \
            "completed" \
            "$LOG"

    else

        echo "[FAIL] $STEP rc=$RC"

        record \
            "$STEP" \
            "FAIL" \
            "rc=$RC" \
            "$LOG"

    fi

    return 0
}


# ============================================================
# 01. SAMPLE SIZE AUDIT
# ============================================================

sample_size_audit() {

    FILE="$OUTROOT/GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"

    if [ ! -s "$FILE" ]; then
        echo "SAMPLE_SIZE_MODEL_MISSING"
        return 11
    fi

    echo "===== SAMPLE SIZE MODEL ====="

    column -t -s $'\t' "$FILE"

    echo
    echo "===== EXPECTED ====="

    python3 - <<'PY'
import csv
from pathlib import Path

p = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/phase8_qc/"
    "GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"
)

expected = {
    "AS": 264655,
    "AIS": 256274,
    "CES": 238168,
    "LAS": 238977,
    "SVS": 242774,
}

seen = {}

with p.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        try:
            seen[
                r["phenotype"]
            ] = int(
                float(
                    r["n_used"]
                )
            )
        except Exception:
            pass


for phenotype, expected_n in expected.items():

    actual = seen.get(
        phenotype
    )

    print(
        phenotype,
        "EXPECTED=",
        expected_n,
        "ACTUAL=",
        actual,
        "PASS="
        + str(
            actual == expected_n
        )
    )


bad = [
    p
    for p, n in expected.items()
    if seen.get(p) != n
]

raise SystemExit(
    0 if not bad else 1
)
PY

    return $?
}


# ============================================================
# 02. BBJ LD DIAGNOSTICS
# ============================================================

bbj_ld_diagnostics() {

Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)


ROOT <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

OUTDIR <- paste0(
    "/srv/is-analysis/results/is/",
    "stage4_cross_eas/phase8_qc/",
    "bbj_ld_diagnostics"
)

dir.create(
    OUTDIR,
    recursive=TRUE,
    showWarnings=FALSE
)


loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)


N_CASE <- 22664
N_CTRL <- 152022

N_TOTAL <- N_CASE + N_CTRL

N_EFF <- round(
    4 /
    (
        1/N_CASE +
        1/N_CTRL
    )
)


read_ld <- function(path, n) {

    expected <- n*n*8

    actual <- file.info(path)$size

    if (
        is.na(actual)
        || actual != expected
    ) {
        stop(
            paste0(
                "LD size mismatch ",
                actual,
                " != ",
                expected
            )
        )
    }

    con <- file(
        path,
        "rb"
    )

    on.exit(
        close(con)
    )

    x <- readBin(
        con,
        what="numeric",
        n=n*n,
        size=8,
        endian="little"
    )

    matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )
}


master <- list()


for (locus in loci) {

    cat(
        "\n=============================\n",
        locus,
        "\n=============================\n",
        sep=""
    )


    result <- tryCatch({

        ssfile <- file.path(
            ROOT,
            "susie_inputs_v3",
            paste0(
                locus,
                ".tsv"
            )
        )

        varsfile <- file.path(
            ROOT,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin.vars"
            )
        )

        ldfile <- file.path(
            ROOT,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin"
            )
        )


        ss <- fread(
            ssfile
        )

        vars <- readLines(
            varsfile
        )

        nvar <- length(vars)


        if (
            nrow(ss) != nvar
            || !all(
                ss$variant_id == vars
            )
        ) {
            stop(
                "summary/LD order mismatch"
            )
        }


        R <- read_ld(
            ldfile,
            nvar
        )

        R <- (
            R +
            t(R)
        ) / 2

        diag(R) <- 1


        z <- (
            as.numeric(
                ss$beta
            ) /
            as.numeric(
                ss$se
            )
        )


        s_neff <- estimate_s_rss(
            z=z,
            R=R,
            n=N_EFF,
            method="null-mle"
        )

        s_total <- estimate_s_rss(
            z=z,
            R=R,
            n=N_TOTAL,
            method="null-mle"
        )

        s_no_n <- estimate_s_rss(
            z=z,
            R=R,
            method="null-mle"
        )


        kr <- kriging_rss(
            z=z,
            R=R,
            n=N_EFF,
            s=s_neff
        )


        tab <- as.data.table(
            kr$table
        )

        tab[
            ,
            variant_id := vars
        ]


        fwrite(
            tab,
            file.path(
                OUTDIR,
                paste0(
                    locus,
                    ".KRIGING.tsv"
                )
            ),
            sep="\t"
        )


        lr_cols <- grep(
            "log.*lr|loglr",
            names(tab),
            ignore.case=TRUE,
            value=TRUE
        )


        n_switch <- NA_integer_


        if (
            length(lr_cols) >= 1
        ) {

            lr <- as.numeric(
                tab[[lr_cols[1]]]
            )

            n_switch <- sum(
                is.finite(lr)
                & lr > 2
                & abs(z) > 2,
                na.rm=TRUE
            )
        }


        cat(
            "s_neff=",
            s_neff,
            "\n",
            sep=""
        )

        cat(
            "s_total=",
            s_total,
            "\n",
            sep=""
        )

        cat(
            "s_no_n=",
            s_no_n,
            "\n",
            sep=""
        )

        cat(
            "possible_switch_n=",
            n_switch,
            "\n",
            sep=""
        )


        list(
            locus=locus,
            status="PASS",
            n_variants=nvar,
            n_eff=N_EFF,
            n_total=N_TOTAL,
            s_neff=s_neff,
            s_total=s_total,
            s_no_n=s_no_n,
            possible_switch_n=n_switch,
            error=""
        )

    }, error=function(e) {

        cat(
            "ERROR=",
            conditionMessage(e),
            "\n",
            sep=""
        )

        list(
            locus=locus,
            status="FAIL",
            n_variants=NA_integer_,
            n_eff=N_EFF,
            n_total=N_TOTAL,
            s_neff=NA_real_,
            s_total=NA_real_,
            s_no_n=NA_real_,
            possible_switch_n=NA_integer_,
            error=conditionMessage(e)
        )
    })


    master[[length(master) + 1]] <- result

}


out <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)


fwrite(
    out,
    file.path(
        OUTDIR,
        "BBJ_LD_DIAGNOSTIC_SUMMARY.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== BBJ LD DIAGNOSTIC =====\n"
)

print(out)

RS

    return $?
}


# ============================================================
# 03. GIGASTROKE LD DIAGNOSTICS
# ============================================================

gigastroke_ld_diagnostics() {

Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)


ROOT <- "/srv/is-analysis/results/is"

BBJ_ROOT <- file.path(
    ROOT,
    "stage3_finemap/japan/bbj"
)

HARM_ROOT <- file.path(
    ROOT,
    "stage4_cross_eas/harmonized"
)

META_FILE <- file.path(
    ROOT,
    "stage4_cross_eas/phase8_qc/",
    "GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"
)

OUTDIR <- file.path(
    ROOT,
    "stage4_cross_eas/phase8_qc/",
    "gigastroke_ld_diagnostics"
)

dir.create(
    OUTDIR,
    recursive=TRUE,
    showWarnings=FALSE
)


meta <- as.data.frame(
    fread(
        META_FILE
    )
)


phenotypes <- c(
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
)

loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)


valid_harm <- c(
    "MATCH",
    "SWAP",
    "COMPLEMENT",
    "COMPLEMENT_SWAP"
)


read_ld <- function(path, n) {

    con <- file(
        path,
        "rb"
    )

    on.exit(
        close(con)
    )

    x <- readBin(
        con,
        what="numeric",
        n=n*n,
        size=8,
        endian="little"
    )

    if (
        length(x) != n*n
    ) {
        stop(
            "LD size mismatch"
        )
    }

    matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )
}


master <- list()


for (locus in loci) {

    varsfile <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    ldfile <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin"
        )
    )


    full_vars <- readLines(
        varsfile
    )

    R_full <- read_ld(
        ldfile,
        length(
            full_vars
        )
    )

    R_full <- (
        R_full +
        t(R_full)
    ) / 2

    diag(R_full) <- 1


    index_map <- setNames(
        seq_along(
            full_vars
        ),
        full_vars
    )


    for (phenotype_name in phenotypes) {

        cat(
            "\n=============================\n",
            phenotype_name,
            " / ",
            locus,
            "\n=============================\n",
            sep=""
        )


        result <- tryCatch({

            m <- meta[
                meta$phenotype == phenotype_name,
                ,
                drop=FALSE
            ]


            if (nrow(m) != 1) {
                stop(
                    paste0(
                        "metadata row count=",
                        nrow(m)
                    )
                )
            }


            n_used <- suppressWarnings(
                as.numeric(
                    m$n_used[1]
                )
            )

            n_mode <- as.character(
                m$n_mode[1]
            )


            if (
                !is.finite(n_used)
                || n_used <= 0
            ) {
                stop(
                    "sample size missing"
                )
            }


            infile <- file.path(
                HARM_ROOT,
                paste0(
                    phenotype_name,
                    ".",
                    locus,
                    ".harmonized.tsv"
                )
            )


            d <- fread(
                infile
            )

            d <- d[
                harmonization %in%
                valid_harm
            ]

            d <- d[
                variant_id %in%
                full_vars
            ]

            d <- d[
                !duplicated(
                    variant_id
                )
            ]


            idx <- unname(
                index_map[
                    d$variant_id
                ]
            )

            keep <- !is.na(
                idx
            )

            d <- d[
                keep
            ]

            idx <- as.integer(
                idx[keep]
            )


            ord <- order(
                idx
            )

            d <- d[
                ord
            ]

            idx <- idx[
                ord
            ]


            R <- R_full[
                idx,
                idx,
                drop=FALSE
            ]


            beta <- as.numeric(
                d$giga_beta_aligned
            )

            p <- as.numeric(
                d$giga_p
            )


            ok <- (
                is.finite(beta)
                & is.finite(p)
                & p > 0
                & p <= 1
            )


            d <- d[
                ok
            ]

            R <- R[
                ok,
                ok,
                drop=FALSE
            ]

            beta <- beta[
                ok
            ]

            p <- p[
                ok
            ]


            z <- sign(
                beta
            ) * qnorm(
                pmax(
                    p,
                    .Machine$double.xmin
                ) / 2,
                lower.tail=FALSE
            )


            s_with_n <- estimate_s_rss(
                z=z,
                R=R,
                n=n_used,
                method="null-mle"
            )

            s_without_n <- estimate_s_rss(
                z=z,
                R=R,
                method="null-mle"
            )


            kr <- kriging_rss(
                z=z,
                R=R,
                n=n_used,
                s=s_with_n
            )


            tab <- as.data.table(
                kr$table
            )

            tab[
                ,
                variant_id := d$variant_id
            ]


            KRDIR <- file.path(
                OUTDIR,
                "kriging"
            )

            dir.create(
                KRDIR,
                recursive=TRUE,
                showWarnings=FALSE
            )


            fwrite(
                tab,
                file.path(
                    KRDIR,
                    paste0(
                        phenotype_name,
                        ".",
                        locus,
                        ".tsv"
                    )
                ),
                sep="\t"
            )


            lr_cols <- grep(
                "log.*lr|loglr",
                names(tab),
                ignore.case=TRUE,
                value=TRUE
            )


            n_switch <- NA_integer_


            if (
                length(lr_cols) >= 1
            ) {

                lr <- as.numeric(
                    tab[[lr_cols[1]]]
                )

                n_switch <- sum(
                    is.finite(lr)
                    & lr > 2
                    & abs(z) > 2,
                    na.rm=TRUE
                )
            }


            cat(
                "n_used=",
                n_used,
                "\n",
                sep=""
            )

            cat(
                "n_mode=",
                n_mode,
                "\n",
                sep=""
            )

            cat(
                "s_with_n=",
                s_with_n,
                "\n",
                sep=""
            )

            cat(
                "s_without_n=",
                s_without_n,
                "\n",
                sep=""
            )

            cat(
                "possible_switch_n=",
                n_switch,
                "\n",
                sep=""
            )


            list(
                phenotype=
                    phenotype_name,

                locus=
                    locus,

                status=
                    "PASS",

                n_variants=
                    length(z),

                n_used=
                    n_used,

                n_mode=
                    n_mode,

                s_with_n=
                    s_with_n,

                s_without_n=
                    s_without_n,

                possible_switch_n=
                    n_switch,

                error=
                    ""
            )

        }, error=function(e) {

            cat(
                "ERROR=",
                conditionMessage(e),
                "\n",
                sep=""
            )

            list(
                phenotype=
                    phenotype_name,

                locus=
                    locus,

                status=
                    "FAIL",

                n_variants=
                    NA_integer_,

                n_used=
                    NA_real_,

                n_mode=
                    "",

                s_with_n=
                    NA_real_,

                s_without_n=
                    NA_real_,

                possible_switch_n=
                    NA_integer_,

                error=
                    conditionMessage(e)
            )
        })


        master[[length(master) + 1]] <- result

    }
}


out <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)


fwrite(
    out,
    file.path(
        OUTDIR,
        "GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== GIGASTROKE LD DIAGNOSTIC =====\n"
)

print(out)

RS

    return $?
}


# ============================================================
# 04. CORRECT N-AWARE SUSIE
# ============================================================

gigastroke_naware_susie() {

Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)


ROOT <- "/srv/is-analysis/results/is"

BBJ_ROOT <- file.path(
    ROOT,
    "stage3_finemap/japan/bbj"
)

HARM_ROOT <- file.path(
    ROOT,
    "stage4_cross_eas/harmonized"
)

META_FILE <- file.path(
    ROOT,
    "stage4_cross_eas/phase8_qc/",
    "GIGASTROKE_SAMPLE_SIZE_MODEL.tsv"
)

OUT_ROOT <- file.path(
    ROOT,
    "stage4_cross_eas/phase8_qc/",
    "susie_naware_corrected"
)


dir.create(
    OUT_ROOT,
    recursive=TRUE,
    showWarnings=FALSE
)


meta <- as.data.frame(
    fread(
        META_FILE
    )
)


phenotypes <- c(
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
)

loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)


valid_harm <- c(
    "MATCH",
    "SWAP",
    "COMPLEMENT",
    "COMPLEMENT_SWAP"
)


read_ld <- function(path, n) {

    con <- file(
        path,
        "rb"
    )

    on.exit(
        close(con)
    )

    x <- readBin(
        con,
        what="numeric",
        n=n*n,
        size=8,
        endian="little"
    )

    matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )
}


master <- list()


for (locus in loci) {

    varsfile <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    ldfile <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin"
        )
    )


    full_vars <- readLines(
        varsfile
    )


    R_full <- read_ld(
        ldfile,
        length(
            full_vars
        )
    )

    R_full <- (
        R_full +
        t(R_full)
    ) / 2

    diag(R_full) <- 1


    index_map <- setNames(
        seq_along(
            full_vars
        ),
        full_vars
    )


    for (phenotype_name in phenotypes) {

        cat(
            "\n=============================\n",
            phenotype_name,
            " / ",
            locus,
            "\n=============================\n",
            sep=""
        )


        result <- tryCatch({

            m <- meta[
                meta$phenotype == phenotype_name,
                ,
                drop=FALSE
            ]


            if (nrow(m) != 1) {
                stop(
                    paste0(
                        "metadata rows=",
                        nrow(m)
                    )
                )
            }


            n_used <- suppressWarnings(
                as.numeric(
                    m$n_used[1]
                )
            )

            n_mode <- as.character(
                m$n_mode[1]
            )


            if (
                !is.finite(n_used)
                || n_used <= 0
            ) {
                stop(
                    "sample size unavailable"
                )
            }


            infile <- file.path(
                HARM_ROOT,
                paste0(
                    phenotype_name,
                    ".",
                    locus,
                    ".harmonized.tsv"
                )
            )


            d <- fread(
                infile
            )

            d <- d[
                harmonization %in%
                valid_harm
            ]

            d <- d[
                variant_id %in%
                full_vars
            ]

            d <- d[
                !duplicated(
                    variant_id
                )
            ]


            idx <- unname(
                index_map[
                    d$variant_id
                ]
            )

            keep <- !is.na(
                idx
            )

            d <- d[
                keep
            ]

            idx <- as.integer(
                idx[keep]
            )


            ord <- order(
                idx
            )

            d <- d[
                ord
            ]

            idx <- idx[
                ord
            ]


            R <- R_full[
                idx,
                idx,
                drop=FALSE
            ]


            beta <- as.numeric(
                d$giga_beta_aligned
            )

            p <- as.numeric(
                d$giga_p
            )


            ok <- (
                is.finite(beta)
                & is.finite(p)
                & p > 0
                & p <= 1
            )


            d <- d[
                ok
            ]

            R <- R[
                ok,
                ok,
                drop=FALSE
            ]

            beta <- beta[
                ok
            ]

            p <- p[
                ok
            ]


            z <- sign(
                beta
            ) * qnorm(
                pmax(
                    p,
                    .Machine$double.xmin
                ) / 2,
                lower.tail=FALSE
            )


            fit <- susie_rss(
                z=z,
                R=R,
                n=n_used,
                L=10,
                estimate_residual_variance=FALSE,
                coverage=0.95,
                min_abs_corr=0.5,
                max_iter=1000,
                verbose=FALSE
            )


            outdir <- file.path(
                OUT_ROOT,
                phenotype_name,
                locus
            )


            dir.create(
                outdir,
                recursive=TRUE,
                showWarnings=FALSE
            )


            pip <- data.table(
                variant_id=
                    d$variant_id,

                pip=
                    fit$pip
            )


            setorder(
                pip,
                -pip
            )


            fwrite(
                pip,
                file.path(
                    outdir,
                    "PIP.tsv"
                ),
                sep="\t"
            )


            cs <- fit$sets$cs

            cs_table <- data.table(
                credible_set=integer(),
                variant_id=character(),
                pip=numeric()
            )


            if (
                !is.null(cs)
                && length(cs) > 0
            ) {

                for (
                    i in seq_along(cs)
                ) {

                    ci <- cs[[i]]

                    cs_table <- rbind(
                        cs_table,
                        data.table(
                            credible_set=
                                i,

                            variant_id=
                                d$variant_id[ci],

                            pip=
                                fit$pip[ci]
                        )
                    )
                }
            }


            fwrite(
                cs_table,
                file.path(
                    outdir,
                    "CREDIBLE_SETS.tsv"
                ),
                sep="\t"
            )


            top <- which.max(
                fit$pip
            )


            ncs <- if (
                is.null(cs)
            ) {
                0L
            } else {
                length(cs)
            }


            cat(
                "n_used=",
                n_used,
                "\n",
                sep=""
            )

            cat(
                "n_mode=",
                n_mode,
                "\n",
                sep=""
            )

            cat(
                "top_variant=",
                d$variant_id[top],
                "\n",
                sep=""
            )

            cat(
                "top_pip=",
                fit$pip[top],
                "\n",
                sep=""
            )

            cat(
                "n_cs=",
                ncs,
                "\n",
                sep=""
            )


            list(
                phenotype=
                    phenotype_name,

                locus=
                    locus,

                status=
                    "PASS",

                converged=
                    fit$converged,

                n_used=
                    n_used,

                n_mode=
                    n_mode,

                n_variants=
                    length(z),

                n_cs=
                    ncs,

                top_variant=
                    d$variant_id[top],

                top_pip=
                    fit$pip[top],

                error=
                    ""
            )

        }, error=function(e) {

            cat(
                "ERROR=",
                conditionMessage(e),
                "\n",
                sep=""
            )

            list(
                phenotype=
                    phenotype_name,

                locus=
                    locus,

                status=
                    "FAIL",

                converged=
                    FALSE,

                n_used=
                    NA_real_,

                n_mode=
                    "",

                n_variants=
                    NA_integer_,

                n_cs=
                    NA_integer_,

                top_variant=
                    NA_character_,

                top_pip=
                    NA_real_,

                error=
                    conditionMessage(e)
            )
        })


        master[[length(master) + 1]] <- result

    }
}


out <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)


fwrite(
    out,
    file.path(
        OUTROOT,
        "GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== CORRECTED N-AWARE SUSIE =====\n"
)

print(out)

RS

    return $?
}


# ============================================================
# 05. Z-ONLY VS CORRECTED N-AWARE
# ============================================================

compare_sensitivity() {

python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas"
)

zroot = (
    root /
    "fine_mapping/susie_zonly"
)

nroot = (
    root /
    "phase8_qc/"
    "susie_naware_corrected"
)

out = (
    root /
    "phase8_qc/"
    "GIGASTROKE_ZONLY_VS_NAWARE_CORRECTED.tsv"
)


phenotypes = [
    "AS",
    "AIS",
    "CES",
    "LAS",
    "SVS"
]

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
]


def read_pip(path):

    d = {}

    if not path.exists():
        return d

    with path.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            try:
                d[
                    r["variant_id"]
                ] = float(
                    r["pip"]
                )
            except Exception:
                pass

    return d


def read_cs(path):

    s = set()

    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        return s

    with path.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            v = r.get(
                "variant_id",
                ""
            )

            if v:
                s.add(v)

    return s


rows = []


for phenotype in phenotypes:

    for locus in loci:

        zp = read_pip(
            zroot /
            phenotype /
            locus /
            "PIP.tsv"
        )

        np = read_pip(
            nroot /
            phenotype /
            locus /
            "PIP.tsv"
        )


        zcs = read_cs(
            zroot /
            phenotype /
            locus /
            "CREDIBLE_SETS.tsv"
        )

        ncs = read_cs(
            nroot /
            phenotype /
            locus /
            "CREDIBLE_SETS.tsv"
        )


        ztop = (
            max(
                zp,
                key=zp.get
            )
            if zp
            else ""
        )

        ntop = (
            max(
                np,
                key=np.get
            )
            if np
            else ""
        )


        inter = (
            zcs &
            ncs
        )

        union = (
            zcs |
            ncs
        )


        common = (
            set(zp)
            &
            set(np)
        )


        maxdiff = (
            max(
                abs(
                    zp[v] -
                    np[v]
                )
                for v in common
            )
            if common
            else ""
        )


        rows.append({
            "phenotype":
                phenotype,

            "locus":
                locus,

            "zonly_top":
                ztop,

            "zonly_top_pip":
                zp.get(
                    ztop,
                    ""
                ),

            "naware_top":
                ntop,

            "naware_top_pip":
                np.get(
                    ntop,
                    ""
                ),

            "same_top":
                (
                    "YES"
                    if (
                        ztop
                        and ntop
                        and ztop == ntop
                    )
                    else "NO"
                ),

            "zonly_cs_size":
                len(zcs),

            "naware_cs_size":
                len(ncs),

            "cs_intersection":
                len(inter),

            "cs_union":
                len(union),

            "cs_jaccard":
                (
                    len(inter) /
                    len(union)
                    if union
                    else 1.0
                ),

            "max_abs_pip_diff":
                maxdiff
        })


fields = list(
    rows[0].keys()
)


with out.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(
        rows
    )


print(
    "\t".join(fields)
)

for r in rows:

    print(
        "\t".join(
            str(r[x])
            for x in fields
        )
    )
PY

    return $?
}


# ============================================================
# 06. REBUILD QC SUMMARY
# ============================================================

build_qc_summary() {

python3 - <<'PY'
import csv
import math
from pathlib import Path

base = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/phase8_qc"
)

gdiag = (
    base /
    "gigastroke_ld_diagnostics/"
    "GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv"
)

sens = (
    base /
    "GIGASTROKE_ZONLY_VS_NAWARE_CORRECTED.tsv"
)

out = (
    base /
    "PHASE8B_INTERPRETATION_QC.tsv"
)


D = {}

with gdiag.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        D[
            (
                r["phenotype"],
                r["locus"]
            )
        ] = r


S = {}

with sens.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        S[
            (
                r["phenotype"],
                r["locus"]
            )
        ] = r


rows = []


for key in sorted(
    set(D) | set(S)
):

    phenotype, locus = key

    d = D.get(
        key,
        {}
    )

    s = S.get(
        key,
        {}
    )


    try:
        sval = float(
            d.get(
                "s_with_n",
                "nan"
            )
        )
    except Exception:
        sval = float("nan")


    try:
        jac = float(
            s.get(
                "cs_jaccard",
                "nan"
            )
        )
    except Exception:
        jac = float("nan")


    if math.isfinite(sval):

        if sval < 0.02:
            ld_flag = "LOW_MISMATCH"

        elif sval < 0.05:
            ld_flag = "MILD_MISMATCH"

        elif sval < 0.10:
            ld_flag = "MODERATE_MISMATCH"

        else:
            ld_flag = "HIGH_MISMATCH"

    else:
        ld_flag = "NO_DIAGNOSTIC"


    if math.isfinite(jac):

        if jac >= 0.8:
            stability = "STABLE"

        elif jac >= 0.5:
            stability = "MODERATELY_STABLE"

        else:
            stability = "SENSITIVE"

    else:
        stability = "NOT_AVAILABLE"


    rows.append({
        "phenotype":
            phenotype,

        "locus":
            locus,

        "n_used":
            d.get(
                "n_used",
                ""
            ),

        "n_mode":
            d.get(
                "n_mode",
                ""
            ),

        "s_with_n":
            d.get(
                "s_with_n",
                ""
            ),

        "s_without_n":
            d.get(
                "s_without_n",
                ""
            ),

        "ld_reference_flag":
            ld_flag,

        "possible_switch_n":
            d.get(
                "possible_switch_n",
                ""
            ),

        "same_top":
            s.get(
                "same_top",
                ""
            ),

        "cs_jaccard":
            s.get(
                "cs_jaccard",
                ""
            ),

        "sample_size_sensitivity":
            stability
    })


fields = list(
    rows[0].keys()
)


with out.open(
    "w",
    newline=""
) as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(
        rows
    )


print(
    "\t".join(fields)
)

for r in rows:

    print(
        "\t".join(
            str(r[x])
            for x in fields
        )
    )
PY

    return $?
}


# ============================================================
# 07. FINAL STATUS
# ============================================================

final_status() {

    BBJ="$OUTROOT/bbj_ld_diagnostics/BBJ_LD_DIAGNOSTIC_SUMMARY.tsv"

    GIGA="$OUTROOT/gigastroke_ld_diagnostics/GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv"

    NAWARE="$OUTROOT/GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv"

    COMP="$OUTROOT/GIGASTROKE_ZONLY_VS_NAWARE_CORRECTED.tsv"

    QC="$OUTROOT/PHASE8B_INTERPRETATION_QC.tsv"

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE8B_STATUS.tsv"


    BBJ_N=0
    GIGA_N=0
    NAWARE_N=0
    COMP_N=0


    if [ -s "$BBJ" ]; then

        BBJ_N="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$BBJ"
        )"
    fi


    if [ -s "$GIGA" ]; then

        GIGA_N="$(
            awk -F '\t' '
                NR>1 && $3=="PASS" {n++}
                END {print n+0}
            ' "$GIGA"
        )"
    fi


    if [ -s "$NAWARE" ]; then

        NAWARE_N="$(
            awk -F '\t' '
                NR>1 && $3=="PASS" {n++}
                END {print n+0}
            ' "$NAWARE"
        )"
    fi


    if [ -s "$COMP" ]; then

        COMP_N="$(
            awk '
                NR>1 {n++}
                END {print n+0}
            ' "$COMP"
        )"
    fi


    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"


    if [ "$BBJ_N" -eq 4 ]; then
        printf "BBJ_LD_DIAGNOSTIC\tREADY\t\n" >> "$OUT"
    else
        printf "BBJ_LD_DIAGNOSTIC\tPARTIAL\t%s/4\n" "$BBJ_N" >> "$OUT"
    fi


    if [ "$GIGA_N" -eq 20 ]; then
        printf "GIGASTROKE_LD_DIAGNOSTIC\tREADY\t\n" >> "$OUT"
    else
        printf "GIGASTROKE_LD_DIAGNOSTIC\tPARTIAL\t%s/20\n" "$GIGA_N" >> "$OUT"
    fi


    if [ "$NAWARE_N" -eq 20 ]; then
        printf "GIGASTROKE_NAWARE_CORRECTED\tREADY\t\n" >> "$OUT"
    else
        printf "GIGASTROKE_NAWARE_CORRECTED\tPARTIAL\t%s/20\n" "$NAWARE_N" >> "$OUT"
    fi


    if [ "$COMP_N" -eq 20 ]; then
        printf "N_SENSITIVITY_CORRECTED\tREADY\t\n" >> "$OUT"
    else
        printf "N_SENSITIVITY_CORRECTED\tPARTIAL\t%s/20\n" "$COMP_N" >> "$OUT"
    fi


    if [ -s "$QC" ]; then
        printf "PHASE8B_INTERPRETATION_QC\tREADY\t\n" >> "$OUT"
    else
        printf "PHASE8B_INTERPRETATION_QC\tBLOCKED\tmissing\n" >> "$OUT"
    fi


    printf "GIGASTROKE_ANCESTRY\tUNVERIFIED\tmetadata ancestry blank\n" >> "$OUT"
    printf "PHASE7_SHARED_SIGNAL_BENCHMARK\tREADY\t\n" >> "$OUT"
    printf "TPMI\tWAITING\trate_limit\n" >> "$OUT"
    printf "CKB\tWAITING\tdecryption_key\n" >> "$OUT"


    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 8B REPAIR"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_sample_size_audit" \
    sample_size_audit

run_step \
    "02_bbj_ld_diagnostics" \
    bbj_ld_diagnostics

run_step \
    "03_gigastroke_ld_diagnostics" \
    gigastroke_ld_diagnostics

run_step \
    "04_correct_naware_susie" \
    gigastroke_naware_susie

run_step \
    "05_compare_sensitivity" \
    compare_sensitivity

run_step \
    "06_build_qc_summary" \
    build_qc_summary

run_step \
    "07_final_status" \
    final_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' \
    "$STATUS"


echo
echo "===== BBJ LD DIAGNOSTIC ====="

column -t -s $'\t' \
    "$OUTROOT/bbj_ld_diagnostics/BBJ_LD_DIAGNOSTIC_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== GIGASTROKE LD DIAGNOSTIC ====="

column -t -s $'\t' \
    "$OUTROOT/gigastroke_ld_diagnostics/GIGASTROKE_LD_DIAGNOSTIC_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== CORRECTED N-AWARE ====="

column -t -s $'\t' \
    "$OUTROOT/GIGASTROKE_SUSIE_NAWARE_CORRECTED_MASTER.tsv" \
    2>/dev/null


echo
echo "===== Z-ONLY VS CORRECTED N-AWARE ====="

column -t -s $'\t' \
    "$OUTROOT/GIGASTROKE_ZONLY_VS_NAWARE_CORRECTED.tsv" \
    2>/dev/null


echo
echo "===== INTERPRETATION QC ====="

column -t -s $'\t' \
    "$OUTROOT/PHASE8B_INTERPRETATION_QC.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE8B_STATUS.tsv" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

