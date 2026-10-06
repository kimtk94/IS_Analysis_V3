#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase7b/$RUN_ID"
STATUS="$LOGROOT/STATUS.tsv"

mkdir -p "$LOGROOT"
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
        record "$STEP" "PASS" "completed" "$LOG"
    else
        echo "[FAIL] $STEP rc=$RC"
        record "$STEP" "FAIL" "rc=$RC" "$LOG"
    fi

    return 0
}


# ============================================================
# 01. CORRECT INPUT AUDIT
# ============================================================

input_audit() {

    BASE="$ROOT/results/is"

    LDROOT="$BASE/stage3_finemap/japan/bbj/ld_v3"
    HROOT="$BASE/stage4_cross_eas/harmonized"

    FAIL=0

    echo "===== LD ====="

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        MATRIX="$LDROOT/${LOCUS}.unphased.vcor1.bin"
        VARS="$LDROOT/${LOCUS}.unphased.vcor1.bin.vars"

        if [ -s "$MATRIX" ] && [ -s "$VARS" ]; then

            N="$(wc -l < "$VARS")"
            BYTES="$(stat -c%s "$MATRIX")"
            EXPECTED=$((N * N * 8))

            echo \
                "$LOCUS N=$N bytes=$BYTES expected=$EXPECTED"

            if [ "$BYTES" -ne "$EXPECTED" ]; then
                echo "SIZE_FAIL $LOCUS"
                FAIL=$((FAIL + 1))
            fi

        else

            echo "MISSING $LOCUS"
            FAIL=$((FAIL + 1))

        fi

    done


    echo
    echo "===== HARMONIZED INPUTS ====="

    for P in AS AIS CES LAS SVS
    do

        for LOCUS in \
            BBJ_IS_L001 \
            BBJ_IS_L002 \
            BBJ_IS_L003 \
            BBJ_IS_L004
        do

            F="$HROOT/${P}.${LOCUS}.harmonized.tsv"

            if [ -s "$F" ]; then
                echo \
                    "OK $P $LOCUS rows=$(awk 'NR>1{n++} END{print n+0}' "$F")"
            else
                echo "MISSING $P $LOCUS"
                FAIL=$((FAIL + 1))
            fi

        done

    done


    echo
    echo "===== SUSIER ====="

    Rscript - <<'RS'
cat(
    "susieR_available=",
    requireNamespace("susieR", quietly=TRUE),
    "\n",
    sep=""
)

if (requireNamespace("susieR", quietly=TRUE)) {
    cat(
        "susieR_version=",
        as.character(packageVersion("susieR")),
        "\n",
        sep=""
    )
}
RS

    if [ "$FAIL" -gt 0 ]; then
        return 11
    fi

    return 0
}


# ============================================================
# 02. REMOVE INVALID PHASE7 DERIVED RESULTS
# ============================================================

clean_invalid_phase7() {

    BASE="$ROOT/results/is/stage4_cross_eas/fine_mapping"

    echo "Removing only failed/empty Phase7 derived outputs."

    rm -rf \
        "$BASE/susie_zonly"

    rm -f \
        "$BASE/GIGASTROKE_SUSIE_ZONLY_MASTER.tsv" \
        "$BASE/BBJ_GIGASTROKE_CS_COMPARISON.tsv" \
        "$BASE/CS_JACCARD_MATRIX.tsv" \
        "$BASE/GIGASTROKE_TOP_PIP_MATRIX.tsv" \
        "$BASE/BBJ_CS_MAX_GIGA_PIP_MATRIX.tsv"

    mkdir -p \
        "$BASE/susie_zonly"

    return 0
}


# ============================================================
# 03. GIGASTROKE Z-ONLY SUSIE
# ============================================================

run_gigastroke_susie() {

Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

if (!requireNamespace(
    "susieR",
    quietly=TRUE
)) {
    stop("susieR unavailable")
}

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

OUT_ROOT <- file.path(
    ROOT,
    "stage4_cross_eas/fine_mapping/susie_zonly"
)

dir.create(
    OUT_ROOT,
    recursive=TRUE,
    showWarnings=FALSE
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

    expected <- n*n*8
    actual <- file.info(path)$size

    if (
        is.na(actual) ||
        actual != expected
    ) {
        stop(
            paste(
                "LD byte mismatch",
                actual,
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

    if (length(x) != n*n) {
        stop(
            "LD element count mismatch"
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

    cat(
        "\n########################################\n",
        locus,
        "\n########################################\n",
        sep=""
    )


    vars_file <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    ld_file <- file.path(
        BBJ_ROOT,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin"
        )
    )


    if (
        !file.exists(vars_file) ||
        !file.exists(ld_file)
    ) {

        for (phenotype in phenotypes) {

            master[[length(master) + 1]] <- list(
                phenotype=phenotype,
                locus=locus,
                analysis="SUSIE_RSS_Z_ONLY",
                status="FAIL",
                converged=FALSE,
                n_variants=NA_integer_,
                n_cs=NA_integer_,
                top_variant=NA_character_,
                top_pip=NA_real_,
                max_abs_z=NA_real_,
                error="locus LD missing"
            )

        }

        next
    }


    full_vars <- readLines(
        vars_file
    )

    n_full <- length(
        full_vars
    )

    R_full <- read_ld(
        ld_file,
        n_full
    )

    R_full <- (
        R_full +
        t(R_full)
    ) / 2

    diag(R_full) <- 1


    index_map <- setNames(
        seq_along(full_vars),
        full_vars
    )


    for (phenotype in phenotypes) {

        cat(
            "\n===== ",
            phenotype,
            " / ",
            locus,
            " =====\n",
            sep=""
        )


        result <- tryCatch({

            infile <- file.path(
                HARM_ROOT,
                paste0(
                    phenotype,
                    ".",
                    locus,
                    ".harmonized.tsv"
                )
            )

            if (!file.exists(infile)) {
                stop(
                    "harmonized input missing"
                )
            }


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


            if (nrow(d) < 10) {
                stop(
                    "too few overlapping variants"
                )
            }


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
                is.finite(beta) &
                is.finite(p) &
                p > 0 &
                p <= 1
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


            if (length(beta) < 10) {
                stop(
                    "too few finite variants"
                )
            }


            # Signed z from two-sided p-value.
            p_safe <- pmax(
                p,
                .Machine$double.xmin
            )

            z <- sign(
                beta
            ) * qnorm(
                p_safe / 2,
                lower.tail=FALSE
            )


            if (
                any(!is.finite(z))
            ) {
                stop(
                    "nonfinite z"
                )
            }


            if (
                any(!is.finite(R))
            ) {
                stop(
                    "nonfinite LD"
                )
            }


            cat(
                "n_variants=",
                length(z),
                "\n",
                sep=""
            )

            cat(
                "max_abs_z=",
                max(abs(z)),
                "\n",
                sep=""
            )


            fit <- susie_rss(
                z=z,
                R=R,
                L=10,
                estimate_residual_variance=FALSE,
                coverage=0.95,
                min_abs_corr=0.5,
                max_iter=1000,
                verbose=FALSE
            )


            outdir <- file.path(
                OUT_ROOT,
                phenotype,
                locus
            )

            dir.create(
                outdir,
                recursive=TRUE,
                showWarnings=FALSE
            )


            saveRDS(
                fit,
                file.path(
                    outdir,
                    "susie_fit.rds"
                )
            )


            pip <- data.table(
                phenotype=phenotype,
                locus=locus,
                variant_id=d$variant_id,
                beta_aligned=beta,
                p=p,
                z=z,
                pip=fit$pip
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

            cs_table <- data.table()


            if (
                !is.null(cs) &&
                length(cs) > 0
            ) {

                for (
                    i in seq_along(cs)
                ) {

                    ci <- cs[[i]]

                    cs_table <- rbind(
                        cs_table,
                        data.table(
                            phenotype=phenotype,
                            locus=locus,
                            credible_set=i,
                            variant_id=d$variant_id[ci],
                            pip=fit$pip[ci]
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
                "converged=",
                fit$converged,
                "\n",
                sep=""
            )

            cat(
                "n_cs=",
                ncs,
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


            list(
                phenotype=phenotype,
                locus=locus,
                analysis="SUSIE_RSS_Z_ONLY",
                status="PASS",
                converged=fit$converged,
                n_variants=length(z),
                n_cs=ncs,
                top_variant=d$variant_id[top],
                top_pip=fit$pip[top],
                max_abs_z=max(abs(z)),
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
                phenotype=phenotype,
                locus=locus,
                analysis="SUSIE_RSS_Z_ONLY",
                status="FAIL",
                converged=FALSE,
                n_variants=NA_integer_,
                n_cs=NA_integer_,
                top_variant=NA_character_,
                top_pip=NA_real_,
                max_abs_z=NA_real_,
                error=conditionMessage(e)
            )
        })


        master[[length(master) + 1]] <- result

    }
}


summary <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)


outfile <- file.path(
    dirname(OUT_ROOT),
    "GIGASTROKE_SUSIE_ZONLY_MASTER.tsv"
)

fwrite(
    summary,
    outfile,
    sep="\t"
)


cat(
    "\n===== GIGASTROKE SUSIE MASTER =====\n"
)

print(
    summary
)


cat(
    "\nPASS_COUNT=",
    sum(
        summary$status == "PASS"
    ),
    "/",
    nrow(summary),
    "\n",
    sep=""
)

RS

    return $?
}


# ============================================================
# 04. VERIFY 20/20 BEFORE COMPARISON
# ============================================================

verify_susie() {

    FILE="$ROOT/results/is/stage4_cross_eas/fine_mapping/GIGASTROKE_SUSIE_ZONLY_MASTER.tsv"

    if [ ! -s "$FILE" ]; then
        echo "MASTER_MISSING"
        return 31
    fi

    TOTAL="$(
        awk '
            NR>1 {n++}
            END {print n+0}
        ' "$FILE"
    )"

    PASS="$(
        awk -F '\t' '
            NR>1 && $4=="PASS" {
                n++
            }
            END {
                print n+0
            }
        ' "$FILE"
    )"

    echo "TOTAL=$TOTAL"
    echo "PASS=$PASS"

    if [ "$TOTAL" -ne 20 ] || [ "$PASS" -ne 20 ]; then
        echo "SUSIE_NOT_COMPLETE"
        return 32
    fi

    return 0
}


# ============================================================
# 05. CREDIBLE SET COMPARISON
# ============================================================

compare_cs() {

python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is"
)

bbj_root = (
    root /
    "stage3_finemap/japan/bbj"
)

giga_root = (
    root /
    "stage4_cross_eas/"
    "fine_mapping/susie_zonly"
)

phase6_file = (
    root /
    "stage4_cross_eas/"
    "GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"
)

out = (
    root /
    "stage4_cross_eas/fine_mapping/"
    "BBJ_GIGASTROKE_CS_COMPARISON.tsv"
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


phase6 = {}

with phase6_file.open() as f:

    for r in csv.DictReader(
        f,
        delimiter="\t"
    ):

        phase6[
            (
                r["phenotype"],
                r["locus"]
            )
        ] = r


rows = []


for locus in loci:

    bpip = read_pip(
        bbj_root /
        "susie_v3" /
        locus /
        "PIP.tsv"
    )

    bcs = read_cs(
        bbj_root /
        "susie_v3" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    btop = max(
        bpip,
        key=bpip.get
    )


    for phenotype in phenotypes:

        gpip = read_pip(
            giga_root /
            phenotype /
            locus /
            "PIP.tsv"
        )

        gcs = read_cs(
            giga_root /
            phenotype /
            locus /
            "CREDIBLE_SETS.tsv"
        )

        gtop = (
            max(
                gpip,
                key=gpip.get
            )
            if gpip
            else ""
        )


        inter = (
            bcs &
            gcs
        )

        union = (
            bcs |
            gcs
        )

        jac = (
            len(inter) /
            len(union)
            if union
            else ""
        )


        phase = phase6[
            (
                phenotype,
                locus
            )
        ]


        rows.append({
            "phenotype":
                phenotype,

            "locus":
                locus,

            "bbj_top_variant":
                btop,

            "bbj_top_pip":
                bpip[btop],

            "giga_top_variant":
                gtop,

            "giga_top_pip":
                gpip.get(
                    gtop,
                    ""
                ),

            "same_top_variant":
                (
                    "YES"
                    if btop == gtop
                    else "NO"
                ),

            "bbj_cs_size":
                len(bcs),

            "giga_cs_size":
                len(gcs),

            "cs_intersection_n":
                len(inter),

            "cs_union_n":
                len(union),

            "cs_jaccard":
                jac,

            "bbj_top_in_giga_cs":
                (
                    "YES"
                    if btop in gcs
                    else "NO"
                ),

            "giga_top_in_bbj_cs":
                (
                    "YES"
                    if gtop in bcs
                    else "NO"
                ),

            "bbj_cs_max_giga_pip":
                (
                    max(
                        [
                            gpip[x]
                            for x in bcs
                            if x in gpip
                        ],
                        default=""
                    )
                ),

            "giga_cs_max_bbj_pip":
                (
                    max(
                        [
                            bpip[x]
                            for x in gcs
                            if x in bpip
                        ],
                        default=""
                    )
                ),

            "bbj_top_giga_p":
                phase[
                    "bbj_top_giga_p"
                ],

            "bbj_top_direction_concordant":
                phase[
                    "bbj_top_direction_concordant"
                ],

            "bbj_cs_best_giga_p":
                phase[
                    "bbj_cs_best_giga_p"
                ],

            "bbj_cs_direction_concordance":
                phase[
                    "bbj_cs_direction_concordance"
                ]
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
    "N=",
    len(rows)
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
# 06. BUILD MATRICES
# ============================================================

build_matrices() {

python3 - <<'PY'
import csv
from pathlib import Path

inp = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/fine_mapping/"
    "BBJ_GIGASTROKE_CS_COMPARISON.tsv"
)

outdir = inp.parent

with inp.open() as f:

    rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
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

index = {
    (
        r["locus"],
        r["phenotype"]
    ):
        r
    for r in rows
}


metrics = [
    (
        "cs_jaccard",
        "CS_JACCARD_MATRIX.tsv"
    ),
    (
        "giga_top_pip",
        "GIGASTROKE_TOP_PIP_MATRIX.tsv"
    ),
    (
        "bbj_cs_max_giga_pip",
        "BBJ_CS_MAX_GIGA_PIP_MATRIX.tsv"
    )
]


for metric, filename in metrics:

    out = outdir / filename

    with out.open(
        "w",
        newline=""
    ) as f:

        fields = [
            "locus"
        ] + phenotypes

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=fields
        )

        w.writeheader()

        for locus in loci:

            row = {
                "locus":
                    locus
            }

            for phenotype in phenotypes:

                row[
                    phenotype
                ] = index[
                    (
                        locus,
                        phenotype
                    )
                ][
                    metric
                ]

            w.writerow(
                row
            )

    print()
    print(
        "=====",
        filename,
        "====="
    )

    print(
        out.read_text()
    )
PY

    return $?
}


# ============================================================
# 07. FINAL STATUS
# ============================================================

final_status() {

    BASE="$ROOT/results/is/stage4_cross_eas/fine_mapping"

    SUS="$BASE/GIGASTROKE_SUSIE_ZONLY_MASTER.tsv"
    COMP="$BASE/BBJ_GIGASTROKE_CS_COMPARISON.tsv"

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE7B_STATUS.tsv"

    SPASS=0
    COMP_N=0

    if [ -s "$SUS" ]; then

        SPASS="$(
            awk -F '\t' '
                NR>1 && $4=="PASS" {
                    n++
                }
                END {
                    print n+0
                }
            ' "$SUS"
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


    if [ "$SPASS" -eq 20 ]; then

        printf \
            "GIGASTROKE_ZONLY_SUSIE\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "GIGASTROKE_ZONLY_SUSIE\tPARTIAL\t%s/20\n" \
            "$SPASS" \
            >> "$OUT"

    fi


    if [ "$COMP_N" -eq 20 ] && \
       [ "$SPASS" -eq 20 ]
    then

        printf \
            "BBJ_GIGASTROKE_CS_COMPARE\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_GIGASTROKE_CS_COMPARE\tBLOCKED\tSuSiE incomplete\n" \
            >> "$OUT"

    fi


    printf "BBJ_PRIMARY_FINEMAP\tREADY\t\n" >> "$OUT"
    printf "BBJ_N_SENSITIVITY\tREADY\t\n" >> "$OUT"
    printf "PHASE6_20_REGION_BENCHMARK\tREADY\t\n" >> "$OUT"
    printf "TPMI\tWAITING\trate_limit\n" >> "$OUT"
    printf "CKB\tWAITING\tdecryption_key\n" >> "$OUT"


    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 7B REPAIR"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_input_audit" \
    input_audit

run_step \
    "02_clean_invalid_phase7" \
    clean_invalid_phase7

run_step \
    "03_run_gigastroke_susie" \
    run_gigastroke_susie

run_step \
    "04_verify_susie" \
    verify_susie

run_step \
    "05_compare_cs" \
    compare_cs

run_step \
    "06_build_matrices" \
    build_matrices

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
echo "===== GIGASTROKE SUSIE ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/GIGASTROKE_SUSIE_ZONLY_MASTER.tsv" \
    2>/dev/null


echo
echo "===== CS COMPARISON ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv" \
    2>/dev/null


echo
echo "===== CS JACCARD ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/CS_JACCARD_MATRIX.tsv" \
    2>/dev/null


echo
echo "===== TOP PIP ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/GIGASTROKE_TOP_PIP_MATRIX.tsv" \
    2>/dev/null


echo
echo "===== N_TOTAL VS N_EFF ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_NTOTAL_VS_NEFF_COMPARISON.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE7B_STATUS.tsv" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="
