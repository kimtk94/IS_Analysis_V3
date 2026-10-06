#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase7/$RUN_ID"
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
# 01. PHASE 6 + INPUT AUDIT
# ============================================================

input_audit() {

    BASE="$ROOT/results/is"

    BENCH="$BASE/stage4_cross_eas/GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"

    echo "===== PHASE6 BENCHMARK ====="

    if [ ! -s "$BENCH" ]; then
        echo "BENCHMARK_MISSING"
        return 11
    fi

    N="$(
        awk '
            NR > 1 { n++ }
            END { print n+0 }
        ' "$BENCH"
    )"

    U="$(
        awk -F '\t' '
            NR > 1 {
                x[$1 FS $2]=1
            }
            END {
                print length(x)
            }
        ' "$BENCH"
    )"

    echo "ROWS=$N"
    echo "UNIQUE_COMPARISONS=$U"

    if [ "$N" -ne 20 ] || [ "$U" -ne 20 ]; then
        echo "EXPECTED_20_COMPARISONS"
        return 12
    fi

    echo
    echo "===== HARMONIZED INPUTS ====="

    HROOT="$BASE/stage4_cross_eas/harmonized"

    FAIL=0

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
                NROW="$(
                    awk '
                        NR>1 {n++}
                        END {print n+0}
                    ' "$F"
                )"

                echo "OK $P $LOCUS rows=$NROW"
            else
                echo "MISSING $P $LOCUS"
                FAIL=$((FAIL + 1))
            fi

        done
    done

    echo
    echo "===== EAS LD ====="

    LDROOT="$BASE/stage3_finemap/japan/bbj/ld_v3"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        LD="$LDROOT/${LOCUS}.unphased.vcor1.bin"
        VARS="$LD.bin.vars"

        if [ -s "$LD" ] && [ -s "$VARS" ]; then
            echo \
                "OK $LOCUS variants=$(wc -l < "$VARS")"
        else
            echo "LD_MISSING $LOCUS"
            FAIL=$((FAIL + 1))
        fi
    done

    echo
    echo "===== susieR ====="

    Rscript - <<'RS'
cat(
    "susieR=",
    requireNamespace(
        "susieR",
        quietly=TRUE
    ),
    "\n",
    sep=""
)

if (
    requireNamespace(
        "susieR",
        quietly=TRUE
    )
) {
    cat(
        "version=",
        as.character(
            packageVersion("susieR")
        ),
        "\n",
        sep=""
    )
}
RS

    if [ "$FAIL" -gt 0 ]; then
        return 13
    fi

    return 0
}


# ============================================================
# 02. HARMONIZATION QC
# ============================================================

harmonization_qc() {

    HROOT="$ROOT/results/is/stage4_cross_eas/harmonized"
    OUTDIR="$ROOT/results/is/stage4_cross_eas/fine_mapping"

    mkdir -p "$OUTDIR"

    OUT="$OUTDIR/HARMONIZATION_QC.tsv"

    python3 - <<'PY'
import csv
from collections import Counter
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/harmonized"
)

out = Path(
    "/srv/is-analysis/results/is/"
    "stage4_cross_eas/fine_mapping/"
    "HARMONIZATION_QC.tsv"
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

valid = {
    "MATCH",
    "SWAP",
    "COMPLEMENT",
    "COMPLEMENT_SWAP"
}

rows = []

for phenotype in phenotypes:

    for locus in loci:

        f = (
            root /
            f"{phenotype}.{locus}.harmonized.tsv"
        )

        if not f.exists():

            rows.append({
                "phenotype": phenotype,
                "locus": locus,
                "n_total": 0,
                "n_valid": 0,
                "valid_rate": 0,
                "n_match": 0,
                "n_swap": 0,
                "n_complement": 0,
                "n_complement_swap": 0,
                "n_ambiguous_palindrome": 0,
                "n_mismatch": 0
            })

            continue

        with f.open() as fh:

            data = list(
                csv.DictReader(
                    fh,
                    delimiter="\t"
                )
            )

        c = Counter(
            x["harmonization"]
            for x in data
        )

        n_valid = sum(
            c[x]
            for x in valid
        )

        rows.append({
            "phenotype":
                phenotype,

            "locus":
                locus,

            "n_total":
                len(data),

            "n_valid":
                n_valid,

            "valid_rate":
                (
                    n_valid / len(data)
                    if data else 0
                ),

            "n_match":
                c["MATCH"],

            "n_swap":
                c["SWAP"],

            "n_complement":
                c["COMPLEMENT"],

            "n_complement_swap":
                c["COMPLEMENT_SWAP"],

            "n_ambiguous_palindrome":
                c["AMBIGUOUS_PALINDROME"],

            "n_mismatch":
                c["ALLELE_MISMATCH"]
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
    w.writerows(rows)


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
# 03. GIGASTROKE Z-ONLY SUSIE-RSS
#
# IMPORTANT:
# - same 1000G EAS LD as BBJ
# - no unverified GIGASTROKE N inserted
# - z reconstructed from aligned beta sign + two-sided p
# ============================================================

run_gigastroke_susie() {

    Rscript - <<'RS'

suppressPackageStartupMessages(
    library(data.table)
)

if (
    !requireNamespace(
        "susieR",
        quietly=TRUE
    )
) {
    cat(
        "ERROR: susieR missing\n"
    )

    quit(status=21)
}

suppressPackageStartupMessages(
    library(susieR)
)


root <- paste0(
    "/srv/is-analysis/results/is/"
)

bbj_root <- file.path(
    root,
    "stage3_finemap/japan/bbj"
)

harm_root <- file.path(
    root,
    "stage4_cross_eas/harmonized"
)

out_root <- file.path(
    root,
    "stage4_cross_eas/fine_mapping/susie_zonly"
)

dir.create(
    out_root,
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


read_ld <- function(
    file,
    n
) {

    expected <- n*n*8

    actual <- file.info(
        file
    )$size

    if (
        is.na(actual) ||
        actual != expected
    ) {

        stop(
            paste0(
                "LD byte mismatch: ",
                actual,
                " vs ",
                expected
            )
        )
    }

    con <- file(
        file,
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
        bbj_root,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    ld_file <- file.path(
        bbj_root,
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

        cat(
            "LOCUS LD MISSING\n"
        )

        next
    }


    full_vars <- readLines(
        vars_file
    )

    n_full <- length(
        full_vars
    )

    full_R <- read_ld(
        ld_file,
        n_full
    )

    full_R <- (
        full_R +
        t(full_R)
    ) / 2

    diag(
        full_R
    ) <- 1


    index_map <- setNames(
        seq_along(
            full_vars
        ),
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

            input <- file.path(
                harm_root,
                paste0(
                    phenotype,
                    ".",
                    locus,
                    ".harmonized.tsv"
                )
            )


            if (
                !file.exists(input)
            ) {
                stop(
                    "harmonized input missing"
                )
            }


            d <- fread(
                input
            )


            d <- d[
                harmonization %in%
                valid_harm
            ]


            d <- d[
                variant_id %in%
                full_vars
            ]


            if (
                nrow(d) < 10
            ) {
                stop(
                    "too few overlapping variants"
                )
            }


            # one row per variant
            d <- d[
                !duplicated(
                    variant_id
                )
            ]


            # Exact LD ordering
            idx <- index_map[
                d$variant_id
            ]


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


            R <- full_R[
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
                p >= 0 &
                p <= 1 &
                beta != 0
            )


            if (
                sum(ok) < 10
            ) {
                stop(
                    "too few finite beta/p values"
                )
            }


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


            # Prevent qnorm(0) = Inf.
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
                any(
                    !is.finite(z)
                )
            ) {
                stop(
                    "nonfinite z reconstructed"
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


            # Z-only RSS benchmark.
            # No unverified total/effective N is inserted.
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
                out_root,
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
                phenotype=
                    phenotype,

                locus=
                    locus,

                variant_id=
                    d$variant_id,

                beta_aligned=
                    beta,

                p=
                    p,

                z=
                    z,

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

            cs_table <- data.table()


            if (
                !is.null(cs) &&
                length(cs) > 0
            ) {

                for (
                    i in seq_along(cs)
                ) {

                    cs_idx <- cs[[i]]

                    cs_table <- rbind(
                        cs_table,
                        data.table(
                            phenotype=
                                phenotype,

                            locus=
                                locus,

                            credible_set=
                                i,

                            variant_id=
                                d$variant_id[
                                    cs_idx
                                ],

                            pip=
                                fit$pip[
                                    cs_idx
                                ]
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
                0
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
                phenotype=
                    phenotype,

                locus=
                    locus,

                analysis=
                    "SUSIE_RSS_Z_ONLY",

                status=
                    "PASS",

                converged=
                    fit$converged,

                n_variants=
                    length(z),

                n_cs=
                    ncs,

                top_variant=
                    d$variant_id[top],

                top_pip=
                    fit$pip[top],

                max_abs_z=
                    max(abs(z)),

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
                    phenotype,

                locus=
                    locus,

                analysis=
                    "SUSIE_RSS_Z_ONLY",

                status=
                    "FAIL",

                converged=
                    FALSE,

                n_variants=
                    NA_integer_,

                n_cs=
                    NA_integer_,

                top_variant=
                    NA_character_,

                top_pip=
                    NA_real_,

                max_abs_z=
                    NA_real_,

                error=
                    conditionMessage(e)
            )
        })


        master[
            [length(master) + 1]
        ] <- result
    }
}


summary <- rbindlist(
    lapply(
        master,
        as.data.table
    ),
    fill=TRUE
)


fwrite(
    summary,
    file.path(
        dirname(out_root),
        "GIGASTROKE_SUSIE_ZONLY_MASTER.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== GIGASTROKE SUSIE MASTER =====\n"
)

print(
    summary
)

RS

    return $?
}


# ============================================================
# 04. BBJ ↔ GIGASTROKE CREDIBLE SET COMPARISON
# ============================================================

compare_credible_sets() {

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

bench_file = (
    root /
    "stage4_cross_eas/"
    "GIGASTROKE_20_REGION_BENCHMARK_SUMMARY.tsv"
)

outdir = (
    root /
    "stage4_cross_eas/fine_mapping"
)

outdir.mkdir(
    parents=True,
    exist_ok=True
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

    out = {}

    if not path.exists():
        return out

    with path.open() as f:

        for r in csv.DictReader(
            f,
            delimiter="\t"
        ):

            try:
                out[
                    r["variant_id"]
                ] = float(
                    r["pip"]
                )
            except Exception:
                pass

    return out


def read_cs(path):

    out = set()

    if (
        not path.exists()
        or path.stat().st_size == 0
    ):
        return out

    with path.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for r in reader:

            v = r.get(
                "variant_id",
                ""
            )

            if v:
                out.add(v)

    return out


# Phase6 directional / p benchmark
phase6 = {}

if bench_file.exists():

    with bench_file.open() as f:

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

    bbj_pip = read_pip(
        bbj_root /
        "susie_v3" /
        locus /
        "PIP.tsv"
    )

    bbj_cs = read_cs(
        bbj_root /
        "susie_v3" /
        locus /
        "CREDIBLE_SETS.tsv"
    )

    bbj_top = (
        max(
            bbj_pip,
            key=bbj_pip.get
        )
        if bbj_pip
        else ""
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


        intersect = (
            bbj_cs &
            gcs
        )

        union = (
            bbj_cs |
            gcs
        )

        jaccard = (
            len(intersect) /
            len(union)
            if union
            else None
        )


        bbj_cs_giga_pip = [
            gpip[x]
            for x in bbj_cs
            if x in gpip
        ]


        giga_cs_bbj_pip = [
            bbj_pip[x]
            for x in gcs
            if x in bbj_pip
        ]


        ph6 = phase6.get(
            (
                phenotype,
                locus
            ),
            {}
        )


        rows.append({
            "phenotype":
                phenotype,

            "locus":
                locus,

            "analysis_role":
                "EAS_CONSORTIUM_BENCHMARK_ZONLY",

            "bbj_top_variant":
                bbj_top,

            "bbj_top_pip":
                bbj_pip.get(
                    bbj_top,
                    ""
                ),

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
                    if (
                        bbj_top
                        and gtop
                        and bbj_top == gtop
                    )
                    else "NO"
                ),

            "bbj_cs_size":
                len(
                    bbj_cs
                ),

            "giga_cs_size":
                len(
                    gcs
                ),

            "cs_intersection_n":
                len(
                    intersect
                ),

            "cs_union_n":
                len(
                    union
                ),

            "cs_jaccard":
                (
                    jaccard
                    if jaccard is not None
                    else ""
                ),

            "bbj_cs_max_giga_pip":
                (
                    max(
                        bbj_cs_giga_pip
                    )
                    if bbj_cs_giga_pip
                    else ""
                ),

            "giga_cs_max_bbj_pip":
                (
                    max(
                        giga_cs_bbj_pip
                    )
                    if giga_cs_bbj_pip
                    else ""
                ),

            "bbj_top_giga_p":
                ph6.get(
                    "bbj_top_giga_p",
                    ""
                ),

            "bbj_top_direction_concordant":
                ph6.get(
                    "bbj_top_direction_concordant",
                    ""
                ),

            "bbj_cs_best_giga_p":
                ph6.get(
                    "bbj_cs_best_giga_p",
                    ""
                ),

            "bbj_cs_direction_concordance":
                ph6.get(
                    "bbj_cs_direction_concordance",
                    ""
                )
        })


out = (
    outdir /
    "BBJ_GIGASTROKE_CS_COMPARISON.tsv"
)

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
    "N_COMPARISONS=",
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
# 05. CREATE MATRICES FOR QUICK REVIEW
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
        x["locus"],
        x["phenotype"]
    ):
        x
    for x in rows
}


for metric, filename in [
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
]:

    out = (
        outdir /
        filename
    )

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

            r = {
                "locus":
                    locus
            }

            for phenotype in phenotypes:

                r[
                    phenotype
                ] = index.get(
                    (
                        locus,
                        phenotype
                    ),
                    {}
                ).get(
                    metric,
                    ""
                )

            w.writerow(
                r
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
# 06. PHASE7 READINESS
# ============================================================

final_status() {

    BASE="$ROOT/results/is/stage4_cross_eas/fine_mapping"

    SUS="$BASE/GIGASTROKE_SUSIE_ZONLY_MASTER.tsv"
    COMP="$BASE/BBJ_GIGASTROKE_CS_COMPARISON.tsv"

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE7_STATUS.tsv"

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
                NR>1 {
                    n++
                }
                END {
                    print n+0
                }
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


    if [ "$COMP_N" -eq 20 ]; then

        printf \
            "BBJ_GIGASTROKE_CS_COMPARE\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_GIGASTROKE_CS_COMPARE\tPARTIAL\t%s/20\n" \
            "$COMP_N" \
            >> "$OUT"

    fi


    printf \
        "BBJ_PRIMARY_FINEMAP\tREADY\t\n" \
        >> "$OUT"

    printf \
        "BBJ_N_SENSITIVITY\tREADY\t\n" \
        >> "$OUT"

    printf \
        "PHASE6_20_REGION_BENCHMARK\tREADY\t\n" \
        >> "$OUT"

    printf \
        "TPMI\tWAITING\trate_limit\n" \
        >> "$OUT"

    printf \
        "CKB\tWAITING\tdecryption_key\n" \
        >> "$OUT"


    column -t -s $'\t' \
        "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 7"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="


run_step \
    "01_input_audit" \
    input_audit


run_step \
    "02_harmonization_qc" \
    harmonization_qc


run_step \
    "03_gigastroke_susie" \
    run_gigastroke_susie


run_step \
    "04_compare_credible_sets" \
    compare_credible_sets


run_step \
    "05_build_matrices" \
    build_matrices


run_step \
    "06_final_status" \
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
echo "===== BBJ vs GIGASTROKE CS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/BBJ_GIGASTROKE_CS_COMPARISON.tsv" \
    2>/dev/null


echo
echo "===== CS JACCARD MATRIX ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/CS_JACCARD_MATRIX.tsv" \
    2>/dev/null


echo
echo "===== GIGASTROKE TOP PIP MATRIX ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage4_cross_eas/fine_mapping/GIGASTROKE_TOP_PIP_MATRIX.tsv" \
    2>/dev/null


echo
echo "===== N_TOTAL vs N_EFF ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_NTOTAL_VS_NEFF_COMPARISON.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE7_STATUS.tsv" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="
