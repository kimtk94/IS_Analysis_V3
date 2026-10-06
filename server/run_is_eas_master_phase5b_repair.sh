#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase5b/$RUN_ID"
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

    RC=${PIPESTATUS[0]}

    if [ "$RC" -eq 0 ]; then
        echo "[PASS] $STEP"
        record "$STEP" "PASS" "completed" "$LOG"
    else
        echo "[FAIL] $STEP rc=$RC"
        record "$STEP" "FAIL" "rc=$RC" "$LOG"
    fi

    # master pipeline continues
    return 0
}


# ============================================================
# 01. INPUT AUDIT
# ============================================================

input_audit() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    OUT="$BASE/PHASE5B_INPUT_AUDIT.tsv"

    printf \
        "locus\tmatched_pgen\tvariants\tbbj_summary\tstatus\n" \
        > "$OUT"

    FAIL=0

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PGEN="$BASE/ld/${LOCUS}.matched.pgen"
        PVAR="$BASE/ld/${LOCUS}.matched.pvar"
        SUM="$BASE/intersection/${LOCUS}.summary.tsv"

        PG="NO"
        NV=0
        SM="NO"
        ST="READY"

        [ -s "$PGEN" ] && PG="YES"
        [ -s "$SUM" ] && SM="YES"

        if [ -s "$PVAR" ]; then
            NV="$(
                awk '
                    $0 !~ /^#/ {n++}
                    END {print n+0}
                ' "$PVAR"
            )"
        fi

        if [ "$PG" != "YES" ] || \
           [ "$SM" != "YES" ] || \
           [ "$NV" -lt 10 ]
        then
            ST="BLOCKED"
            FAIL=$((FAIL + 1))
        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" "$PG" "$NV" "$SM" "$ST" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    if [ "$FAIL" -gt 0 ]; then
        return 11
    fi

    return 0
}


# ============================================================
# 02. CREATE CORRECT SIGNED LD MATRICES
#
# CRITICAL:
# --r-unphased = signed genotype correlation
# ref-based     = deterministic REF allele orientation
# bin           = FLOAT64 / 8 bytes
# ============================================================

create_signed_ld() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    LDDIR="$BASE/ld"

    FAIL=0

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PREFIX="$LDDIR/${LOCUS}.matched"
        OUT="$LDDIR/$LOCUS"

        echo
        echo "===================================================="
        echo "$LOCUS"
        echo "===================================================="

        if [ ! -s "$PREFIX.pgen" ]; then
            echo "MATCHED_PGEN_MISSING"
            FAIL=$((FAIL + 1))
            continue
        fi

        rm -f \
            "$OUT.unphased.vcor1.bin" \
            "$OUT.unphased.vcor1.bin.vars"

        plink2 \
            --pfile "$PREFIX" \
            --r-unphased square bin ref-based \
            --out "$OUT"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "LD_FAIL=$LOCUS rc=$RC"
            FAIL=$((FAIL + 1))
            continue
        fi

        MATRIX="$OUT.unphased.vcor1.bin"
        VARS="$OUT.unphased.vcor1.bin.vars"

        if [ ! -s "$MATRIX" ] || [ ! -s "$VARS" ]; then
            echo "EXPECTED_LD_OUTPUT_MISSING"
            FAIL=$((FAIL + 1))

            echo "Produced files:"
            ls -lh "$OUT"* 2>/dev/null
            continue
        fi

        N="$(wc -l < "$VARS")"
        BYTES="$(stat -c%s "$MATRIX")"

        # bin = float64
        EXPECTED=$((N * N * 8))

        echo "N_VARIANTS=$N"
        echo "MATRIX_BYTES=$BYTES"
        echo "EXPECTED_BYTES=$EXPECTED"

        if [ "$BYTES" -eq "$EXPECTED" ]; then
            echo "LD_MATRIX_SIZE=PASS"
        else
            echo "LD_MATRIX_SIZE=FAIL"
            FAIL=$((FAIL + 1))
        fi

    done

    if [ "$FAIL" -gt 0 ]; then
        return 21
    fi

    return 0
}


# ============================================================
# 03. LD MATRIX AUDIT
# ============================================================

ld_audit() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    LDDIR="$BASE/ld"
    OUT="$BASE/BBJ_LD_MATRIX_AUDIT_V2.tsv"

    printf \
        "locus\tn_variants\tmatrix_bytes\texpected_bytes\tvars_match_summary\tstatus\n" \
        > "$OUT"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        MATRIX="$LDDIR/$LOCUS.unphased.vcor1.bin"
        VARS="$LDDIR/$LOCUS.unphased.vcor1.bin.vars"
        SUM="$BASE/intersection/${LOCUS}.summary.tsv"

        N=0
        BYTES=0
        EXPECTED=0
        ORDER="NO"
        ST="MISSING"

        if [ -s "$VARS" ]; then
            N="$(wc -l < "$VARS")"
        fi

        if [ -s "$MATRIX" ]; then
            BYTES="$(stat -c%s "$MATRIX")"
        fi

        EXPECTED=$((N * N * 8))

        if [ -s "$VARS" ] && [ -s "$SUM" ]; then

            python3 - "$VARS" "$SUM" <<'PY'
import csv
import sys

varfile = sys.argv[1]
sumfile = sys.argv[2]

vars_ = [
    x.strip()
    for x in open(varfile)
    if x.strip()
]

summary = {}

with open(sumfile) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        summary[r["variant_id"]] = 1

missing = [
    x for x in vars_
    if x not in summary
]

raise SystemExit(
    0 if not missing else 1
)
PY

            if [ "$?" -eq 0 ]; then
                ORDER="YES"
            fi

        fi

        if [ "$N" -gt 0 ] && \
           [ "$BYTES" -eq "$EXPECTED" ] && \
           [ "$ORDER" = "YES" ]
        then
            ST="PASS"
        elif [ "$BYTES" -gt 0 ]; then
            ST="CHECK"
        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" \
            "$N" \
            "$BYTES" \
            "$EXPECTED" \
            "$ORDER" \
            "$ST" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    PASS_N="$(
        awk -F '\t' '
            NR>1 && $6=="PASS" {n++}
            END {print n+0}
        ' "$OUT"
    )"

    if [ "$PASS_N" -ne 4 ]; then
        return 31
    fi

    return 0
}


# ============================================================
# 04. INSTALL susieR IF NEEDED
# ============================================================

install_susie() {

    echo "===== CURRENT R ====="
    Rscript --version

    echo
    echo "===== PACKAGE CHECK ====="

    HAVE="$(
        Rscript -e '
        cat(
            requireNamespace(
                "susieR",
                quietly=TRUE
            )
        )
        '
    )"

    echo "susieR_before=$HAVE"

    if [ "$HAVE" = "TRUE" ]; then

        Rscript -e '
        cat(
            "susieR_version=",
            as.character(
                packageVersion("susieR")
            ),
            "\n",
            sep=""
        )
        '

        return 0
    fi

    echo
    echo "===== INSTALL susieR ====="

    Rscript - <<'RS'
options(
    repos=c(
        CRAN="https://cloud.r-project.org"
    )
)

install.packages(
    "susieR",
    dependencies=TRUE
)

if (!requireNamespace(
    "susieR",
    quietly=TRUE
)) {
    quit(status=41)
}

cat(
    "susieR_version=",
    as.character(
        packageVersion("susieR")
    ),
    "\n",
    sep=""
)
RS

    return $?
}


# ============================================================
# 05. LD / SUMMARY ORDER PREPARATION
#
# Matrix .vars order is canonical.
# Summary statistics are reordered to that sequence.
# ============================================================

prepare_susie_inputs() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    OUTDIR="$BASE/susie_inputs"

    mkdir -p "$OUTDIR"

    python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/stage3_finemap/"
    "japan/bbj"
)

outdir = root / "susie_inputs"
outdir.mkdir(parents=True, exist_ok=True)

for locus in [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]:

    varfile = (
        root /
        "ld" /
        f"{locus}.unphased.vcor1.bin.vars"
    )

    sumfile = (
        root /
        "intersection" /
        f"{locus}.summary.tsv"
    )

    outfile = (
        outdir /
        f"{locus}.susie_summary.tsv"
    )

    if not varfile.exists():
        print(locus, "VARS_MISSING")
        continue

    if not sumfile.exists():
        print(locus, "SUMMARY_MISSING")
        continue

    order = [
        x.strip()
        for x in varfile.read_text().splitlines()
        if x.strip()
    ]

    with sumfile.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        rows = {
            r["variant_id"]: r
            for r in reader
        }

    missing = [
        x for x in order
        if x not in rows
    ]

    if missing:

        print(
            locus,
            "MISSING_SUMMARY_VARIANTS=",
            len(missing)
        )

        continue

    ordered = [
        rows[x]
        for x in order
    ]

    fields = list(ordered[0])

    with outfile.open(
        "w",
        newline=""
    ) as f:

        w = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=fields
        )

        w.writeheader()
        w.writerows(ordered)

    print(
        locus,
        "N=",
        len(ordered),
        "OUTPUT=",
        outfile
    )
PY

    return $?
}


# ============================================================
# 06. RUN SUSIE RSS
# ============================================================

run_susie() {

    Rscript - <<'RS'
suppressPackageStartupMessages(
    library(data.table)
)

if (!requireNamespace(
    "susieR",
    quietly=TRUE
)) {

    cat(
        "ERROR: susieR not installed\n"
    )

    quit(status=51)
}

suppressPackageStartupMessages(
    library(susieR)
)

root <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

loci <- c(
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004"
)


read_float64_square <- function(
    file,
    n
) {

    expected <- n * n * 8

    actual <- file.info(file)$size

    if (
        is.na(actual) ||
        actual != expected
    ) {

        stop(
            paste0(
                "matrix byte mismatch: expected=",
                expected,
                " actual=",
                actual
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
            "incorrect matrix element count"
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
        "\n========================================\n",
        locus,
        "\n========================================\n",
        sep=""
    )

    res <- tryCatch({

        sumfile <- file.path(
            root,
            "susie_inputs",
            paste0(
                locus,
                ".susie_summary.tsv"
            )
        )

        varsfile <- file.path(
            root,
            "ld",
            paste0(
                locus,
                ".unphased.vcor1.bin.vars"
            )
        )

        matrixfile <- file.path(
            root,
            "ld",
            paste0(
                locus,
                ".unphased.vcor1.bin"
            )
        )

        if (
            !file.exists(sumfile) ||
            !file.exists(varsfile) ||
            !file.exists(matrixfile)
        ) {

            stop(
                "required input missing"
            )
        }

        ss <- fread(sumfile)

        vars <- fread(
            varsfile,
            header=FALSE
        )[[1]]

        nvar <- length(vars)

        if (
            nrow(ss) != nvar
        ) {

            stop(
                "summary/matrix dimension mismatch"
            )
        }

        if (
            !all(
                ss$variant_id == vars
            )
        ) {

            stop(
                "summary and LD variant order mismatch"
            )
        }

        R <- read_float64_square(
            matrixfile,
            nvar
        )

        # numerical symmetry guard
        R <- (
            R +
            t(R)
        ) / 2

        diag(R) <- 1

        if (
            any(!is.finite(R))
        ) {

            stop(
                "nonfinite LD matrix"
            )
        }

        z <- (
            as.numeric(ss$beta) /
            as.numeric(ss$se)
        )

        if (
            any(!is.finite(z))
        ) {

            stop(
                "nonfinite z score"
            )
        }

        cat(
            "variants=",
            nvar,
            "\n",
            sep=""
        )

        cat(
            "max_abs_z=",
            max(abs(z)),
            "\n",
            sep=""
        )

        # BBJ total N
        N <- 174686

        fit <- susie_rss(
            z=z,
            R=R,
            n=N,
            L=10,
            estimate_residual_variance=FALSE,
            max_iter=1000,
            coverage=0.95,
            min_abs_corr=0.5,
            verbose=FALSE
        )

        outdir <- file.path(
            root,
            "susie",
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
            variant_id=vars,
            beta=ss$beta,
            se=ss$se,
            z=z,
            p=ss$p,
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

        csout <- data.table()

        cs <- fit$sets$cs

        if (
            !is.null(cs) &&
            length(cs) > 0
        ) {

            for (
                i in seq_along(cs)
            ) {

                idx <- cs[[i]]

                tmp <- data.table(
                    credible_set=i,
                    variant_index=idx,
                    variant_id=vars[idx],
                    pip=fit$pip[idx]
                )

                csout <- rbind(
                    csout,
                    tmp
                )
            }
        }

        fwrite(
            csout,
            file.path(
                outdir,
                "CREDIBLE_SETS.tsv"
            ),
            sep="\t"
        )

        top_idx <- which.max(
            fit$pip
        )

        cat(
            "converged=",
            fit$converged,
            "\n",
            sep=""
        )

        cat(
            "credible_sets=",
            ifelse(
                is.null(cs),
                0,
                length(cs)
            ),
            "\n",
            sep=""
        )

        cat(
            "top_variant=",
            vars[top_idx],
            "\n",
            sep=""
        )

        cat(
            "top_pip=",
            fit$pip[top_idx],
            "\n",
            sep=""
        )

        list(
            locus=locus,
            status="PASS",
            converged=fit$converged,
            nvar=nvar,
            ncs=ifelse(
                is.null(cs),
                0,
                length(cs)
            ),
            top_variant=vars[top_idx],
            top_pip=fit$pip[top_idx],
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
            converged=FALSE,
            nvar=NA,
            ncs=NA,
            top_variant=NA,
            top_pip=NA,
            error=conditionMessage(e)
        )
    })

    master[[length(master)+1]] <- res
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
        root,
        "SUSIE_MASTER_SUMMARY_V2.tsv"
    ),
    sep="\t"
)

cat(
    "\n===== MASTER =====\n"
)

print(summary)
RS

    return $?
}


# ============================================================
# 07. OUTPUT SANITY CHECK
# ============================================================

susie_audit() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    SUM="$BASE/SUSIE_MASTER_SUMMARY_V2.tsv"

    if [ ! -s "$SUM" ]; then
        echo "SUSIE_MASTER_SUMMARY_MISSING"
        return 71
    fi

    echo "===== SUSIE MASTER ====="

    column -t -s $'\t' "$SUM"

    echo
    echo "===== TOP PIP PER LOCUS ====="

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PIP="$BASE/susie/$LOCUS/PIP.tsv"

        echo
        echo "--- $LOCUS ---"

        if [ -s "$PIP" ]; then
            column -t -s $'\t' "$PIP" \
                | head -11
        else
            echo "NO_PIP_FILE"
        fi

    done

    return 0
}


# ============================================================
# 08. FINAL READINESS
# ============================================================

final_status() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE5B_STATUS.tsv"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    LD="$BASE/BBJ_LD_MATRIX_AUDIT_V2.tsv"

    NLPASS=0

    if [ -s "$LD" ]; then

        NLPASS="$(
            awk -F '\t' '
                NR>1 && $6=="PASS" {n++}
                END {print n+0}
            ' "$LD"
        )"

    fi

    if [ "$NLPASS" -eq 4 ]; then
        printf "BBJ_SIGNED_LD\tREADY\t\n" >> "$OUT"
    else
        printf \
            "BBJ_SIGNED_LD\tPARTIAL\t%s/4\n" \
            "$NLPASS" \
            >> "$OUT"
    fi


    SUS="$BASE/SUSIE_MASTER_SUMMARY_V2.tsv"

    SPASS=0

    if [ -s "$SUS" ]; then

        SPASS="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$SUS"
        )"

    fi

    if [ "$SPASS" -eq 4 ]; then
        printf "BBJ_SUSIE\tREADY\t\n" >> "$OUT"
    else
        printf \
            "BBJ_SUSIE\tPARTIAL\t%s/4\n" \
            "$SPASS" \
            >> "$OUT"
    fi

    printf \
        "GIGASTROKE_EAS\tCANONICAL_READY\t\n" \
        >> "$OUT"

    printf \
        "TPMI\tWAITING\trate_limit\n" \
        >> "$OUT"

    printf \
        "CKB\tWAITING\tdecryption_key\n" \
        >> "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 5B REPAIR"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_input_audit" \
    input_audit

run_step \
    "02_create_signed_ld" \
    create_signed_ld

run_step \
    "03_ld_audit" \
    ld_audit

run_step \
    "04_install_susie" \
    install_susie

run_step \
    "05_prepare_susie_inputs" \
    prepare_susie_inputs

run_step \
    "06_run_susie" \
    run_susie

run_step \
    "07_susie_audit" \
    susie_audit

run_step \
    "08_final_status" \
    final_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== LD AUDIT ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_LD_MATRIX_AUDIT_V2.tsv" \
    2>/dev/null

echo
echo "===== SUSIE ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY_V2.tsv" \
    2>/dev/null

echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE5B_STATUS.tsv" \
    2>/dev/null

echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "STATUS=$STATUS"
echo "===================================================="
