#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

TOOLDIR="$REPO/tools/plink2_portable"
PLINK2="$TOOLDIR/plink2"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase5d/$RUN_ID"
STATUS="$LOGROOT/STATUS.tsv"

mkdir -p "$LOGROOT" "$TOOLDIR"
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

    return 0
}


# ============================================================
# 01. INSTALL CPU-COMPATIBLE MODERN PLINK2
# ============================================================

install_plink2() {

    ZIP="$TOOLDIR/plink2.zip"

    echo "===== CPU FLAGS ====="

    grep -m1 '^flags' /proc/cpuinfo \
        | fold -w 120 \
        || true

    echo
    echo "===== SYSTEM PLINK2 ====="

    /usr/bin/plink2 --version 2>/dev/null || true

    echo
    echo "===== CPU MODE ====="

    if grep -m1 '^flags' /proc/cpuinfo \
        | grep -qw avx2
    then

        MODE="AVX2"

        URL="https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20260929.zip"

    else

        MODE="PORTABLE_X86_64"

        URL="https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260929.zip"

    fi

    echo "CPU_MODE=$MODE"
    echo "URL=$URL"

    # Existing binary can be reused only when it actually runs
    # and contains modern LD functionality.
    if [ -x "$PLINK2" ]; then

        echo
        echo "===== EXISTING PROJECT PLINK2 ====="

        "$PLINK2" --version

        RC=$?

        if [ "$RC" -eq 0 ]; then

            if "$PLINK2" --help r-unphased 2>&1 \
                | grep -q 'ref-based'
            then

                echo "MODERN_PLINK2_ALREADY_READY"
                return 0

            fi

        fi

    fi

    echo
    echo "===== DOWNLOAD ====="

    rm -f \
        "$ZIP" \
        "$PLINK2"

    wget \
        --tries=5 \
        --timeout=60 \
        -O "$ZIP" \
        "$URL"

    RC=$?

    if [ "$RC" -ne 0 ] || [ ! -s "$ZIP" ]; then
        echo "PLINK2_DOWNLOAD_FAILED"
        return 11
    fi

    unzip -o \
        "$ZIP" \
        -d "$TOOLDIR"

    chmod +x "$PLINK2"

    echo
    echo "===== VERSION ====="

    "$PLINK2" --version

    RC=$?

    if [ "$RC" -ne 0 ]; then
        echo "PLINK2_RUNTIME_FAILED"
        return 12
    fi

    echo
    echo "===== FEATURE CHECK ====="

    "$PLINK2" --help r-unphased 2>&1 \
        | grep -E \
          'r-unphased|square|triangle|bin4|bin|ref-based' \
        | head -50

    if ! "$PLINK2" --help r-unphased 2>&1 \
        | grep -q 'ref-based'
    then
        echo "REF_BASED_NOT_SUPPORTED"
        return 13
    fi

    echo
    echo "PLINK2_READY=YES"

    return 0
}


# ============================================================
# 02. CREATE SIGNED EAS LD
# ============================================================

create_ld() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    INPUT="$BASE/ld"
    OUTDIR="$BASE/ld_v3"

    mkdir -p "$OUTDIR"

    SUMMARY="$OUTDIR/LD_BUILD_STATUS.tsv"

    printf \
        "locus\tstatus\tn_variants\tmatrix_bytes\texpected_bytes\n" \
        > "$SUMMARY"

    FAIL=0

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PREFIX="$INPUT/${LOCUS}.matched"
        OUT="$OUTDIR/$LOCUS"

        echo
        echo "===================================================="
        echo "$LOCUS"
        echo "===================================================="

        if [ ! -s "$PREFIX.pgen" ]; then

            echo "MATCHED_PGEN_MISSING"

            printf \
                "%s\tMISSING_INPUT\t0\t0\t0\n" \
                "$LOCUS" \
                >> "$SUMMARY"

            FAIL=$((FAIL + 1))
            continue

        fi

        rm -f \
            "$OUT.unphased.vcor1.bin" \
            "$OUT.unphased.vcor1.bin.vars" \
            "$OUT.log"

        "$PLINK2" \
            --pfile "$PREFIX" \
            --r-unphased \
              square \
              bin \
              ref-based \
            --out "$OUT"

        RC=$?

        if [ "$RC" -ne 0 ]; then

            echo "LD_COMMAND_FAIL rc=$RC"

            printf \
                "%s\tPLINK_FAIL\t0\t0\t0\n" \
                "$LOCUS" \
                >> "$SUMMARY"

            FAIL=$((FAIL + 1))
            continue

        fi

        MATRIX="$OUT.unphased.vcor1.bin"
        VARS="$MATRIX.vars"

        echo
        echo "===== GENERATED FILES ====="

        ls -lh "$OUT"* 2>/dev/null || true

        if [ ! -s "$MATRIX" ]; then

            echo "MATRIX_MISSING"

            printf \
                "%s\tMATRIX_MISSING\t0\t0\t0\n" \
                "$LOCUS" \
                >> "$SUMMARY"

            FAIL=$((FAIL + 1))
            continue

        fi

        if [ ! -s "$VARS" ]; then

            echo "VARS_MISSING"

            printf \
                "%s\tVARS_MISSING\t0\t0\t0\n" \
                "$LOCUS" \
                >> "$SUMMARY"

            FAIL=$((FAIL + 1))
            continue

        fi

        N="$(wc -l < "$VARS")"
        BYTES="$(stat -c%s "$MATRIX")"

        # PLINK 'bin' = double precision
        EXPECTED=$((N * N * 8))

        echo "N_VARIANTS=$N"
        echo "MATRIX_BYTES=$BYTES"
        echo "EXPECTED_BYTES=$EXPECTED"

        if [ "$BYTES" -eq "$EXPECTED" ]; then

            ST="PASS"

        else

            ST="SIZE_MISMATCH"
            FAIL=$((FAIL + 1))

        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" \
            "$ST" \
            "$N" \
            "$BYTES" \
            "$EXPECTED" \
            >> "$SUMMARY"

    done

    echo
    echo "===== LD BUILD SUMMARY ====="

    column -t -s $'\t' "$SUMMARY"

    if [ "$FAIL" -gt 0 ]; then
        return 21
    fi

    return 0
}


# ============================================================
# 03. PREPARE SUSIE INPUTS USING EXACT .vars ORDER
# ============================================================

prepare_susie_inputs() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    LD="$BASE/ld_v3"
    SUM="$BASE/intersection"

    OUTDIR="$BASE/susie_inputs_v3"

    mkdir -p "$OUTDIR"

    python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/"
    "stage3_finemap/japan/bbj"
)

loci = [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]

for locus in loci:

    varfile = (
        root /
        "ld_v3" /
        f"{locus}.unphased.vcor1.bin.vars"
    )

    sumfile = (
        root /
        "intersection" /
        f"{locus}.summary.tsv"
    )

    outfile = (
        root /
        "susie_inputs_v3" /
        f"{locus}.tsv"
    )

    if not varfile.exists():

        print(
            locus,
            "VARS_MISSING"
        )
        continue

    if not sumfile.exists():

        print(
            locus,
            "SUMMARY_MISSING"
        )
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
        x
        for x in order
        if x not in rows
    ]

    if missing:

        print(
            locus,
            "MISSING_SUMMARY_VARIANTS=",
            len(missing)
        )

        print(
            "EXAMPLE=",
            missing[:5]
        )

        continue

    ordered = [
        rows[x]
        for x in order
    ]

    if not ordered:

        print(
            locus,
            "NO_ORDERED_VARIANTS"
        )
        continue

    with outfile.open(
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=list(
                ordered[0].keys()
            ),
            delimiter="\t"
        )

        writer.writeheader()
        writer.writerows(
            ordered
        )

    print(
        locus,
        "READY",
        "N=",
        len(ordered),
        "OUT=",
        outfile
    )
PY

    return $?
}


# ============================================================
# 04. LD NUMERICAL SANITY CHECK
# ============================================================

ld_sanity() {

    Rscript - <<'RS'

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

for (locus in loci) {

    cat(
        "\n==============================\n",
        locus,
        "\n==============================\n",
        sep=""
    )

    varsfile <- file.path(
        root,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    matrixfile <- file.path(
        root,
        "ld_v3",
        paste0(
            locus,
            ".unphased.vcor1.bin"
        )
    )

    if (
        !file.exists(varsfile) ||
        !file.exists(matrixfile)
    ) {

        cat(
            "STATUS=MISSING\n"
        )

        next
    }

    vars <- readLines(
        varsfile
    )

    n <- length(vars)

    expected <- n * n * 8

    actual <- file.info(
        matrixfile
    )$size

    cat(
        "N=",
        n,
        "\n",
        sep=""
    )

    cat(
        "BYTES=",
        actual,
        "\n",
        sep=""
    )

    cat(
        "EXPECTED=",
        expected,
        "\n",
        sep=""
    )

    if (actual != expected) {

        cat(
            "STATUS=SIZE_FAIL\n"
        )

        next
    }

    con <- file(
        matrixfile,
        "rb"
    )

    x <- readBin(
        con,
        what="numeric",
        n=n*n,
        size=8,
        endian="little"
    )

    close(con)

    R <- matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )

    asym <- max(
        abs(
            R - t(R)
        )
    )

    dr <- range(
        diag(R)
    )

    rr <- range(
        R
    )

    nf <- sum(
        !is.finite(R)
    )

    cat(
        "MAX_ASYMMETRY=",
        asym,
        "\n",
        sep=""
    )

    cat(
        "DIAG_MIN=",
        dr[1],
        "\n",
        sep=""
    )

    cat(
        "DIAG_MAX=",
        dr[2],
        "\n",
        sep=""
    )

    cat(
        "R_MIN=",
        rr[1],
        "\n",
        sep=""
    )

    cat(
        "R_MAX=",
        rr[2],
        "\n",
        sep=""
    )

    cat(
        "NONFINITE=",
        nf,
        "\n",
        sep=""
    )

    if (
        nf == 0 &&
        asym < 1e-6 &&
        min(dr) > 0.999 &&
        max(dr) < 1.001
    ) {

        cat(
            "STATUS=PASS\n"
        )

    } else {

        cat(
            "STATUS=CHECK\n"
        )
    }
}
RS

    return $?
}


# ============================================================
# 05. RUN SUSIE RSS
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
        "susieR unavailable\n"
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


read_ld <- function(
    path,
    n
) {

    expected <- n * n * 8

    actual <- file.info(
        path
    )$size

    if (
        is.na(actual) ||
        actual != expected
    ) {

        stop(
            paste0(
                "LD byte mismatch ",
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
        "\n========================================\n",
        locus,
        "\n========================================\n",
        sep=""
    )

    result <- tryCatch({

        ssfile <- file.path(
            root,
            "susie_inputs_v3",
            paste0(
                locus,
                ".tsv"
            )
        )

        varsfile <- file.path(
            root,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin.vars"
            )
        )

        matrixfile <- file.path(
            root,
            "ld_v3",
            paste0(
                locus,
                ".unphased.vcor1.bin"
            )
        )

        if (
            !file.exists(ssfile) ||
            !file.exists(varsfile) ||
            !file.exists(matrixfile)
        ) {

            stop(
                "required input missing"
            )
        }

        ss <- fread(
            ssfile
        )

        vars <- readLines(
            varsfile
        )

        nvar <- length(
            vars
        )

        if (
            nrow(ss) != nvar
        ) {

            stop(
                "summary/LD dimensions differ"
            )
        }

        if (
            !all(
                ss$variant_id == vars
            )
        ) {

            stop(
                "summary/LD order differs"
            )
        }

        R <- read_ld(
            matrixfile,
            nvar
        )

        # numerical guard only
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

        if (
            any(
                !is.finite(z)
            )
        ) {

            stop(
                "nonfinite z"
            )
        }

        cat(
            "nvar=",
            nvar,
            "\n",
            sep=""
        )

        cat(
            "max_abs_z=",
            max(
                abs(z)
            ),
            "\n",
            sep=""
        )

        fit <- susie_rss(
            z=z,
            R=R,
            n=174686,
            L=10,
            estimate_residual_variance=FALSE,
            coverage=0.95,
            min_abs_corr=0.5,
            max_iter=1000,
            verbose=FALSE
        )

        outdir <- file.path(
            root,
            "susie_v3",
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
            p=ss$p,
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

        csout <- data.table()

        if (
            !is.null(cs) &&
            length(cs) > 0
        ) {

            for (
                i in seq_along(cs)
            ) {

                idx <- cs[[i]]

                csout <- rbind(
                    csout,
                    data.table(
                        credible_set=i,
                        variant_index=idx,
                        variant_id=vars[idx],
                        pip=fit$pip[idx]
                    )
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
            "credible_sets=",
            ncs,
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
            n_cs=ncs,
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
            nvar=NA_integer_,
            n_cs=NA_integer_,
            top_variant=NA_character_,
            top_pip=NA_real_,
            error=conditionMessage(e)
        )
    })


    master[[length(master) + 1]] <- result
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
        "SUSIE_MASTER_SUMMARY_V4.tsv"
    ),
    sep="\t"
)


cat(
    "\n===== SUSIE MASTER =====\n"
)

print(
    summary
)

RS

    return $?
}


# ============================================================
# 06. CREATE RESEARCH SUMMARY
# ============================================================

research_summary() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    SUS="$BASE/SUSIE_MASTER_SUMMARY_V4.tsv"

    OUT="$BASE/BBJ_FINEMAP_RESEARCH_SUMMARY.tsv"

    if [ ! -s "$SUS" ]; then
        echo "SUSIE_SUMMARY_MISSING"
        return 61
    fi

    python3 - <<'PY'
import csv
from pathlib import Path

base = Path(
    "/srv/is-analysis/results/is/"
    "stage3_finemap/japan/bbj"
)

inp = (
    base /
    "SUSIE_MASTER_SUMMARY_V4.tsv"
)

out = (
    base /
    "BBJ_FINEMAP_RESEARCH_SUMMARY.tsv"
)

rows = []

with inp.open() as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for r in reader:

        locus = r["locus"]

        pipfile = (
            base /
            "susie_v3" /
            locus /
            "PIP.tsv"
        )

        csfile = (
            base /
            "susie_v3" /
            locus /
            "CREDIBLE_SETS.tsv"
        )

        n_pip_01 = 0
        n_pip_05 = 0

        if pipfile.exists():

            with pipfile.open() as pf:

                pr = csv.DictReader(
                    pf,
                    delimiter="\t"
                )

                for x in pr:

                    try:
                        p = float(
                            x["pip"]
                        )
                    except Exception:
                        continue

                    if p >= 0.1:
                        n_pip_01 += 1

                    if p >= 0.5:
                        n_pip_05 += 1

        cs_sizes = {}

        if csfile.exists() and csfile.stat().st_size > 0:

            with csfile.open() as cf:

                cr = csv.DictReader(
                    cf,
                    delimiter="\t"
                )

                for x in cr:

                    cs = x.get(
                        "credible_set",
                        ""
                    )

                    if cs:
                        cs_sizes[cs] = \
                            cs_sizes.get(
                                cs,
                                0
                            ) + 1

        rows.append({
            "locus":
                locus,

            "status":
                r["status"],

            "converged":
                r["converged"],

            "n_variants":
                r["nvar"],

            "n_credible_sets":
                r["n_cs"],

            "top_variant":
                r["top_variant"],

            "top_pip":
                r["top_pip"],

            "n_variants_pip_ge_0.1":
                n_pip_01,

            "n_variants_pip_ge_0.5":
                n_pip_05,

            "credible_set_sizes":
                ";".join(
                    f"CS{k}:{v}"
                    for k, v
                    in sorted(
                        cs_sizes.items(),
                        key=lambda x:
                        int(x[0])
                    )
                )
        })


fields = list(
    rows[0].keys()
)

with out.open(
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(
        rows
    )


print(
    "\t".join(fields)
)

for r in rows:

    print(
        "\t".join(
            str(
                r[x]
            )
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

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    SUS="$BASE/SUSIE_MASTER_SUMMARY_V4.tsv"
    LD="$BASE/ld_v3/LD_BUILD_STATUS.tsv"

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE5D_STATUS.tsv"

    LDS=0
    SUSN=0

    if [ -s "$LD" ]; then

        LDS="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$LD"
        )"

    fi

    if [ -s "$SUS" ]; then

        SUSN="$(
            awk -F '\t' '
                NR>1 && $2=="PASS" {n++}
                END {print n+0}
            ' "$SUS"
        )"

    fi

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    if [ "$LDS" -eq 4 ]; then

        printf \
            "BBJ_SIGNED_LD\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_SIGNED_LD\tPARTIAL\t%s/4\n" \
            "$LDS" \
            >> "$OUT"

    fi

    if [ "$SUSN" -eq 4 ]; then

        printf \
            "BBJ_SUSIE\tREADY\t\n" \
            >> "$OUT"

    else

        printf \
            "BBJ_SUSIE\tPARTIAL\t%s/4\n" \
            "$SUSN" \
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


    echo "===== READINESS ====="

    column -t -s $'\t' \
        "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 5D"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="


run_step \
    "01_install_plink2" \
    install_plink2


run_step \
    "02_create_ld" \
    create_ld


run_step \
    "03_prepare_susie_inputs" \
    prepare_susie_inputs


run_step \
    "04_ld_sanity" \
    ld_sanity


run_step \
    "05_run_susie" \
    run_susie


run_step \
    "06_research_summary" \
    research_summary


run_step \
    "07_final_status" \
    final_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"


echo
echo "===== LD BUILD ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/LD_BUILD_STATUS.tsv" \
    2>/dev/null


echo
echo "===== SUSIE MASTER ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY_V4.tsv" \
    2>/dev/null


echo
echo "===== RESEARCH SUMMARY ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_FINEMAP_RESEARCH_SUMMARY.tsv" \
    2>/dev/null


echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE5D_STATUS.tsv" \
    2>/dev/null


echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

