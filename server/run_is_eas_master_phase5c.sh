#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"
TOOLDIR="$REPO/tools/plink2_latest"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase5c/$RUN_ID"
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
# 01. INSTALL PROJECT-LOCAL CURRENT PLINK2
# ============================================================

install_plink2() {

    ZIP="$TOOLDIR/plink2.zip"
    BIN="$TOOLDIR/plink2"

    echo "===== SYSTEM PLINK2 ====="
    /usr/bin/plink2 --version 2>/dev/null || true

    if [ -x "$BIN" ]; then

        echo
        echo "===== EXISTING PROJECT PLINK2 ====="
        "$BIN" --version

        HELP="$("$BIN" --help r-unphased 2>&1)"

        if printf "%s" "$HELP" | grep -q "ref-based"; then
            echo "PROJECT_PLINK2_MODERN=YES"
            return 0
        fi

    fi

    echo
    echo "===== DOWNLOAD CURRENT PLINK2 ====="

    rm -f "$ZIP"

    # Official PLINK2 Linux AVX2 Intel current/development binary.
    # If AVX2 download fails, fallback to generic 64-bit.
    URL1="https://s3.amazonaws.com/plink2-assets/plink2_linux_avx2_20260928.zip"
    URL2="https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20260928.zip"

    wget \
        --tries=3 \
        --timeout=60 \
        -O "$ZIP" \
        "$URL1"

    RC=$?

    if [ "$RC" -ne 0 ] || [ ! -s "$ZIP" ]; then

        echo "AVX2_DOWNLOAD_FAILED; TRY_GENERIC"

        rm -f "$ZIP"

        wget \
            --tries=3 \
            --timeout=60 \
            -O "$ZIP" \
            "$URL2"
    fi

    if [ ! -s "$ZIP" ]; then
        echo "PLINK2_DOWNLOAD_FAILED"
        return 11
    fi

    unzip -o "$ZIP" -d "$TOOLDIR"

    chmod +x "$BIN"

    echo
    echo "===== PROJECT PLINK2 ====="

    "$BIN" --version

    echo
    echo "===== FEATURE CHECK ====="

    "$BIN" --help r-unphased \
        | grep -E \
          'r-unphased|ref-based|square|bin4|bin' \
        | head -30

    if ! "$BIN" --help r-unphased 2>&1 \
        | grep -q "ref-based"
    then
        echo "MODERN_LD_FEATURES_MISSING"
        return 12
    fi

    return 0
}


# ============================================================
# 02. CREATE SIGNED LD WITH MODERN PLINK2
# ============================================================

create_ld() {

    PLINK2="$TOOLDIR/plink2"

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    LDDIR="$BASE/ld_v2"

    mkdir -p "$LDDIR"

    FAIL=0

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        IN="$BASE/ld/${LOCUS}.matched"
        OUT="$LDDIR/$LOCUS"

        echo
        echo "===================================================="
        echo "$LOCUS"
        echo "===================================================="

        if [ ! -s "$IN.pgen" ]; then
            echo "MATCHED_PGEN_MISSING"
            FAIL=$((FAIL+1))
            continue
        fi

        rm -f "$OUT".*

        "$PLINK2" \
            --pfile "$IN" \
            --r-unphased \
              square \
              bin \
              ref-based \
            --out "$OUT"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "PLINK2_LD_FAIL=$RC"
            FAIL=$((FAIL+1))
            continue
        fi

        echo
        echo "===== PRODUCED ====="

        ls -lh "$OUT"* || true

        MATRIX="$OUT.unphased.vcor1.bin"
        VARS="$MATRIX.vars"

        if [ ! -s "$MATRIX" ]; then
            echo "MATRIX_MISSING"
            FAIL=$((FAIL+1))
            continue
        fi

        if [ ! -s "$VARS" ]; then
            echo "VARS_MISSING"
            FAIL=$((FAIL+1))
            continue
        fi

        N="$(wc -l < "$VARS")"
        BYTES="$(stat -c%s "$MATRIX")"
        EXPECTED=$((N * N * 8))

        echo "N=$N"
        echo "BYTES=$BYTES"
        echo "EXPECTED=$EXPECTED"

        if [ "$BYTES" -ne "$EXPECTED" ]; then
            echo "MATRIX_SIZE_MISMATCH"
            FAIL=$((FAIL+1))
            continue
        fi

        echo "LD_MATRIX_OK"

    done

    [ "$FAIL" -eq 0 ]
}


# ============================================================
# 03. PREPARE SUSIE SUMMARY IN EXACT LD ORDER
# ============================================================

prepare_inputs() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"

    mkdir -p "$BASE/susie_inputs_v2"

    python3 - <<'PY'
import csv
from pathlib import Path

root = Path(
    "/srv/is-analysis/results/is/stage3_finemap/"
    "japan/bbj"
)

for locus in [
    "BBJ_IS_L001",
    "BBJ_IS_L002",
    "BBJ_IS_L003",
    "BBJ_IS_L004",
]:

    vars_file = (
        root /
        "ld_v2" /
        f"{locus}.unphased.vcor1.bin.vars"
    )

    summary_file = (
        root /
        "intersection" /
        f"{locus}.summary.tsv"
    )

    out = (
        root /
        "susie_inputs_v2" /
        f"{locus}.tsv"
    )

    if not vars_file.exists():
        print(locus, "VARS_MISSING")
        continue

    if not summary_file.exists():
        print(locus, "SUMMARY_MISSING")
        continue

    order = [
        x.strip()
        for x in vars_file.read_text().splitlines()
        if x.strip()
    ]

    with summary_file.open() as f:

        rows = {
            r["variant_id"]: r
            for r in csv.DictReader(
                f,
                delimiter="\t"
            )
        }

    missing = [
        x for x in order
        if x not in rows
    ]

    if missing:
        print(
            locus,
            "SUMMARY_MISSING_FOR_LD=",
            len(missing)
        )
        continue

    ordered = [
        rows[x]
        for x in order
    ]

    with out.open("w", newline="") as f:

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=list(ordered[0].keys())
        )

        writer.writeheader()
        writer.writerows(ordered)

    print(
        locus,
        "READY",
        len(ordered)
    )
PY

    return $?
}


# ============================================================
# 04. LD SANITY CHECK IN R
# ============================================================

ld_sanity() {

    Rscript - <<'RS'
root <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

loci <- sprintf(
    "BBJ_IS_L%03d",
    1:4
)

for (locus in loci) {

    cat(
        "\n===== ",
        locus,
        " =====\n",
        sep=""
    )

    varsfile <- file.path(
        root,
        "ld_v2",
        paste0(
            locus,
            ".unphased.vcor1.bin.vars"
        )
    )

    binfile <- file.path(
        root,
        "ld_v2",
        paste0(
            locus,
            ".unphased.vcor1.bin"
        )
    )

    if (
        !file.exists(varsfile) ||
        !file.exists(binfile)
    ) {
        cat("MISSING\n")
        next
    }

    vars <- readLines(varsfile)
    n <- length(vars)

    con <- file(binfile, "rb")

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
        byrow=TRUE
    )

    cat(
        "N=",
        n,
        "\n",
        sep=""
    )

    cat(
        "MAX_ASYMMETRY=",
        max(abs(R-t(R))),
        "\n",
        sep=""
    )

    cat(
        "DIAG_RANGE=",
        paste(
            range(diag(R)),
            collapse=","
        ),
        "\n",
        sep=""
    )

    cat(
        "R_RANGE=",
        paste(
            range(R),
            collapse=","
        ),
        "\n",
        sep=""
    )

    if (any(!is.finite(R))) {
        cat("NONFINITE=YES\n")
    } else {
        cat("NONFINITE=NO\n")
    }
}
RS

    return $?
}


# ============================================================
# 05. RUN SUSIE-RSS
# ============================================================

run_susie() {

    Rscript - <<'RS'
suppressPackageStartupMessages(
    library(data.table)
)

suppressPackageStartupMessages(
    library(susieR)
)

root <- paste0(
    "/srv/is-analysis/results/is/",
    "stage3_finemap/japan/bbj"
)

loci <- sprintf(
    "BBJ_IS_L%03d",
    1:4
)

master <- list()


read_matrix <- function(
    file,
    n
) {

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
            "LD matrix size mismatch"
        )
    }

    matrix(
        x,
        nrow=n,
        ncol=n,
        byrow=TRUE
    )
}


for (locus in loci) {

    cat(
        "\n====================================\n",
        locus,
        "\n====================================\n",
        sep=""
    )

    ans <- tryCatch({

        ssfile <- file.path(
            root,
            "susie_inputs_v2",
            paste0(
                locus,
                ".tsv"
            )
        )

        varsfile <- file.path(
            root,
            "ld_v2",
            paste0(
                locus,
                ".unphased.vcor1.bin.vars"
            )
        )

        matrixfile <- file.path(
            root,
            "ld_v2",
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

        ss <- fread(ssfile)

        vars <- readLines(
            varsfile
        )

        nvar <- length(vars)

        if (
            nrow(ss) != nvar
        ) {
            stop(
                "summary/LD dimension mismatch"
            )
        }

        if (
            !all(
                ss$variant_id == vars
            )
        ) {
            stop(
                "variant order mismatch"
            )
        }

        R <- read_matrix(
            matrixfile,
            nvar
        )

        # guard tiny numerical asymmetry
        R <- (R+t(R))/2
        diag(R) <- 1

        z <- (
            as.numeric(ss$beta) /
            as.numeric(ss$se)
        )

        if (
            any(!is.finite(z))
        ) {
            stop(
                "invalid z values"
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
            max(abs(z)),
            "\n",
            sep=""
        )

        fit <- susie_rss(
            z=z,
            R=R,
            n=174686,
            L=10,
            estimate_residual_variance=FALSE,
            max_iter=1000,
            coverage=0.95,
            min_abs_corr=0.5,
            verbose=FALSE
        )

        outdir <- file.path(
            root,
            "susie_v2",
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

        cs_table <- data.table()

        if (
            !is.null(cs) &&
            length(cs) > 0
        ) {

            for (
                i in seq_along(cs)
            ) {

                idx <- cs[[i]]

                cs_table <- rbind(
                    cs_table,
                    data.table(
                        credible_set=i,
                        variant_id=vars[idx],
                        pip=fit$pip[idx]
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

        cat(
            "converged=",
            fit$converged,
            "\n",
            sep=""
        )

        cat(
            "n_cs=",
            ifelse(
                is.null(cs),
                0,
                length(cs)
            ),
            "\n",
            sep=""
        )

        cat(
            "top=",
            vars[top],
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
            locus=locus,
            status="PASS",
            converged=fit$converged,
            nvar=nvar,
            n_cs=ifelse(
                is.null(cs),
                0,
                length(cs)
            ),
            top_variant=vars[top],
            top_pip=fit$pip[top],
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
            n_cs=NA,
            top_variant=NA,
            top_pip=NA,
            error=conditionMessage(e)
        )
    })

    master[
        [length(master)+1]
    ] <- ans
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
        "SUSIE_MASTER_SUMMARY_V3.tsv"
    ),
    sep="\t"
)

cat(
    "\n===== FINAL SUSIE =====\n"
)

print(summary)
RS

    return $?
}


# ============================================================
# 06. FINAL AUDIT
# ============================================================

final_audit() {

    BASE="$ROOT/results/is/stage3_finemap/japan/bbj"
    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE5C_STATUS.tsv"

    echo "===== SUSIE SUMMARY ====="

    column -t -s $'\t' \
        "$BASE/SUSIE_MASTER_SUMMARY_V3.tsv" \
        2>/dev/null

    PASS="$(
        awk -F '\t' '
            NR>1 && $2=="PASS" {n++}
            END {print n+0}
        ' "$BASE/SUSIE_MASTER_SUMMARY_V3.tsv" \
        2>/dev/null
    )"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    if [ "$PASS" -eq 4 ]; then
        printf "BBJ_SUSIE\tREADY\t\n" >> "$OUT"
    else
        printf \
            "BBJ_SUSIE\tPARTIAL\t%s/4\n" \
            "$PASS" \
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

    echo
    column -t -s $'\t' "$OUT"

    echo
    echo "===== TOP PIPS ====="

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        FILE="$BASE/susie_v2/$LOCUS/PIP.tsv"

        echo
        echo "--- $LOCUS ---"

        if [ -s "$FILE" ]; then
            column -t -s $'\t' "$FILE" \
                | head -11
        else
            echo "NO_RESULT"
        fi

    done

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 5C"
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
    "03_prepare_inputs" \
    prepare_inputs

run_step \
    "04_ld_sanity" \
    ld_sanity

run_step \
    "05_run_susie" \
    run_susie

run_step \
    "06_final_audit" \
    final_audit


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== SUSIE MASTER ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY_V3.tsv" \
    2>/dev/null

echo
echo "===== READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE5C_STATUS.tsv" \
    2>/dev/null

echo
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

