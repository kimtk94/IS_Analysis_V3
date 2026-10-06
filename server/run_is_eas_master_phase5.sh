#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase5/$RUN_ID"
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

    return 0
}


# ============================================================
# 01. REBUILD 1KG EAS PGEN FROM GT EXPLICITLY
# ============================================================

rebuild_gt_pgen() {

    VCFDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/regions"
    OUTDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/pgen_gt"

    mkdir -p "$OUTDIR"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        VCF="$VCFDIR/${LOCUS}.1KG_EAS.GRCh37.vcf.gz"
        PREFIX="$OUTDIR/${LOCUS}.1KG_EAS.GRCh37"

        echo
        echo "===== $LOCUS ====="

        if [ ! -s "$VCF" ]; then
            echo "VCF_MISSING"
            continue
        fi

        rm -f \
            "$PREFIX.pgen" \
            "$PREFIX.pvar" \
            "$PREFIX.psam"

        plink2 \
            --vcf "$VCF" \
            --double-id \
            --max-alleles 2 \
            --set-all-var-ids '@:#:$r:$a' \
            --make-pgen \
            --out "$PREFIX"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "PGEN_BUILD_FAIL=$LOCUS"
            continue
        fi

        echo "PGEN_OK"

        ls -lh \
            "$PREFIX.pgen" \
            "$PREFIX.pvar" \
            "$PREFIX.psam"

        echo
        echo "SAMPLES:"
        awk 'NR>1{n++} END{print n+0}' "$PREFIX.psam"

        echo
        echo "VARIANTS:"
        awk 'NR>1{n++} END{print n+0}' "$PREFIX.pvar"

    done

    return 0
}


# ============================================================
# 02. GENOTYPE QC
# ============================================================

genotype_qc() {

    PDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/pgen_gt"
    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj/ld_qc"

    mkdir -p "$OUTDIR"

    SUMMARY="$OUTDIR/GENOTYPE_QC_SUMMARY.tsv"

    printf \
        "locus\tn_samples\tn_variants_raw\tn_variants_maf01\tgenotyping_rate\n" \
        > "$SUMMARY"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PREFIX="$PDIR/${LOCUS}.1KG_EAS.GRCh37"

        echo
        echo "===== $LOCUS ====="

        if [ ! -s "$PREFIX.pgen" ]; then
            echo "PGEN_MISSING"
            continue
        fi

        plink2 \
            --pfile "$PREFIX" \
            --freq \
            --missing \
            --out "$OUTDIR/$LOCUS"

        RAW="$(awk 'NR>1{n++} END{print n+0}' "$PREFIX.pvar")"
        NS="$(awk 'NR>1{n++} END{print n+0}' "$PREFIX.psam")"

        plink2 \
            --pfile "$PREFIX" \
            --maf 0.01 \
            --geno 0.05 \
            --make-pgen \
            --out "$OUTDIR/${LOCUS}.QC"

        QC="$(awk 'NR>1{n++} END{print n+0}' \
            "$OUTDIR/${LOCUS}.QC.pvar" 2>/dev/null)"

        if [ -s "$OUTDIR/${LOCUS}.smiss" ]; then

            RATE="$(awk '
                NR>1 {
                    total++;
                    sum += (1-$5)
                }
                END {
                    if(total>0) print sum/total;
                    else print "NA"
                }
            ' "$OUTDIR/${LOCUS}.smiss")"

        else
            RATE="NA"
        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" \
            "$NS" \
            "$RAW" \
            "$QC" \
            "$RATE" \
            >> "$SUMMARY"

    done

    column -t -s $'\t' "$SUMMARY"

    return 0
}


# ============================================================
# 03. BUILD EXACT BBJ/1KG INTERSECTED DATA
# ============================================================

build_intersection() {

    OVERDIR="$ROOT/results/is/stage3_finemap/japan/bbj/overlap"
    QCDIR="$ROOT/results/is/stage3_finemap/japan/bbj/ld_qc"

    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj/intersection"

    mkdir -p "$OUTDIR"

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

    overlap = (
        root /
        "overlap" /
        f"{locus}.OVERLAP.tsv"
    )

    pvar = (
        root /
        "ld_qc" /
        f"{locus}.QC.pvar"
    )

    outdir = root / "intersection"
    outdir.mkdir(parents=True, exist_ok=True)

    if not overlap.exists():
        print(locus, "OVERLAP_MISSING")
        continue

    if not pvar.exists():
        print(locus, "PVAR_MISSING")
        continue

    kg = set()

    with pvar.open() as f:

        header = f.readline()

        for line in f:

            c = line.rstrip("\n").split("\t")

            if len(c) < 5:
                continue

            chrom = c[0]
            pos = c[1]
            ref = c[3]
            alt = c[4]

            kg.add(
                f"{chrom}:{pos}:{ref}:{alt}"
            )

    kept = []

    with overlap.open() as f:

        reader = csv.DictReader(
            f,
            delimiter="\t"
        )

        for r in reader:

            if r["status"] != "MATCH":
                continue

            try:
                eaf = float(r["eaf"])
                info = float(r["info"])
            except Exception:
                continue

            maf = min(eaf, 1-eaf)

            # Primary fine-mapping QC
            if maf < 0.01:
                continue

            if info < 0.8:
                continue

            key = (
                f'{r["chr"]}:'
                f'{r["pos"]}:'
                f'{r["bbj_ref"]}:'
                f'{r["bbj_alt"]}'
            )

            if key not in kg:
                continue

            x = dict(r)
            x["variant_id"] = key
            x["maf"] = maf

            kept.append(x)

    kept.sort(
        key=lambda x: int(x["pos"])
    )

    extract = (
        outdir /
        f"{locus}.variants.txt"
    )

    extract.write_text(
        "\n".join(
            x["variant_id"]
            for x in kept
        ) + "\n"
    )

    table = (
        outdir /
        f"{locus}.summary.tsv"
    )

    if kept:

        fields = list(kept[0])

        with table.open(
            "w",
            newline=""
        ) as fo:

            w = csv.DictWriter(
                fo,
                fieldnames=fields,
                delimiter="\t"
            )

            w.writeheader()
            w.writerows(kept)

    print(
        locus,
        "N_INTERSECT=",
        len(kept)
    )
PY

    return $?
}


# ============================================================
# 04. CREATE MATCHED PGEN + LD MATRIX
# ============================================================

create_ld() {

    QCDIR="$ROOT/results/is/stage3_finemap/japan/bbj/ld_qc"
    INTDIR="$ROOT/results/is/stage3_finemap/japan/bbj/intersection"
    LDDIR="$ROOT/results/is/stage3_finemap/japan/bbj/ld"

    mkdir -p "$LDDIR"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PREFIX="$QCDIR/${LOCUS}.QC"
        VARS="$INTDIR/${LOCUS}.variants.txt"
        OUT="$LDDIR/$LOCUS"

        echo
        echo "===== $LOCUS ====="

        if [ ! -s "$PREFIX.pgen" ]; then
            echo "QC_PGEN_MISSING"
            continue
        fi

        if [ ! -s "$VARS" ]; then
            echo "VARIANT_LIST_MISSING"
            continue
        fi

        N="$(wc -l < "$VARS")"

        echo "N_VARIANTS=$N"

        if [ "$N" -lt 10 ]; then
            echo "TOO_FEW_VARIANTS"
            continue
        fi

        # matched pgen
        plink2 \
            --pfile "$PREFIX" \
            --extract "$VARS" \
            --make-pgen \
            --out "$OUT.matched"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "MATCHED_PGEN_FAIL"
            continue
        fi

        # LD matrix, binary double
        plink2 \
            --pfile "$OUT.matched" \
            --r square bin \
            --out "$OUT"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "LD_FAIL"
            continue
        fi

        echo "LD_OK"

        ls -lh "$OUT"* | head -20

    done

    return 0
}


# ============================================================
# 05. LD OUTPUT AUDIT
# ============================================================

ld_audit() {

    LDDIR="$ROOT/results/is/stage3_finemap/japan/bbj/ld"
    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj"
    OUT="$OUTDIR/BBJ_LD_MATRIX_AUDIT.tsv"

    printf \
        "locus\tn_variants\tmatrix_file\tmatrix_bytes\texpected_float_bytes\tstatus\n" \
        > "$OUT"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        PVAR="$LDDIR/${LOCUS}.matched.pvar"

        MATRIX=""

        for CAND in \
            "$LDDIR/${LOCUS}.vcor1.bin" \
            "$LDDIR/${LOCUS}.vcor.bin" \
            "$LDDIR/${LOCUS}.unphased.vcor1.bin"
        do
            if [ -s "$CAND" ]; then
                MATRIX="$CAND"
                break
            fi
        done

        if [ -s "$PVAR" ]; then
            N="$(awk 'NR>1{n++} END{print n+0}' "$PVAR")"
        else
            N=0
        fi

        if [ -n "$MATRIX" ]; then
            BYTES="$(stat -c%s "$MATRIX")"
        else
            BYTES=0
        fi

        # PLINK2 binary correlations are float32
        EXPECTED=$(( N * N * 4 ))

        if [ "$N" -gt 0 ] && [ "$BYTES" -eq "$EXPECTED" ]; then
            ST="PASS"
        elif [ "$BYTES" -gt 0 ]; then
            ST="CHECK_FORMAT"
        else
            ST="MISSING"
        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" \
            "$N" \
            "$MATRIX" \
            "$BYTES" \
            "$EXPECTED" \
            "$ST" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 06. R / SUSIER CHECK
# ============================================================

susie_check() {

    echo "===== R ====="

    Rscript --version

    echo
    echo "===== susieR ====="

    Rscript - <<'RS'
ok <- requireNamespace(
    "susieR",
    quietly = TRUE
)

cat(
    "susieR_available=",
    ok,
    "\n",
    sep=""
)

if (ok) {
    cat(
        "susieR_version=",
        as.character(
            packageVersion("susieR")
        ),
        "\n",
        sep=""
    )
}
RS

    return 0
}


# ============================================================
# 07. RUN SUSIE-RSS
# ============================================================

run_susie() {

    Rscript - <<'RS'
library(data.table)

if (!requireNamespace(
    "susieR",
    quietly = TRUE
)) {

    cat("susieR missing; skipping\n")
    quit(status=0)
}

library(susieR)

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

read_matrix <- function(path, n) {

    con <- file(
        path,
        "rb"
    )

    on.exit(
        close(con)
    )

    x <- readBin(
        con,
        what=numeric(),
        n=n*n,
        size=4,
        endian="little"
    )

    if (length(x) != n*n) {
        stop(
            "matrix size mismatch"
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
        "\n==============================\n",
        locus,
        "\n==============================\n",
        sep=""
    )

    tryCatch({

        sumfile <- file.path(
            root,
            "intersection",
            paste0(
                locus,
                ".summary.tsv"
            )
        )

        pvar <- file.path(
            root,
            "ld",
            paste0(
                locus,
                ".matched.pvar"
            )
        )

        candidates <- c(
            file.path(
                root,
                "ld",
                paste0(
                    locus,
                    ".vcor1.bin"
                )
            ),
            file.path(
                root,
                "ld",
                paste0(
                    locus,
                    ".vcor.bin"
                )
            ),
            file.path(
                root,
                "ld",
                paste0(
                    locus,
                    ".unphased.vcor1.bin"
                )
            )
        )

        matrixfile <- candidates[
            file.exists(candidates)
        ][1]

        if (
            !file.exists(sumfile) ||
            !file.exists(pvar) ||
            is.na(matrixfile)
        ) {

            stop(
                "required input missing"
            )
        }

        ss <- fread(sumfile)

        pv <- fread(
            pvar,
            skip=0
        )

        # PLINK pvar names
        if ("#CHROM" %in% names(pv)) {
            setnames(
                pv,
                "#CHROM",
                "CHROM"
            )
        }

        pv[, variant_id :=
            paste(
                CHROM,
                POS,
                REF,
                ALT,
                sep=":"
            )
        ]

        ss <- ss[
            match(
                pv$variant_id,
                ss$variant_id
            )
        ]

        if (
            anyNA(ss$variant_id)
        ) {
            stop(
                "summary/pvar order mismatch"
            )
        }

        nvar <- nrow(pv)

        cat(
            "nvar=",
            nvar,
            "\n",
            sep=""
        )

        R <- read_matrix(
            matrixfile,
            nvar
        )

        R <- (
            R +
            t(R)
        ) / 2

        diag(R) <- 1

        z <- (
            as.numeric(ss$beta) /
            as.numeric(ss$se)
        )

        fit <- susie_rss(
            z=z,
            R=R,
            n=174686,
            L=min(
                10,
                max(
                    1,
                    floor(nvar/100)
                )
            ),
            estimate_residual_variance=FALSE,
            max_iter=1000,
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
            variant_id=pv$variant_id,
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

        if (!is.null(cs)) {

            for (
                i in seq_along(cs)
            ) {

                idx <- cs[[i]]

                x <- data.table(
                    credible_set=i,
                    variant_index=idx,
                    variant_id=
                        pv$variant_id[idx],
                    pip=
                        fit$pip[idx]
                )

                csout <- rbind(
                    csout,
                    x
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

        cat(
            "converged=",
            fit$converged,
            "\n",
            sep=""
        )

        cat(
            "credible_sets=",
            length(cs),
            "\n",
            sep=""
        )

        cat(
            "top_pip=",
            max(fit$pip),
            "\n",
            sep=""
        )

    }, error=function(e) {

        cat(
            "ERROR=",
            conditionMessage(e),
            "\n",
            sep=""
        )

    })
}
RS

    return 0
}


# ============================================================
# 08. SUSIE SUMMARY
# ============================================================

susie_summary() {

    ROOTSUS="$ROOT/results/is/stage3_finemap/japan/bbj/susie"
    OUT="$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY.tsv"

    printf \
        "locus\tfit_exists\tn_credible_sets\ttop_variant\ttop_pip\n" \
        > "$OUT"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        DIR="$ROOTSUS/$LOCUS"
        PIP="$DIR/PIP.tsv"
        CS="$DIR/CREDIBLE_SETS.tsv"

        FIT="NO"
        NCS=0
        TOP="NA"
        PIPV="NA"

        if [ -s "$DIR/susie_fit.rds" ]; then
            FIT="YES"
        fi

        if [ -s "$CS" ]; then
            NCS="$(awk -F '\t' '
                NR>1 {x[$1]=1}
                END {print length(x)}
            ' "$CS")"
        fi

        if [ -s "$PIP" ]; then
            TOP="$(awk -F '\t' 'NR==2{print $1}' "$PIP")"
            PIPV="$(awk -F '\t' 'NR==2{print $2}' "$PIP")"
        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\n" \
            "$LOCUS" \
            "$FIT" \
            "$NCS" \
            "$TOP" \
            "$PIPV" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 09. PHASE5 READINESS
# ============================================================

phase5_status() {

    OUT="$ROOT/results/is/stage0_registry/IS_EAS_PHASE5_STATUS.tsv"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    AUDIT="$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_LD_MATRIX_AUDIT.tsv"

    PASS_N=0

    if [ -s "$AUDIT" ]; then
        PASS_N="$(
            awk -F '\t' '
                NR>1 && $6=="PASS" {n++}
                END {print n+0}
            ' "$AUDIT"
        )"
    fi

    if [ "$PASS_N" -eq 4 ]; then
        printf "BBJ_LD\tREADY\t\n" >> "$OUT"
    else
        printf "BBJ_LD\tPARTIAL\t%s/4\n" "$PASS_N" >> "$OUT"
    fi

    SUS="$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY.tsv"

    NFIT=0

    if [ -s "$SUS" ]; then
        NFIT="$(
            awk -F '\t' '
                NR>1 && $2=="YES" {n++}
                END {print n+0}
            ' "$SUS"
        )"
    fi

    if [ "$NFIT" -eq 4 ]; then
        printf "BBJ_SUSIE\tREADY\t\n" >> "$OUT"
    else
        printf "BBJ_SUSIE\tPARTIAL\t%s/4\n" "$NFIT" >> "$OUT"
    fi

    printf "GIGASTROKE_EAS\tCANONICAL_READY\t\n" >> "$OUT"
    printf "TPMI\tWAITING\trate_limit\n" >> "$OUT"
    printf "CKB\tWAITING\tdecryption_key\n" >> "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 5"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_rebuild_gt_pgen" \
    rebuild_gt_pgen

run_step \
    "02_genotype_qc" \
    genotype_qc

run_step \
    "03_build_intersection" \
    build_intersection

run_step \
    "04_create_ld" \
    create_ld

run_step \
    "05_ld_audit" \
    ld_audit

run_step \
    "06_susie_check" \
    susie_check

run_step \
    "07_run_susie" \
    run_susie

run_step \
    "08_susie_summary" \
    susie_summary

run_step \
    "09_phase5_status" \
    phase5_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== GENOTYPE QC ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/ld_qc/GENOTYPE_QC_SUMMARY.tsv" \
    2>/dev/null

echo
echo "===== LD AUDIT ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_LD_MATRIX_AUDIT.tsv" \
    2>/dev/null

echo
echo "===== SUSIE SUMMARY ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/SUSIE_MASTER_SUMMARY.tsv" \
    2>/dev/null

echo
echo "===== PHASE5 STATUS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE5_STATUS.tsv" \
    2>/dev/null

echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "STATUS=$STATUS"
echo "===================================================="

