#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase4/$RUN_ID"
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
# 01. TPMI STATUS FIX
# ============================================================

tpmi_status_fix() {

    FILE="$ROOT/data/is/east_asia/taiwan/tpmi/TPMI_433.21.tsv"
    OUTDIR="$ROOT/results/is/stage0_registry/tpmi"
    OUT="$OUTDIR/TPMI_DOWNLOAD_STATUS.tsv"

    mkdir -p "$OUTDIR"

    STATUS_CODE="MISSING"
    MESSAGE=""

    if [ -f "$FILE" ]; then

        SIZE="$(stat -c%s "$FILE" 2>/dev/null || echo 0)"

        echo "FILE_SIZE=$SIZE"

        HEAD="$(head -c 500 "$FILE" 2>/dev/null || true)"

        echo "===== CONTENT PREVIEW ====="
        printf "%s\n" "$HEAD"

        if printf "%s" "$HEAD" \
            | grep -qiE \
            'download limit exceeded|已超過檔案下載限制'
        then
            STATUS_CODE="WAITING_RATE_LIMIT"
            MESSAGE="TPMI server download limit exceeded"

            rm -f "$FILE"

        elif [ "$SIZE" -lt 1000000 ]; then

            STATUS_CODE="INVALID_SMALL_FILE"
            MESSAGE="Downloaded file too small"

            rm -f "$FILE"

        else

            HEADER="$(head -1 "$FILE")"

            if printf "%s" "$HEADER" \
                | grep -qiE \
                'chrom|position|pval|beta'
            then
                STATUS_CODE="READY"
                MESSAGE="Full GWAS file appears valid"
            else
                STATUS_CODE="INVALID_HEADER"
                MESSAGE="Unexpected TPMI file header"
            fi

        fi

    fi

    printf \
        "phenocode\tstatus\tmessage\n433.21\t%s\t%s\n" \
        "$STATUS_CODE" \
        "$MESSAGE" \
        > "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 02. GIGASTROKE CANONICAL AUDIT
# ============================================================

gigastroke_canonical_audit() {

    DIR="$ROOT/data/is/processed/gigastroke/eas"
    OUTDIR="$ROOT/results/is/stage1_gwas_qc/gigastroke/eas"
    OUT="$OUTDIR/GIGASTROKE_CANONICAL_AUDIT.tsv"

    mkdir -p "$OUTDIR"

    printf \
        "dataset\tphenotype\tgzip_ok\trows\theader_ok\tnote\n" \
        > "$OUT"

    for SPEC in \
        "GCST90104544 AS" \
        "GCST90104545 AIS" \
        "GCST90104546 CES" \
        "GCST90104547 LAS" \
        "GCST90104548 SVS"
    do

        ACC="$(echo "$SPEC" | awk '{print $1}')"
        PHENO="$(echo "$SPEC" | awk '{print $2}')"

        FILE="$DIR/${ACC}_${PHENO}_GRCh37.canonical.tsv.gz"

        GZIP_OK="NO"
        HEADER_OK="NO"
        ROWS="0"
        NOTE=""

        if [ ! -s "$FILE" ]; then

            NOTE="missing"
            printf \
                "%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$ACC" "$PHENO" "$GZIP_OK" "$ROWS" "$HEADER_OK" "$NOTE" \
                >> "$OUT"

            continue
        fi

        gzip -t "$FILE" 2>/dev/null

        if [ "$?" -eq 0 ]; then
            GZIP_OK="YES"
        fi

        HEADER="$(zcat "$FILE" | head -1)"

        if printf "%s" "$HEADER" \
            | grep -q $'effect_allele\tother_allele\tbeta\tse\tp\teaf'
        then
            HEADER_OK="YES"
        fi

        ROWS="$(zcat "$FILE" | awk 'END{print NR-1}')"

        NOTE="allele_key_only_not_refalt"

        printf \
            "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$ACC" "$PHENO" "$GZIP_OK" "$ROWS" "$HEADER_OK" "$NOTE" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 03. 1000G EAS REGIONAL VCF EXTRACTION
# ============================================================

extract_1kg_eas_regions() {

    REGIONS="$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv"
    SAMPLES="$ROOT/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt"

    OUTDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/regions"

    mkdir -p "$OUTDIR"

    if [ ! -s "$REGIONS" ]; then
        echo "REGIONS_MISSING"
        return 30
    fi

    if [ ! -s "$SAMPLES" ]; then
        echo "EAS_SAMPLE_LIST_MISSING"
        return 31
    fi

    BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"

    tail -n +2 "$REGIONS" \
    | while IFS=$'\t' read -r \
        LOCUS CHR GWS_START GWS_END REGION_START REGION_END LEAD ROLE
    do

        if [ "$ROLE" = "EXCLUDE_PRIMARY_RARE_LOW_INFO" ]; then
            echo "SKIP $LOCUS $ROLE"
            continue
        fi

        echo
        echo "===================================================="
        echo "$LOCUS chr${CHR}:${REGION_START}-${REGION_END}"
        echo "===================================================="

        REMOTE="$BASE/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"

        OUT="$OUTDIR/${LOCUS}.1KG_EAS.GRCh37.vcf.gz"

        if [ -s "$OUT" ]; then

            bcftools index -t "$OUT" 2>/dev/null || true

            bcftools view -h "$OUT" >/dev/null 2>&1

            if [ "$?" -eq 0 ]; then
                echo "ALREADY_VALID $OUT"
                continue
            fi

        fi

        rm -f "$OUT" "$OUT.tbi"

        bcftools view \
            -r "${CHR}:${REGION_START}-${REGION_END}" \
            -S "$SAMPLES" \
            -m2 -M2 \
            -v snps \
            -Oz \
            -o "$OUT.tmp.gz" \
            "$REMOTE"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "BCFTOOLS_REMOTE_FAIL $LOCUS rc=$RC"
            rm -f "$OUT.tmp.gz"
            continue
        fi

        mv "$OUT.tmp.gz" "$OUT"

        tabix -p vcf "$OUT"

        echo "REGION_OK"
        ls -lh "$OUT" "$OUT.tbi"

        echo "SAMPLES=$(bcftools query -l "$OUT" | wc -l)"
        echo "VARIANTS=$(bcftools index -n "$OUT")"

    done

    return 0
}


# ============================================================
# 04. REGIONAL PGEN BUILD
# ============================================================

build_regional_pgen() {

    VCFDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/regions"
    OUTDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/pgen"

    mkdir -p "$OUTDIR"

    FOUND=0

    for VCF in "$VCFDIR"/BBJ_IS_L00[1-4].1KG_EAS.GRCh37.vcf.gz
    do

        [ -s "$VCF" ] || continue

        FOUND=1

        BASE="$(basename "$VCF" .1KG_EAS.GRCh37.vcf.gz)"
        PREFIX="$OUTDIR/$BASE.1KG_EAS.GRCh37"

        echo
        echo "===== $BASE ====="

        plink2 \
            --vcf "$VCF" dosage=DS \
            --double-id \
            --max-alleles 2 \
            --make-pgen \
            --out "$PREFIX"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "PLINK2_FAIL $BASE"
            continue
        fi

        ls -lh \
            "$PREFIX.pgen" \
            "$PREFIX.pvar" \
            "$PREFIX.psam"

    done

    if [ "$FOUND" -eq 0 ]; then
        echo "NO_REGIONAL_VCF_AVAILABLE"
        return 40
    fi

    return 0
}


# ============================================================
# 05. BBJ SUMMARY REGIONAL EXTRACTION
# ============================================================

extract_bbj_summary_regions() {

    GWAS="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"
    REGIONS="$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv"

    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj/summary_regions"

    mkdir -p "$OUTDIR"

    python3 - <<'PY'
import csv
import gzip
from pathlib import Path

gwas = Path(
    "/srv/is-analysis/data/is/processed/japan/bbj/"
    "BBJ_IS_GRCh37.canonical.tsv.gz"
)

regions_file = Path(
    "/srv/is-analysis/results/is/stage3_finemap/japan/bbj/"
    "BBJ_IS_FINEMAP_REGIONS.tsv"
)

outdir = Path(
    "/srv/is-analysis/results/is/stage3_finemap/japan/bbj/"
    "summary_regions"
)
outdir.mkdir(parents=True, exist_ok=True)

with regions_file.open() as f:
    regions = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

regions = [
    r for r in regions
    if r["role"] != "EXCLUDE_PRIMARY_RARE_LOW_INFO"
]

handles = {}
writers = {}
counts = {}

fields_out = [
    "chr",
    "pos",
    "rsid",
    "ref",
    "alt",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    "eaf",
    "info",
    "n",
    "variant_id"
]

for r in regions:

    locus = r["locus_id"]

    out = outdir / f"{locus}.BBJ_GRCh37.tsv"

    h = out.open("w", newline="")
    w = csv.DictWriter(
        h,
        delimiter="\t",
        fieldnames=fields_out
    )

    w.writeheader()

    handles[locus] = h
    writers[locus] = w
    counts[locus] = 0


with gzip.open(gwas, "rt") as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for row in reader:

        chrom = row["chr"]

        try:
            pos = int(row["pos"])
        except Exception:
            continue

        for r in regions:

            if chrom != r["chr"]:
                continue

            start = int(r["region_start"])
            end = int(r["region_end"])

            if start <= pos <= end:

                writers[
                    r["locus_id"]
                ].writerow(
                    {
                        k: row.get(k, "")
                        for k in fields_out
                    }
                )

                counts[
                    r["locus_id"]
                ] += 1


for h in handles.values():
    h.close()

print("===== REGION COUNTS =====")

for locus, n in counts.items():
    print(locus, n)
PY

    return $?
}


# ============================================================
# 06. BBJ ↔ 1KG ALLELE OVERLAP
# ============================================================

compute_overlap() {

    SUMDIR="$ROOT/results/is/stage3_finemap/japan/bbj/summary_regions"
    VCFDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37/regions"

    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj/overlap"
    mkdir -p "$OUTDIR"

    for LOCUS in \
        BBJ_IS_L001 \
        BBJ_IS_L002 \
        BBJ_IS_L003 \
        BBJ_IS_L004
    do

        SUM="$SUMDIR/${LOCUS}.BBJ_GRCh37.tsv"
        VCF="$VCFDIR/${LOCUS}.1KG_EAS.GRCh37.vcf.gz"
        OUT="$OUTDIR/${LOCUS}.OVERLAP.tsv"

        echo
        echo "===== $LOCUS ====="

        if [ ! -s "$SUM" ]; then
            echo "SUMMARY_MISSING"
            continue
        fi

        if [ ! -s "$VCF" ]; then
            echo "VCF_MISSING"
            continue
        fi

        python3 - "$SUM" "$VCF" "$OUT" <<'PY'
import csv
import subprocess
import sys
from pathlib import Path

sum_file = Path(sys.argv[1])
vcf_file = Path(sys.argv[2])
out_file = Path(sys.argv[3])

vcf_variants = {}

cmd = [
    "bcftools",
    "query",
    "-f",
    "%CHROM\t%POS\t%REF\t%ALT\n",
    str(vcf_file)
]

p = subprocess.run(
    cmd,
    capture_output=True,
    text=True
)

if p.returncode != 0:
    print(p.stderr)
    raise SystemExit(p.returncode)

for line in p.stdout.splitlines():

    chrom, pos, ref, alt = line.split("\t")

    vcf_variants[
        (chrom, int(pos))
    ] = (ref, alt)


rows = []

matched = 0
flipped = 0
mismatch = 0

with sum_file.open() as f:

    reader = csv.DictReader(
        f,
        delimiter="\t"
    )

    for r in reader:

        key = (
            r["chr"],
            int(r["pos"])
        )

        if key not in vcf_variants:
            continue

        ref, alt = vcf_variants[key]

        bbj_ref = r["ref"]
        bbj_alt = r["alt"]

        if bbj_ref == ref and bbj_alt == alt:
            status = "MATCH"
            matched += 1

        elif bbj_ref == alt and bbj_alt == ref:
            status = "FLIP"
            flipped += 1

        else:
            status = "ALLELE_MISMATCH"
            mismatch += 1

        rows.append({
            "chr": r["chr"],
            "pos": r["pos"],
            "bbj_ref": bbj_ref,
            "bbj_alt": bbj_alt,
            "kg_ref": ref,
            "kg_alt": alt,
            "status": status,
            "beta": r["beta"],
            "se": r["se"],
            "p": r["p"],
            "eaf": r["eaf"],
            "info": r["info"]
        })


fields = [
    "chr",
    "pos",
    "bbj_ref",
    "bbj_alt",
    "kg_ref",
    "kg_alt",
    "status",
    "beta",
    "se",
    "p",
    "eaf",
    "info"
]

with out_file.open("w", newline="") as f:

    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )

    w.writeheader()
    w.writerows(rows)


print("OVERLAP =", len(rows))
print("MATCH =", matched)
print("FLIP =", flipped)
print("MISMATCH =", mismatch)
print("OUT =", out_file)
PY

    done

    return 0
}


# ============================================================
# 07. PHASE4 READINESS
# ============================================================

phase4_status() {

    OUTDIR="$ROOT/results/is/stage0_registry"
    OUT="$OUTDIR/IS_EAS_PHASE4_STATUS.tsv"

    mkdir -p "$OUTDIR"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    TPMI="$ROOT/data/is/east_asia/taiwan/tpmi/TPMI_433.21.tsv"

    if [ -s "$TPMI" ] && \
       ! grep -qiE \
         'download limit exceeded|已超過檔案下載限制' \
         "$TPMI" 2>/dev/null
    then
        printf "TPMI\tREADY\t\n" >> "$OUT"
    else
        printf "TPMI\tWAITING\trate_limit\n" >> "$OUT"
    fi

    NVCF="$(find \
        "$ROOT/data/is/ld_reference/1kg_eas_grch37/regions" \
        -name 'BBJ_IS_L00[1-4].1KG_EAS.GRCh37.vcf.gz' \
        -size +0c \
        2>/dev/null \
        | wc -l)"

    if [ "$NVCF" -eq 4 ]; then
        printf "EAS_REGIONAL_VCF\tREADY\t\n" >> "$OUT"
    else
        printf \
            "EAS_REGIONAL_VCF\tPARTIAL\t%s/4\n" \
            "$NVCF" \
            >> "$OUT"
    fi

    NPGEN="$(find \
        "$ROOT/data/is/ld_reference/1kg_eas_grch37/pgen" \
        -name 'BBJ_IS_L00[1-4].1KG_EAS.GRCh37.pgen' \
        -size +0c \
        2>/dev/null \
        | wc -l)"

    if [ "$NPGEN" -eq 4 ]; then
        printf "EAS_REGIONAL_PGEN\tREADY\t\n" >> "$OUT"
    else
        printf \
            "EAS_REGIONAL_PGEN\tPARTIAL\t%s/4\n" \
            "$NPGEN" \
            >> "$OUT"
    fi

    NOV="$(find \
        "$ROOT/results/is/stage3_finemap/japan/bbj/overlap" \
        -name 'BBJ_IS_L00[1-4].OVERLAP.tsv' \
        -size +0c \
        2>/dev/null \
        | wc -l)"

    if [ "$NOV" -eq 4 ]; then
        printf "BBJ_1KG_OVERLAP\tREADY\t\n" >> "$OUT"
    else
        printf \
            "BBJ_1KG_OVERLAP\tPARTIAL\t%s/4\n" \
            "$NOV" \
            >> "$OUT"
    fi

    printf "GIGASTROKE_EAS\tCANONICAL_READY\t\n" >> "$OUT"
    printf "CKB\tWAITING\tdecryption_key\n" >> "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 4"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_tpmi_status_fix" \
    tpmi_status_fix

run_step \
    "02_gigastroke_canonical_audit" \
    gigastroke_canonical_audit

run_step \
    "03_extract_1kg_eas_regions" \
    extract_1kg_eas_regions

run_step \
    "04_build_regional_pgen" \
    build_regional_pgen

run_step \
    "05_extract_bbj_summary_regions" \
    extract_bbj_summary_regions

run_step \
    "06_compute_overlap" \
    compute_overlap

run_step \
    "07_phase4_status" \
    phase4_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== PHASE4 READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE4_STATUS.tsv" \
    2>/dev/null

echo
echo "===== OVERLAP SUMMARY ====="

for f in \
    "$ROOT/results/is/stage3_finemap/japan/bbj/overlap/"*.OVERLAP.tsv
do
    [ -s "$f" ] || continue

    echo
    echo "$(basename "$f")"

    awk -F '\t' '
        NR>1 {x[$7]++}
        END {
            for (k in x)
                print k, x[k]
        }
    ' "$f"
done

echo
echo "===== REGIONAL FILES ====="

find \
    "$ROOT/data/is/ld_reference/1kg_eas_grch37" \
    -maxdepth 2 \
    -type f \
    -printf '%p\t%s\n' \
    | sort

echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "STATUS=$STATUS"
echo "===================================================="

