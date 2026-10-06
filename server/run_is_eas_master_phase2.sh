#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase2/$RUN_ID"
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
        record "$STEP" "PASS" "completed" "$LOG"
        echo "[PASS] $STEP"
    else
        record "$STEP" "FAIL" "rc=$RC" "$LOG"
        echo "[FAIL] $STEP rc=$RC"
    fi

    return 0
}


# ============================================================
# 01. TPMI FULL DOWNLOAD
# ============================================================

tpmi_download() {

    DIR="$ROOT/data/is/east_asia/taiwan/tpmi"
    OUT="$DIR/TPMI_433.21.tsv"
    URL="https://pheweb.ibms.sinica.edu.tw/download/433.21"

    mkdir -p "$DIR"

    echo "===== TPMI 433.21 ====="

    if [ -s "$OUT" ]; then
        echo "EXISTING_FILE"
        ls -lh "$OUT"
        return 0
    fi

    wget \
        --continue \
        --tries=10 \
        --timeout=60 \
        -O "$OUT" \
        "$URL"

    RC=$?

    echo
    echo "wget_rc=$RC"

    if [ "$RC" -ne 0 ]; then
        return "$RC"
    fi

    echo
    ls -lh "$OUT"

    echo
    echo "===== HEADER ====="
    head -3 "$OUT"

    echo
    echo "===== LINE COUNT ====="
    wc -l "$OUT"

    return 0
}


# ============================================================
# 02. TPMI SCHEMA AUDIT
# ============================================================

tpmi_schema() {

    IN="$ROOT/data/is/east_asia/taiwan/tpmi/TPMI_433.21.tsv"
    OUTDIR="$ROOT/results/is/stage1_gwas_qc/taiwan/tpmi"

    mkdir -p "$OUTDIR"

    if [ ! -s "$IN" ]; then
        echo "TPMI_FULL_FILE_MISSING"
        return 20
    fi

    python3 - <<'PY'
import csv
import json
from pathlib import Path

inp = Path(
    "/srv/is-analysis/data/is/east_asia/"
    "taiwan/tpmi/TPMI_433.21.tsv"
)

outdir = Path(
    "/srv/is-analysis/results/is/"
    "stage1_gwas_qc/taiwan/tpmi"
)

with inp.open(errors="replace") as f:
    reader = csv.DictReader(f, delimiter="\t")

    fields = reader.fieldnames or []

    rows = []

    for i, r in enumerate(reader):
        rows.append(r)

        if i >= 4:
            break

print("N_COLUMNS =", len(fields))
print("COLUMNS =", fields)

print()
print("===== FIRST ROWS =====")

for r in rows:
    print(r)

lower = {
    x.lower(): x
    for x in fields
}

groups = {
    "chr": [
        "chrom", "chr", "chromosome"
    ],
    "pos": [
        "pos", "position", "base_pair_location"
    ],
    "ref": [
        "ref", "reference_allele"
    ],
    "alt": [
        "alt", "alternate_allele",
        "effect_allele"
    ],
    "beta": [
        "beta", "effect", "effect_size"
    ],
    "se": [
        "se", "stderr",
        "standard_error", "sebeta"
    ],
    "p": [
        "pval", "p", "pvalue", "p_value"
    ],
    "eaf": [
        "af", "eaf",
        "effect_allele_frequency"
    ]
}

mapping = {}

for canonical, candidates in groups.items():

    mapping[canonical] = ""

    for c in candidates:

        if c in lower:
            mapping[canonical] = lower[c]
            break

result = {
    "n_columns": len(fields),
    "columns": fields,
    "detected_mapping": mapping
}

with (
    outdir /
    "TPMI_43321_SCHEMA.json"
).open("w") as f:

    json.dump(
        result,
        f,
        indent=2
    )

print()
print("===== MAPPING =====")
print(json.dumps(mapping, indent=2))
PY

    return $?
}


# ============================================================
# 03. GIGASTROKE SCHEMA AUDIT
# ============================================================

gigastroke_schema() {

    DIR="$ROOT/data/is/reference/gigastroke/eas"
    OUTDIR="$ROOT/results/is/stage1_gwas_qc/gigastroke/eas"

    mkdir -p "$OUTDIR"

    OUT="$OUTDIR/GIGASTROKE_EAS_SCHEMA_AUDIT.tsv"

    printf \
        "accession\tphenotype\tn_columns\theader\n" \
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

        FILE="$DIR/${ACC}_buildGRCh37.tsv.gz"

        echo
        echo "===== $ACC / $PHENO ====="

        if [ ! -s "$FILE" ]; then
            echo "MISSING"
            continue
        fi

        gzip -t "$FILE" 2>/dev/null

        if [ "$?" -ne 0 ]; then
            echo "INVALID_GZIP"
            continue
        fi

        HEADER="$(zcat "$FILE" | head -1)"
        NCOL="$(printf "%s\n" "$HEADER" | awk -F '\t' '{print NF}')"

        echo "$HEADER"
        echo "N_COLUMNS=$NCOL"

        SAFE_HEADER="$(printf "%s" "$HEADER" | tr '\t' '|')"

        printf \
            "%s\t%s\t%s\t%s\n" \
            "$ACC" \
            "$PHENO" \
            "$NCOL" \
            "$SAFE_HEADER" \
            >> "$OUT"

    done

    echo
    echo "===== SUMMARY ====="
    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 04. BBJ LOCUS SANITY AUDIT
# ============================================================

bbj_locus_audit() {

    GWAS="$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_GWS_VARIANTS.tsv"
    LOCI="$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_PROVISIONAL_LOCI.tsv"

    OUTDIR="$ROOT/results/is/stage1_gwas_qc/japan/bbj"
    OUT="$OUTDIR/BBJ_IS_LOCUS_SANITY_AUDIT.tsv"

    if [ ! -s "$GWAS" ] || [ ! -s "$LOCI" ]; then
        echo "BBJ_QC_OUTPUT_MISSING"
        return 40
    fi

    python3 - <<'PY'
import csv
from pathlib import Path

gwas = Path(
    "/srv/is-analysis/results/is/stage1_gwas_qc/"
    "japan/bbj/BBJ_IS_GWS_VARIANTS.tsv"
)

loci = Path(
    "/srv/is-analysis/results/is/stage1_gwas_qc/"
    "japan/bbj/BBJ_IS_PROVISIONAL_LOCI.tsv"
)

out = Path(
    "/srv/is-analysis/results/is/stage1_gwas_qc/"
    "japan/bbj/BBJ_IS_LOCUS_SANITY_AUDIT.tsv"
)

with gwas.open() as f:
    variants = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

with loci.open() as f:
    locus_rows = list(
        csv.DictReader(
            f,
            delimiter="\t"
        )
    )

by_chr = {}

for v in variants:
    by_chr.setdefault(
        v["chr"],
        []
    ).append(v)


audit = []

for L in locus_rows:

    chrom = L["chr"]
    start = int(L["start"])
    end = int(L["end"])

    vv = [
        x for x in by_chr.get(chrom, [])
        if start <= int(x["pos"]) <= end
    ]

    vv.sort(
        key=lambda x: int(x["pos"])
    )

    positions = [
        int(x["pos"])
        for x in vv
    ]

    gaps = [
        positions[i] - positions[i-1]
        for i in range(1, len(positions))
    ]

    max_gap = max(gaps) if gaps else 0

    lead_id = L["lead_variant"]

    lead = next(
        (
            x for x in vv
            if x["variant_id"] == lead_id
        ),
        None
    )

    lead_eaf = (
        float(lead["eaf"])
        if lead else None
    )

    lead_maf = (
        min(
            lead_eaf,
            1-lead_eaf
        )
        if lead_eaf is not None
        else None
    )

    lead_info = (
        float(lead["info"])
        if lead else None
    )

    lead_beta = (
        float(lead["beta"])
        if lead else None
    )

    width = end - start

    flags = []

    if width > 1_000_000:
        flags.append(
            "LOCUS_WIDTH_GT_1MB"
        )

    if max_gap > 500_000:
        flags.append(
            "LARGE_INTERNAL_GAP"
        )

    if (
        lead_maf is not None
        and lead_maf < 0.01
    ):
        flags.append(
            "LEAD_MAF_LT_1PCT"
        )

    if (
        lead_info is not None
        and lead_info < 0.9
    ):
        flags.append(
            "LEAD_INFO_LT_0.9"
        )

    if (
        lead_beta is not None
        and abs(lead_beta) > 0.5
    ):
        flags.append(
            "LARGE_ABS_BETA"
        )

    audit.append({
        "locus_id":
            L["locus_id"],

        "chr":
            chrom,

        "start":
            start,

        "end":
            end,

        "width_bp":
            width,

        "n_gws":
            len(vv),

        "max_internal_gap_bp":
            max_gap,

        "lead_variant":
            lead_id,

        "lead_maf":
            lead_maf,

        "lead_info":
            lead_info,

        "lead_beta":
            lead_beta,

        "flags":
            ";".join(flags)
    })


fields = list(audit[0])

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
    w.writerows(audit)

for x in audit:
    print(
        x["locus_id"],
        x["width_bp"],
        x["lead_maf"],
        x["lead_info"],
        x["lead_beta"],
        x["flags"],
        sep="\t"
    )
PY

    echo
    echo "===== AUDIT ====="

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 05. EXISTING LD REFERENCE DISCOVERY
# ============================================================

ld_discovery() {

    OUTDIR="$ROOT/results/is/stage0_registry/ld"
    OUT="$OUTDIR/EAS_LD_DISCOVERY.txt"

    mkdir -p "$OUTDIR"

    {
        echo "===== CANDIDATE LD DIRECTORIES ====="

        find "$ROOT/data" \
            -maxdepth 6 \
            \( \
              -iname '*1000g*' -o \
              -iname '*1kg*' -o \
              -iname '*eas*' -o \
              -iname '*japan*' -o \
              -iname '*japanese*' -o \
              -iname '*bbj*' \
            \) \
            -print \
            2>/dev/null \
            | head -300

        echo
        echo "===== PLINK/BGEN/PGEN FILES ====="

        find "$ROOT/data" \
            -maxdepth 7 \
            \( \
              -name '*.bed' -o \
              -name '*.bim' -o \
              -name '*.fam' -o \
              -name '*.pgen' -o \
              -name '*.pvar' -o \
              -name '*.psam' -o \
              -name '*.bgen' \
            \) \
            -print \
            2>/dev/null \
            | head -500

    } | tee "$OUT"

    return 0
}


# ============================================================
# 06. MASTER READINESS
# ============================================================

phase2_status() {

    OUTDIR="$ROOT/results/is/stage0_registry"
    OUT="$OUTDIR/IS_EAS_PHASE2_STATUS.tsv"

    mkdir -p "$OUTDIR"

    TPMI="$ROOT/data/is/east_asia/taiwan/tpmi/TPMI_433.21.tsv"

    GS="$ROOT/results/is/stage1_gwas_qc/gigastroke/eas/GIGASTROKE_EAS_SCHEMA_AUDIT.tsv"

    BBJ="$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_LOCUS_SANITY_AUDIT.tsv"

    LD="$ROOT/results/is/stage0_registry/ld/EAS_LD_DISCOVERY.txt"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    if [ -s "$TPMI" ]; then
        printf "TPMI_43321\tREADY\t\n" >> "$OUT"
    else
        printf "TPMI_43321\tBLOCKED\tdownload\n" >> "$OUT"
    fi

    if [ -s "$GS" ]; then
        printf "GIGASTROKE_EAS\tREADY\t\n" >> "$OUT"
    else
        printf "GIGASTROKE_EAS\tBLOCKED\tschema_audit\n" >> "$OUT"
    fi

    if [ -s "$BBJ" ]; then
        printf "BBJ_LOCI\tREADY\t\n" >> "$OUT"
    else
        printf "BBJ_LOCI\tBLOCKED\tlocus_audit\n" >> "$OUT"
    fi

    if [ -s "$LD" ]; then
        printf "LD_REFERENCE_DISCOVERY\tDONE\tmanual_selection_pending\n" >> "$OUT"
    else
        printf "LD_REFERENCE_DISCOVERY\tBLOCKED\tdiscovery\n" >> "$OUT"
    fi

    printf "CKB\tBLOCKED\tdecryption_key\n" >> "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# RUN
# ============================================================

echo "===================================================="
echo "IS EAS MASTER PHASE 2"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_tpmi_download" \
    tpmi_download

run_step \
    "02_tpmi_schema" \
    tpmi_schema

run_step \
    "03_gigastroke_schema" \
    gigastroke_schema

run_step \
    "04_bbj_locus_audit" \
    bbj_locus_audit

run_step \
    "05_ld_discovery" \
    ld_discovery

run_step \
    "06_phase2_status" \
    phase2_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== PHASE2 READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE2_STATUS.tsv" \
    2>/dev/null

echo
echo "===== BBJ LOCUS AUDIT ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_LOCUS_SANITY_AUDIT.tsv" \
    2>/dev/null

echo
echo "===== TPMI SCHEMA ====="

cat \
    "$ROOT/results/is/stage1_gwas_qc/taiwan/tpmi/TPMI_43321_SCHEMA.json" \
    2>/dev/null

echo
echo "===================================================="
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "STATUS=$STATUS"
echo "===================================================="
