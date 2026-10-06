#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_phase3/$RUN_ID"
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
# 01. TPMI RETRY CHECK
# ============================================================

tpmi_retry() {

    DIR="$ROOT/data/is/east_asia/taiwan/tpmi"
    OUT="$DIR/TPMI_433.21.tsv"
    URL="https://pheweb.ibms.sinica.edu.tw/download/433.21"

    mkdir -p "$DIR"

    if [ -s "$OUT" ]; then
        echo "TPMI_FULL_ALREADY_PRESENT"
        ls -lh "$OUT"
        return 0
    fi

    echo "===== TPMI RETRY ====="

    curl \
        -L \
        --retry 2 \
        --connect-timeout 30 \
        --max-time 180 \
        -A "Mozilla/5.0" \
        -D "$DIR/TPMI_43321_headers.txt" \
        "$URL" \
        -o "$OUT.tmp"

    RC=$?

    echo "curl_rc=$RC"

    if [ "$RC" -ne 0 ]; then
        rm -f "$OUT.tmp"
        return 10
    fi

    TYPE="$(file -b "$OUT.tmp")"

    echo "FILE_TYPE=$TYPE"

    if grep -qi "403" "$DIR/TPMI_43321_headers.txt"; then
        echo "TPMI_TEMPORARILY_RATE_LIMITED"
        rm -f "$OUT.tmp"
        return 11
    fi

    if [ ! -s "$OUT.tmp" ]; then
        echo "TPMI_EMPTY"
        rm -f "$OUT.tmp"
        return 12
    fi

    mv "$OUT.tmp" "$OUT"

    ls -lh "$OUT"
    head -3 "$OUT"

    return 0
}


# ============================================================
# 02. GIGASTROKE CANONICALIZATION
# ============================================================

gigastroke_canonicalize() {

    INDIR="$ROOT/data/is/reference/gigastroke/eas"
    OUTDIR="$ROOT/data/is/processed/gigastroke/eas"

    mkdir -p "$OUTDIR"

    python3 - <<'PY'
import csv
import gzip
from pathlib import Path

indir = Path(
    "/srv/is-analysis/data/is/reference/gigastroke/eas"
)

outdir = Path(
    "/srv/is-analysis/data/is/processed/gigastroke/eas"
)
outdir.mkdir(parents=True, exist_ok=True)

specs = {
    "GCST90104544": "AS",
    "GCST90104545": "AIS",
    "GCST90104546": "CES",
    "GCST90104547": "LAS",
    "GCST90104548": "SVS",
}

fields = [
    "dataset",
    "phenotype",
    "ancestry",
    "build",
    "chr",
    "pos",
    "effect_allele",
    "other_allele",
    "beta",
    "se",
    "p",
    "eaf",
    "or",
    "variant_id"
]

for acc, pheno in specs.items():

    inp = indir / f"{acc}_buildGRCh37.tsv.gz"
    out = outdir / f"{acc}_{pheno}_GRCh37.canonical.tsv.gz"

    if not inp.exists():
        print("MISSING", acc)
        continue

    n = 0
    bad = 0

    with gzip.open(inp, "rt") as fi, \
         gzip.open(out, "wt") as fo:

        reader = csv.DictReader(fi, delimiter="\t")

        writer = csv.DictWriter(
            fo,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()

        for r in reader:

            try:
                chrom = str(r["chromosome"])
                pos = int(r["base_pair_location"])

                ea = r["effect_allele"]
                oa = r["other_allele"]

                beta = float(r["beta"])
                se = float(r["standard_error"])
                p = float(r["p_value"])
                eaf = float(r["effect_allele_frequency"])

                OR = r.get("odds_ratio", "")

            except Exception:
                bad += 1
                continue

            writer.writerow({
                "dataset": acc,
                "phenotype": pheno,
                "ancestry": "EAS",
                "build": "GRCh37",
                "chr": chrom,
                "pos": pos,
                "effect_allele": ea,
                "other_allele": oa,
                "beta": beta,
                "se": se,
                "p": p,
                "eaf": eaf,
                "or": OR,
                "variant_id":
                    f"{chrom}:{pos}:{oa}:{ea}"
            })

            n += 1

    print(
        acc,
        pheno,
        "rows=",
        n,
        "bad=",
        bad,
        "out=",
        out
    )
PY

    return $?
}


# ============================================================
# 03. BBJ FINEMAP LOCUS PREP
# ============================================================

bbj_finemap_regions() {

    LOCI="$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_LOCUS_SANITY_AUDIT.tsv"
    OUTDIR="$ROOT/results/is/stage3_finemap/japan/bbj"

    mkdir -p "$OUTDIR"

    OUT="$OUTDIR/BBJ_IS_FINEMAP_REGIONS.tsv"

    python3 - <<'PY'
import csv
from pathlib import Path

inp = Path(
    "/srv/is-analysis/results/is/stage1_gwas_qc/"
    "japan/bbj/BBJ_IS_LOCUS_SANITY_AUDIT.tsv"
)

out = Path(
    "/srv/is-analysis/results/is/stage3_finemap/"
    "japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv"
)

rows = []

with inp.open() as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:

        locus = r["locus_id"]

        # rare low-info signal excluded from primary finemap
        if locus == "BBJ_IS_L005":
            role = "EXCLUDE_PRIMARY_RARE_LOW_INFO"
        elif locus == "BBJ_IS_L003":
            role = "PRIMARY_WIDE_REQUIRES_LD_SPLIT"
        else:
            role = "PRIMARY"

        start = int(r["start"])
        end = int(r["end"])

        # +/-500 kb context
        region_start = max(
            1,
            start - 500_000
        )

        region_end = end + 500_000

        rows.append({
            "locus_id": locus,
            "chr": r["chr"],
            "gws_start": start,
            "gws_end": end,
            "region_start": region_start,
            "region_end": region_end,
            "lead_variant": r["lead_variant"],
            "role": role
        })

fields = list(rows[0])

with out.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )
    w.writeheader()
    w.writerows(rows)

for r in rows:
    print(
        r["locus_id"],
        r["chr"],
        r["region_start"],
        r["region_end"],
        r["role"],
        sep="\t"
    )
PY

    column -t -s $'\t' "$OUT"

    return 0
}


# ============================================================
# 04. PREPARE 1000G EAS SAMPLE PANEL
# ============================================================

prepare_1kg_eas_panel() {

    OUTDIR="$ROOT/data/is/ld_reference/1kg_eas_grch37"
    mkdir -p "$OUTDIR"

    PANEL="$OUTDIR/integrated_call_samples_v3.20130502.ALL.panel"

    URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

    if [ ! -s "$PANEL" ]; then

        wget \
            --continue \
            --tries=10 \
            --timeout=60 \
            -O "$PANEL" \
            "$URL"

    fi

    echo "===== PANEL ====="
    head "$PANEL"

    python3 - <<'PY'
from pathlib import Path

panel = Path(
    "/srv/is-analysis/data/is/ld_reference/"
    "1kg_eas_grch37/"
    "integrated_call_samples_v3.20130502.ALL.panel"
)

out = panel.parent / "EAS.samples.txt"

eas_pops = {
    "CHB",
    "JPT",
    "CHS",
    "CDX",
    "KHV"
}

samples = []

for line in panel.read_text().splitlines():

    if not line.strip():
        continue

    cols = line.split("\t")

    if cols[0].lower() == "sample":
        continue

    if len(cols) < 2:
        continue

    sample = cols[0]
    pop = cols[1]

    if pop in eas_pops:
        samples.append(sample)

out.write_text(
    "\n".join(samples) + "\n"
)

print("EAS_SAMPLES =", len(samples))
print("OUTPUT =", out)
PY

    return $?
}


# ============================================================
# 05. LD TOOL CHECK
# ============================================================

ld_tool_check() {

    echo "===== TOOLS ====="

    for x in \
        bcftools \
        tabix \
        bgzip \
        plink2 \
        plink
    do

        printf "%-10s " "$x"

        if command -v "$x" >/dev/null 2>&1; then
            command -v "$x"
        else
            echo "MISSING"
        fi

    done

    return 0
}


# ============================================================
# 06. PHASE3 STATUS
# ============================================================

phase3_status() {

    OUTDIR="$ROOT/results/is/stage0_registry"
    OUT="$OUTDIR/IS_EAS_PHASE3_STATUS.tsv"

    mkdir -p "$OUTDIR"

    printf \
        "component\tstatus\tblocker\n" \
        > "$OUT"

    TPMI="$ROOT/data/is/east_asia/taiwan/tpmi/TPMI_433.21.tsv"

    if [ -s "$TPMI" ]; then
        printf "TPMI\tREADY\t\n" >> "$OUT"
    else
        printf "TPMI\tWAITING\ttemporary_download_restriction\n" >> "$OUT"
    fi

    GIGA="$ROOT/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz"

    if [ -s "$GIGA" ]; then
        printf "GIGASTROKE_EAS\tCANONICAL_READY\t\n" >> "$OUT"
    else
        printf "GIGASTROKE_EAS\tBLOCKED\tcanonicalization\n" >> "$OUT"
    fi

    REGIONS="$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv"

    if [ -s "$REGIONS" ]; then
        printf "BBJ_FINEMAP_REGIONS\tREADY\t\n" >> "$OUT"
    else
        printf "BBJ_FINEMAP_REGIONS\tBLOCKED\tregion_prep\n" >> "$OUT"
    fi

    SAMPLES="$ROOT/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt"

    if [ -s "$SAMPLES" ]; then
        printf "EAS_LD_PANEL\tREADY\tregional_vcf_pending\n" >> "$OUT"
    else
        printf "EAS_LD_PANEL\tBLOCKED\tpanel_download\n" >> "$OUT"
    fi

    printf "CKB\tWAITING\tdecryption_key\n" >> "$OUT"

    column -t -s $'\t' "$OUT"

    return 0
}


echo "===================================================="
echo "IS EAS MASTER PHASE 3"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "===================================================="

run_step \
    "01_tpmi_retry" \
    tpmi_retry

run_step \
    "02_gigastroke_canonicalize" \
    gigastroke_canonicalize

run_step \
    "03_bbj_finemap_regions" \
    bbj_finemap_regions

run_step \
    "04_prepare_1kg_eas_panel" \
    prepare_1kg_eas_panel

run_step \
    "05_ld_tool_check" \
    ld_tool_check

run_step \
    "06_phase3_status" \
    phase3_status


echo
echo "===================================================="
echo "FINAL STATUS"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== PHASE3 READINESS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_PHASE3_STATUS.tsv" \
    2>/dev/null

echo
echo "===== BBJ FINEMAP REGIONS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage3_finemap/japan/bbj/BBJ_IS_FINEMAP_REGIONS.tsv" \
    2>/dev/null

echo
echo "===== EAS SAMPLE COUNT ====="

wc -l \
    "$ROOT/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt" \
    2>/dev/null

echo
echo "END=$(date)"
echo "LOGROOT=$LOGROOT"
echo "STATUS=$STATUS"
