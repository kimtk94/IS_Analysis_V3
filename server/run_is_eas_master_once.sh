#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="$ROOT/IS_Analysis_V3"

RUN_ID="$(date +%Y%m%d_%H%M%S)"
LOGROOT="$ROOT/logs/is_eas_master/$RUN_ID"
STATUS="$LOGROOT/STATUS.tsv"

mkdir -p "$LOGROOT"

cd "$REPO" || exit 1

printf "step\tstatus\tmessage\tlog\n" > "$STATUS"

echo "===================================================="
echo "IS EAS MASTER ORCHESTRATOR"
echo "RUN_ID=$RUN_ID"
echo "START=$(date)"
echo "LOGROOT=$LOGROOT"
echo "===================================================="


record_status() {
    STEP="$1"
    STATUS_CODE="$2"
    MESSAGE="$3"
    LOG="$4"

    printf "%s\t%s\t%s\t%s\n" \
        "$STEP" \
        "$STATUS_CODE" \
        "$MESSAGE" \
        "$LOG" \
        >> "$STATUS"
}


run_step() {
    STEP="$1"
    shift

    LOG="$LOGROOT/${STEP}.log"

    echo
    echo "===================================================="
    echo "STEP: $STEP"
    echo "START: $(date)"
    echo "LOG: $LOG"
    echo "===================================================="

    (
        "$@"
    ) > >(tee "$LOG") 2>&1

    RC=${PIPESTATUS[0]}

    if [ "$RC" -eq 0 ]; then
        echo "[PASS] $STEP"
        record_status "$STEP" "PASS" "completed" "$LOG"
    else
        echo "[FAIL] $STEP rc=$RC"
        record_status "$STEP" "FAIL" "rc=$RC" "$LOG"
    fi

    return 0
}


skip_step() {
    STEP="$1"
    MESSAGE="$2"
    LOG="$LOGROOT/${STEP}.log"

    echo
    echo "[SKIP] $STEP : $MESSAGE"

    echo "$MESSAGE" > "$LOG"

    record_status \
        "$STEP" \
        "SKIP" \
        "$MESSAGE" \
        "$LOG"
}


# ==========================================================
# STEP 00
# MASTER INVENTORY
# ==========================================================

step00_inventory() {

    echo "===== SYSTEM ====="
    date
    uptime
    df -h "$ROOT" || true

    echo
    echo "===== BBJ ====="

    find \
        "$ROOT/data/is/east_asia/japan/bbj" \
        -maxdepth 4 \
        -type f \
        -printf '%p\t%s\n' \
        2>/dev/null \
        | head -50

    echo
    echo "===== CKB ====="

    find \
        "$ROOT/data/is/east_asia/china/ckb" \
        -maxdepth 2 \
        -type f \
        -printf '%p\t%s\n' \
        2>/dev/null

    echo
    echo "===== TPMI ====="

    find \
        "$ROOT/data/is/east_asia/taiwan/tpmi" \
        -maxdepth 2 \
        -type f \
        -printf '%p\t%s\n' \
        2>/dev/null

    echo
    echo "===== GIGASTROKE ====="

    find \
        "$ROOT/data/is/reference/gigastroke" \
        -maxdepth 3 \
        -type f \
        -printf '%p\t%s\n' \
        2>/dev/null \
        | head -100

    return 0
}


# ==========================================================
# STEP 01
# BBJ CANONICAL INTEGRITY
# ==========================================================

step01_bbj_integrity() {

    FILE="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"

    if [ ! -f "$FILE" ]; then
        echo "BBJ_CANONICAL_MISSING"
        return 10
    fi

    echo "===== FILE ====="
    ls -lh "$FILE"

    echo
    echo "===== GZIP TEST ====="

    gzip -t "$FILE"

    RC=$?

    if [ "$RC" -ne 0 ]; then
        echo "BBJ_CANONICAL_GZIP_INVALID"
        return 11
    fi

    echo "BBJ_CANONICAL_GZIP_OK"

    echo
    echo "===== HEADER ====="

    zcat "$FILE" | head -2

    return 0
}


# ==========================================================
# STEP 02
# BBJ QC + PROVISIONAL LOCI
# ==========================================================

step02_bbj_qc() {

    IN="$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz"
    OUTDIR="$ROOT/results/is/stage1_gwas_qc/japan/bbj"

    mkdir -p "$OUTDIR"

    if [ ! -f "$IN" ]; then
        echo "BBJ_INPUT_MISSING"
        return 20
    fi

    gzip -t "$IN"

    if [ "$?" -ne 0 ]; then
        echo "BBJ_INPUT_INVALID_GZIP"
        return 21
    fi

    python3 - <<'PY'
import csv
import gzip
import json
from pathlib import Path
from collections import Counter, defaultdict

inp = Path(
    "/srv/is-analysis/data/is/processed/japan/bbj/"
    "BBJ_IS_GRCh37.canonical.tsv.gz"
)

outdir = Path(
    "/srv/is-analysis/results/is/stage1_gwas_qc/japan/bbj"
)
outdir.mkdir(parents=True, exist_ok=True)

summary_file = outdir / "BBJ_IS_QC_SUMMARY.json"
loci_file = outdir / "BBJ_IS_PROVISIONAL_LOCI.tsv"
gws_file = outdir / "BBJ_IS_GWS_VARIANTS.tsv"

counts = Counter()
gws = []

with gzip.open(inp, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for r in reader:
        counts["total"] += 1

        try:
            p = float(r["p"])
            eaf = float(r["eaf"])
            info = float(r["info"])
            beta = float(r["beta"])
            se = float(r["se"])
            pos = int(r["pos"])
        except Exception:
            counts["numeric_fail"] += 1
            continue

        if info < 0.7:
            counts["info_lt_0.7"] += 1
            continue

        counts["info_ge_0.7"] += 1

        if not 0 < eaf < 1:
            counts["invalid_eaf"] += 1
            continue

        counts["valid_eaf"] += 1

        maf = min(eaf, 1-eaf)

        if maf >= 0.001:
            counts["maf_ge_0.001"] += 1

        if maf >= 0.01:
            counts["maf_ge_0.01"] += 1

        if p <= 5e-8:
            counts["gws"] += 1

            gws.append({
                "chr": r["chr"],
                "pos": pos,
                "rsid": r["rsid"],
                "ref": r["ref"],
                "alt": r["alt"],
                "effect_allele": r["effect_allele"],
                "other_allele": r["other_allele"],
                "beta": beta,
                "se": se,
                "p": p,
                "eaf": eaf,
                "info": info,
                "variant_id": r["variant_id"]
            })


gws.sort(key=lambda x: x["p"])

fields = [
    "chr","pos","rsid","ref","alt",
    "effect_allele","other_allele",
    "beta","se","p","eaf","info",
    "variant_id"
]

with gws_file.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )
    w.writeheader()
    w.writerows(gws)


by_chr = defaultdict(list)

for x in gws:
    by_chr[x["chr"]].append(x)

loci = []

for chrom, variants in by_chr.items():

    variants.sort(key=lambda x: x["pos"])

    block = []

    for v in variants:

        if not block:
            block = [v]

        elif v["pos"] - block[-1]["pos"] <= 1_000_000:
            block.append(v)

        else:
            loci.append((chrom, block))
            block = [v]

    if block:
        loci.append((chrom, block))


def chromosome_key(chrom):
    try:
        return int(chrom)
    except Exception:
        return 99


loci.sort(
    key=lambda x: (
        chromosome_key(x[0]),
        min(v["pos"] for v in x[1])
    )
)

locus_rows = []

for i, (chrom, variants) in enumerate(loci, 1):

    lead = min(
        variants,
        key=lambda x: x["p"]
    )

    locus_rows.append({
        "locus_id": f"BBJ_IS_L{i:03d}",
        "chr": chrom,
        "start": min(v["pos"] for v in variants),
        "end": max(v["pos"] for v in variants),
        "n_gws_variants": len(variants),
        "lead_variant": lead["variant_id"],
        "lead_rsid": lead["rsid"],
        "lead_p": lead["p"],
        "lead_beta": lead["beta"]
    })


fields2 = [
    "locus_id",
    "chr",
    "start",
    "end",
    "n_gws_variants",
    "lead_variant",
    "lead_rsid",
    "lead_p",
    "lead_beta"
]

with loci_file.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields2
    )
    w.writeheader()
    w.writerows(locus_rows)


summary = {
    "dataset": "BBJ",
    "phenotype": "ischemic_stroke",
    "build": "GRCh37",
    "counts": dict(counts),
    "gws_variants": len(gws),
    "provisional_loci": len(locus_rows)
}

with summary_file.open("w") as f:
    json.dump(
        summary,
        f,
        indent=2
    )

print(json.dumps(summary, indent=2))

print()
print("===== TOP 20 GWS =====")

for x in gws[:20]:
    print(
        x["variant_id"],
        x["rsid"],
        x["beta"],
        x["p"]
    )
PY

    return $?
}


# ==========================================================
# STEP 03
# TPMI PHENOTYPE SELECTION
# ==========================================================

step03_tpmi_select() {

    IN="$ROOT/data/is/east_asia/taiwan/tpmi/probe_download_phenotypes.tsv"
    OUTDIR="$ROOT/results/is/stage0_registry/tpmi"

    mkdir -p "$OUTDIR"

    if [ ! -f "$IN" ]; then
        echo "TPMI_PHENOTYPE_TABLE_MISSING"
        return 30
    fi

    python3 - <<'PY'
import csv
from pathlib import Path

inp = Path(
    "/srv/is-analysis/data/is/east_asia/taiwan/tpmi/"
    "probe_download_phenotypes.tsv"
)

outdir = Path(
    "/srv/is-analysis/results/is/stage0_registry/tpmi"
)

targets = {
    "433.21": "PRIMARY_IS",
    "433.2": "SECONDARY_OCCLUSION",
    "433": "BROAD_CEREBROVASCULAR"
}

with inp.open() as f:
    rows = list(
        csv.DictReader(f, delimiter="\t")
    )

selected = []

for r in rows:
    code = r.get("phenocode", "")

    if code in targets:
        r["analysis_role"] = targets[code]
        selected.append(r)

if not selected:
    raise SystemExit(
        "No target TPMI stroke phenotypes found"
    )

out = outdir / "TPMI_IS_PHENOTYPE_SET.tsv"

fields = list(selected[0].keys())

with out.open("w", newline="") as f:
    w = csv.DictWriter(
        f,
        delimiter="\t",
        fieldnames=fields
    )
    w.writeheader()
    w.writerows(selected)

print(
    "phenocode\trole\tphenotype\tcases\tcontrols"
)

for r in selected:
    print(
        r["phenocode"],
        r["analysis_role"],
        r["phenostring"],
        r["num_cases"],
        r["num_controls"],
        sep="\t"
    )
PY

    return $?
}


# ==========================================================
# STEP 04
# TPMI 433.21 DOWNLOAD DISCOVERY
# ==========================================================

step04_tpmi_download_probe() {

    OUTDIR="$ROOT/data/is/east_asia/taiwan/tpmi/download_probe"
    mkdir -p "$OUTDIR"

    CODE="433.21"
    BASE="https://pheweb.ibms.sinica.edu.tw"

    echo "===== TPMI DOWNLOAD ENDPOINT PROBE ====="

    python3 - <<'PY'
import urllib.request
import urllib.error
from pathlib import Path

code = "433.21"

base = "https://pheweb.ibms.sinica.edu.tw"

targets = [
    f"/download/{code}",
    f"/download/{code}.tsv",
    f"/download/{code}.tsv.gz",
    f"/pheno/{code}",
    f"/api/pheno/{code}",
]

outdir = Path(
    "/srv/is-analysis/data/is/east_asia/"
    "taiwan/tpmi/download_probe"
)

for path in targets:

    url = base + path

    print()
    print("TRY", url)

    req = urllib.request.Request(
        url,
        headers={
            "User-Agent":
            "Mozilla/5.0 IS-MASTER-V2"
        }
    )

    try:

        with urllib.request.urlopen(
            req,
            timeout=60
        ) as r:

            data = r.read(2_000_000)

            ctype = r.headers.get(
                "Content-Type",
                ""
            )

            print(
                "STATUS",
                r.status
            )

            print(
                "CONTENT_TYPE",
                ctype
            )

            print(
                "CONTENT_LENGTH",
                r.headers.get(
                    "Content-Length"
                )
            )

            safe = (
                path.strip("/")
                .replace("/", "_")
            )

            outfile = (
                outdir /
                f"{safe}.probe"
            )

            outfile.write_bytes(data)

            print(
                "SAVED",
                outfile,
                len(data)
            )

    except Exception as e:

        print(
            "ERROR",
            type(e).__name__,
            str(e)
        )
PY

    return 0
}


# ==========================================================
# STEP 05
# GIGASTROKE EAS DOWNLOAD
# ==========================================================

step05_gigastroke_download() {

    OUTDIR="$ROOT/data/is/reference/gigastroke/eas"

    mkdir -p "$OUTDIR"

    BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000"

    for ACC in \
        GCST90104544 \
        GCST90104545 \
        GCST90104546 \
        GCST90104547 \
        GCST90104548
    do

        FILE="${ACC}_buildGRCh37.tsv.gz"
        URL="$BASE/$ACC/$FILE"
        OUT="$OUTDIR/$FILE"

        echo
        echo "===== $ACC ====="

        if [ -s "$OUT" ]; then

            gzip -t "$OUT" 2>/dev/null

            if [ "$?" -eq 0 ]; then
                echo "ALREADY_COMPLETE"
                continue
            fi

            echo "EXISTING_FILE_INVALID_RESUME"
        fi

        wget \
            --continue \
            --tries=10 \
            --timeout=60 \
            -O "$OUT" \
            "$URL"

        RC=$?

        if [ "$RC" -ne 0 ]; then
            echo "DOWNLOAD_FAIL $ACC"
            continue
        fi

        gzip -t "$OUT"

        if [ "$?" -eq 0 ]; then
            echo "DOWNLOAD_OK $ACC"
        else
            echo "GZIP_FAIL $ACC"
        fi

    done

    echo
    echo "===== INVENTORY ====="

    ls -lh "$OUTDIR" 2>/dev/null

    return 0
}


# ==========================================================
# STEP 06
# CKB STATUS
# ==========================================================

step06_ckb_status() {

    ZIP="$ROOT/data/is/east_asia/china/ckb/CKB_i63_IS.zip"

    if [ ! -f "$ZIP" ]; then
        echo "CKB_ZIP_MISSING"
        return 60
    fi

    echo "===== CKB ARCHIVE ====="

    ls -lh "$ZIP"

    echo
    echo "===== CONTENT ====="

    unzip -l "$ZIP"

    echo
    echo "===== SECURITY ====="

    zipinfo -v "$ZIP" \
        | grep -Ei \
          'file security status|encryption|encrypted|compression' \
        | head -30

    echo
    echo "CKB_STATUS=WAITING_DECRYPTION_KEY"

    return 0
}


# ==========================================================
# STEP 07
# GIGASTROKE INVENTORY QC
# ==========================================================

step07_gigastroke_qc() {

    DIR="$ROOT/data/is/reference/gigastroke/eas"
    OUTDIR="$ROOT/results/is/stage1_gwas_qc/gigastroke/eas"

    mkdir -p "$OUTDIR"

    OUT="$OUTDIR/GIGASTROKE_EAS_DOWNLOAD_QC.tsv"

    printf \
        "accession\tphenotype\tfile_exists\tgzip_ok\tsize_bytes\n" \
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

        EXISTS="NO"
        GZIP="NO"
        SIZE="0"

        if [ -f "$FILE" ]; then

            EXISTS="YES"

            SIZE="$(stat -c%s "$FILE" 2>/dev/null || echo 0)"

            gzip -t "$FILE" 2>/dev/null

            if [ "$?" -eq 0 ]; then
                GZIP="YES"
            fi

        fi

        printf \
            "%s\t%s\t%s\t%s\t%s\n" \
            "$ACC" \
            "$PHENO" \
            "$EXISTS" \
            "$GZIP" \
            "$SIZE" \
            >> "$OUT"

    done

    column -t -s $'\t' "$OUT"

    return 0
}


# ==========================================================
# STEP 08
# REGISTRY REFRESH
# ==========================================================

step08_registry_refresh() {

    OUTDIR="$ROOT/results/is/stage0_registry"
    OUT="$OUTDIR/IS_EAS_MASTER_STATUS.tsv"

    mkdir -p "$OUTDIR"

    python3 - <<'PY'
from pathlib import Path
import csv
import gzip
import json
import zipfile

root = Path("/srv/is-analysis")

rows = []

# BBJ
bbj = (
    root /
    "data/is/processed/japan/bbj/"
    "BBJ_IS_GRCh37.canonical.tsv.gz"
)

bbj_qc = (
    root /
    "results/is/stage1_gwas_qc/"
    "japan/bbj/BBJ_IS_QC_SUMMARY.json"
)

rows.append({
    "dataset": "BBJ",
    "country": "Japan",
    "phenotype": "ischemic_stroke",
    "build": "GRCh37",
    "raw_status": "AVAILABLE",
    "canonical_status":
        "DONE" if bbj.exists() else "MISSING",
    "qc_status":
        "DONE" if bbj_qc.exists() else "PENDING",
    "blocker": ""
})


# TPMI
tpmi_reg = (
    root /
    "results/is/stage0_registry/tpmi/"
    "TPMI_IS_PHENOTYPE_SET.tsv"
)

rows.append({
    "dataset": "TPMI",
    "country": "Taiwan",
    "phenotype": "433.21",
    "build": "GRCh38",
    "raw_status":
        "PHENOTYPE_REGISTRY_AVAILABLE",
    "canonical_status": "PENDING",
    "qc_status": "PENDING",
    "blocker":
        "summary_stats_download"
})


# CKB
ckb = (
    root /
    "data/is/east_asia/china/ckb/"
    "CKB_i63_IS.zip"
)

rows.append({
    "dataset": "CKB",
    "country": "China",
    "phenotype": "I63",
    "build": "TBD",
    "raw_status":
        "ENCRYPTED_AVAILABLE"
        if ckb.exists()
        else "MISSING",
    "canonical_status": "BLOCKED",
    "qc_status": "BLOCKED",
    "blocker":
        "decryption_key"
})


# GIGASTROKE
gd = (
    root /
    "data/is/reference/gigastroke/eas"
)

files = list(
    gd.glob(
        "GCST9010454*_buildGRCh37.tsv.gz"
    )
)

rows.append({
    "dataset": "GIGASTROKE_EAS",
    "country": "EastAsia",
    "phenotype": "AS/AIS/CES/LAS/SVS",
    "build": "GRCh37",
    "raw_status":
        f"{len(files)}/5_FILES_PRESENT",
    "canonical_status": "PENDING",
    "qc_status": "DOWNLOAD_QC",
    "blocker":
        "" if len(files) == 5
        else "download_incomplete"
})


out = (
    root /
    "results/is/stage0_registry/"
    "IS_EAS_MASTER_STATUS.tsv"
)

fields = [
    "dataset",
    "country",
    "phenotype",
    "build",
    "raw_status",
    "canonical_status",
    "qc_status",
    "blocker"
]

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

print(out)
PY

    echo
    column -t -s $'\t' "$OUT"

    return 0
}


# ==========================================================
# RUN
# ==========================================================

run_step \
    "00_inventory" \
    step00_inventory


run_step \
    "01_bbj_integrity" \
    step01_bbj_integrity


if gzip -t \
    "$ROOT/data/is/processed/japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz" \
    2>/dev/null
then

    run_step \
        "02_bbj_qc" \
        step02_bbj_qc

else

    skip_step \
        "02_bbj_qc" \
        "BBJ canonical gzip unavailable or invalid"

fi


run_step \
    "03_tpmi_select" \
    step03_tpmi_select


run_step \
    "04_tpmi_download_probe" \
    step04_tpmi_download_probe


run_step \
    "05_gigastroke_download" \
    step05_gigastroke_download


run_step \
    "06_ckb_status" \
    step06_ckb_status


run_step \
    "07_gigastroke_qc" \
    step07_gigastroke_qc


run_step \
    "08_registry_refresh" \
    step08_registry_refresh


# ==========================================================
# FINAL SUMMARY
# ==========================================================

echo
echo
echo "===================================================="
echo "IS EAS MASTER FINAL SUMMARY"
echo "===================================================="

column -t -s $'\t' "$STATUS"

echo
echo "===== MASTER DATA STATUS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage0_registry/IS_EAS_MASTER_STATUS.tsv" \
    2>/dev/null

echo
echo "===== BBJ QC ====="

cat \
    "$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_QC_SUMMARY.json" \
    2>/dev/null

echo
echo
echo "===== BBJ LOCI ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage1_gwas_qc/japan/bbj/BBJ_IS_PROVISIONAL_LOCI.tsv" \
    2>/dev/null \
    | head -40

echo
echo "===== GIGASTROKE EAS ====="

column -t -s $'\t' \
    "$ROOT/results/is/stage1_gwas_qc/gigastroke/eas/GIGASTROKE_EAS_DOWNLOAD_QC.tsv" \
    2>/dev/null

echo
echo "===================================================="
echo "END=$(date)"
echo "STATUS=$STATUS"
echo "LOGROOT=$LOGROOT"
echo "===================================================="

