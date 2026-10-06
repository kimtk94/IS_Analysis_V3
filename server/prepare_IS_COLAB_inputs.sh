#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REMOTE="gdrive:MASTER_DEGREE/IS_COLAB/input"

FAIL=0
PASS=0

echo "===================================================="
echo "IS COLAB INPUT EXPORT"
echo "START=$(date)"
echo "REMOTE=$REMOTE"
echo "===================================================="


copy_one() {

    SRC="$1"
    DST="$2"

    echo
    echo "----------------------------------------------------"
    echo "SRC=$SRC"
    echo "DST=$REMOTE/$DST"
    echo "----------------------------------------------------"

    if [ ! -s "$SRC" ]; then

        echo "[FAIL] SOURCE_MISSING_OR_EMPTY"

        FAIL=$((FAIL + 1))

        return 0
    fi


    rclone copyto \
        "$SRC" \
        "$REMOTE/$DST" \
        --progress


    RC=$?


    if [ "$RC" -eq 0 ]; then

        echo "[PASS] $DST"

        PASS=$((PASS + 1))

    else

        echo "[FAIL] rc=$RC $DST"

        FAIL=$((FAIL + 1))
    fi


    return 0
}


# ============================================================
# 1. PHASE 9 TABLES
# ============================================================

copy_one \
    "$ROOT/results/is/stage5_functional/phase9/BBJ_CREDIBLE_SET_MASTER.tsv" \
    "BBJ_CREDIBLE_SET_MASTER.tsv"


copy_one \
    "$ROOT/results/is/stage5_functional/phase9/phase9b/PHASE9B_LOCUS_MASTER.tsv" \
    "PHASE9B_LOCUS_MASTER.tsv"


copy_one \
    "$ROOT/results/is/stage5_functional/phase9/GRCH37_PROTEIN_CODING_GENES.tsv" \
    "GRCH37_PROTEIN_CODING_GENES.tsv"


# ============================================================
# 2. FOUR BBJ REGIONAL INPUTS
# ============================================================

for LOCUS in \
    BBJ_IS_L001 \
    BBJ_IS_L002 \
    BBJ_IS_L003 \
    BBJ_IS_L004
do

    copy_one \
        "$ROOT/results/is/stage3_finemap/japan/bbj/susie_inputs_v3/${LOCUS}.tsv" \
        "bbj_locus/${LOCUS}.tsv"

done


# ============================================================
# 3. FOUR LD MATRICES + FOUR VAR LISTS
# ============================================================

for LOCUS in \
    BBJ_IS_L001 \
    BBJ_IS_L002 \
    BBJ_IS_L003 \
    BBJ_IS_L004
do

    copy_one \
        "$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/${LOCUS}.unphased.vcor1.bin" \
        "ld/${LOCUS}.unphased.vcor1.bin"


    copy_one \
        "$ROOT/results/is/stage3_finemap/japan/bbj/ld_v3/${LOCUS}.unphased.vcor1.bin.vars" \
        "ld/${LOCUS}.unphased.vcor1.bin.vars"

done


# ============================================================
# 4. OPTIONAL L004 KRIGING SWITCH FILE
# ============================================================

SWITCH="$ROOT/results/is/stage4_cross_eas/phase8_qc/kriging_fixed/BBJ_IS_L004.possible_switch.tsv"

if [ -s "$SWITCH" ]; then

    copy_one \
        "$SWITCH" \
        "BBJ_IS_L004.possible_switch.tsv"

else

    echo
    echo "[WARN] OPTIONAL SWITCH FILE MISSING"
    echo "$SWITCH"

fi


# ============================================================
# 5. REMOTE AUDIT
# ============================================================

echo
echo "===================================================="
echo "REMOTE FILE AUDIT"
echo "===================================================="

rclone lsl \
    "$REMOTE" \
    | sort


echo
echo "===================================================="
echo "REQUIRED FILE CHECK"
echo "===================================================="


python3 - <<'PY'
import subprocess

remote = "gdrive:MASTER_DEGREE/IS_COLAB/input"

required = [
    "BBJ_CREDIBLE_SET_MASTER.tsv",
    "PHASE9B_LOCUS_MASTER.tsv",
    "GRCH37_PROTEIN_CODING_GENES.tsv",

    "bbj_locus/BBJ_IS_L001.tsv",
    "bbj_locus/BBJ_IS_L002.tsv",
    "bbj_locus/BBJ_IS_L003.tsv",
    "bbj_locus/BBJ_IS_L004.tsv",

    "ld/BBJ_IS_L001.unphased.vcor1.bin",
    "ld/BBJ_IS_L001.unphased.vcor1.bin.vars",

    "ld/BBJ_IS_L002.unphased.vcor1.bin",
    "ld/BBJ_IS_L002.unphased.vcor1.bin.vars",

    "ld/BBJ_IS_L003.unphased.vcor1.bin",
    "ld/BBJ_IS_L003.unphased.vcor1.bin.vars",

    "ld/BBJ_IS_L004.unphased.vcor1.bin",
    "ld/BBJ_IS_L004.unphased.vcor1.bin.vars",
]

p = subprocess.run(
    [
        "rclone",
        "lsf",
        remote,
        "--recursive"
    ],
    capture_output=True,
    text=True
)

if p.returncode != 0:

    print(
        "REMOTE_LIST_FAILED"
    )

    print(
        p.stderr
    )

    raise SystemExit(1)


found = {
    x.strip()
    for x in p.stdout.splitlines()
    if x.strip()
}


missing = [
    x
    for x in required
    if x not in found
]


print(
    "REQUIRED=",
    len(required)
)

print(
    "FOUND=",
    len(required) - len(missing)
)

print(
    "MISSING=",
    len(missing)
)


for x in missing:
    print(
        "MISSING",
        x
    )


if missing:

    print(
        "EXPORT_VALIDATION=FAIL"
    )

    raise SystemExit(1)


print(
    "EXPORT_VALIDATION=PASS"
)
PY

VALIDATE_RC=$?


echo
echo "===================================================="
echo "EXPORT SUMMARY"
echo "===================================================="

echo "COPY_PASS=$PASS"
echo "COPY_FAIL=$FAIL"
echo "VALIDATE_RC=$VALIDATE_RC"


if [ "$FAIL" -eq 0 ] && \
   [ "$VALIDATE_RC" -eq 0 ]
then

    echo "EXPORT_STATUS=PASS"

else

    echo "EXPORT_STATUS=PARTIAL"

fi


echo "END=$(date)"
echo "===================================================="
