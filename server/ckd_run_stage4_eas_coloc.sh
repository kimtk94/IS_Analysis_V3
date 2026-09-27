#!/usr/bin/env bash

ROOT="${ROOT:-/srv/is-analysis}"
REPO="${REPO:-$ROOT/IS_Analysis_V3}"

ANCHOR="$ROOT/results/ckd/stage4_koges/STAGE4_KOGES_ANCHOR_PANEL.tsv"

PQTL="$ROOT/results/ckd/stage2b_coloc/pqtl_locus"

EGFR="$ROOT/data/ckd/eas/outcome/EAS/eas_chen2024_egfr_meta/TWB2_BBJ_eGFR_hg19_METAL_FUMA_noNA.gz"

BUN="$ROOT/data/ckd/eas/support/EAS/eas_chen2024_bun_meta/TWB2_BBJ_BUN_hg19_METAL_FUMA_noNA.gz"

OUT="$ROOT/results/ckd/stage4_eas/regional"

R_LIB="${CKD_R_LIB:-/srv/is-analysis/.Rlib}"

mkdir -p \
  "$OUT/eGFR" \
  "$OUT/BUN" \
  "$R_LIB"

echo "============================================================"
echo "CKD STAGE 4 EAS REGIONAL COLOCALIZATION"
echo "============================================================"

echo
echo "===== INPUT ====="

for f in \
  "$ANCHOR" \
  "$EGFR" \
  "$BUN"
do
  if [ -s "$f" ]; then
    echo "PASS $f"
  else
    echo "MISSING $f"
    exit 2
  fi
done

if [ ! -d "$PQTL" ]; then
  echo "MISSING $PQTL"
  exit 2
fi

echo
echo "===== EAS eGFR REGIONAL PREP ====="

python3 \
  "$REPO/scripts/prepare_ckd_stage4_eas_regional.py" \
  --anchor "$ANCHOR" \
  --pqtl-locus-dir "$PQTL" \
  --outcome "$EGFR" \
  --phenotype eGFRcrea \
  --outcome-n 244952 \
  --output-root "$OUT/eGFR"

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "EAS eGFR preparation failed: $RC"
  exit "$RC"
fi

echo
echo "===== EAS BUN REGIONAL PREP ====="

python3 \
  "$REPO/scripts/prepare_ckd_stage4_eas_regional.py" \
  --anchor "$ANCHOR" \
  --pqtl-locus-dir "$PQTL" \
  --outcome "$BUN" \
  --phenotype BUN \
  --outcome-n 241112 \
  --output-root "$OUT/BUN"

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "EAS BUN preparation failed: $RC"
  exit "$RC"
fi

echo
echo "===== INSTALL/VERIFY COLOC ====="

export R_LIBS_USER="$R_LIB"

Rscript -e '
if (!requireNamespace("coloc", quietly=TRUE)) {
  install.packages(
    "coloc",
    repos="https://cloud.r-project.org",
    lib=Sys.getenv("R_LIBS_USER")
  )
}
'

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "R coloc package setup failed: $RC"
  exit "$RC"
fi

echo
echo "===== EAS eGFR COLOC ====="

mkdir -p "$OUT/eGFR/coloc_results"

Rscript \
  "$REPO/scripts/run_ckd_stage4_eas_coloc.R" \
  "$OUT/eGFR/coloc_input" \
  eGFRcrea \
  "$OUT/eGFR/coloc_results"

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "EAS eGFR coloc failed: $RC"
  exit "$RC"
fi

echo
echo "===== EAS BUN COLOC ====="

mkdir -p "$OUT/BUN/coloc_results"

Rscript \
  "$REPO/scripts/run_ckd_stage4_eas_coloc.R" \
  "$OUT/BUN/coloc_input" \
  BUN \
  "$OUT/BUN/coloc_results"

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "EAS BUN coloc failed: $RC"
  exit "$RC"
fi

echo
echo "===== eGFR DEFAULT ====="

column -t -s $'\t' \
  "$OUT/eGFR/coloc_results/EAS_COLOC_DEFAULT.tsv"

echo
echo "===== BUN DEFAULT ====="

column -t -s $'\t' \
  "$OUT/BUN/coloc_results/EAS_COLOC_DEFAULT.tsv"

echo
echo "===== SHA256 ====="

find "$OUT" \
  -type f \
  ! -name SHA256SUMS.txt \
  -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > "$OUT/SHA256SUMS.txt"

echo
echo "CKD_STAGE4_EAS_REGIONAL_COLOC_PASS"
echo "RESULT=$OUT"
