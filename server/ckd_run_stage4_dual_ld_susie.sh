#!/usr/bin/env bash

ROOT="${ROOT:-/srv/is-analysis}"
REPO="${REPO:-$ROOT/IS_Analysis_V3}"

ANCHOR="$ROOT/results/ckd/stage4_koges/STAGE4_KOGES_ANCHOR_PANEL.tsv"

REGIONAL="$ROOT/results/ckd/stage4_eas/regional"

EGFR="$REGIONAL/eGFR/coloc_input"
BUN="$REGIONAL/BUN/coloc_input"

WORK="$ROOT/data/ckd/stage4_eas_dual_ld"

OUT="$ROOT/results/ckd/stage4_eas/dual_ld_susie"

R_LIB="${CKD_R_LIB:-/srv/is-analysis/.Rlib}"
PLINK2="${PLINK2:-plink2}"

mkdir -p \
  "$WORK" \
  "$OUT" \
  "$R_LIB"

echo "============================================================"
echo "CKD STAGE 4 DUAL-LD SUSIE"
echo "============================================================"

echo
echo "===== INPUT CHECK ====="

for f in "$ANCHOR"; do
  if [ -s "$f" ]; then
    echo "PASS $f"
  else
    echo "MISSING $f"
    exit 2
  fi
done

for d in "$EGFR" "$BUN"; do
  if [ -d "$d" ]; then
    echo "PASS $d"
  else
    echo "MISSING $d"
    exit 2
  fi
done

for exe in python3 bcftools "$PLINK2" Rscript curl; do
  if command -v "$exe" >/dev/null 2>&1; then
    echo "PASS executable: $exe"
  else
    echo "MISSING executable: $exe"
    exit 3
  fi
done

echo
echo "===== PREPARE DUAL LD ====="

python3 \
  "$REPO/scripts/prepare_ckd_stage4_dual_ld.py" \
  --anchor "$ANCHOR" \
  --egfr-input "$EGFR" \
  --bun-input "$BUN" \
  --work-root "$WORK" \
  --output-root "$OUT" \
  --plink2 "$PLINK2" \
  --window-bp 500000 \
  --min-reference-maf 0.01 \
  --min-ld-snps 100 \
  --threads "${DUAL_LD_THREADS:-2}" \
  --memory-mb "${DUAL_LD_MEMORY_MB:-3500}"

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "DUAL_LD_PREP_FAILED rc=$RC"
  exit "$RC"
fi

echo
echo "===== DUAL LD QC ====="

column -t -s $'\t' \
  "$OUT/DUAL_LD_QC.tsv"

echo
echo "===== R PACKAGE CHECK ====="

export R_LIBS_USER="$R_LIB"

Rscript -e '
pkgs <- c("coloc","susieR")
miss <- pkgs[
  !vapply(
    pkgs,
    requireNamespace,
    logical(1),
    quietly=TRUE
  )
]
if (length(miss)) {
  install.packages(
    miss,
    repos="https://cloud.r-project.org",
    lib=Sys.getenv("R_LIBS_USER")
  )
}
'

RC=$?

if [ "$RC" -ne 0 ]; then
  echo "R_PACKAGE_SETUP_FAILED rc=$RC"
  exit "$RC"
fi

echo
echo "===== RUN DUAL-LD SUSIE ====="

mkdir -p \
  "$OUT/susie_results"

Rscript \
  "$REPO/scripts/run_ckd_stage4_dual_ld_susie.R" \
  "$OUT/susie_input" \
  "$OUT/ld" \
  "$REGIONAL" \
  "$OUT/susie_results"

RC=$?

echo
echo "===== DEFAULT RESULTS ====="

if [ -s "$OUT/susie_results/DUAL_LD_SUSIE_DEFAULT.tsv" ]; then
  column -t -s $'\t' \
    "$OUT/susie_results/DUAL_LD_SUSIE_DEFAULT.tsv"
else
  echo "DEFAULT RESULT NOT FOUND"
fi

echo
echo "===== FAILURES ====="

if [ -s "$OUT/susie_results/DUAL_LD_SUSIE_FAILURES.tsv" ]; then
  column -t -s $'\t' \
    "$OUT/susie_results/DUAL_LD_SUSIE_FAILURES.tsv"
else
  echo "NONE"
fi

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
echo "DUAL_LD_SUSIE_RC=$RC"

if [ "$RC" -eq 0 ]; then
  echo "CKD_STAGE4_DUAL_LD_PIPELINE_PASS"
fi

exit "$RC"
