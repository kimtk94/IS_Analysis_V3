#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
cd "$ROOT" || exit 1

set +e
set +u
set +o pipefail 2>/dev/null || true

CODE="${CODE_DIR:-$ROOT/code/metabolic_resilience}"
OUT="$ROOT/results/metabolic_resilience/stage3_confirmatory_mr"
HARM="$OUT/harmonized/STAGE3D2_HARMONIZED.tsv.gz"

mkdir -p "$OUT"

echo "===================================================="
echo "STAGE 3-D : CONFIRMATORY MULTI-SNP MR"
echo "===================================================="

python3 "$CODE/40_stage3d0_discover_outcomes.py"
RC0=$?
echo "D0 rc=$RC0"

if [ "$RC0" != "0" ]; then
  echo "[WARN] review outcome registry before continuing"
fi

python3 "$CODE/41_stage3d1_extract_outcome_instruments.py"
RC1=$?
echo "D1 rc=$RC1"

python3 "$CODE/42_stage3d2_harmonize.py"
RC2=$?
echo "D2 rc=$RC2"

if [ "$RC2" = "0" ] && [ -s "$HARM" ]; then
  Rscript "$CODE/43_stage3d3_mr.R" "$HARM" "$OUT"
  RC3=$?
else
  RC3=99
fi

echo "D3 rc=$RC3"

if [ "$RC3" = "0" ]; then
  python3 "$CODE/45_stage3d5_summarize_mr.py"
  RC4=$?
else
  RC4=99
fi

echo "D5 rc=$RC4"

echo "===================================================="
echo "STAGE 3-D END"
echo "===================================================="
