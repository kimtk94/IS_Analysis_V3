#!/usr/bin/env bash

ROOT="${ROOT:-/srv/is-analysis/IS_Analysis_V3}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$ROOT/results/skin_beauty/stage1_mr}"
EXPOSURE="${EXPOSURE:-/srv/is-analysis/results/metabolic_resilience/stage1_pqtl/instruments/UKBPPP_CIS_INSTRUMENTS_PRIMARY_NON_MHC.tsv.gz}"
MODE="${1:-preflight}"

mkdir -p "$OUT"

if [ -f "$HOME/.config/is-analysis/opengwas.env" ]; then
  source "$HOME/.config/is-analysis/opengwas.env"
fi

echo "===== SKIN BEAUTY STAGE 1 ====="
echo "root     = $ROOT"
echo "exposure = $EXPOSURE"
echo "out      = $OUT"
echo "mode     = $MODE"

if [ "$MODE" = "preflight" ]; then
  TARGET="$OUT/preflight"
  python3 "$SCRIPT_DIR/run_skin_stage1_mr.py" \
    --root "$ROOT" \
    --exposure "$EXPOSURE" \
    --out-dir "$TARGET" \
    --exposure-only
  STATUS=$?
elif [ "$MODE" = "pilot" ]; then
  TARGET="$OUT/pilot"
  if [ -z "${OPENGWAS_JWT:-}" ]; then
    echo "ERROR: OPENGWAS_JWT missing"
    exit 2
  fi
  python3 "$SCRIPT_DIR/run_skin_stage1_mr.py" \
    --root "$ROOT" \
    --exposure "$EXPOSURE" \
    --out-dir "$TARGET" \
    --pilot-variants 512 \
    --batch-size 32
  STATUS=$?
elif [ "$MODE" = "full" ]; then
  TARGET="$OUT/full"
  if [ -z "${OPENGWAS_JWT:-}" ]; then
    echo "ERROR: OPENGWAS_JWT missing"
    exit 2
  fi
  python3 "$SCRIPT_DIR/run_skin_stage1_mr.py" \
    --root "$ROOT" \
    --exposure "$EXPOSURE" \
    --out-dir "$TARGET" \
    --batch-size 32
  STATUS=$?
else
  echo "ERROR: mode must be one of: preflight, pilot, full"
  exit 2
fi

echo
echo "===== STAGE 1 STATUS ====="
echo "status=$STATUS"
echo "result_dir=$TARGET"

echo
echo "===== EXPOSURE QC ====="
if [ -f "$TARGET/STAGE1_EXPOSURE_QC.json" ]; then
  cat "$TARGET/STAGE1_EXPOSURE_QC.json"
else
  echo "STAGE1_EXPOSURE_QC.json not found"
fi

if [ "$MODE" != "preflight" ]; then
  echo
  echo "===== STAGE 1 SUMMARY ====="
  if [ -f "$TARGET/STAGE1_SUMMARY.json" ]; then
    cat "$TARGET/STAGE1_SUMMARY.json"
  else
    echo "STAGE1_SUMMARY.json not found"
  fi
fi

exit "$STATUS"
