#!/usr/bin/env bash

ROOT="${ROOT:-/srv/is-analysis/IS_Analysis_V3}"
OUT="$ROOT/results/skin_beauty/stage0_audit"

mkdir -p "$OUT"

echo "===== SKIN BEAUTY STAGE 0A: EXPOSURE AUDIT ====="
python3 "$ROOT/scripts/skin_beauty/prepare_skin_stage0.py" \
  --root "$ROOT" \
  --search-root "/srv/is-analysis" \
  --data-dir "$ROOT/data/skin_beauty/stage0_gwas" \
  --out-dir "$OUT"
AUDIT_STATUS=$?

echo
echo "===== SKIN BEAUTY STAGE 0B: OPENGWAS API ====="
python3 "$ROOT/scripts/skin_beauty/query_skin_wrinkle_opengwas.py"
API_STATUS=$?

echo
echo "===== STAGE 0 STATUS ====="
echo "python_audit_status=$AUDIT_STATUS"
echo "opengwas_api_status=$API_STATUS"
echo "result_dir=$OUT"

echo
echo "===== OUTPUTS ====="
ls -1 "$OUT" 2>/dev/null

if [ "$AUDIT_STATUS" -eq 0 ] && [ "$API_STATUS" -eq 0 ]; then
  echo
  echo "STAGE0_GATE=PASS"
else
  echo
  echo "STAGE0_GATE=REVIEW"
fi

exit 0
