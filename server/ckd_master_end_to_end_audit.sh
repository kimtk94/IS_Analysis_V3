#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
MAP="${CKD_MASTER_STAGE_MAP:-$ROOT/config/ckd/master_stage_map.tsv}"
OUTROOT="${CKD_MASTER_ROOT:-/srv/is-analysis/results/ckd/master}"
PYTHON="${CKD_PYTHON:-python3}"

mkdir -p "$OUTROOT"

echo "============================================================"
echo " CKD MASTER — EXISTING OUTPUT AUDIT / END-TO-END PLAN"
echo "============================================================"
echo "ROOT=$ROOT"
echo "MAP=$MAP"
echo "OUTROOT=$OUTROOT"

"$PYTHON" "$ROOT/scripts/audit_ckd_master_existing_outputs.py" \
  --map "$MAP" \
  --output "$OUTROOT/CKD_MASTER_STAGE_READINESS.tsv" \
  --summary "$OUTROOT/CKD_MASTER_STAGE_READINESS.json"

echo
echo "===== READINESS ====="
column -ts $'\t' "$OUTROOT/CKD_MASTER_STAGE_READINESS.tsv" 2>/dev/null || \
  cat "$OUTROOT/CKD_MASTER_STAGE_READINESS.tsv"

echo
echo "===== GENERIC MASTER DRY RUN ====="
"$PYTHON" "$ROOT/workflow/run_master.py" --disease ckd --stages 0-19 --dry-run

echo
echo "IMPORTANT:"
echo "This wrapper is read-only. It does not execute missing Stage 4/7-12/16-19 analyses."
echo "Use CKD_MASTER_STAGE_READINESS.tsv to bridge existing outputs first, then run only stages with explicit validated inputs."
echo "CKD_MASTER_AUDIT_PASS"
