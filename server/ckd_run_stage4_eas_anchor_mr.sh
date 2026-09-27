#!/usr/bin/env bash

ROOT="${ROOT:-/srv/is-analysis}"
REPO="${REPO:-$ROOT/IS_Analysis_V3}"

ANCHOR="$ROOT/results/ckd/stage4_koges/STAGE4_KOGES_ANCHOR_PANEL.tsv"

EGFR="$ROOT/results/ckd/stage4_eas/anchor_screen/EAS_eGFR_ANCHOR_MATCH.tsv.gz"
BUN="$ROOT/results/ckd/stage4_eas/anchor_screen/EAS_BUN_ANCHOR_MATCH.tsv.gz"

OUT="$ROOT/results/ckd/stage4_eas/anchor_mr"

mkdir -p "$OUT"

echo "============================================================"
echo "CKD STAGE 4 EAS ANCHOR MR"
echo "============================================================"

for f in "$ANCHOR" "$EGFR" "$BUN"; do
    if [ -s "$f" ]; then
        echo "PASS $f"
    else
        echo "MISSING $f"
    fi
done

echo
python3 \
  "$REPO/scripts/run_ckd_stage4_eas_anchor_mr.py" \
  --anchor "$ANCHOR" \
  --egfr "$EGFR" \
  --bun "$BUN" \
  --outdir "$OUT"

RC=$?

echo
echo "===== OUTPUT ====="

find "$OUT" \
  -maxdepth 1 \
  -type f \
  -printf '%f\t%s bytes\n' \
  | sort

echo
echo "===== INTEGRATED EVIDENCE ====="

if [ -s "$OUT/EAS_ANCHOR_INTEGRATED_EVIDENCE.tsv" ]; then
    column -t -s $'\t' \
      "$OUT/EAS_ANCHOR_INTEGRATED_EVIDENCE.tsv"
fi

echo
echo "===== SUMMARY ====="

if [ -s "$OUT/EAS_ANCHOR_MR_SUMMARY.json" ]; then
    cat "$OUT/EAS_ANCHOR_MR_SUMMARY.json"
fi

exit "$RC"
