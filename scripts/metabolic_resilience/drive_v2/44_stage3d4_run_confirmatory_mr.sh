#!/usr/bin/env bash
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CODE="${CODE_DIR:-$SCRIPT_DIR}"
OUT="$ROOT/results/metabolic_resilience/stage3_confirmatory_mr"
HARM="$OUT/harmonized/STAGE3D2_HARMONIZED.tsv.gz"

cd "$ROOT" || {
  echo "[FAIL] cannot cd to $ROOT"
}

set +e
set +u
set +o pipefail 2>/dev/null || true

mkdir -p "$OUT"

echo "===================================================="
echo "STAGE 3-D : CONFIRMATORY MULTI-SNP MR"
echo "===================================================="
echo "CODE=$CODE"

READY=1

python3 "$CODE/35_stage3c4_apply_rsid_overlay.py"
RC_OVERLAY=$?
echo "C4 overlay rc=$RC_OVERLAY"
[ "$RC_OVERLAY" = "0" ] || READY=0

if [ "$READY" = "1" ]; then
  python3 "$CODE/40_stage3d0_discover_outcomes.py"
  RC0=$?
else
  RC0=90
fi
echo "D0 rc=$RC0"
[ "$RC0" = "0" ] || READY=0

if [ "$READY" = "1" ]; then
  python3 "$CODE/41_stage3d1_extract_outcome_instruments.py"
  RC1=$?
else
  RC1=90
fi
echo "D1 rc=$RC1"
[ "$RC1" = "0" ] || READY=0

if [ "$READY" = "1" ]; then
  python3 "$CODE/42_stage3d2_harmonize.py"
  RC2=$?
else
  RC2=90
fi
echo "D2 rc=$RC2"
[ "$RC2" = "0" ] || READY=0

if [ "$READY" = "1" ] && [ -s "$HARM" ]; then
  Rscript "$CODE/43_stage3d3_mr.R" "$HARM" "$OUT"
  RC3=$?
else
  RC3=90
fi
echo "D3 rc=$RC3"
[ "$RC3" = "0" ] || READY=0

if [ "$READY" = "1" ]; then
  python3 "$CODE/45_stage3d5_summarize_mr.py"
  RC5=$?
else
  RC5=90
fi
echo "D5 rc=$RC5"

echo "===================================================="
echo "STAGE 3-D END"
echo "===================================================="
echo "C4_OVERLAY_RC=$RC_OVERLAY"
echo "D0_RC=$RC0"
echo "D1_RC=$RC1"
echo "D2_RC=$RC2"
echo "D3_RC=$RC3"
echo "D5_RC=$RC5"

if [ "$RC_OVERLAY" = "0" ] && [ "$RC0" = "0" ] && [ "$RC1" = "0" ] && [ "$RC2" = "0" ] && [ "$RC3" = "0" ] && [ "$RC5" = "0" ]; then
  echo "[PASS] STAGE 3-D COMPLETE"
else
  echo "[STOP / REVIEW] STAGE 3-D incomplete"
fi


FINAL_RC=0

for RC in "$RC_OVERLAY" "$RC0" "$RC1" "$RC2" "$RC3" "$RC5"
do
  if [ "$RC" != "0" ]; then
    FINAL_RC=2
  fi
done

echo "FINAL_RC=$FINAL_RC"

# Fail closed so parent wrappers cannot mistake an incomplete pipeline for success.
if [ "$FINAL_RC" = "0" ]; then
  true
else
  false
fi
