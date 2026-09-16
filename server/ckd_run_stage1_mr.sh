#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
READY="${CKD_READY_ROOT:-/srv/is-analysis/data/ckd/analysis_ready}"
RESULT="${CKD_STAGE1_ROOT:-/srv/is-analysis/results/ckd/stage1}"
F_THRESHOLD="${CKD_F_THRESHOLD:-10}"
PAL_MAF="${CKD_PAL_MAF_THRESHOLD:-0.42}"
FREQ_TOL="${CKD_FREQ_TOLERANCE:-0.10}"

if [[ -x "$ROOT/.venv-ckd/bin/python" ]]; then
  PYTHON="${CKD_PYTHON:-$ROOT/.venv-ckd/bin/python}"
else
  PYTHON="${CKD_PYTHON:-python3}"
fi

test -f "$READY/EUR/eGFRcrea.tsv.gz"
test -f "$READY/EAS/eGFRcrea.tsv.gz"

mkdir -p "$RESULT"

echo "[stage1] ready=$READY"
echo "[stage1] result=$RESULT"
echo "[stage1] python=$PYTHON"

"$PYTHON" "$ROOT/scripts/run_ckd_stage1_mr.py" \
  --ready-root "$READY" \
  --output-root "$RESULT" \
  --f-threshold "$F_THRESHOLD" \
  --palindrome-maf-threshold "$PAL_MAF" \
  --freq-tolerance "$FREQ_TOL"

(
  cd "$RESULT"
  find . -type f ! -name 'SHA256SUMS.txt' -print0 \
    | sort -z \
    | xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)

echo "CKD_STAGE1_LOCAL_PASS"
echo "RESULT=$RESULT"

if [[ "${SYNC_DRIVE:-0}" == "1" ]]; then
  RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
  DRIVE_RESULTS_BASE="${DRIVE_RESULTS_BASE:-IS_Analysis_V3/results/ckd/stage1}"
  command -v rclone >/dev/null
  rclone copy "$RESULT" "${RCLONE_REMOTE}:${DRIVE_RESULTS_BASE}" \
    --create-empty-src-dirs \
    --transfers 4 \
    --checkers 8 \
    --checksum \
    --progress
  echo "CKD_STAGE1_DRIVE_SYNC_PASS"
  echo "REMOTE=${RCLONE_REMOTE}:${DRIVE_RESULTS_BASE}"
fi
