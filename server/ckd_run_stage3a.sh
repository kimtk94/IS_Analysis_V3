#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAGE2C="${CKD_STAGE2C_ROOT:-/srv/is-analysis/results/ckd/stage2c_susie}"
RAW="${CKD_STAGE3A_RAW_ROOT:-/srv/is-analysis/data/ckd/stage3a_kidney}"
OUT="${CKD_STAGE3A_ROOT:-/srv/is-analysis/results/ckd/stage3a_kidney}"
PYTHON="${CKD_PYTHON:-$ROOT/.venv-ckd/bin/python}"

DEFAULT="$STAGE2C/susie_results/STAGE2C_SUSIE_DEFAULT.tsv"
test -x "$PYTHON"
test -f "$DEFAULT"
command -v curl >/dev/null || { echo "curl is required" >&2; exit 7; }

if [[ "${ACCEPT_SUSZTAK_TERMS:-0}" != "1" ]]; then
  cat >&2 <<'EOF'
CKD Stage 3A uses kidney QTL files linked from the Susztak Kidney Biobank.
Read the current agreement first:
  https://susztaklab.com/agree.php
If you accept it, rerun with ACCEPT_SUSZTAK_TERMS=1.
EOF
  exit 9
fi

mkdir -p "$RAW" "$OUT"

if ! "$PYTHON" -c 'import openpyxl' >/dev/null 2>&1; then
  "$PYTHON" -m pip install --disable-pip-version-check openpyxl
fi

"$PYTHON" "$ROOT/scripts/run_ckd_stage3a_kidney_evidence.py" \
  --stage2c-default "$DEFAULT" \
  --raw-root "$RAW" \
  --output-root "$OUT" \
  --accept-susztak-terms

(
  cd "$OUT"
  find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)

if [[ "${SYNC_DRIVE:-0}" == "1" ]]; then
  RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
  DEST="${RCLONE_REMOTE}:IS_Analysis_V3/results/ckd/stage3a_kidney"
  rclone copy "$OUT" "$DEST" \
    --checksum --transfers 1 --checkers 2 --progress \
    --filter '+ /STAGE3A_*.tsv' \
    --filter '+ /STAGE3A_*.json' \
    --filter '+ /SHA256SUMS.txt' \
    --filter '- *'
  echo "CKD_STAGE3A_DRIVE_SYNC_PASS"
fi

echo "CKD_STAGE3A_PASS"
echo "RESULT=$OUT"
echo "RAW=$RAW"
