#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAGE1="${CKD_STAGE1_ROOT:-/srv/is-analysis/results/ckd/stage1}"
STAGE2="${CKD_STAGE2_ROOT:-/srv/is-analysis/results/ckd/stage2}"
PQTL="${CKD_STAGE2_PQTL_ROOT:-/srv/is-analysis/data/ckd/stage2_pqtl}"
MANIFEST="${UKB_PPP_MANIFEST:-$ROOT/data/metadata/ukb_ppp_download_manifest.tsv}"
PYTHON="${CKD_PYTHON:-$ROOT/.venv-ckd/bin/python}"

test -x "$PYTHON"
test -f "$STAGE1/stage1_protein_summary.tsv"
test -f "$STAGE1/instrument_qc.tsv.gz"
test -f "$MANIFEST"
mkdir -p "$STAGE2" "$PQTL"

"$PYTHON" "$ROOT/scripts/run_ckd_stage2_candidates.py" \
  --stage1-root "$STAGE1" \
  --download-manifest "$MANIFEST" \
  --output-root "$STAGE2"

echo "CKD_STAGE2_PLAN_LOCAL_PASS"

if [[ "${STAGE2_DOWNLOAD:-0}" == "1" ]]; then
  if ! "$PYTHON" -c 'import synapseclient' 2>/dev/null; then
    echo "[env] installing synapseclient into CKD venv"
    "$PYTHON" -m pip install --disable-pip-version-check synapseclient
  fi
  export PATH="$(dirname "$PYTHON"):$PATH"
  if ! "$PYTHON" -c 'import synapseclient; synapseclient.login(silent=True)' >/dev/null 2>&1; then
    cat >&2 <<'EOF'
Synapse authentication is not available to non-interactive Python/CLI calls.
Configure one persistent method, then rerun:
  1) ~/.synapseConfig via: synapse config
  2) or environment variable: SYNAPSE_AUTH_TOKEN
Verification:
  python -c 'import synapseclient; synapseclient.login()'
EOF
    exit 4
  fi

  "$PYTHON" "$ROOT/scripts/download_ckd_stage2_pqtl.py" \
    --manifest "$STAGE2/stage2_pqtl_download_manifest.tsv" \
    --dest "$PQTL" \
    --ancestry EUR \
    --progress "$STAGE2/stage2_pqtl_download_progress.tsv"

  "$PYTHON" "$ROOT/scripts/inspect_ckd_stage2_pqtl.py" \
    --root "$PQTL" \
    --output "$STAGE2/stage2_pqtl_schema.tsv"
fi

(
  cd "$STAGE2"
  find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)

if [[ "${SYNC_DRIVE:-0}" == "1" ]]; then
  RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
  rclone copy "$STAGE2" "${RCLONE_REMOTE}:IS_Analysis_V3/results/ckd/stage2" --checksum --progress
  if [[ "${STAGE2_DOWNLOAD:-0}" == "1" ]]; then
    rclone copy "$PQTL" "${RCLONE_REMOTE}:IS_Analysis_V3/data/rawdata/ckd/stage2_pqtl" --checksum --progress
  fi
  echo "CKD_STAGE2_DRIVE_SYNC_PASS"
fi

echo "CKD_STAGE2_PASS"
echo "RESULT=$STAGE2"
echo "PQTL=$PQTL"
