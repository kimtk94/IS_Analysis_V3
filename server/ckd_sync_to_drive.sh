#!/usr/bin/env bash
set -euo pipefail

# Sync verified CKD public data from server-local storage to Google Drive.
# Requires rclone and an existing Google Drive remote (default: gdrive).
#
# Example:
#   RCLONE_REMOTE=gdrive \
#   LOCAL_CKD_ROOT=/srv/is-analysis/data/ckd \
#   DRIVE_BASE='IS_Analysis_V3/data' \
#   server/ckd_sync_to_drive.sh

RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
LOCAL_CKD_ROOT="${LOCAL_CKD_ROOT:-/srv/is-analysis/data/ckd}"
DRIVE_BASE="${DRIVE_BASE:-IS_Analysis_V3/data}"
RAW_LOCAL="${RAW_LOCAL:-$LOCAL_CKD_ROOT/rawdata}"
READY_LOCAL="${READY_LOCAL:-$LOCAL_CKD_ROOT/analysis_ready}"
RAW_REMOTE="${RCLONE_REMOTE}:${DRIVE_BASE}/rawdata/ckd"
READY_REMOTE="${RCLONE_REMOTE}:${DRIVE_BASE}/analysis_ready/ckd"

command -v rclone >/dev/null 2>&1 || {
  echo "rclone is required. Install/configure it first." >&2
  exit 2
}

if ! rclone listremotes | grep -qx "${RCLONE_REMOTE}:"; then
  echo "rclone remote '${RCLONE_REMOTE}' is not configured." >&2
  echo "Run: rclone config" >&2
  exit 3
fi

test -d "$RAW_LOCAL"
test -d "$READY_LOCAL"
test -f "$READY_LOCAL/SHA256SUMS.txt"

echo "[verify] local SHA256"
(
  cd "$READY_LOCAL"
  sha256sum -c SHA256SUMS.txt
)

echo "[sync] raw public data -> $RAW_REMOTE"
rclone copy "$RAW_LOCAL" "$RAW_REMOTE" \
  --create-empty-src-dirs \
  --transfers 4 \
  --checkers 8 \
  --checksum \
  --progress

echo "[sync] analysis-ready -> $READY_REMOTE"
rclone copy "$READY_LOCAL" "$READY_REMOTE" \
  --create-empty-src-dirs \
  --transfers 4 \
  --checkers 8 \
  --checksum \
  --progress

echo "[verify] remote analysis-ready listing"
rclone lsf "$READY_REMOTE" --recursive

echo "CKD_DRIVE_SYNC_PASS"
echo "RAW_REMOTE=$RAW_REMOTE"
echo "READY_REMOTE=$READY_REMOTE"
