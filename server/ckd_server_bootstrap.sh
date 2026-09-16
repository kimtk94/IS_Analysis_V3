#!/usr/bin/env bash
set -euo pipefail

# Bootstrap CKD thesis repository on a Linux server.
# Code lives locally; large data are staged locally then synced to Google Drive.

REPO_URL="${REPO_URL:-https://github.com/kimtk94/IS_Analysis_V3.git}"
APP_ROOT="${APP_ROOT:-/srv/is-analysis}"
REPO_DIR="${REPO_DIR:-$APP_ROOT/IS_Analysis_V3}"
SCRATCH_ROOT="${SCRATCH_ROOT:-$APP_ROOT/data/ckd}"

if ! command -v git >/dev/null 2>&1; then
  echo "git is required" >&2
  exit 2
fi

sudo mkdir -p "$APP_ROOT"
sudo chown -R "${USER}:${USER}" "$APP_ROOT"

if [[ -d "$REPO_DIR/.git" ]]; then
  echo "[git] updating $REPO_DIR"
  git -C "$REPO_DIR" fetch --prune origin
  git -C "$REPO_DIR" checkout main
  git -C "$REPO_DIR" pull --ff-only origin main
else
  echo "[git] cloning $REPO_URL -> $REPO_DIR"
  git clone --branch main "$REPO_URL" "$REPO_DIR"
fi

mkdir -p "$SCRATCH_ROOT/rawdata" "$SCRATCH_ROOT/analysis_ready"

echo "SERVER_BOOTSTRAP_PASS"
echo "REPO_DIR=$REPO_DIR"
echo "SCRATCH_ROOT=$SCRATCH_ROOT"
echo
echo "Next:"
echo "  cd $REPO_DIR"
echo "  CKD_DATA_ROOT=$SCRATCH_ROOT/rawdata CKD_READY_ROOT=$SCRATCH_ROOT/analysis_ready server/ckd_acquire_public_data.sh"
