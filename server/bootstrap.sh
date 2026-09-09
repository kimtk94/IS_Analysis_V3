#!/usr/bin/env bash
set -euo pipefail

ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
CODE_ROOT="${CODE_ROOT:-$ROOT/app}"
OWNER="${SUDO_USER:-${USER}}"

if [[ $EUID -ne 0 ]]; then
  echo "Run with sudo: sudo bash server/bootstrap.sh" >&2
  exit 1
fi

mkdir -p \
  "$ROOT/data/staging" \
  "$ROOT/data/rawdata" \
  "$ROOT/data/reference" \
  "$ROOT/data/standardized/exposure" \
  "$ROOT/data/standardized/outcome" \
  "$ROOT/results/instruments" \
  "$ROOT/results/harmonized" \
  "$ROOT/results/mr" \
  "$ROOT/results/coloc" \
  "$ROOT/results/replication" \
  "$ROOT/results/sc_context" \
  "$ROOT/reports/latest" \
  "$ROOT/reports/archive" \
  "$ROOT/manifests" \
  "$ROOT/logs" \
  "$ROOT/state" \
  "$ROOT/backup"

chown -R "$OWNER":"$OWNER" "$ROOT"
chmod 750 "$ROOT"

if [[ ! -e "$CODE_ROOT" ]]; then
  ln -s "$(pwd)" "$CODE_ROOT"
fi

python3 scripts/storage_guard.py \
  --path "$ROOT" \
  --staging "$ROOT/data/staging" \
  --min-free-gib "${MIN_FREE_GIB:-40}" \
  --max-staging-gib "${MAX_STAGING_GIB:-45}"

echo
echo "Bootstrap complete: $ROOT"
echo "Next: copy reviewed V2/V3 metadata/reference files, then run the smoke test."
