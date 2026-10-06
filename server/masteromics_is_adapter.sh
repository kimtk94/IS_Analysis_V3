#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="/srv/is-analysis/MasterOmics"
OUT="${1:-/srv/is-analysis/results/masteromics/is_baseline}"
HASH_MODE="${MASTEROMICS_IS_HASH_LARGE:-0}"

cd "$REPO" || exit 1

ARGS=("$ROOT" "$OUT")
if [ "$HASH_MODE" = "1" ]; then
  ARGS+=("--hash-large")
fi

python3 -m masteromics is-adapter "${ARGS[@]}"
RC=$?

echo "MASTEROMICS_IS_ADAPTER_RC=$RC"
echo "OUTPUT=$OUT"

exit "$RC"
