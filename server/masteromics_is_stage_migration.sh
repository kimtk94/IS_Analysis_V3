#!/usr/bin/env bash

ROOT="/srv/is-analysis"
REPO="/srv/is-analysis/MasterOmics"
BASELINE="${MASTEROMICS_IS_BASELINE:-/srv/is-analysis/results/masteromics/is_baseline}"
OUT="${1:-/srv/is-analysis/results/masteromics/is_stage_migration}"

cd "$REPO" || exit 1

python3 -m masteromics is-stage-migration "$ROOT" "$BASELINE" "$OUT"
RC=$?

echo "MASTEROMICS_IS_STAGE_MIGRATION_RC=$RC"
echo "BASELINE=$BASELINE"
echo "OUTPUT=$OUT"

exit "$RC"
