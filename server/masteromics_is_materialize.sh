#!/usr/bin/env bash

REPO="/srv/is-analysis/MasterOmics"
MIGRATION="${MASTEROMICS_IS_MIGRATION:-/srv/is-analysis/results/masteromics/is_stage_migration}"
OUT="${1:-/srv/is-analysis/results/masteromics/is_frozen_prefix_verified}"

cd "$REPO" || exit 1

python3 -m masteromics is-bind materialize "$MIGRATION" "$OUT"
RC=$?

echo "MASTEROMICS_IS_MATERIALIZE_RC=$RC"
echo "MIGRATION=$MIGRATION"
echo "OUTPUT=$OUT"

exit "$RC"
