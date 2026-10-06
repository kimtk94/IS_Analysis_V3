#!/usr/bin/env bash

REPO="/srv/is-analysis/MasterOmics"
PROJECT="${MASTEROMICS_IS_BLUEPRINT:-/srv/is-analysis/masteromics_workspace/ischemic_stroke/project.json}"
MIGRATION="${MASTEROMICS_IS_MIGRATION:-/srv/is-analysis/results/masteromics/is_stage_migration}"
OUT="${1:-/srv/is-analysis/masteromics_workspace/ischemic_stroke/config/project.frozen.json}"

cd "$REPO" || exit 1

python3 -m masteromics is-bind bind "$PROJECT" "$MIGRATION" "$OUT"
RC=$?

echo "MASTEROMICS_IS_BIND_RC=$RC"
echo "PROJECT=$PROJECT"
echo "MIGRATION=$MIGRATION"
echo "OUTPUT=$OUT"

exit "$RC"
