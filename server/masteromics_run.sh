#!/usr/bin/env bash
# Explicit status handling; no set -e.
REPO_ROOT="${MASTEROMICS_CODE_ROOT:-/srv/is-analysis/IS_Analysis_V3}"
PYTHON_BIN="${MASTEROMICS_PYTHON:-python3}"
cd "$REPO_ROOT" || exit 2
if [ "$#" -eq 0 ]; then
  echo 'Usage: bash server/masteromics_run.sh projects/ckd.json projects/ischemic_stroke.json --registry projects/datasets.json --jobs 4'
  exit 2
fi
"$PYTHON_BIN" -m masteromics run "$@"
run_status=$?
if [ "$run_status" -ne 0 ]; then
  echo "MASTEROMICS_FAILED status=$run_status; inspect run_summary.json and logs"
  exit "$run_status"
fi
echo 'MASTEROMICS_SUCCESS'
