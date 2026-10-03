#!/usr/bin/env bash
# No set -e: capture failures explicitly and retain the comparison report.
set -u
set -o pipefail
code_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 1
cd "$code_root" || exit 1
python_bin="${MASTEROMICS_PYTHON:-python3}"
if [[ -z "${MASTEROMICS_PYTHON:-}" && -x /srv/is-analysis/.venv-ckd/bin/python ]]; then
  python_bin=/srv/is-analysis/.venv-ckd/bin/python
fi
"$python_bin" -c 'import numpy,pandas,scipy' || { echo 'Python numpy/pandas/scipy environment required; set MASTEROMICS_PYTHON.' >&2; exit 2; }
# The adapter records R dependency failures per locus; it never installs or
# modifies the existing scientific environment automatically.
output="${MASTEROMICS_REGRESSION_OUT:-/srv/is-analysis/results/masteromics/regression/ckd_$(date -u +%Y%m%dT%H%M%SZ)_$$}"
"$python_bin" -m masteromics.regression --out "$output" "$@"
result=$?
echo "Regression exit: $result; report: $output/REGRESSION_SUMMARY.json"
exit "$result"
