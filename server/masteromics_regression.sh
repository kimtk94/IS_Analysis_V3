#!/usr/bin/env bash
# No set -e: capture failures explicitly and retain the comparison report.
set -u
set -o pipefail
code_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 1
cd "$code_root" || exit 1
python_bin="${MASTEROMICS_PYTHON:-python3}"
if [[ -z "${MASTEROMICS_PYTHON:-}" ]]; then
  for candidate in "$code_root/.venv-ckd/bin/python" /srv/is-analysis/IS_Analysis_V3/.venv-ckd/bin/python /srv/is-analysis/.venv-ckd/bin/python; do
    if [[ -x "$candidate" ]]; then python_bin="$candidate"; break; fi
  done
fi
if [[ -n "${CKD_R_LIB:-}" ]]; then
  export R_LIBS_USER="$CKD_R_LIB"
elif [[ -z "${R_LIBS_USER:-}" && -d /srv/is-analysis/.Rlib ]]; then
  export R_LIBS_USER=/srv/is-analysis/.Rlib
fi
"$python_bin" -c 'import numpy,pandas,scipy' || { echo 'Python numpy/pandas/scipy environment required; set MASTEROMICS_PYTHON.' >&2; exit 2; }
# The adapter records R dependency failures per locus; it never installs or
# modifies the existing scientific environment automatically.
output="${MASTEROMICS_REGRESSION_OUT:-/srv/is-analysis/results/masteromics/regression/ckd_$(date -u +%Y%m%dT%H%M%SZ)_$$}"
"$python_bin" -m masteromics.regression --out "$output" "$@"
result=$?
echo "Regression exit: $result; report: $output/REGRESSION_SUMMARY.json"
exit "$result"
