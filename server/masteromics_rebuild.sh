#!/usr/bin/env bash
# Explicit exit handling; no set -e and no changes to legacy result roots.
set -u
code_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 1
cd "$code_root" || exit 1
python_bin="${MASTEROMICS_PYTHON:-$code_root/.venv-masteromics/bin/python}"
if [[ -n "${CKD_R_LIB:-}" ]]; then export R_LIBS_USER="$CKD_R_LIB"; fi
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
"$python_bin" -m masteromics.rebuild "$@"
result=$?
exit "$result"
