#!/usr/bin/env bash
# Explicit selection and pinning only; does not install packages or log in.
set -u
code_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 1
cd "$code_root" || exit 1
python_bin="${MASTEROMICS_PYTHON:-python3}"
"$python_bin" -S -m masteromics resources collect "$@"
result=$?
exit "$result"
