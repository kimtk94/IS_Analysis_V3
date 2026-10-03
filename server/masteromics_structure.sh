#!/usr/bin/env bash
# Creates structure only; never reads production datasets or runs analysis.
set -u
code_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" || exit 1
cd "$code_root" || exit 1
python_bin="${MASTEROMICS_PYTHON:-python3}"
"$python_bin" -m masteromics blueprint "$@"
result=$?
exit "$result"
