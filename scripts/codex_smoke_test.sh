#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
case "${1:-}" in
  --fixture-only) suite=fixtures ;;
  --with-legacy) suite=all ;;
  "")
    # Existing source-audit workflows prepare this pinned checkout. Preserve
    # their full integration gate while allowing a fresh clone to run offline.
    if [[ -d legacy ]]; then suite=all; else suite=fixtures; fi
    ;;
  *) echo "Usage: $0 [--fixture-only|--with-legacy]" >&2; exit 2 ;;
esac
if [[ $# -gt 1 ]]; then
  echo "Usage: $0 [--fixture-only|--with-legacy]" >&2
  exit 2
fi
python3 scripts/run_repository_tests.py --suite "$suite"
python3 -m py_compile scripts/ukb_ppp_batch_manifest_runner_fast.py workflow/*.py workflow/annotation/*.py
