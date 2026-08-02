#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -m unittest discover -s tests -p 'test_*.py' -v
python3 -m py_compile scripts/ukb_ppp_batch_manifest_runner_fast.py workflow/*.py workflow/annotation/*.py
