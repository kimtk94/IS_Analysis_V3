#!/usr/bin/env bash
set +e
set +u
set +o pipefail 2>/dev/null || true
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
BASE="$ROOT/results/metabolic_resilience"
OUT="$BASE/reproducibility"
mkdir -p "$OUT"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
find "$BASE" -type f ! -path "$OUT/*" -print0 2>/dev/null | sort -z | xargs -0 sha256sum > "$OUT/RESULT_SHA256_$STAMP.txt"
echo "[PASS] checksums -> $OUT/RESULT_SHA256_$STAMP.txt"
