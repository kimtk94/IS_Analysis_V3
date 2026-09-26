#!/usr/bin/env bash
cd "${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
set +e
set +u
set +o pipefail 2>/dev/null || true
ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
OUT="$ROOT/results/metabolic_resilience/reproducibility"
mkdir -p "$OUT"
STAMP="$(date -u +%Y%m%dT%H%M%SZ)"
{
  echo "utc=$STAMP"
  echo "hostname=$(hostname)"
  echo "kernel=$(uname -a)"
  echo "python=$(python3 --version 2>&1)"
  echo "r=$(Rscript --version 2>&1 | head -n1)"
  echo "plink2=$(plink2 --version 2>&1 | head -n1)"
  echo "bcftools=$(bcftools --version 2>&1 | head -n1)"
  echo "git_head=$(git -C "$ROOT/IS_Analysis_V3" rev-parse HEAD 2>/dev/null)"
  echo "git_branch=$(git -C "$ROOT/IS_Analysis_V3" branch --show-current 2>/dev/null)"
} > "$OUT/ENVIRONMENT_$STAMP.txt"
python3 -m pip freeze > "$OUT/PIP_FREEZE_$STAMP.txt" 2>/dev/null
Rscript -e 'writeLines(capture.output(sessionInfo()))' > "$OUT/R_SESSIONINFO_$STAMP.txt" 2>/dev/null
find "$ROOT/results/metabolic_resilience" -maxdepth 3 -type f -printf '%p\t%s\n' 2>/dev/null | sort > "$OUT/RESULT_FILE_INVENTORY_$STAMP.tsv"
echo "[PASS] reproducibility snapshot -> $OUT"
