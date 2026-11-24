#!/usr/bin/env bash
set -euo pipefail
CODE_ROOT="${CODE_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
RESULTS_ROOT="${DASHBOARD_RESULTS_ROOT:-${IS_ANALYSIS_ROOT:-/srv/is-analysis}/results}"
DASH="${DASHBOARD_REPO:-/srv/is-analysis/brain_research_mr_scrna_seq}"

echo "===== MASTER DEGREE DASHBOARD SYNC ====="
echo "RESULTS=$RESULTS_ROOT"
echo "DASHBOARD=$DASH"

if [ ! -d "$DASH/.git" ]; then
  echo "Dashboard repository clone is missing."
  echo "Clone the private repo first, for example:"
  echo "  cd /srv/is-analysis"
  echo "  git clone git@github.com:kimtk94/brain_research_mr_scrna_seq.git"
  exit 2
fi

python3 "$CODE_ROOT/server/publish_research_dashboard.py" \
  --results-root "$RESULTS_ROOT" \
  --dashboard-repo "$DASH" \
  --git-push
