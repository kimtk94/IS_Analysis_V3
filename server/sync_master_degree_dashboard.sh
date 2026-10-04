#!/usr/bin/env bash
ROOT="/srv/is-analysis/IS_Analysis_V3"
DASH="${DASHBOARD_REPO:-/srv/is-analysis/brain_research_mr_scrna_seq}"

echo "===== MASTER DEGREE DASHBOARD SYNC ====="
echo "SOURCE=$ROOT"
echo "DASHBOARD=$DASH"

if [ ! -d "$DASH/.git" ]; then
  echo "Dashboard repository clone is missing."
  echo "Clone the private repo first, for example:"
  echo "  cd /srv/is-analysis"
  echo "  git clone git@github.com:kimtk94/brain_research_mr_scrna_seq.git"
  exit 2
fi

python3 "$ROOT/server/publish_research_dashboard.py" \
  --source "$ROOT" \
  --dashboard-repo "$DASH" \
  --git-push

RC=$?
echo "SYNC_RC=$RC"
exit $RC
