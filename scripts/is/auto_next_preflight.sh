#!/usr/bin/env bash
# Read-only preflight: does not install units or start AUTO_NEXT.
# Deliberately no set -e, safe on servers.
ROOT="${IS_REPO_ROOT:-/srv/is-analysis/IS_Analysis_V3}"
FAILED=0
echo "AUTO_NEXT PREFLIGHT (read only)"
echo "utc=$(date -u +%FT%TZ)"
if [ ! -f "$ROOT/scripts/is/auto_next.py" ]; then
  echo "FAIL: controller missing at $ROOT/scripts/is/auto_next.py"
  FAILED=1
else
  python3 -m py_compile "$ROOT/scripts/is/auto_next.py" || FAILED=1
fi
if [ ! -f "$ROOT/scripts/is/github_issue_queue_worker.py" ]; then
  echo "FAIL: existing worker unavailable"
  FAILED=1
fi
echo "Existing related user timers:"
systemctl --user list-timers --all --no-pager 2>&1 | grep -Ei 'is|harness|queue|auto.next' || true
echo "Existing related processes:"
pgrep -af '[g]ithub_issue_queue_worker.py|[a]uto_next.py|[h]arness_research_jobs.py' || true
echo "AUTO_NEXT file status:"
git -C "$ROOT" status --short -- scripts/is/auto_next.py scripts/is/github_issue_queue_worker.py 2>&1 || FAILED=1
if [ "$FAILED" -eq 0 ]; then
  echo "PREFLIGHT_STATIC_PASS: runtime conflict clearance remains MANUAL"
else
  echo "PREFLIGHT_FAILED"
fi
exit "$FAILED"
