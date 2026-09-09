#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

ENV_FILE="${SERVER_ENV_FILE:-$REPO_ROOT/server/.env}"
if [[ ! -f "$ENV_FILE" ]]; then
  echo "Missing $ENV_FILE. Create it with: cp server/.env.example server/.env" >&2
  exit 1
fi

set -a
# shellcheck disable=SC1090
source "$ENV_FILE"
set +a

ROOT="${IS_ANALYSIS_ROOT:-/srv/is-analysis}"
SYNAPSE_ENV_FILE="${SYNAPSE_ENV_FILE:-$HOME/.config/is-analysis/synapse.env}"

if [[ "$POSTGRES_PASSWORD" == "CHANGE_ME_TO_A_LONG_RANDOM_PASSWORD" || -z "${POSTGRES_PASSWORD:-}" ]]; then
  echo "POSTGRES_PASSWORD is not configured in server/.env" >&2
  exit 1
fi

if [[ ! -f "$SYNAPSE_ENV_FILE" ]]; then
  echo "Missing Synapse secret file: $SYNAPSE_ENV_FILE" >&2
  exit 1
fi

# Reject group/world-readable secret files.
perm="$(stat -c '%a' "$SYNAPSE_ENV_FILE")"
if (( 10#$perm % 100 != 0 )); then
  echo "Synapse secret file must not be group/world accessible (current mode: $perm)" >&2
  exit 1
fi

set -a
# shellcheck disable=SC1090
source "$SYNAPSE_ENV_FILE"
set +a

if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]]; then
  echo "SYNAPSE_AUTH_TOKEN is empty" >&2
  exit 1
fi

command -v docker >/dev/null || { echo "docker is not installed" >&2; exit 1; }
docker compose version >/dev/null

echo "[1/5] Storage guard"
python3 scripts/storage_guard.py \
  --path "$ROOT" \
  --staging "$ROOT/data/staging" \
  --min-free-gib "${MIN_FREE_GIB:-40}" \
  --max-staging-gib "${MAX_STAGING_GIB:-35}"

echo "[2/5] Compose validation"
docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml config -q

echo "[3/5] PostgreSQL"
docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml up -d postgres
for _ in $(seq 1 30); do
  if docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml exec -T postgres \
      pg_isready -U "${POSTGRES_USER:-is_analysis}" -d "${POSTGRES_DB:-is_analysis}" >/dev/null 2>&1; then
    break
  fi
  sleep 2
done
docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml exec -T postgres \
  psql -v ON_ERROR_STOP=1 -U "${POSTGRES_USER:-is_analysis}" -d "${POSTGRES_DB:-is_analysis}" \
  -c 'SELECT 1 AS db_ok;'

echo "[4/5] Pipeline image"
docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml build pipeline

echo "[5/5] Synapse authentication + metadata access"
docker compose --env-file "$ENV_FILE" -f server/docker-compose.yml run --rm \
  -e SYNAPSE_AUTH_TOKEN pipeline \
  python3 -c 'import synapseclient; s=synapseclient.login(silent=True); e=s.get("syn52363617", downloadFile=False); print("SYNAPSE_OK", e.id, getattr(e, "name", ""))'

echo
echo "SERVER_SMOKE_PASS"
echo "No production raw data was downloaded."
