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
NETWORK="${IS_ANALYSIS_DOCKER_NETWORK:-is-analysis-net}"
POSTGRES_CONTAINER="${IS_ANALYSIS_POSTGRES_CONTAINER:-is-analysis-postgres}"
PIPELINE_IMAGE="${IS_ANALYSIS_PIPELINE_IMAGE:-is-analysis-pipeline:local}"

if [[ "$POSTGRES_PASSWORD" == "CHANGE_ME_TO_A_LONG_RANDOM_PASSWORD" || -z "${POSTGRES_PASSWORD:-}" ]]; then
  echo "POSTGRES_PASSWORD is not configured in server/.env" >&2
  exit 1
fi

if [[ ! -f "$SYNAPSE_ENV_FILE" ]]; then
  echo "Missing Synapse secret file: $SYNAPSE_ENV_FILE" >&2
  exit 1
fi

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
docker info >/dev/null 2>&1 || { echo "Docker daemon is not accessible by this user" >&2; exit 1; }

echo "[1/4] Storage guard"
python3 scripts/storage_guard.py \
  --path "$ROOT" \
  --staging "$ROOT/data/staging" \
  --min-free-gib "${MIN_FREE_GIB:-40}" \
  --max-staging-gib "${MAX_STAGING_GIB:-35}"

echo "[2/4] Plain Docker stack"
bash server/start_stack.sh

echo "[3/4] PostgreSQL query"
docker exec "$POSTGRES_CONTAINER" \
  psql -v ON_ERROR_STOP=1 \
  -U "${POSTGRES_USER:-is_analysis}" \
  -d "${POSTGRES_DB:-is_analysis}" \
  -c 'SELECT 1 AS db_ok;'

echo "[4/4] Synapse authentication + metadata access"
docker run --rm \
  --network "$NETWORK" \
  --cpus "2.00" \
  --memory "3g" \
  -e SYNAPSE_AUTH_TOKEN \
  -e WORK_ROOT="$ROOT" \
  -e DATABASE_URL="postgresql://${POSTGRES_USER:-is_analysis}:${POSTGRES_PASSWORD}@${POSTGRES_CONTAINER}:5432/${POSTGRES_DB:-is_analysis}" \
  -e OMP_NUM_THREADS=2 \
  -e OPENBLAS_NUM_THREADS=2 \
  -e MKL_NUM_THREADS=2 \
  -v "$ROOT:$ROOT" \
  "$PIPELINE_IMAGE" \
  python3 -c 'import synapseclient; s=synapseclient.login(silent=True); e=s.get("syn52363617", downloadFile=False); print("SYNAPSE_OK", e.id, getattr(e, "name", ""))'

echo
echo "SERVER_SMOKE_PASS"
echo "No production raw data was downloaded."
