#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
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
NETWORK="${IS_ANALYSIS_DOCKER_NETWORK:-is-analysis-net}"
POSTGRES_CONTAINER="${IS_ANALYSIS_POSTGRES_CONTAINER:-is-analysis-postgres}"
POSTGRES_IMAGE="${IS_ANALYSIS_POSTGRES_IMAGE:-postgres:16-alpine}"
PIPELINE_IMAGE="${IS_ANALYSIS_PIPELINE_IMAGE:-is-analysis-pipeline:local}"

if [[ -z "${POSTGRES_PASSWORD:-}" || "$POSTGRES_PASSWORD" == "CHANGE_ME_TO_A_LONG_RANDOM_PASSWORD" ]]; then
  echo "POSTGRES_PASSWORD is not configured in $ENV_FILE" >&2
  exit 1
fi

command -v docker >/dev/null || { echo "docker is not installed" >&2; exit 1; }
docker info >/dev/null 2>&1 || { echo "Docker daemon is not accessible by this user" >&2; exit 1; }

mkdir -p "$ROOT/state/postgres"

if ! docker network inspect "$NETWORK" >/dev/null 2>&1; then
  docker network create "$NETWORK" >/dev/null
fi

if docker container inspect "$POSTGRES_CONTAINER" >/dev/null 2>&1; then
  if [[ "$(docker inspect -f '{{.State.Running}}' "$POSTGRES_CONTAINER")" != "true" ]]; then
    docker start "$POSTGRES_CONTAINER" >/dev/null
  fi
else
  docker run -d \
    --name "$POSTGRES_CONTAINER" \
    --restart unless-stopped \
    --network "$NETWORK" \
    --cpus "0.50" \
    --memory "768m" \
    -e POSTGRES_DB="${POSTGRES_DB:-is_analysis}" \
    -e POSTGRES_USER="${POSTGRES_USER:-is_analysis}" \
    -e POSTGRES_PASSWORD="$POSTGRES_PASSWORD" \
    -v "$ROOT/state/postgres:/var/lib/postgresql/data" \
    -v "$REPO_ROOT/server/sql:/docker-entrypoint-initdb.d:ro" \
    "$POSTGRES_IMAGE" >/dev/null
fi

echo "Waiting for PostgreSQL..."
for _ in $(seq 1 45); do
  if docker exec "$POSTGRES_CONTAINER" \
      pg_isready -U "${POSTGRES_USER:-is_analysis}" -d "${POSTGRES_DB:-is_analysis}" >/dev/null 2>&1; then
    echo "PostgreSQL ready"
    break
  fi
  sleep 2
done

if ! docker exec "$POSTGRES_CONTAINER" \
    pg_isready -U "${POSTGRES_USER:-is_analysis}" -d "${POSTGRES_DB:-is_analysis}" >/dev/null 2>&1; then
  echo "PostgreSQL did not become ready" >&2
  docker logs --tail 80 "$POSTGRES_CONTAINER" >&2 || true
  exit 1
fi

if ! docker image inspect "$PIPELINE_IMAGE" >/dev/null 2>&1; then
  echo "Building pipeline image: $PIPELINE_IMAGE"
  docker build -t "$PIPELINE_IMAGE" -f "$REPO_ROOT/server/Dockerfile" "$REPO_ROOT"
fi

echo "PLAIN_DOCKER_STACK_READY"
echo "network=$NETWORK postgres=$POSTGRES_CONTAINER pipeline_image=$PIPELINE_IMAGE"
