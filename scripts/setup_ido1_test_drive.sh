#!/usr/bin/env bash
set -euo pipefail

: "${CODE_ROOT:=/content/IS_Analysis_V3}"
: "${WORK_ROOT:=/content/drive/MyDrive/IS_Analysis_V3}"
: "${SOURCE_WORK_ROOT:=/content/drive/MyDrive/IS_Analysis_V2}"

case "$WORK_ROOT" in
  /content/drive/*) ;;
  *) echo "ERROR: WORK_ROOT must be under /content/drive: $WORK_ROOT" >&2; exit 2 ;;
esac

mkdir -p \
  "$WORK_ROOT/data/metadata" "$WORK_ROOT/data/reference" \
  "$WORK_ROOT/data/rawdata/pqtl/selected_targets/EUR" \
  "$WORK_ROOT/data/rawdata/pqtl/selected_targets/EAS" \
  "$WORK_ROOT/results/qc/ido1_test_pipeline" \
  "$WORK_ROOT/results/test/ido1/exposure_batches" \
  "$WORK_ROOT/results/test/ido1/standardized/pqtl" \
  "$WORK_ROOT/results/test/ido1/instrument_candidates"

copy_if_missing() {
  local relative=$1 source="$SOURCE_WORK_ROOT/$1" target="$WORK_ROOT/$1"
  if [[ -e "$target" ]]; then
    printf 'Preserving existing %s\n' "$target"
  elif [[ -f "$source" ]]; then
    cp "$source" "$target"
    printf 'Copied %s\n' "$relative"
  else
    printf 'WARNING: source not found: %s\n' "$source" >&2
  fi
}

copy_if_missing data/metadata/ukb_ppp_download_manifest.tsv
copy_if_missing data/reference/gene_coordinates_hg38.tsv

env_file="$WORK_ROOT/ido1_test.env"
cat >"$env_file" <<EOF
export CODE_ROOT=$(printf '%q' "$CODE_ROOT")
export WORK_ROOT=$(printf '%q' "$WORK_ROOT")
export SOURCE_WORK_ROOT=$(printf '%q' "$SOURCE_WORK_ROOT")
export IDO1_TEST_QC_DIR=$(printf '%q' "$WORK_ROOT/results/qc/ido1_test_pipeline")
export IDO1_TEST_OUTDIR=$(printf '%q' "$WORK_ROOT/results/test/ido1/exposure_batches")
export IDO1_TEST_STANDARDIZED_DIR=$(printf '%q' "$WORK_ROOT/results/test/ido1/standardized/pqtl")
export IDO1_TEST_INSTRUMENT_DIR=$(printf '%q' "$WORK_ROOT/results/test/ido1/instrument_candidates")
EOF
echo "Created $env_file"
