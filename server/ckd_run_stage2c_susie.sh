#!/usr/bin/env bash
set -euo pipefail

ROOT="${ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
STAGE2="${CKD_STAGE2_ROOT:-/srv/is-analysis/results/ckd/stage2}"
STAGE2B="${CKD_STAGE2B_ROOT:-/srv/is-analysis/results/ckd/stage2b_coloc}"
OUT="${CKD_STAGE2C_ROOT:-/srv/is-analysis/results/ckd/stage2c_susie}"
WORK="${CKD_STAGE2C_WORK_ROOT:-/srv/is-analysis/data/ckd/stage2c_ld}"
PYTHON="${CKD_PYTHON:-$ROOT/.venv-ckd/bin/python}"
R_LIB="${CKD_R_LIB:-/srv/is-analysis/.Rlib}"
PLINK2="${PLINK2:-plink2}"

test -x "$PYTHON"
test -f "$STAGE2/stage2_candidates.tsv"
test -d "$STAGE2B/coloc_input"
test -f "$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv"
command -v curl >/dev/null || { echo "curl is required" >&2; exit 7; }
command -v Rscript >/dev/null || { echo "Rscript is required" >&2; exit 7; }
command -v bcftools >/dev/null || {
  echo "bcftools is required. Ubuntu: sudo apt update && sudo apt install -y bcftools" >&2
  exit 8
}
command -v "$PLINK2" >/dev/null || {
  echo "plink2 is required. Ubuntu 24.04: sudo apt update && sudo apt install -y plink2" >&2
  exit 8
}

mkdir -p "$OUT" "$WORK" "$R_LIB"

cleanup_arg=()
if [[ "${STAGE2C_CLEANUP_CHR_CACHE:-1}" == "1" ]]; then
  cleanup_arg+=(--cleanup-chromosome-cache)
fi

"$PYTHON" "$ROOT/scripts/prepare_ckd_stage2c_ld.py" \
  --candidates "$STAGE2/stage2_candidates.tsv" \
  --stage2b-input "$STAGE2B/coloc_input" \
  --work-root "$WORK" \
  --output-root "$OUT" \
  --plink2 "$PLINK2" \
  --threads "${STAGE2C_THREADS:-2}" \
  --memory-mb "${STAGE2C_MEMORY_MB:-2500}" \
  "${cleanup_arg[@]}"

export R_LIBS_USER="$R_LIB"
Rscript -e 'pkgs<-c("coloc","susieR"); miss<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]; if(length(miss)) install.packages(miss,repos="https://cloud.r-project.org",lib=Sys.getenv("R_LIBS_USER"))'

mkdir -p "$OUT/susie_results"
Rscript "$ROOT/scripts/run_ckd_stage2c_susie.R" \
  "$OUT/susie_input" \
  "$WORK/ld" \
  "$STAGE2B/coloc_results/STAGE2B_COLOC_DEFAULT.tsv" \
  "$OUT/susie_results"

(
  cd "$OUT"
  find . -type f ! -name SHA256SUMS.txt -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS.txt
  sha256sum -c SHA256SUMS.txt
)

if [[ "${SYNC_DRIVE:-0}" == "1" ]]; then
  RCLONE_REMOTE="${RCLONE_REMOTE:-gdrive}"
  rclone copy "$OUT" "${RCLONE_REMOTE}:IS_Analysis_V3/results/ckd/stage2c_susie" \
    --checksum --transfers 1 --checkers 2 --progress
  echo "CKD_STAGE2C_DRIVE_SYNC_PASS"
fi

echo "CKD_STAGE2C_PASS"
echo "RESULT=$OUT"
echo "LD_WORK=$WORK"
